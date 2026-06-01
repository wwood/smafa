use needletail::parse_fastx_file;
use rayon::prelude::*;

use std::collections::HashMap;
use std::collections::HashSet;
use std::error::Error;
use std::hash::{BuildHasherDefault, Hasher};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use log::info;

use crate::{SeqEncoding, SeqEncodingLength, WindowSet};

/// Number of unique input sequences processed per parallel block. Larger blocks
/// expose more parallelism in phase 1, but increase the serial phase-2
/// reconciliation cost, which is bounded by O(BLOCK_SIZE^2) per block.
const BLOCK_SIZE: usize = 8192;

/// Hamming distance between two equal-length packed encodings, in number of
/// differing nucleotide positions. Each mismatch flips two bits in the 5-bit
/// one-hot encoding, hence the `/ 2`.
#[inline]
fn hamming(a: &SeqEncoding, b: &SeqEncoding) -> usize {
    a.0.iter()
        .zip(b.0.iter())
        .map(|(x, y)| (x ^ y).count_ones() as usize)
        .sum::<usize>()
        / 2
}

/// FNV-1a hash of the (canonical) 5-bit codes covering nucleotide positions
/// `[start, end)` of an encoding. Identical band content always yields an
/// identical key; hash collisions only ever add extra candidates (never drop
/// true ones), so they cost time, not correctness.
#[inline]
fn band_hash(enc: &SeqEncoding, start: usize, end: usize) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for i in start..end {
        let code = (enc.0[i / 12] >> (5 * (i % 12))) & 0x1f;
        h ^= code;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// Identity hasher for band keys: the keys are already well-mixed FNV hashes,
/// so we avoid SipHash and bucket directly on the u64.
#[derive(Default)]
struct U64Hasher(u64);
impl Hasher for U64Hasher {
    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }
    #[inline]
    fn write_u64(&mut self, i: u64) {
        self.0 = i;
    }
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = self.0.rotate_left(8) ^ b as u64;
        }
    }
}
type U64Map = HashMap<u64, Vec<u32>, BuildHasherDefault<U64Hasher>>;

/// Pigeonhole / banding index over centroids.
///
/// For a maximum Hamming divergence `d`, the sequence is split into `d + 1`
/// disjoint contiguous bands. Any two sequences within distance `d` must have
/// at least one band that is identical (d mismatches cannot touch all d+1
/// bands), so a centroid within `d` of a query is guaranteed to share a band
/// with it. The index therefore lets us restrict the (expensive) full-distance
/// comparisons to centroids that share at least one band, turning the greedy
/// scan from O(K) per query towards near-constant when clusters are dense.
///
/// This is exact for Hamming distance: it never misses a centroid within `d`.
/// (Used only when `d + 1 < len`; otherwise a single all-differing pair could
/// share no band, so the caller falls back to a full scan.)
struct BandIndex {
    bounds: Vec<(usize, usize)>,
    maps: Vec<U64Map>,
}

impl BandIndex {
    fn new(max_divergence: usize, len: usize) -> Self {
        let num_bands = (max_divergence + 1).min(len.max(1));
        let bounds = (0..num_bands)
            .map(|i| (i * len / num_bands, (i + 1) * len / num_bands))
            .collect();
        let maps = (0..num_bands).map(|_| U64Map::default()).collect();
        BandIndex { bounds, maps }
    }

    fn keys(&self, enc: &SeqEncoding) -> Vec<u64> {
        self.bounds
            .iter()
            .map(|&(s, e)| band_hash(enc, s, e))
            .collect()
    }

    /// Centroid indices that share at least one band with the given keys,
    /// sorted ascending and deduplicated. Ascending order means a later
    /// strict-`<` scan keeps the lowest-index centroid on distance ties,
    /// matching the original greedy behaviour.
    fn candidates_sorted(&self, keys: &[u64]) -> Vec<u32> {
        let mut v = Vec::new();
        for (b, &k) in keys.iter().enumerate() {
            if let Some(list) = self.maps[b].get(&k) {
                v.extend_from_slice(list);
            }
        }
        v.sort_unstable();
        v.dedup();
        v
    }

    fn insert(&mut self, idx: u32, keys: &[u64]) {
        for (b, &k) in keys.iter().enumerate() {
            self.maps[b].entry(k).or_default().push(idx);
        }
    }
}

/// Greedy single-linkage-style centroid clustering, parallelised by block
/// (phase 1 across sequences) and, when divergence is small relative to the
/// sequence length, accelerated with a pigeonhole/banding index so that each
/// sequence is only compared against plausible centroids rather than all of
/// them.
///
/// Output is identical to the original serial implementation: the same
/// lowest-index-at-minimum-distance assignment, the same new-centroid ordering
/// and the same exact-duplicate filtering.
pub fn cluster(
    input_fasta: &Path,
    max_divergence: u32,
    print_stream: &mut dyn std::io::Write,
) -> Result<(), Box<dyn Error>> {
    let start = Instant::now();
    let max_divergence_usize = max_divergence as usize;

    // All accepted centroids, in creation (= input) order.
    let mut centroids = WindowSet::new(0);

    // Exact-duplicate filter: identical encodings are collapsed and not emitted.
    let mut seen_sequences = HashSet::<Vec<u64>>::new();

    let mut query_reader = parse_fastx_file(input_fasta).expect("valid path/file of input fasta");

    // Buffer output: the original wrote each line straight to an unbuffered
    // stdout, a large syscall overhead at tens of millions of lines.
    let mut out = BufWriter::new(print_stream);

    info!("Clustering ..");
    let mut query_number: u64 = 0; // total records read (including duplicates)
    let mut n_unique: u64 = 0;

    let mut block: Vec<(Vec<u8>, SeqEncodingLength)> = Vec::with_capacity(BLOCK_SIZE);

    // Banding strategy is decided once the sequence length is known (all
    // sequences must be the same length). `band_index` is `None` when banding
    // is not used, in which case phase 1 falls back to a full centroid scan.
    let mut initialized = false;
    let mut band_index: Option<BandIndex> = None;

    let mut reader_done = false;
    while !reader_done {
        // ---- Fill a block with up to BLOCK_SIZE unique sequences ----
        block.clear();
        while block.len() < BLOCK_SIZE {
            match query_reader.next() {
                Some(record) => {
                    query_number += 1;
                    let record = record.expect("Failed to parse input sequence");
                    let seq = record.seq();
                    let enc = SeqEncodingLength::from_bytes(record.id(), &seq);
                    if !seen_sequences.insert(enc.encoding.0.clone()) {
                        continue;
                    }
                    block.push((seq.to_vec(), enc));
                }
                None => {
                    reader_done = true;
                    break;
                }
            }
        }
        if block.is_empty() {
            break;
        }
        n_unique += block.len() as u64;

        if !initialized {
            let len = block[0].1.len;
            // Banding needs d + 1 < len to guarantee a shared band for every
            // within-d pair; otherwise fall back to full scans.
            if max_divergence_usize + 1 < len {
                band_index = Some(BandIndex::new(max_divergence_usize, len));
            }
            initialized = true;
        }

        // ---- Phase 1 (parallel): nearest pre-block centroid for each member ----
        // `centroids` and `band_index` are immutable here, so each member is
        // scored independently across threads. Returns the member's band keys
        // (for later insertion) alongside its best pre-block hit.
        let phase1: Vec<(Vec<u64>, (usize, usize))> = {
            let centroids_ref = &centroids;
            let bi_ref = band_index.as_ref();
            block
                .par_iter()
                .with_min_len(64)
                .map(|(_, enc)| match bi_ref {
                    Some(bi) => {
                        let keys = bi.keys(&enc.encoding);
                        let mut best = (usize::MAX, usize::MAX);
                        for idx in bi.candidates_sorted(&keys) {
                            let d = hamming(&centroids_ref.windows[idx as usize], &enc.encoding);
                            if d < best.0 {
                                best = (d, idx as usize);
                            }
                        }
                        (keys, best)
                    }
                    None => (Vec::new(), centroids_ref.nearest(enc)),
                })
                .collect()
        };

        // ---- Phase 2 (serial): reconcile within the block, in input order ----
        // Members may match a centroid minted earlier in this same block; that
        // set (`block_new`) is small, so a direct scan is fine.
        let base = centroids.windows.len();
        let mut block_new = WindowSet::new(0);
        let mut assignments: Vec<usize> = Vec::with_capacity(block.len());
        let mut minted: Vec<usize> = Vec::new(); // member indices that became centroids

        for (j, (_, enc)) in block.iter().enumerate() {
            let cand_existing = phase1[j].1; // (dist, global index < base) or (MAX, MAX)

            let cand_new = if block_new.windows.is_empty() {
                (usize::MAX, usize::MAX)
            } else {
                let (d, i) = block_new.nearest(enc);
                (d, base + i)
            };

            // Pre-block indices < base <= any in-block index, so a lexicographic
            // min reproduces the lowest-index tie-break across both sets.
            let (best_d, best_idx) = std::cmp::min(cand_existing, cand_new);

            if best_d <= max_divergence_usize {
                assignments.push(best_idx);
            } else {
                let local = block_new.windows.len();
                block_new.push_encoding(enc.clone());
                assignments.push(base + local);
                minted.push(j);
            }
        }

        // Add this block's new centroids to the band index (global indices are
        // assigned in mint order: the k-th minted centroid is `base + k`).
        if let Some(bi) = band_index.as_mut() {
            for (k, &j) in minted.iter().enumerate() {
                bi.insert((base + k) as u32, &phase1[j].0);
            }
        }

        // Make the new centroids globally visible before rendering output.
        centroids.merge(block_new);

        // ---- Emit output for the block, in input order ----
        for ((seq_bytes, _), assigned) in block.iter().zip(assignments.iter()) {
            writeln!(
                out,
                "{}\t{}",
                std::str::from_utf8(seq_bytes).unwrap(),
                centroids.get_as_string(*assigned)
            )?;
        }
    }

    out.flush()?;

    info!(
        "Clustering complete, took {} seconds. Clustered {} sequences ({} unique) into {} clusters.",
        start.elapsed().as_secs(),
        query_number,
        n_unique,
        centroids.windows.len()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_simple() {
        let mut stream = Cursor::new(Vec::new());
        cluster(Path::new("tests/data/cluster_dummy1.fna"), 1, &mut stream).unwrap();
        assert_eq!(
            "ATGC\tATGC
ATGG\tATGC
AAAA\tAAAA
",
            std::str::from_utf8(stream.get_ref()).unwrap()
        )
    }

    #[test]
    fn test_bug1() {
        let mut stream = Cursor::new(Vec::new());
        cluster(Path::new("tests/data/cluster_bug1.fna"), 2, &mut stream).unwrap();
        assert_eq!(
            "ATGCAAAAA\tATGCAAAAA\n\
             ATAAAAAAA\tATGCAAAAA\n\
             TTAAAAAAA\tTTAAAAAAA\n",
            std::str::from_utf8(stream.get_ref()).unwrap()
        )
    }

    #[test]
    fn test_best_hit_changes_bug() {
        // seq4 in the file shouldn't be reported otherwise there are two
        // sequences that are the same but are given different centroids.
        let mut stream = Cursor::new(Vec::new());
        cluster(
            Path::new("tests/data/cluster_best_hit_changes.fna"),
            2,
            &mut stream,
        )
        .unwrap();
        assert_eq!(
            "ATGCAAAAA\tATGCAAAAA\n\
             ATAAAAAAA\tATGCAAAAA\n\
             TTAAAAAAA\tTTAAAAAAA\n",
            std::str::from_utf8(stream.get_ref()).unwrap()
        )
    }
}
