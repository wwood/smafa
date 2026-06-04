use needletail::parse_fastx_file;
use serde::{Deserialize, Serialize};

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::io::{Read, Write};
use std::num::{NonZeroU8, NonZeroUsize};
use std::path::{Path, PathBuf};
use std::time::Instant;
use std::{error::Error, fs::File};

use log::{debug, info};

mod cluster;
pub use cluster::cluster;

pub const AUTHOR_AND_EMAIL: &str =
    "Ben J. Woodcroft, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology <benjwoodcroft near gmail.com>";

pub const CURRENT_DB_VERSION: u32 = 2;
// NOTE: when this version is next bumped, cache the per-column conservation
// weights (see `column_weights`) in the new DB format. `query` and `cluster`
// currently recompute them at load to build the balanced band partition (~10 ms
// for ~137k subjects); storing them would make that free. See BANDING_BENCHMARK.md.

#[derive(Serialize, Deserialize, Debug, Clone)]
struct SeqEncoding(Vec<u64>);

#[derive(Clone)]
struct SeqEncodingLength {
    encoding: SeqEncoding,
    len: usize,
}

impl SeqEncodingLength {
    fn from_bytes(identifier: &[u8], seq: &[u8]) -> Self {
        // We encode chunks of 12 nucleotides to a u64
        let encoding = seq.chunks(12).enumerate().map(|(chunk_num, chunk)| {
            // The index i here is just to throw a good error message
            chunk.iter().enumerate().fold(0u64, |acc, (i, &byte)| {
                let b = encode_single(byte)
                    .unwrap_or_else(|| {
                        let seqname = String::from_utf8_lossy(identifier);
                        panic!(
                            "Byte {} cannot be interpreted as nucleotide, in sequence \"{}\" at position {}",
                            byte, seqname, 12 * chunk_num + i
                        )
                    })
                    .get();
                acc | ((b as u64) << (5 * i))
            })
        }).collect();
        Self {
            encoding: SeqEncoding(encoding),
            len: seq.len(),
        }
    }
}

#[derive(Serialize, Deserialize, Debug)]
struct WindowSet {
    version: u32,
    windows: Vec<SeqEncoding>,
    // None if there are no windows. Else, they must all be the same size
    len: Option<NonZeroUsize>,
}

impl WindowSet {
    fn new(version: u32) -> Self {
        WindowSet {
            version,
            windows: Vec::new(),
            len: None,
        }
    }

    fn get_distances(&self, seq: &SeqEncodingLength, distances: &mut [usize]) {
        if let Some(n) = self.len {
            if n.get() != seq.len {
                panic!(
                    "{}",
                    &format!("Cannot compute distances between seq of length {} and windows of lengths {}", seq.len, n.get())
                )
            }
        }
        for (window, distance) in self.windows.iter().zip(distances.iter_mut()) {
            *distance = window
                .0
                .iter()
                .zip(seq.encoding.0.iter())
                .map(|(a, b)| (a ^ b).count_ones() as usize)
                .sum::<usize>()
                / 2
        }
    }

    /// Serial nearest-centroid search. Returns `(distance, index)` of the
    /// lowest-indexed window with the smallest Hamming distance to `seq`, or
    /// `(usize::MAX, usize::MAX)` when there are no windows. The strict `<`
    /// comparison keeps the first (lowest-index) window on ties, matching the
    /// original greedy clustering behaviour.
    fn nearest(&self, seq: &SeqEncodingLength) -> (usize, usize) {
        if let Some(n) = self.len {
            if n.get() != seq.len {
                panic!(
                    "{}",
                    &format!("Cannot compute distances between seq of length {} and windows of lengths {}", seq.len, n.get())
                )
            }
        }
        let mut best = (usize::MAX, usize::MAX);
        for (i, window) in self.windows.iter().enumerate() {
            let d = window
                .0
                .iter()
                .zip(seq.encoding.0.iter())
                .map(|(a, b)| (a ^ b).count_ones() as usize)
                .sum::<usize>()
                / 2;
            if d < best.0 {
                best = (d, i);
            }
        }
        best
    }

    /// Move all windows from `other` into `self`, preserving order.
    fn merge(&mut self, mut other: WindowSet) {
        if other.windows.is_empty() {
            return;
        }
        match self.len {
            Some(n) => {
                if let Some(m) = other.len {
                    if n.get() != m.get() {
                        panic!(
                            "Cannot merge WindowSets with differing sequence lengths {} and {}",
                            n.get(),
                            m.get()
                        );
                    }
                }
            }
            None => self.len = other.len,
        }
        self.windows.append(&mut other.windows);
    }

    fn push_encoding(&mut self, encoding: SeqEncodingLength) {
        if let Some(n) = self.len {
            if n.get() != encoding.len {
                panic!(
                    "{}",
                    &format!(
                        "WindowSet seq length is {}, got a new sequence of length {}",
                        n, encoding.len
                    )
                )
            }
        } else {
            self.len = Some(
                encoding
                    .len
                    .try_into()
                    .expect("Cannot add empty sequence to WindowSet"),
            );
        }
        self.windows.push(encoding.encoding)
    }

    fn get_as_string(&self, index: usize) -> String {
        let uints = &self.windows[index].0;
        let v = (0..self.len.map(NonZeroUsize::get).unwrap_or(0))
            .map(|i| {
                let d = i / 12;
                let r = i % 12;
                let b = ((uints[d] >> (5 * r)) & 31) as u8;
                match b {
                    0b10000 => b'A',
                    0b01000 => b'C',
                    0b00100 => b'G',
                    0b00010 => b'T',
                    0b00001 => b'N',
                    _ => {
                        panic!("Invalid character in query sequence: {b}")
                    }
                }
            })
            .collect();
        // Safety: All the bytes above are ASCII, so it will never fail
        unsafe { String::from_utf8_unchecked(v) }
    }
}

/// Hamming distance between two equal-length packed encodings, in number of
/// differing nucleotide positions. Each mismatch flips two bits in the 5-bit
/// one-hot encoding, hence the `/ 2`.
#[inline]
pub(crate) fn hamming(a: &SeqEncoding, b: &SeqEncoding) -> usize {
    a.0.iter()
        .zip(b.0.iter())
        .map(|(x, y)| (x ^ y).count_ones() as usize)
        .sum::<usize>()
        / 2
}

/// A band is the set of nucleotide column indices it covers. For the original
/// contiguous bands this is just a range; allowing an arbitrary column list lets
/// a band be entropy-balanced (non-contiguous) while the pigeonhole guarantee —
/// which only needs the bands to be a disjoint cover — still holds.
type Band = Vec<usize>;

/// FNV-1a hash of the (canonical) 5-bit codes at the given nucleotide column
/// positions of an encoding. Identical band content (same columns, same order)
/// always yields an identical key; hash collisions only ever add extra
/// candidates (never drop true ones), so they cost time, not correctness.
///
/// FNV is order-dependent, so the same `cols` ordering must be used to index a
/// subject and to probe with a query — guaranteed here because both go through
/// the same `Partition`.
#[inline]
fn band_hash(enc: &SeqEncoding, cols: &[usize]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &i in cols {
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

/// The number of bands in one pigeonhole partition for divergence `d`: `d + 1`,
/// capped at `len` so every band is non-empty.
fn num_bands(d: usize, len: usize) -> usize {
    (d + 1).min(len.max(1))
}

/// The original contiguous partition: `d + 1` equal-width contiguous bands.
fn contiguous_bands(d: usize, len: usize) -> Vec<Band> {
    let nb = num_bands(d, len);
    (0..nb)
        .map(|i| (i * len / nb..(i + 1) * len / nb).collect())
        .collect()
}

/// Per-column collision weight `w = -ln(sum_i p_i^2)` over the subjects' allele
/// frequencies. A fully conserved column has `w ≈ 0`; a uniform/wobble column
/// has `w ≈ ln(4)`. Bands with equal summed weight have equal expected
/// background collision probability, which (by convexity of `exp(-x)`) minimises
/// the total candidate count.
fn column_weights(windows: &[SeqEncoding], len: usize) -> Vec<f64> {
    (0..len)
        .map(|c| {
            let mut counts = [0u64; 32];
            for w in windows {
                let code = ((w.0[c / 12] >> (5 * (c % 12))) & 0x1f) as usize;
                counts[code] += 1;
            }
            let total: u64 = counts.iter().sum();
            if total == 0 {
                return 0.0;
            }
            let m: f64 = counts
                .iter()
                .map(|&ct| {
                    let p = ct as f64 / total as f64;
                    p * p
                })
                .sum();
            -m.max(1e-12).ln()
        })
        .collect()
}

/// Entropy-balanced band assignment. Columns are ranked by weight (descending);
/// the column at rank `r` is round-robined to band `r % nb`, so every band draws
/// one column from each weight tier and ends up with near-equal summed weight —
/// no band is left all-conserved (which is what creates a giant candidate
/// bucket on coding data). Returns the band index for each column.
fn balanced_assignment(d: usize, len: usize, weights: &[f64]) -> Vec<usize> {
    let nb = num_bands(d, len);
    let mut order: Vec<usize> = (0..len).collect();
    order.sort_by(|&a, &b| {
        weights[b]
            .partial_cmp(&weights[a])
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.cmp(&b))
    });
    let mut slot = vec![0usize; len];
    for (rank, &col) in order.iter().enumerate() {
        slot[col] = rank % nb;
    }
    slot
}

fn slots_to_bands(slots: &[usize], nb: usize) -> Vec<Band> {
    let mut bands = vec![Vec::new(); nb];
    for (col, &s) in slots.iter().enumerate() {
        bands[s].push(col);
    }
    bands
}

/// Pigeonhole / banding index.
///
/// For a maximum Hamming divergence `d`, the sequence is split into `d + 1`
/// disjoint bands. Any two sequences within distance `d` must share at least one
/// identical band (d mismatches cannot touch all d+1 bands), so the candidate set
/// is the union of the query's band buckets — every within-`d` subject is in it.
/// The bands may be contiguous (`new`) or entropy-balanced (`single_balanced`);
/// balancing spreads conserved columns across bands so no single band's bucket
/// holds most of the database, which is the dominant cost on coding data.
///
/// Exact for Hamming distance: it never misses a within-`d` subject. Valid only
/// when `d + 1 <= len` (otherwise an all-differing pair could share no band);
/// callers use it only when `d + 1 < len` and fall back to a full scan
/// otherwise.
pub(crate) struct BandIndex {
    bands: Vec<Band>,
    maps: Vec<U64Map>,
}

impl BandIndex {
    fn from_bands(bands: Vec<Band>) -> Self {
        let maps = (0..bands.len()).map(|_| U64Map::default()).collect();
        BandIndex { bands, maps }
    }

    /// `d + 1` equal contiguous bands — the original scheme.
    pub(crate) fn new(max_divergence: usize, len: usize) -> Self {
        Self::from_bands(contiguous_bands(max_divergence, len))
    }

    /// A single entropy-balanced partition. Same number of bands and lookups as
    /// `new`, but the columns are balanced by conservation (computed from
    /// `windows`), which shrinks the candidate sets on coding data.
    pub(crate) fn single_balanced(
        max_divergence: usize,
        len: usize,
        windows: &[SeqEncoding],
    ) -> Self {
        let nb = num_bands(max_divergence, len);
        let w = column_weights(windows, len);
        let a = balanced_assignment(max_divergence, len, &w);
        Self::from_bands(slots_to_bands(&a, nb))
    }

    pub(crate) fn keys(&self, enc: &SeqEncoding) -> Vec<u64> {
        self.bands.iter().map(|b| band_hash(enc, b)).collect()
    }

    /// Subject indices sharing at least one band with the query, sorted ascending
    /// and deduplicated. Ascending order means a later strict-`<` scan keeps the
    /// lowest-index subject on distance ties.
    pub(crate) fn candidates_sorted(&self, keys: &[u64]) -> Vec<u32> {
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

    pub(crate) fn insert(&mut self, idx: u32, keys: &[u64]) {
        for (b, &k) in keys.iter().enumerate() {
            self.maps[b].entry(k).or_default().push(idx);
        }
    }
}

pub fn makedb(subject_fasta: &Path, db_path: &Path) -> Result<(), Box<dyn Error>> {
    // Iterate over lines, creating a vector of u8, where the lowest 5 bits of the u8
    // are the input nucleotides, one-hot encoded

    // Open the query file as a fasta file.
    debug!("Opening subject fasta file: {:?}", subject_fasta);
    let mut subject_reader =
        parse_fastx_file(subject_fasta).expect("valid path/file of subject fasta");

    info!("Encoding subject sequences ..");
    let mut windows = WindowSet::new(CURRENT_DB_VERSION);
    while let Some(record) = subject_reader.next() {
        let record = record.expect("valid record");
        let encoded = SeqEncodingLength::from_bytes(record.id(), &record.seq());
        windows.push_encoding(encoded);
    }

    info!(
        "Encoding of {} sequences complete, writing db file {}",
        windows.windows.len(),
        db_path.to_string_lossy()
    );

    // Encode
    let mut ferris_file = File::create(db_path)?;
    ferris_file.write_all(&postcard::to_allocvec(&windows).unwrap())?;
    info!("DB file written");
    Ok(())
}

const fn create_lut() -> [u8; 256] {
    let mut lut = [0; 256];
    let mut i = 0;
    while i < 256 {
        let b = match i as u8 {
            b'A' | b'a' => 0b10000,
            b'C' | b'c' => 0b01000,
            b'G' | b'g' => 0b00100,
            b'T' | b't' | b'U' | b'u' => 0b00010,
            b'N' | b'W' | b'S' | b'M' | b'K' | b'R' | b'Y' | b'B' | b'D' | b'H' | b'V' | b'-'
            | b'n' | b'w' | b's' | b'm' | b'k' | b'r' | b'y' | b'b' | b'd' | b'h' | b'v' => 0b00001,
            _ => 0,
        };
        lut[i] = b;
        i += 1;
    }
    lut
}

const BYTE_LUT: [u8; 256] = create_lut();

// inline this function, performance affects untested, guessing it's better
#[inline(always)]
fn encode_single(c: u8) -> Option<NonZeroU8> {
    let lut: [u8; 256] = BYTE_LUT; // statically verify lut has 256 elements

    // Safety: We just verified it has indices 0-255, so a u8 can't be out of bounds
    let encoding = unsafe { *lut.get_unchecked(c as usize) };
    NonZeroU8::new(encoding)
}

pub fn query(
    db_path: &Path,
    query_fasta: &Path,
    max_divergence: Option<u32>,
    max_num_hits: Option<u32>,
    limit_per_sequence: Option<u32>,
    no_banding: bool,
) -> Result<(), Box<dyn Error>> {
    // Decode
    info!("Decoding db file {:?}", db_path);
    let start = Instant::now();
    let mut ferris_file = File::open(db_path)?;
    let mut buffer = Vec::new();
    ferris_file.read_to_end(&mut buffer)?;

    // Check that the version of the db file is the most recent. We do not
    // support backwards compatibility.
    let version: u32 = postcard::from_bytes(&buffer[0..4])?;
    if version != CURRENT_DB_VERSION {
        panic!("Unsupported db file version: {}. This version of smafa only works with version {} databases. The last version to support version 1 databases was v0.7.1.", version, CURRENT_DB_VERSION);
    }
    let windows: WindowSet = postcard::from_bytes(&buffer)?;

    // Open the query file as a fasta file.
    let mut query_reader = parse_fastx_file(query_fasta).expect("valid path/file of query fasta");

    // 1 is a special case, it is equivalent to None.
    let max_divergence_for_match = max_num_hits.filter(|&max_num_hits| max_num_hits != 1);

    // Banding (the pigeonhole prefilter) can only prune when there is a
    // divergence bound: without one, the nearest hit could be arbitrarily far
    // and share no band, so we must scan everything. It is also only used when
    // d + 1 < window length (a tighter bound makes the bands degenerate and the
    // candidate sets explode), and can be disabled explicitly for large d.
    let window_len = windows.len.map(NonZeroUsize::get).unwrap_or(0);
    let band_index: Option<BandIndex> = match (no_banding, max_divergence) {
        (false, Some(d)) if !windows.windows.is_empty() && (d as usize) + 1 < window_len => {
            let d = d as usize;
            info!(
                "Building balanced band index over {} subjects ({} bands) ..",
                windows.windows.len(),
                d + 1
            );
            // Entropy-balanced partition: column conservation is recomputed from
            // the subjects at load (~10 ms for ~137k subjects), so no DB-format
            // change is needed. Exact like any d+1-band partition.
            let mut bi = BandIndex::single_balanced(d, window_len, &windows.windows);
            for (i, w) in windows.windows.iter().enumerate() {
                let keys = bi.keys(w);
                bi.insert(i as u32, &keys);
            }
            Some(bi)
        }
        (false, Some(d)) if (d as usize) + 1 >= window_len && window_len > 0 => {
            info!("Banding disabled (max-divergence too large for the window length); scanning all subjects.");
            None
        }
        _ => None,
    };

    // Pre-initialise the distances vector so don't have to continually
    // reallocate. Only needed for the full-scan path.
    let mut distances = if band_index.is_none() {
        vec![0; windows.windows.len()]
    } else {
        Vec::new()
    };

    // Iterate over the query file.
    info!("Querying ..");
    let mut query_number: u32 = 0;
    while let Some(record) = query_reader.next() {
        // encode a line from stdin as a vector of bools
        let record = record.expect("Failed to parse query sequence");
        let query_vec = SeqEncodingLength::from_bytes(record.id(), &record.seq());

        if let Some(bi) = band_index.as_ref() {
            // ---- Banded path (only reached when max_divergence is Some) ----
            let max_div = max_divergence.unwrap() as usize;
            let keys = bi.keys(&query_vec.encoding);
            let mut cand: Vec<(usize, usize)> = bi
                .candidates_sorted(&keys)
                .into_iter()
                .map(|i| {
                    (
                        hamming(&windows.windows[i as usize], &query_vec.encoding),
                        i as usize,
                    )
                })
                .collect();
            // Sort by (distance, index): the within-d subset and its ordering
            // then match the full-scan path exactly, including the lowest-index
            // tie-break.
            cand.sort_unstable();

            match max_divergence_for_match {
                Some(max_num_hits) => {
                    // k-th smallest distance among candidates. Output is gated by
                    // `<= max_div`, and every within-d subject is a candidate, so
                    // this reproduces the full-scan top-k-within-d result.
                    let max_distance = if max_num_hits > cand.len() as u32 {
                        cand.iter().map(|(d, _)| *d).max().unwrap_or(0)
                    } else {
                        cand[(max_num_hits - 1) as usize].0
                    };

                    let mut last_sequence: Option<(String, u32)> = None;
                    let mut new_last_sequence: Option<(String, u32)>;
                    for (distance, i) in cand.iter() {
                        if *distance <= max_distance && *distance <= max_div {
                            let s = windows.get_as_string(*i);
                            debug!("Found hit sequence {} at distance {}", s, distance);

                            if let Some(limit_per_sequence_unwrapped) = limit_per_sequence {
                                match &last_sequence {
                                    Some((last_seq, last_seq_count)) if last_seq == &s => {
                                        if last_seq_count >= &limit_per_sequence_unwrapped {
                                            continue;
                                        } else {
                                            new_last_sequence =
                                                Some((s.clone(), last_seq_count + 1));
                                        }
                                    }
                                    _ => {
                                        new_last_sequence = Some((s.clone(), 1));
                                    }
                                }
                                last_sequence = new_last_sequence;
                            }

                            println!("{}\t{}\t{}\t{}", query_number, i, distance, s);
                        }
                    }
                }
                None => {
                    if limit_per_sequence.is_some() {
                        panic!("limit_per_sequence is implemented unless max_num_hits > 1. It can be implemented by analogy, just haven't gotten around to it.");
                    }
                    // Closest candidate; print all at that distance if within d.
                    if let Some(min_distance) = cand.iter().map(|(d, _)| *d).min() {
                        if min_distance <= max_div {
                            for (distance, i) in cand.iter() {
                                if *distance == min_distance {
                                    let s = windows.get_as_string(*i);
                                    println!("{}\t{}\t{}\t{}", query_number, i, distance, s);
                                }
                            }
                        }
                    }
                }
            }

            query_number += 1;
            continue;
        }

        // ---- Full-scan path (no divergence bound, or banding disabled) ----
        // Get the minimum distance between the query and each window using xor.
        windows.get_distances(&query_vec, &mut distances);

        // Find the max_num_hits'th minimum distance.
        match max_divergence_for_match {
            Some(max_num_hits) => {
                let mut min_distances = distances
                    .iter()
                    .enumerate()
                    .map(|(i, d)| (*d, i))
                    .collect::<Vec<_>>();
                // There might be a faster way of doing this using a priority
                // queue, but eh for now unless it really is slow.
                min_distances.sort();

                // If max num hits is greater than the number of windows, just print them all.
                let max_distance = match max_num_hits > min_distances.len() as u32 {
                    true => *distances.iter().max().unwrap(),
                    false => min_distances[(max_num_hits - 1) as usize].0,
                };

                // Print out the windows that qualify in order of increasing distance.
                let mut last_sequence: Option<(String, u32)> = None;
                let mut new_last_sequence: Option<(String, u32)>; // to get around borrow checker
                for (distance, i) in min_distances.iter() {
                    if *distance <= max_distance
                        && (max_divergence.is_none()
                            || *distance <= max_divergence.unwrap() as usize)
                    {
                        let s = windows.get_as_string(*i);
                        debug!("Found hit sequence {} at distance {}", s, distance);

                        if let Some(limit_per_sequence_unwrapped) = limit_per_sequence {
                            // limit per sequence
                            match &last_sequence {
                                Some((last_seq, last_seq_count)) if last_seq == &s => {
                                    if last_seq_count >= &limit_per_sequence_unwrapped {
                                        continue;
                                    } else {
                                        new_last_sequence = Some((s.clone(), last_seq_count + 1));
                                    }
                                }
                                _ => {
                                    new_last_sequence = Some((s.clone(), 1));
                                }
                            }
                            last_sequence = new_last_sequence;
                        }

                        // Print the window if we make it here.
                        println!("{}\t{}\t{}\t{}", query_number, i, distance, s);
                    }
                }
            }
            None => {
                // Find the minimum distance.
                let min_distance = distances.iter().min().unwrap();
                debug!("Min distance: {}", min_distance);

                if limit_per_sequence.is_some() {
                    panic!("limit_per_sequence is implemented unless max_num_hits > 1. It can be implemented by analogy, just haven't gotten around to it.");
                }

                // Print the windows with the minimum distance.
                if max_divergence.is_none() || *min_distance <= max_divergence.unwrap() as usize {
                    for (i, distance) in distances.iter().enumerate() {
                        if distance == min_distance {
                            let s = windows.get_as_string(i);
                            println!("{}\t{}\t{}\t{}", query_number, i, distance, s);
                        }
                    }
                }
            }
        }

        query_number += 1;
    }

    info!(
        "Querying complete, took {} seconds",
        start.elapsed().as_secs()
    );
    Ok(())
}

/// Load a v2 DB file into a `WindowSet`.
fn load_windows(db_path: &Path) -> Result<WindowSet, Box<dyn Error>> {
    let mut ferris_file = File::open(db_path)?;
    let mut buffer = Vec::new();
    ferris_file.read_to_end(&mut buffer)?;
    let version: u32 = postcard::from_bytes(&buffer[0..4])?;
    if version != CURRENT_DB_VERSION {
        panic!("Unsupported db file version: {version}");
    }
    Ok(postcard::from_bytes(&buffer)?)
}

/// One banding strategy under test in `bench_banding`.
struct BenchMethod {
    label: &'static str,
    index: BandIndex,
    build_secs: f64,
    cand_total: u128,
    cand_max: usize,
    scan_secs: f64,
}

/// Benchmark the query-side banding strategies against each other on a real DB
/// and query set: the single contiguous band (the current default) versus the
/// single entropy-balanced band. For each divergence it reports the candidate-set
/// size (the work that decides query speed) and the wall time, and asserts both
/// strategies return the identical set of within-`d` hits — so
/// the only thing that differs is how many false candidates each scans.
pub fn bench_banding(
    db_path: &Path,
    query_fasta: &Path,
    divergences: &[u32],
    max_queries: Option<usize>,
) -> Result<(), Box<dyn Error>> {
    info!("Loading DB {db_path:?} ..");
    let windows = load_windows(db_path)?;
    let len = windows.len.map(NonZeroUsize::get).unwrap_or(0);
    let subjects = &windows.windows;
    info!("DB has {} subjects of length {len}", subjects.len());

    info!("Loading queries {query_fasta:?} ..");
    let mut query_reader = parse_fastx_file(query_fasta)?;
    let mut queries: Vec<SeqEncoding> = Vec::new();
    while let Some(record) = query_reader.next() {
        let record = record?;
        let enc = SeqEncodingLength::from_bytes(record.id(), &record.seq());
        if enc.len != len {
            panic!("Query length {} != subject length {len}", enc.len);
        }
        queries.push(enc.encoding);
        if let Some(m) = max_queries {
            if queries.len() >= m {
                break;
            }
        }
    }
    let q = queries.len();
    info!("Loaded {q} queries");

    // How long does the entropy (per-column conservation) computation take? It
    // depends only on the subjects, so it is paid once per query invocation if
    // recomputed at load rather than stored in the DB.
    {
        let reps = 20;
        let t = Instant::now();
        let mut sink = 0.0f64;
        for _ in 0..reps {
            let w = column_weights(subjects, len);
            sink += w.iter().sum::<f64>();
        }
        std::hint::black_box(sink);
        println!(
            "# entropy (column_weights) over {} subjects: {:.3} ms/call",
            subjects.len(),
            t.elapsed().as_secs_f64() * 1000.0 / reps as f64
        );
    }

    println!("# subjects={} len={} queries={}", subjects.len(), len, q);
    println!(
        "{:<14} {:>11} {:>11} {:>9} {:>11} {:>10} {:>9}",
        "method", "entries/sub", "build_s", "mean_cand", "total_cand", "max_cand", "scan_s"
    );

    for &dd in divergences {
        let d = dd as usize;
        if d + 1 >= len {
            println!("# d={dd}: skipped (d+1 >= len, banding degenerate)");
            continue;
        }
        let nb = num_bands(d, len);

        // Build each index, timing the build (insert of every subject).
        let build = |label: &'static str, mut index: BandIndex| -> BenchMethod {
            let t = Instant::now();
            for (i, w) in subjects.iter().enumerate() {
                let keys = index.keys(w);
                index.insert(i as u32, &keys);
            }
            BenchMethod {
                label,
                index,
                build_secs: t.elapsed().as_secs_f64(),
                cand_total: 0,
                cand_max: 0,
                scan_secs: 0.0,
            }
        };

        let mut methods = vec![
            build("contiguous", BandIndex::new(d, len)),
            build("balanced", BandIndex::single_balanced(d, len, subjects)),
        ];

        // One pass over the queries; for each query every method computes its
        // candidate set and its within-d hit set. The hit sets must match.
        let mut mismatches = 0u64;
        for query in &queries {
            let mut reference_hits: Option<Vec<u32>> = None;
            for m in methods.iter_mut() {
                let t = Instant::now();
                let keys = m.index.keys(query);
                let cand = m.index.candidates_sorted(&keys);
                let hits: Vec<u32> = cand
                    .iter()
                    .copied()
                    .filter(|&i| hamming(&subjects[i as usize], query) <= d)
                    .collect();
                m.scan_secs += t.elapsed().as_secs_f64();
                m.cand_total += cand.len() as u128;
                m.cand_max = m.cand_max.max(cand.len());
                match &reference_hits {
                    None => reference_hits = Some(hits),
                    Some(r) => {
                        if &hits != r {
                            mismatches += 1;
                        }
                    }
                }
            }
        }

        println!("# d={dd}  bands={nb}");
        for m in &methods {
            println!(
                "{:<14} {:>11} {:>11.3} {:>9.2} {:>11} {:>10} {:>9.3}",
                m.label,
                nb,
                m.build_secs,
                m.cand_total as f64 / q as f64,
                m.cand_total,
                m.cand_max,
                m.scan_secs,
            );
        }
        if mismatches == 0 {
            println!("# d={dd}: all methods returned identical within-d hits (OK)");
        } else {
            println!("# d={dd}: *** {mismatches} queries had differing hit sets ***");
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;

    #[test]
    fn test_makedb() {
        // Create a temporary directory to store the test DB file.
        let temp_dir = tempfile::tempdir().unwrap();
        let db_path = temp_dir.path().join("test.db");
        let subject_fasta = std::path::PathBuf::from_str("tests/data/subjects.fa").unwrap();

        // Call the makedb function with the test subject FASTA file and the path
        // to the test DB file.
        assert!(makedb(&subject_fasta, &db_path).is_ok());

        // Check that the DB file exists.
        assert!(db_path.exists());

        // Open the DB file and decode it to a WindowSet struct.
        let mut ferris_file = File::open(&db_path).unwrap();
        let mut encoded = Vec::new();
        ferris_file.read_to_end(&mut encoded).unwrap();
        let windows = postcard::from_bytes::<WindowSet>(&encoded).unwrap();

        // Check that the WindowSet struct has the expected number of sequences.
        assert_eq!(windows.windows.len(), 5);

        // Check that the first sequence has the expected one-hot encoded values.
        let expected_encoded = vec![
            vec![0b10000],
            vec![0b01000],
            vec![0b00100],
            vec![0b00010],
            vec![0b00001],
        ];
        for (i, j) in expected_encoded.iter().zip(windows.windows.iter()) {
            assert_eq!(i, &j.0)
        }
    }
}

// Derive IntoJson
#[derive(Serialize, Deserialize, Debug)]
struct CountResult {
    path: PathBuf,
    num_reads: usize,
    num_bases: usize,
}

pub fn count<T: Iterator<Item = P>, P: AsRef<Path>>(paths: T) -> Result<(), Box<dyn Error>> {
    let mut results = Vec::new();
    for path in paths {
        let mut reader = parse_fastx_file(&path)?;
        let mut read_count = 0;
        let mut bases_count = 0;
        while let Some(record) = reader.next() {
            let record = record?;
            read_count += 1;
            bases_count += record.seq().len();
        }
        results.push(CountResult {
            path: path.as_ref().to_owned(),
            num_reads: read_count,
            num_bases: bases_count,
        });
    }
    // Print output in JSON format including input path
    println!("{}", serde_json::to_string(&results).unwrap());
    Ok(())
}
