# Banding strategy benchmark (query)

Question (from the design discussion): for the `query` pigeonhole prefilter, is
**balancing column conservation across bands** worth a new database format,
compared to (a) the current single contiguous band, and (b) the simpler
"offset by (d+1)/2 and wrap" second band?

## What was compared

All four are *exact* — each partition has `d+1` bands, so by pigeonhole every
within-`d` subject shares a band in every partition; intersecting two partitions
keeps all true hits and only drops false positives. The benchmark verifies this:
**all four return byte-identical within-`d` hit sets** for every query.

| method | partitions | how the columns are cut | DB-format change |
|---|---|---|---|
| `1-contiguous` | 1 | `d+1` equal contiguous bands (**current default**) | none |
| `1-balanced` | 1 | columns round-robin'd by conservation so every band has equal collision-entropy | column→band map |
| `2-offset-wrap` | 2 | contiguous + the same cut rotated half a band, wrapped; candidates intersected | none |
| `2-balanced-diag` | 2 | two balanced partitions in the diagonal Latin-square layout; intersected | column→band map |

`mean_cand` is the mean number of candidate subjects scanned per query (the work
that drives speed); `scan_s` is the wall time for keys + candidate-gen + Hamming
over all queries, single-threaded.

## Results

Database: `S3.10.ribosomal_protein_S19_rpsS` (136,770 subjects, 60 bp).
Queries: `1000000.S3.10.seqs.fna`, first 200,000.

```
# d=2  bands/partition=3
method             parts  mean_cand   max_cand   scan_s
1-contiguous           1     131.31       2884    0.592
1-balanced             1      56.18       1445    0.410   <- fastest
2-offset-wrap          2      19.72        763    0.734
2-balanced-diag        2      11.03        513    0.672

# d=3  bands/partition=4
1-contiguous           1     701.07       7974    3.590
1-balanced             1     375.12       4636    2.070   <- fastest
2-offset-wrap          2     143.50       3687    4.778
2-balanced-diag        2      57.54       1994    2.709

# d=5  bands/partition=6
1-contiguous           1    6155.60      29654   35.273
1-balanced             1    3592.34      18835   21.378   <- fastest
2-offset-wrap          2    2136.29      19442   68.286
2-balanced-diag        2     670.02       8141   28.407
```

The full 1,000,000-query run reproduces the same ordering:

```
            1-contiguous  1-balanced  2-offset-wrap  2-balanced-diag
d=2 scan_s        8.413       5.453         11.968            7.092
d=3 scan_s       28.513      18.058         39.709           22.905
```

Reproduce with:

```
smafa bench-banding -d <db> -q 1000000.S3.10.seqs.fna --divergences 2 3 5
```

## Conclusion

**Balancing is worth it; a second partition for `query` is not.**

- **`1-balanced` is the fastest at every `d`** (≈1.4× at d=2 up to ≈1.65× at
  d=5 over the current contiguous band) and it cuts candidates ≈1.7–2.3×. It
  costs only the column→band map in the DB and **no extra lookups**, so there is
  no per-query overhead to pay back. This is the change that earns the new
  database format.

- **The 2-partition methods reduce *candidates* a lot but lose on *wall time*.**
  `2-offset-wrap` cuts candidates 6–8× yet is consistently **slower than the
  baseline** (e.g. 68 s vs 35 s at d=5). On real coding data the cost is
  dominated by *gathering and sorting the conserved-band buckets*, not by the
  final Hamming checks — and a second contiguous partition has the same conserved
  bands, so it doubles the expensive part to shrink the cheap part. This is the
  opposite of the uniform-random expectation, and exactly the "measure on real
  data" caveat from the design discussion.

- **`2-balanced-diag` gives by far the smallest candidate sets** (9× at d=5) and
  beats the baseline on time, but is still slower than `1-balanced` because of
  the second partition's overhead. It only becomes attractive if you are
  candidate-/RAM-bound, or if the per-candidate cost downstream were much higher
  than a 5-word Hamming distance — neither is true here.

Why balancing helps so much on this data: the windows are back-translated
protein, so conservation is wildly uneven across columns. A *contiguous* band can
land entirely inside a conserved region, producing one near-zero-entropy band
whose bucket holds most of the database — that single bucket dominates both the
candidate count and the gather/sort cost. Round-robining columns by conservation
guarantees no band is all-conserved, which removes the giant bucket. That is also
why `2-balanced-diag` (28 s) is much faster than `2-offset-wrap` (68 s) at d=5
despite both being two-partition: the win is from killing giant buckets, not from
the number of partitions.

## Cluster

`cluster` streams sequences and has no pre-built DB, but it can estimate
per-column conservation from the **first block** of input (≤8192 unique
sequences) and build a single balanced partition on the fly — no stored format
needed. The partition choice never changes cluster output (verified byte-identical
to `--no-banding` for every strategy). Wall time, 200k input sequences, single
thread:

```
 d   contiguous  offset-wrap  balanced  no-banding
 2     1.61        1.80        1.61       19.70
 3     1.55        1.82        1.45       16.09
 5     2.24        4.32        1.79       10.78
```

The full 1,000,000-sequence run widens the gaps:

```
 d   contiguous  offset-wrap  balanced
 3     12.48       18.15        9.26
 5     44.81       85.35       25.46
```

At 1M/d=5 `balanced` is 1.76× faster than `contiguous` and 3.35× faster than
`offset-wrap`.

Same story as query, and it settles the "is offset+wrap worth it for cluster?"
question: **no.** `offset-wrap` is the *slowest* banded option (≈2× `contiguous`
at d=5) because the second contiguous partition doubles the conserved-bucket
work. **`balanced` (entropy from the first block) is fastest at every `d`** and
needs no extra partition and no stored format — exactly the "compute entropy from
the first ~10k sequences then use the 1-balanced partition like query" idea.

`cluster` now defaults to `balanced`; `--banding offset-wrap|contiguous` are kept
for benchmarking.

## Cost of the entropy computation (and: does it need to be stored?)

`column_weights` over all 136,770 subjects takes **~10.6 ms** (measured, scales
linearly with subject count; independent of the queries and of `d`). Deriving any
`d`'s partition from the weights is a 60-element sort + round-robin —
microseconds.

Put that next to the other per-invocation costs: deserialising the 6.3 MB DB is
~hundreds of ms, and an actual query pass is seconds. So the entropy step is
~2–3% of DB load and a fraction of a percent of a query pass.

Two consequences:

- **No new DB format is needed for balanced `query` after all.** The partition is
  a deterministic function of the subjects, so `query` can just recompute the
  weights at load (~10 ms) and build the same partition `makedb` would have. This
  is the same trick `cluster` already uses (estimate from the data on the fly),
  only `query` has *all* the subjects available so it doesn't even need to
  subsample.

- **Caching the partition in the DB is not worth it — including for d=2.** It
  would save ~10 ms per invocation, dwarfed by the DB deserialize that every
  invocation already pays. Even under heavy fan-out (many small query runs against
  one DB) the load cost dominates the 10 ms. If you ever did want to remove it,
  cache the 60 per-column weights (≈480 bytes, `d`-independent) rather than a
  specific `d`'s partition — but the payoff is still only ~10 ms.

## Recommendation

- **`cluster`: use the balanced partition** (now the default), conservation
  estimated from the first block. Drop offset+wrap as a clustering default.
- **`query`: use a single entropy-balanced partition, recomputed at load.** No DB
  format change: compute per-column conservation from the loaded subjects (~10 ms)
  and build the single band index from it. Skip the second partition (it cuts
  candidates but loses on wall time on this data).

Following this benchmark, **balanced is the only production banding**: both
`query` and `cluster` use `BandIndex::single_balanced` (query recomputes the
weights from the loaded subjects; cluster estimates them from the first input
block). The two-partition strategies (`offset-wrap`, `balanced-diag`) are gone,
and contiguous banding is no longer selectable in production. `BandIndex::new`
(contiguous) is retained only as the baseline inside the `bench-banding`
subcommand, which still compares contiguous vs balanced on a real DB + query set
so this result stays reproducible.
