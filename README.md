# ShorDEC 0.1
## Inference of viral quasispecies tool
### Purpose:
Given a set of reads (only single_end right now. In a close future pair_end, and maybe no so close long reads), return all haplotypes from these reads.
### Steps:
1. Correct the reads following LorDEC's approach
..* Count k-mers (by hand so far, in a close future we are going to use Jellyfish or so).
..* Build/Prune(soft) de Bruijn Graph from solid k-mers (k-mers above a threshold).
..* Find paths using de Bruijn Graph built between reads following solid k-mers.
..* Build a Path Graph for each read.
..* Find the shortest path (using Dijkstra with prior queue so far)
..* Repeat the process "n" times with differents k-mer size (K -> 2*K -> 2*2*K, right now)

2. Find Unitigs in the de Bruijn Graph and report results
..* Still a Naive approach based on finding paths within a graph.
..* Only unitigs available -> we are planning to extend it to omnitig or other longer contigs.

### Use:
..1 make
..2 ./src/test -f [dir or file to read] -o [output] -u [unitigs file output] -k [kmer first size] -h [kmer threshold]
