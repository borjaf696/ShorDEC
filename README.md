# viaDBG 1.0
## Inference of viral quasispecies multi-assembler
### Purpose:
Given a set of paired_end reads from a viral sample, it returns the largest independent contigs that are contained into the sample. The entry is expected to have:
	* High coverage.
	* Genomes should be almost full covered.

The code accepts both fasta and fastq files.
### Overview:
1. Correct the reads following LorDEC's approach (optional)
	* Count k-mers (by hand so far, in a close future we are going to use Jellyfish or so).
	* Build/Prune(soft) de Bruijn Graph from solid k-mers (k-mers above a threshold).
	* Find paths using de Bruijn Graph built between reads following solid k-mers.
	* Build a Path Graph for each read.
	* Find the shortest path (using Dijkstra with prior queue so far)
	* Repeat the process "n" times with differents k-mer size (K -> 2*K -> 2*2*K, right now)

2. Build our own version of the Approximate Paired de Bruijn Graph:
	* Adds paired_end information.
	* Polishs paired_end information (optional) 
	* Calculates cliques for each two adjacent nodes.
	* Splits and merges the nodes of the graph.
	* Returns unitigs from the new representation.

## Tips:
1. For simulated data we discourage to use either polish or reads correction. It is not useful at all.

2. Estimated error is the trickiest part (next versions will avoid the need of using this parameter):
	* If error is much higher the graph will not have enough nodes (solid information) and probably it will erase much more information than needed.
	* If error is much lower the graph will have too much nodes (fake solid information) and the result contigs will not be as good as you expect.

### Use:
	* make clean && make
	* ./bin/viaDBG -s [single_end_reads] -p [paired_end_reads(dir)] -o [output] -u [unitigs file output] -k [kmer size] -h [polish (default No)] -b [if error correction(default No)] -d [if revcomplement paired end reads] -r [estimated error]
