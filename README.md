# viaDBG 1.0
## Inference of viral quasispecies multi-assembler
### Purpose:
Given a set of paired_end reads from a viral sample, it returns the largest independent contigs that are contained into the sample. The entry is expected to have:
	* High coverage.
	* Genomes should be almost full covered.

The code accepts both fasta and fastq files.
### Overview:
1. Correct the reads following LorDEC's approach (optional)
	* Count k-mers (DSK is the default and the recommended way)
	* Build/Prune(soft) de Bruijn Graph from solid k-mers (k-mers above a threshold).
	* Find paths using de Bruijn Graph built between reads following solid k-mers.
	* Build a Path Graph for each read.
	* Find the shortest path (using Dijkstra with prior queue so far)
	* Repeat the process "n" times with differents k-mer size (K -> 2*K -> 2*2*K, right now)

2. Representative nodes.
	* Each unitig, arbritrarily long, is identify only by the three nodes: first, last and middle one.

3. Own version of the Approximate Paired de Bruijn Graph:
	* Removing duplicated reads.
	* Adding paired_end information.
	* Polishing paired_end information (optional) 
	* Building cliques for each pair of adjacent nodes.
	* Split and merging the nodes of the graph.
	* Reporting unitigs from the new representation.

## Tips:
1. Error and polishing are only useful when data is extremely complex. Furthermore, both have a huge efficiency impact, thus they are only recommended when result with the regular execution is far from being good. Otherwise, we encourage to use the regular execution.

2. Estimated error is the trickiest part (next versions will avoid the need of using this parameter):
	* If error is much higher the graph will not have enough nodes (solid information) and probably it will erase much more information than needed.
	* If error is much lower the graph will have too much nodes (fake solid information) and the result contigs will not be as good as you expect.

### Input:
1. If pear:
	* -s single-end-reads.fasta/fastq and -p paired-end-dir/
2. No-pear:
	* -p paired-end-dir/
3. Error-ratio: 
	* Default: 0.0008
	* Recommended: from 0.0008 to 0.003.
4. Full information:
	* If true - both reverse complementary and forward reads are going to be added as paired-end information. Only recommended for complex dataset, more is not always better.
	* If not true - only forward information is added. Meaning forward as the original shape for the reads.
	* Default: False
5. Remove duplicates:
	* Ideally it makes no changes over the results. Unfortunately, it can improve or get worse, as before more is not always better. Nevertheless, speed is improved when removing duplicated reads.
6. K-mer size:
	* With the current dsk compilaton maximum k-mer size is 192.
	* A higher one can be used recompiling DSK and modifying Utils/scripts/dsk_script.
### Use:
	* make clean && make
	* ./bin/viaDBG -s [single_end_reads] -p [paired_end_reads(dir)] -o [output] -u [unitigs file output] -k [kmer size] -h [polish (default No)] -b [if error correction(default No)] -d [if revcomplement paired end reads] -r [estimated error] -n [remove duplicated] -f [add full information] -t [num threads]

#### Example:
	* ./bin/viaDBG -s ../Datasets/Helsinki2.0/definitivo/pear.assembled.fastq -p ../Datasets/Helsinki2.0/definitivo/pair/ -o Output/SequenceContainer.fasta -u ../Output/UnitigsDiscovered_new.gfa -k 120 -r 0.00095 -c dsk -n -t 32
	* ./bin/viaDBG -p ../Datasets/Helsinki2.0/3-strain-ZIKV-20000x/ -o Output/SequenceContainer.fasta -u ../Output/UnitigsDiscovered_new.gfa -k 120 -r 0.0008 -c dsk -n -t 32
