---
title: "Metagenome Assembly"
teaching: 30
exercises: 30
questions:
- "Why do raw reads need to be assembled?"
- "How can we assemble a metagenome?"
objectives:
- "Understand the difference between reads and contigs."  
- "Run a metagenomics assembly workflow."
keypoints:
- "Assembly merges reads into contigs."
- "Flye can be used as a metagenomes assembler."
- "Assemblers take FASTQ files as input and produce a FASTA file as output."
---

## Assembling reads
<img align="right" width="325" height="316" src="{{ page.root }}fig/03_short_analysis_flowchart_short_asm.png" alt="Analysis flow diagram that shows the steps: Sequence reads, Quality control and assembly." />

Now we have put our raw reads through quality control we are going to move onto the next step in the process which is metagenomics assembly.

he assembly strategy differs based on the sequencing technology used to generate the raw reads.
Here we're using raw data from [Nanopore sequencing](https://nanoporetech.com/applications/dna-nanopore-sequencing) as the basis for this
metagenome assembly so we need to use a metagenome assembler appropriate for this problem.
<br clear="right"/>



<span style="color:red"> Introduction to assembly - jigsaw analogy + image
One genome assembly is one jigsaw puzzle but a metagenome is 500+ puzzles, without pictures.
This makes assembly and further analysis a challenge - must first assembly which parts of the puzzle (sequences) we can with assembly. Then separate these larger parts into different "puzzles" / organisms using binning.
</span>

T

[Flye](https://github.com/fenderglass/Flye) is a long-read de novo assembler
for assembling large and complex metagenomics data, and it is one of the
most used and recommended.

Some of the problems faced by metagenomics assembly are: i) the differences in coverage between the genomes, due to the differences in abundance in the sample, ii) the fact that different species often share conserved regions, iii) and the presence of several strains of a single species in the community. How flye deals with this.

How flye deals with long read sequences

> ## Refresher on command line programs
> To fill!
{: .callout}

## Flye is a metagenomics assembler

Let's see what happens if we enter the `flye` command on our terminal.

~~~
 flye
~~~
{: .source}

~~~
$ flye
usage: flye (--pacbio-raw | --pacbio-corr | --pacbio-hifi | --nano-raw |
      --nano-corr | --nano-hq ) file1 [file_2 ...]
      --out-dir PATH

      [--genome-size SIZE] [--threads int] [--iterations int]
      [--meta] [--polish-target] [--min-overlap SIZE]
      [--keep-haplotypes] [--debug] [--version] [--help]
      [--scaffold] [--resume] [--resume-from] [--stop-after]
      [--read-error float] [--extra-params]
flye: error: the following arguments are required: -o/--out-dir
~~~
{: .output}

The reason we're seeing this is because we haven't provided any arguments for Flye to be able to run the assembly.

This output gives us a bit of information about how to run Flye but we can use the help command to show it in full.

~~~
$ flye --help
usage: flye (--pacbio-raw | --pacbio-corr | --pacbio-hifi | --nano-raw |
	     --nano-corr | --nano-hq ) file1 [file_2 ...]
	     --out-dir PATH

	     [--genome-size SIZE] [--threads int] [--iterations int]
	     [--meta] [--polish-target] [--min-overlap SIZE]
	     [--keep-haplotypes] [--debug] [--version] [--help]
	     [--scaffold] [--resume] [--resume-from] [--stop-after]
	     [--read-error float] [--extra-params]

Assembly of long reads with repeat graphs

optional arguments:
 -h, --help            show this help message and exit
 --pacbio-raw path [path ...]
                       PacBio regular CLR reads (<20% error)
 --pacbio-corr path [path ...]
                       PacBio reads that were corrected with other methods (<3% error)
 --pacbio-hifi path [path ...]
                       PacBio HiFi reads (<1% error)
 --nano-raw path [path ...]
                       ONT regular reads, pre-Guppy5 (<20% error)
 --nano-corr path [path ...]
                       ONT reads that were corrected with other methods (<3% error)
 --nano-hq path [path ...]
                       ONT high-quality reads: Guppy5+ or Q20 (<5% error)
 --subassemblies path [path ...]
                       [deprecated] high-quality contigs input
 -g size, --genome-size size
                       estimated genome size (for example, 5m or 2.6g)
 -o path, --out-dir path
                       Output directory
 -t int, --threads int
                       number of parallel threads [1]
 -i int, --iterations int
                       number of polishing iterations [1]
 -m int, --min-overlap int
                       minimum overlap between reads [auto]
 --asm-coverage int    reduced coverage for initial disjointig assembly [not set]
 --hifi-error float    [deprecated] same as --read-error
 --read-error float    adjust parameters for given read error rate (as fraction e.g. 0.03)
 --extra-params extra_params
                       extra configuration parameters list (comma-separated)
 --plasmids            unused (retained for backward compatibility)
 --meta                metagenome / uneven coverage mode
 --keep-haplotypes     do not collapse alternative haplotypes
 --scaffold            enable scaffolding using graph [disabled by default]
 --trestle             [deprecated] enable Trestle [disabled by default]
 --polish-target path  run polisher on the target sequence
 --resume              resume from the last completed stage
 --resume-from stage_name
                       resume from a custom stage
 --stop-after stage_name
                       stop after the specified stage completed
 --debug               enable debug output
 -v, --version         show program's version number and exit

Input reads can be in FASTA or FASTQ format, uncompressed
or compressed with gz. Currently, PacBio (CLR, HiFi, corrected)
and ONT reads (regular, HQ, corrected) are supported. Expected error rates are
<15% for PB CLR/regular ONT; <5% for ONT HQ, <3% for corrected, and <1% for HiFi. Note that Flye
was primarily developed to run on uncorrected reads. You may specify multiple
files with reads (separated by spaces). Mixing different read
types is not yet supported. The --meta option enables the mode
for metagenome/uneven coverage assembly.

To reduce memory consumption for large genome assemblies,
you can use a subset of the longest reads for initial disjointig
assembly by specifying --asm-coverage and --genome-size options. Typically,
40x coverage is enough to produce good disjointigs.

You can run Flye polisher as a standalone tool using
--polish-target option.
~~~
{: .output}

> ## Running commands on the background
> The `&` sign that we are using at the end of the command is for telling
the machine to run the command on the background, this will help us to avoid
the cancelation of the operation in case the connection with the AWS machine is unstable.
{: .callout}

## Updated to here

When the run is finished it shows this message:

~~~
======= SPAdes pipeline finished.

SPAdes log can be found here: /home/dcuser/dc_workshop/results/assembly_JC1A/spades.log

Thank you for using SPAdes!

~~~
{: .bash}

Now we need to press enter to exit from the background, and a message like this will be displayed:
~~~
[1]+  Done                    metaspades.py -1 JC1A_R1.trim.fastq.gz -2 JC1A_R2.trim.fastq.gz -o ../../results/assembly_JC1A
~~~
{: .output}
This is becacause of the use of the `&`. Now, let's go to the files:
~~~
$ cd ../../results/assembly_JC1A
$ ls -F
~~~
{: .bash}

~~~
assembly_graph_after_simplification.gfa
assembly_graph.fastg
assembly_graph_with_scaffolds.gfa
before_rr.fasta
contigs.paths
corrected
dataset.info
first_pe_contigs.fasta
input_dataset.yaml
contigs.fasta
scaffolds.fasta
K21
K33
K55
misc
params.txt
pipeline_state
run_spades.sh
run_spades.yaml
scaffolds.paths
spades.log
strain_graph.gfa
tmp

~~~
{: .output}

As we can see, MetaSPAdes gave us a lot of files. The ones with the assembly are the `contigs.fasta` and the `scaffolds.fasta`.
Also, we found three `K` folders: _K21, K33, and K55_, this contains the individual result files for an assembly
with k-mers equal to those numbers: 21, 33, and 55. The best assembled results are
the ones that are displayed outside this k-folders. The folder `corrected` hold the corrected reads
with the SPAdes algorithm. Moreover, the file
`assembly_graph_with_scaffolds.gfa` have the information needed to visualize
our assembly by different means, like programs as [Bandage](https://rrwick.github.io/Bandage/).

The contigs are just made from assembled reads, but the scaffolds are the result
from a subsequent process in which the contigs are ordered, oriented, and connected with Ns.

We can recognize which sample our assembly outputs corresponds to because they are inside
the assembly results folder: `assembly_JC1A/`. However, the files within it do not have the
sample ID. It is very useful to rename these files, in case we need them out of their folder.


> ## Exercise 1: Rename all files in a folder
>
> Add JC1A (the sample ID) separated by "_"  at the beggining of the names of all the contents in the assembly_JC1A directory. Remember that many solutions are possible.
>
> A)  mv * JC1A_    
> B)  mv * JC1A_*    
> C)  for name in *; do mv $name JC1A_; done     
> D)  for name in *; do mv $name JC1A_$name; done      
>    
>> ## Solution
>> ~~~
>>
>>  A)  No, this option is going to give you as error mv: target 'JC1A_' is not a directory
>>  This is because mv has two options
>>  mv file_1 file_2
>>  mv file_1, file_2, ..... file_n, directory
>>  When a list of files is passed to mv, the mv expects the last parameters to be a directory
>>  Here, * gives you a list of all the files in the directory
>>  The last parameter is JC1A_ (which mv expects to be a directory)  
>>  B)  No, Again every file is send to the same file.
>>  C)  No, every file is sent to the same file JC1A_
>>  D)  Yes, this is one of the possible solutions.
>>
>> ¿Do you have another solution?
>> ~~~
>> {: .bash}
> {: .solution}
{: .challenge}

> ## Exercise 2: Compare two fasta files from the assembly output
> You want to know how many contigs and how many scaffolds results for the assembly. Use `contigs.fasta`  and `scaffolds.fasta ` files and sort the commands to create correct code lines.
>  Do they have the same number of lines? Why?
> Hint: You can use the following commands: grep, | (pipe), -l, “>”, wc, filename.fasta
>
>> ## Solution
>>
>> ~~~
>> $ grep “>” contigs.fasta | wc -l
>> $ grep “>” scaffolds.fasta | wc -l
>> ~~~
>> {: .bash}
>>
>> ~~~
>> A contig is created from reads and then a scaffold from group of cotings so we expect less lines in the `scaffolds.fasta ` .
>> ~~~
> {: .solution}
>
{: .challenge}

While you're waiting for the assembly to finish you might want to read about the
different programs devoted to
[genome](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2874646/) and
metagenome assembly, which can use different assembly strategies such as; Greedy extension, OLC and De Bruijn charts.
