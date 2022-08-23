---
title: "Metagenome Assembly"
teaching: 40
exercises: 40
questions:
- "Why do raw reads need to be assembled?"
- "How does metagenomic assembly differ from genomic assembly?"
- "How can we assemble a metagenome?"
objectives:
- "Run a metagenomic assembly workflow."
- "Assess the quality of an assembly using SeqKit"
keypoints:
- "Assembly merges raw reads into contigs"
- "Flye can be used as a metagenomic assembler."
- "Certain statistics can be used to describe the quality of an assembly"
---

> ## WARNING
> This lesson will take several hours to run and complete! You can find some recommended reading at the end of the page you might want to read whilst you're waiting.
{: .callout}

## Assembling reads
<img align="right" width="325" height="316" src="{{ page.root }}/fig/03_short_analysis_flowchart_short_asm.png" alt="Analysis flow diagram that shows the steps: Sequence reads, Quality control and assembly." />

Now we have put our raw reads through quality control we are going to move onto the next step in the process which is assembly of the metagenome.


### Genomic assembly

Genomic assembly refers to the act of joining smaller fragments of DNA (i.e. reads) to make longer segments to try and reconstruct the original genome.

You can think of this like a jigsaw puzzle, each raw read corresponds to a piece of the puzzle and you're aiming to complete the puzzle by joining these pieces together.

There are two main strategies for genome assembly.
1. a reference-mapping approach when you have a reference genome of what you have sequenced to map your smaller reads onto
2. a _de novo_ approach, this is an assembly approach that doesn't use a reference and instead assembles reads together based on the content of the reads (the specific approach depends on which assembly software you are using)

<br clear="right"/>

<img align="left" width="775" height="717" src="{{ page.root }}/fig/03_genomics_v_metagenomics.png" alt="Metagenomic flow diagram with the steps raw reads, assembly and polishing and binning ." />

Continuing the jigsaw analogy, the reference-mapping approach would be when you have an image of the final puzzle to compare your assembly to. Whereas, a _de novo_ approach you would have no image to look at and have to determine which pieces fit together based on their shape and their content.

### Metagenomic assembly

Metagenomic sequencing adds another layer to the challenge of assembly. Instead of having one organism to assemble you now have multiple! Depending on the complexity of a metagenome you could have anywhere from a handful of organisms in a community to thousands.

This means the single jigsaw puzzle of a genome assembly has now become multiple different jigsaw puzzles in one.

As many of the communities sequenced using metagenomics contain previously uncultured microbes, often known as microbial dark matter, they are unlikely to have a reference genome you can use and often you don't know before sequencing what organisms make up a community.

Note: The data we're using is a mock metagenome so we _do_ actually know what organims make up the community and have reference sequences for them so we could use a reference-mapping approach to assemble this metagenome but as this is unlikely with real-world data we're going to use a _de novo_ approach in this tutorial  

This now means assembling our metaphorical jigsaw will be a challenge! Not only is it now potentially thousands of different jigsaws in one but we also don't have any images to refer back to!

Luckily there are programs, known as assemblers, that will do this for us!

Metagenomic assembly faces additional problems, which means we need an assembler built to handle metagenomes. These additional problems include:
i) the differences in coverage between the genomes, due to the differences in abundance across the sample
ii) the fact that different species often share conserved regions
iii) and the presence of several strains of a single species in the community

The assembly strategy also differs based on the sequencing technology used to generate the raw reads. Here we're using raw data from [Nanopore sequencing](https://nanoporetech.com/applications/dna-nanopore-sequencing) as the basis for this metagenome assembly so we need to use a metagenome assembler appropriate for long-read sequencing.

We will be using [Flye](https://github.com/fenderglass/Flye), which is a **long-read** _de novo_ assembler for assembling large and complex data with a metagenomic mode.

<br clear="left"/>

## Flye is a long-read assembler

Important: Make sure you're still logged into your instance [Logging onto the Cloud](https://cloud-span.github.io/metagenomics01-qc-assembly/01-logging-onto-cloud/index.html)

Flye has been pre-installed onto your instance. First navigate to the `analysis` directory you made in a previous step.
Let's see what happens if we enter the `flye` command on our terminal.

~~~
 cd analysis
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

The above output gives us a bit of information about how to run Flye but we can use the help command to show it in full.
~~~
  flye --help
~~~
{: .bash}

> ## `flye` help documentation
> ~~~
> usage: flye (--pacbio-raw | --pacbio-corr | --pacbio-hifi | --nano-raw |
> 	     --nano-corr | --nano-hq ) file1 [file_2 ...]
> 	     --out-dir PATH
>
> 	     [--genome-size SIZE] [--threads int] [--iterations int]
> 	     [--meta] [--polish-target] [--min-overlap SIZE]
> 	     [--keep-haplotypes] [--debug] [--version] [--help]
> 	     [--scaffold] [--resume] [--resume-from] [--stop-after]
> 	     [--read-error float] [--extra-params]
>
> Assembly of long reads with repeat graphs
>
> optional arguments:
>   -h, --help            show this help message and exit
>   --pacbio-raw path [path ...]
>                         PacBio regular CLR reads (<20% error)
>   --pacbio-corr path [path ...]
>                         PacBio reads that were corrected with other methods (<3% error)
>   --pacbio-hifi path [path ...]
>                         PacBio HiFi reads (<1% error)
>   --nano-raw path [path ...]
>                         ONT regular reads, pre-Guppy5 (<20% error)
>   --nano-corr path [path ...]
>                         ONT reads that were corrected with other methods (<3% error)
>   --nano-hq path [path ...]
>                         ONT high-quality reads: Guppy5+ or Q20 (<5% error)
>   --subassemblies path [path ...]
>                         [deprecated] high-quality contigs input
>   -g size, --genome-size size
>                         estimated genome size (for example, 5m or 2.6g)
>   -o path, --out-dir path
>                         Output directory
>   -t int, --threads int
>                         number of parallel threads [1]
>   -i int, --iterations int
>                         number of polishing iterations [1]
>   -m int, --min-overlap int
>                         minimum overlap between reads [auto]
>   --asm-coverage int    reduced coverage for initial disjointig assembly [not set]
>   --hifi-error float    [deprecated] same as --read-error
>   --read-error float    adjust parameters for given read error rate (as fraction e.g. 0.03)
>   --extra-params extra_params
>                         extra configuration parameters list (comma-separated)
>   --plasmids            unused (retained for backward compatibility)
>   --meta                metagenome / uneven coverage mode
>   --keep-haplotypes     do not collapse alternative haplotypes
>   --scaffold            enable scaffolding using graph [disabled by default]
>   --trestle             [deprecated] enable Trestle [disabled by default]
>   --polish-target path  run polisher on the target sequence
>   --resume              resume from the last completed stage
>   --resume-from stage_name
>                         resume from a custom stage
>   --stop-after stage_name
>                         stop after the specified stage completed
>   --debug               enable debug output
>   -v, --version         show program's version number and exit
>
> Input reads can be in FASTA or FASTQ format, uncompressed
> or compressed with gz. Currently, PacBio (CLR, HiFi, corrected)
> and ONT reads (regular, HQ, corrected) are supported. Expected error rates are
> <15% for PB CLR/regular ONT; <5% for ONT HQ, <3% for corrected, and <1% for HiFi. Note that Flye
> was primarily developed to run on uncorrected reads. You may specify multiple
> files with reads (separated by spaces). Mixing different read
> types is not yet supported. The --meta option enables the mode
> for metagenome/uneven coverage assembly.
>
> To reduce memory consumption for large genome assemblies,
> you can use a subset of the longest reads for initial disjointig
> assembly by specifying --asm-coverage and --genome-size options. Typically,
> 40x coverage is enough to produce good disjointigs.
>
> You can run Flye polisher as a standalone tool using
> --polish-target option.
> ~~~
> {: .output}
{: .solution}


**We can see that Flye has multiple different options available so we need to work out which ones are appropriate for our dataset.**
* We know we have Nanopore raw reads and if we look back at the paper [Nicholls _et al._ 2019](https://academic.oup.com/gigascience/article/8/5/giz043/5486468) we can see that the reads were basecalled with Guppy v2.2.2. Therefore the flag `--nano-raw` for `ONT regular reads, pre-Guppy5 (<20% error)` is most appropriate for this dataset. The `path` after the flag in the help document indicates that the we should put the location of the input file after this flag.
* The next flag we are interested in is `-o` or `--outdir` to specify the location of the flye output. Again, we need to specify a path afterwards.
* We should also make use of the `-t` or `--threads` flag in order to run the assembly on more compute in order to speed it up.
* After making the initial assembly, flye will continue to further improve the assembly using the original raw data using a process called polishing. We can specify the number of times `flye` will polish this data using `-i` or `--iterations` - `number of polishing iterations [1]`. By default the number of iterations is 1 however 3 iterations is often used as standard.
* Finally, as this dataset is a metagenome we need to use the `--meta` option for `metagenome / uneven coverage mode`.

> ## Unused parameters
> There's a lot of parameters that we won't be using; some are deprecated, some are only appropriate for certain types of data (e.g. `--pacbio-raw`) and some are useful to allow tweaking to try further improve an assembly (e.g. `--genome-size` and `--read-error`).  
> Most bioinformatics programs have an associated website (which is often a GitHub page) with a whole manual to use the program.  
> The [Flye Manual](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md) contains a lot of further information about the parameters avaiable. If you're going to try using Flye on your own long-read dataset this is a good place to start.  
{: .callout}

Now we've worked out what parameters are appropriate for our data we can put them all together in one command.
We will be using the filtered Nanopore file we generated in the previous step which should be in the location `~/data/nano_fastq/ERR3152367_sub5_filtered.fastq`
We're also going to get Flye to create the `assembly` directory as its output directory.

~~~
 flye --nano-raw ~/data/nano_fastq/ERR3152367_sub5_filtered.fastq \
     --out-dir assembly \
     --threads 4 \
     --iterations 3 \
     --meta
~~~
{: .bash}

**<span style="color:red"> Don't run this command yet! If you have, you can press <kbd>Ctrl</kbd>+<kbd>z</kbd> to stop the command.</span>**

Now we've built our command we could stop here **but** metagenomic assembly takes a long time.  
If we were to run this command as is we'd have to stay logged into the instance (aka leaving your computer running) for hours.  
Luckily we don't have to do that as we're using a remote computer (as that's what the instance/cloud computing is).

### Running a command in the background

The commands we've previously run in this course have all been run in the foreground - aka they've been run directly in the terminal window we've been using and occupy the window until they've finished.

We can instead run a job in the background so it doesn't take over the terminal window.

To do this we take the command we want to run and then follow it by an ampersand (`&`) symbol.

This puts the job into the background so we can do other things in the terminal but it will still stop running if you log out of the instance.
Finally we need to redirect the output Flye reports to the terminal into a file with `>`.

Once we add these into our command we get the following:
~~~
flye --nano-raw ~/data/nano_fastq/ERR3152367_sub5_filtered.fastq \
     --out-dir assembly \
     --threads 4 \
     --iterations 3 \
     --meta &> flye_output.txt &
~~~
{: .bash}

`&>` redirects any logging information by the program that would originally come to the terminal and save to file. Note the lack of a space between `&>`. The second `&` then runs this command in the background.

We can now press enter to run the command.
Unlike when we have previously run code, your prompt should immediately return. This doesn't mean that the code has finished already, it should now be running in the background.

> ## Running commands on different servers
> There's many different options to run commands in the background in terminal.  
> How you run these commands (also known as jobs) will depend on the computing resources (and their fair use policies) you are running the command on.  
> The main options include:
> * `&` which we've covered here, depending on the infrastructure you're running the command on you may also need to use [`nohup`](https://www.digitalocean.com/community/tutorials/nohup-command-in-linux) to prevent the background job from being killed when you close the terminal.  
> * The command line program [`screen`](https://linuxize.com/post/how-to-use-linux-screen/), which allows you to create a shell session that can be completely detached from a terminal and re-attached when needed.
> * Queuing system - many shared computing resources, like  High Performance Computers (HPC) some Universities have, operate a queuing system, e.g. SLURM or SGE, so each user gets their fair share of computing resources. With these you submit your command / job to the queueing system, which will then handle when to run the job on the resources available.
{: .callout}

As we're running the command in the background we no longer see the output on the terminal. Luckily we have two options available for us to check on the progress of the assembly.

We can view all the currently running background commands with the command `jobs`. The output of this will be similar to:

~~~
jobs
~~~
{: .bash}
~~~
[1]+  Running                 flye --nano-raw ~/data/nano_fastq/ERR3152367_sub5_filtered.fastq --out-dir assembly --threads 4 --iterations 3 --meta &> flye.out &
~~~
{: .output}

Note: if you disconnect from the instance while Flye is running you won't be able to track the jobs progress through this method.

Flye also generates a log file when running within the output folder it has generated.
Using less we can navigate through this file.
~~~
less assembly/flye.log
~~~
{: .bash}

The contents of the file will depend on how far through the assembly Flye is.
At the start of an assembly you'll probably see something like this:
~~~
[2022-08-15 19:58:59] INFO: Starting Flye 2.9-b1778
[2022-08-15 19:58:59] INFO: >>>STAGE: configure
[2022-08-15 19:58:59] INFO: Configuring run
[2022-08-15 19:59:55] INFO: Total read length: 3082258211
[2022-08-15 19:59:55] INFO: Reads N50/N90: 5373 / 2571
[2022-08-15 19:59:55] INFO: Minimum overlap set to 3000
[2022-08-15 19:59:56] INFO: >>>STAGE: assembly
~~~
{: .output}

Note: this log file will contain similar to the `flye_output.txt` file we're generating when redirecting the terminal output. But it's easier to look at the log file as flye will always generate that even if you're running the command differently (e.g. in the foreground).

>## Navigation commands in `less`:
>
> | key     | action |  
> | ------- | ---------- |  
> | <kbd>Space</kbd> | to go forward |  
> |  <kbd>b</kbd>    | to go backward |  
> |  <kbd>g</kbd>    | to go to the beginning |  
> |  <kbd>G</kbd>    | to go to the end |  
> |  <kbd>q</kbd>    | to quit |  
>
> See [Prenomics - Working with Files and Directories](https://cloud-span.github.io/prenomics02-command-line/02-working-with-file/index.html) for a full overview on using less.
{: .callout}


Flye is likely to take a *couple of hours* to finish assembling.
You don't need to remain connected to the instance during this time but once you have disconnected it does make it harder to track the progress of Flye.

In the meantime, if you wanted to read more about assembly and metagenomics there's a few papers and resources below with recommended reading.

### Determining if the assembly has finished

After leaving it at least a couple of hours (or even longer!), Flye should have finished assembling.

If you remained connected to the instance during the process you will get the following output in your terminal when the command has finished.

~~~
[2]+  Done      flye --nano-raw ~/data/nano_fastq/ERR3152367_sub5_filtered.fastq --out-dir assembly --threads 4 --iterations 3 --meta &> flye.out &
~~~
{: .output}

If you disconnected from the instance for whatever reason during the assembly process we need to check the `flye.log` file.  

~~~
cd analysis/assembly/
less flye.log
~~~
{: .bash}

If you navigate to the end of the file you should see something like:
~~~
[2022-08-11 14:03:51] INFO: >>>STAGE: finalize
[2022-08-11 14:03:52] INFO: Assembly statistics:

        Total length:   14953273
        Fragments:      146
        Fragments N50:  2976503
        Largest frg:    6068630
        Scaffolds:      0
        Mean coverage:  181

[2022-08-11 14:03:52] INFO: Final assembly: ~/analysis/assembly/assembly.fasta
~~~
{: .output}

It also contains some basic statistics about the assembly created.

### Assembly Output

If we `ls` in the assembly directory we can see the that Flye has created many different files.

~~~
00-assembly   20-repeat     40-polishing    assembly_graph.gfa  assembly_info.txt  params.json
10-consensus  30-contigger  assembly.fasta  assembly_graph.gv   flye.log
~~~
{: .output}

We've already looked at `flye.log` which contains all the info Flye generates during an assembly.

* Flye generates a directory to contain the output for each step of the assembly process. (These are the `00-assembly`, `10-consensus`, `20-repeat`, `30-contigger` and `40-polishing` directories.)  
* We also have a file containing the parameters we ran the assembly under `params.json` which is useful to keep our analysis reproducible.  
* The assembled contigs are in FASTA format (`assembly.fasta`).  
* There's a text file which contains more information about each contig created (`assembly_info.txt`).
* Finally we have two files for a repeat graph (`assembly_graph.{gfa|gv}`) which is a visual way to view the assembly - see the optional exercise below.    

You can see more about the output for Flye in the [documentation on GitHub](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#output).

> ## Contigs vs. reads
> We have seen reads in the raw sequencing data - these are our individual jigsaw pieces.  
> After assembly we introduce contigs. Contigs (from the word *contiguous*) are fragments of DNA produced after raw reads are joined together by the assembly process. (These are the chunks of the jigsaw puzzle the assembler has managed to complete).  
> As a result contigs are usually much longer than raw reads and also vary in length and number depending on how successful the assembly has been.  
{: .callout}


> ## Optional exercise: Viewing the repeat graph
> If you are interested you can explore the repeat graph created by Flye - this is an entirely optional exercise as not every assembler will generate a repeat graph and usually metagenomic repeat graphs are are large and complicated so don't tell you much. Also it requires some software to be downloaded which can be challenging if you don't have admin rights to your computer.  
> A repeat or assembly graph will show the final contigs of an assembly and how they interact with each other - see [graph theory](https://en.wikipedia.org/wiki/Graph_theory) for more information about what we mean by graph in this context.  
> **If you don't want / aren't able to complete these steps you can see a bandage output of this assembly in the drop down below.**   
> > ## Exercise: Using Bandage to view the repeat graph
> > * In order to view the repeat graph we need to install a program called Bandage onto your local computer. You can download this from the [Bandage Website](https://rrwick.github.io/Bandage/) by selecting your operating system down the side.   
> > * Once the software is installed, we also need to download the graph file generated by the assembler. To do this, adapt the `scp` steps to download the `assembly_graph.gfa` file onto your local computer.  
> > * Once you have the file downloaded and the Bandage software installed, we can view the graph in bandage - we will be following the [Getting Started](https://github.com/rrwick/Bandage/wiki/Getting-started) steps from the Bandage documentation.  
> > * With Bandage open, load the `gfa` graph you have downloaded (follow the steps in the documentation above) and then click the <kbd>draw graph</kbd> button.  
> > * Bandage will take some time to draw the graph, however this depends on the size of the graph.  
> > * You should then see the graph of the assembly and be able to change colours, zoom and save your view to an image (File>Save image)
> {: .solution}
> > ## Solution: Bandage output of this assembly  
> > Your output may look a little different to this as Bandage may have drawn your graph a little differently.  
>> <a href="{{ page.root }}/fig/03_bandage_graph.png">
> > <img align="center" width="713" height="611" src="{{ page.root }}/fig/03_bandage_graph.png" alt="Repeat graph of the assembly showing some contigs joined together in a circle, and some more with small fragments" />
> > </a>
> > In the above graph, there are two large circularised contigs. Which indicates that they're likely complete genomes. The larger of the two has a smaller blue circle attached which could be a plasmid or some form of insertion, though it also could be an artifact of the assembly. (This contig could be run through BLAST to work out its identity however we will be using a different analysis workflow here that's more appropriate for read world metagenomes.)
> Aside from these contigs, the rest of the contigs seem to be a lot shorter with few interactions between them.
> {: .solution}
{: .challenge}


## Assembly Statistics

As we've just seen, Flye has finished the assembly and also given us some basic statistics about the size of the assembly. Not every assembler will give you this information so we will be using assembly FASTA file and the program [Seqkit](https://bioinf.shenwei.me/seqkit/) again, but this time with a different command, `stats`, to generate basic statistics.

We can view the help documentation for this command:
~~~
seqkit stats -h
~~~
{: .bash}

> ## seqkit stats help documentation
> ~~~
> simple statistics of FASTA/Q files
>
> Tips:
>   1. For lots of small files (especially on SDD), use big value of '-j' to
>      parallelize counting.
>
> Usage:
>   seqkit stats [flags]
>
> Aliases:
>   stats, stat
>
> Flags:
>   -a, --all                  all statistics, including quartiles of seq length, sum_gap, N50
>   -b, --basename             only output basename of files
>   -E, --fq-encoding string   fastq quality encoding. available values: 'sanger', 'solexa', 'illumina-1.3+', 'illumina-1.5+', 'illumina-1.8+'. (default "sanger")
>   -G, --gap-letters string   gap letters (default "- .")
>   -h, --help                 help for stats
>   -e, --skip-err             skip error, only show warning message
>   -i, --stdin-label string   label for replacing default "-" for stdin (default "-")
>   -T, --tabular              output in machine-friendly tabular format
>
> Global Flags:
>       --alphabet-guess-seq-length int   length of sequence prefix of the first FASTA record based on which seqkit guesses the sequence type (0 for whole seq) (default 10000)
>       --id-ncbi                         FASTA head is NCBI-style, e.g. >gi|110645304|ref|NC_002516.2| Pseud...
>       --id-regexp string                regular expression for parsing ID (default "^(\\S+)\\s?")
>       --infile-list string              file of input files list (one file per line), if given, they are appended to files from cli arguments
>   -w, --line-width int                  line width when outputting FASTA format (0 for no wrap) (default 60)
>   -o, --out-file string                 out file ("-" for stdout, suffix .gz for gzipped out) (default "-")
>       --quiet                           be quiet and do not show extra information
>   -t, --seq-type string                 sequence type (dna|rna|protein|unlimit|auto) (for auto, it automatically detect by the first sequence) (default "auto")
>   -j, --threads int                     number of CPUs. can also set with environment variable SEQKIT_THREADS) (default 4)
> ~~~
> {: .output}
{: .solution}

We can see there's many different options available, however we don't need many of them to get the basic statistics for the assembly.

First we're going to use the default options to get the statistics.
~~~
seqkit stats assembly.fasta
~~~
{: .bash}

SeqKit is fast so we are just running this directly in the terminal foreground. It should take a couple of seconds to process this assembly (however, it can take longer with more sequencing data).

Once it has finished you should see an output table like this:
~~~
file            format  type  num_seqs     sum_len  min_len    avg_len    max_len
assembly.fasta  FASTA   DNA        146  14,953,273    3,164  102,419.7  6,068,630
~~~
{: .output}

In this table we can see the input file, the format of the file, the type of sequence and other statistics.

Using this table of statistics, answer the questions below.
> ## Exercise X: Looking at basic statistics
> a) How many contigs are in this assembly?  
> b) How many bases in total have been assembled?  
> c) What is the shortest and longest contig produced by this assembly?  
>> ## Solution
>> a) From `num_seqs` we can see that this assembly is made up of 146 contigs  
>> b) Looking at `sum_length` we can see that the assembly is 14,953,273bp in total (almost 15 million bp!)  
>> c) From `min_length` we can see the shortest contig is 3,164bp and from `max_length` the longest contig is 6,068,630bp  
> {: .solution}
{: .challenge}


### The N50 length

As we touched on in the previous lesson, the N50 length is a useful statistic when looking at sequences of varying length as it indicates that 50% of the total sequence is in reads (i.e. chunks) that are that size or larger. See the webpage [What's N50?](https://www.molecularecologist.com/2017/03/29/whats-n50/) for a good explanation.

This is a useful statistic to describe an assembly as it indicates the average size of the contigs the assembly software has produced.

A higher N50 length, means that more of the assembly is in longer fragments. So the chunks of sequence produced by the assembler are, on average, larger.

While it isn't calculated by default, `seqkit stats` has an option to calculate the N50 length. Using the help documentation for seqkit stats answer the exercise below.

> ## Exercise X: Flag to get the N50 length
> a) Using the help documentation, what flag can we add to get the N50 length for this assembly?  
> b) What would the new command be if we added this flag?  
> Bonus exercise: What flag would enable us to save the output table in a tabular (i.e. tsv) format?
>> ## Solution
>> a) We can see from the help documentation that the flag `-a` or `--all` will calculate `all statistics, including quartiles of seq length, sum_gap, N50`.  
>> b) The new command would be `seqkit stats -a assembly.fasta` or `seqkit stats --all assembly.fasta`  
>> Bonus: The flag `-T` allows us to save it in a tabular output - this makes the table easier to use in other command line programs or programming languages such as R and Python. The command could be either `seqkit stats -a -T assembly.fasta` or we can combine the two flags `seqkit stats -aT assembly.fasta`
> {: .solution}
{: .challenge}

Next, run the command to calculate the N50 length and answer the questions below about the output.

> ## Hint: Seqkit stats N50 command
> ~~~
> seqkit stats -a assembly.fasta
> ~~~
{: .solution}


> ## Exercise X: Calculating the N50 length
> a) What is the output if we run the new command from the above exercise?  
> b) What new statistics do we get that we didn't have with the original command?  
> c) What is the N50 length of this assembly?  
> Bonus exercise: Looking at the [information available online for Seqkit stats](https://bioinf.shenwei.me/seqkit/usage/#stats), can you work out what the extra statistics other than N50 tell us?
>> ## Solution
>> a)
>> ~~~
>> file         format  type  num_seqs     sum_len  min_len    avg_len    max_len     Q1        Q2      Q3  sum_gap        N50  Q20(%)  Q30(%)  GC(%)
>> assembly.fasta  FASTA   DNA        146  14,953,273    3,164  102,419.7  6,068,630  7,364  13,415.5  35,259        0  2,976,503       0       0  52.48
>> ~~~
>> {: .output}
>> b) Comparing the header line from this command to the original command we can see we've now got statistics for Q1, Q2, Q3, sum_gap, N50, Q20(%) and Q30(%)  
>> c) The N50 length for this assembly is 2,976,503 bp, this tells us that 50% of the assembly is in fragments that are almost 3m bases long or longer!  
>> Bonus: `Q1`, `Q2`, `Q3` is the quartile range of sequence length, `sum_gap` is the total number of ambiguous bases (N's) in the sequence, N50 we have covered, Q20(%) is the percentage of bases with a PHRED score over 20, Q30(%) is the percentage of bases with a PHRED score over 30. GC(%) is the [guanine-cytosine content](https://en.wikipedia.org/wiki/GC-content) of the sequence.   
> {: .solution}
{: .challenge}


> ## Recommended reading:
> While you're waiting for the assembly to finish here's some things you might want to read about:
> * An overall background to the history of DNA sequencing in [DNA sequencing at 40: past, present and future](https://www.nature.com/articles/nature24286)  
> * An overview of a metagenomics project  [Shotgun metagenomics, from sampling to analysis](https://www.nature.com/articles/nbt.3935) - though note this paper is from 2017 so some techniques and software will be different now.  
> * The challenges of genomic and metagenomic assembly and the algorithms have been built to overcome these in [Assembly Algorithms for Next-Generation Sequencing Data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2874646/)  
> * The approach Flye uses to assemble metagenomes is covered in [metaFlye: scalable long-read metagenome assembly using repeat graphs](https://www.nature.com/articles/s41592-020-00971-x)
{: .callout}  
