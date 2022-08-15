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
<img align="right" width="325" height="316" src="{{ page.root }}/fig/03_short_analysis_flowchart_short_asm.png" alt="Analysis flow diagram that shows the steps: Sequence reads, Quality control and assembly." />

Now we have put our raw reads through quality control we are going to move onto the next step in the process which is assembly of the metagenome.

### Genomic assembly

Genomic assembly refers to the act of joining smaller fragments of DNA (i.e. reads) to make longer segments to try and reconstruct the original genome.

You can think of this like a jigsaw puzzle, each raw read corresponds to a piece of the puzzle and you're aiming to complete the puzzle by joining these pieces together.

There are two main strategies for genome assembly.
1. a reference-mapping approach when you have a reference genome of what you have sequenced to map your smaller reads onto
2. a _de novo_ approach, this is an assembly approach that doesn't use a reference and instead assembles reads together based on the content of the reads (the specific approach depends on which assembly software you are using)

Continuing the jigsaw analogy, the reference-mapping approach would be when you have an image of the final puzzle to compare your assembly to. Whereas, a _de novo_ approach you would have no image to look at and have to determine which pieces fit together based on their shape and their content.

<br clear="right"/>

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

[Flye](https://github.com/fenderglass/Flye) is a **long-read** _de novo_ assembler
for assembling large and complex data with a metagenomic mode.

> ## Choosing an appropriate program
> To fill!
{: .callout}

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


We can see that Flye has multiple different option available so we need to work out which ones are appropriate for our dataset.
We know we have Nanopore raw reads and if we look back at the paper [Nicholls _et al._ 2019](https://academic.oup.com/gigascience/article/8/5/giz043/5486468) we can see that the reads were basecalled with Guppy v2.2.2

Therefore the flag `--nano-raw` for `ONT regular reads, pre-Guppy5 (<20% error)` is most appropriate for this dataset. The `path` after the flag in the help document indicates that the we should put the location of the input file after this flag.

The next flag we are interested in is `-o` or `--outdir` to specify the location of the flye output. Again, we need to specify a path afterwards.

We should also make use of the `-t` or `--threads` flag in order to run the assembly on more compute in order to speed it up.

After making the initial assembly, flye will continue to further improve the assembly using the original raw data using a process called polishing. We can specify the number of times `flye` will polish this data using `-i` or `--iterations` - `number of polishing iterations [1]`. By default the number of iterations is 1 however 3 iterations is often used as standard.

Finally, as this dataset is a metagenome we need to use the `--meta` option for `metagenome / uneven coverage mode`.

> ## Unused parameters
> There's a lot of parameters that we won't be using; some are deprecated, some are only appropriate for certain types of data (e.g. `--pacbio-raw`) and some are useful to allow tweaking to try further improve an assembly (e.g. `--genome-size` and `--read-error`).  
> Most bioinformatics programs have an associated website (which is often a GitHub page) with a whole manual to use the program.  
> The [Flye Manual](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md) contains a lot of further information about the parameters avaiable. If you're going to try using Flye on your own long-read dataset this is a good place to start.  
{: .callout}

Now we've worked out what parameters are appropriate for our data we can put them all together in one command.

~~~
 flye --nano-raw ~/data/nano_fastq/ERR3152367_sub5_filtered.fastq \
     --out-dir assembly \
     --threads 4 \
     --iterations 3 \
     --meta
~~~
{: .bash}

We will be using the filtered Nanopore file we generated in the previous step which should be in the location `~/data/nano_fastq/ERR3152367_sub5_filtered.fastq`
We're also going to get Flye to create the `assembly` directory as its output directory.

**<span style="color:red"> Don't run this command yet!</span>**

Now we've built our command we could stop here **but** metagenomic assembly takes a long time!

If we were to run this command as is we'd have to stay logged into the instance (aka leaving your computer running) for hours.

Luckily we don't have to do that using a remote computer (as that's what the instance/cloud computing is).

## Running a command in the background

The commands we've previously run in this course have all been run in the foreground - aka they've been run directly in the terminal window we've been using and occupy the window until they've finished.

We can instead run a job in the background that doesn't take over the terminal window.

To do this we take the command we want to run and then follow it by an ampersand (`&`) symbol.

This puts the job into the background so we can do other things in the terminal but it will still stop running if you log out of the instance.
So we need to also add no hangup (`nohup`) to prevent the job from terminating when we log off from the instance.
Finally we need to redirect the output Flye reports to the terminal into a file with `>`.

Once we add these three things into our command we get the following:
~~~
nohup flye --nano-raw ~/data/nano_fastq/ERR3152367_sub5_filtered.fastq \
     --out-dir assembly \
     --threads 4 \
     --iterations 3 \
     --meta &> flye.out
~~~
{: .bash}

Note the lack of a space between `&>`.

> ## Running commands on different servers
> Depending on your computing resources -
> Those with shared computing resources often get given the option to run a job in a queue --> slurm etc
> Can run in screen
> The `&` sign that we are using at the end of the command is for telling
the machine to run the command on the background, this will help us to avoid
the cancelation of the operation in case the connection with the AWS machine is unstable.
> TO FILL 
{: .callout}

We can now press enter to run the command.
Unlike when we have previously run code, your prompt should immediately return. This doesn't mean that the code has finished already, it should now be running in the background.
As we have used the `nohup` command to disconnect the command from our terminal this makes it slightly more challenging to keep an eye

We can check this two ways:

1. The file `flye.out` will be created in this location and will be updated by Flye as it continues the assembly. So we can look at the contents of this file to see how far the assembly has got.
Using less we can navigate through this file.
~~~
less flye.out
~~~
{: .bash}

The contents of the file will depend on how far through the assembly flye is.
At the start of an assembly you'll probably see something like this:
~~~
TO FILL
~~~
{: .output}

>## Reminder for how to use less
> `less`
> `g` to get to the end
> `q` to quit etc
> TO FILL
{: .callout}

2. While the job is still running we can also use a unix command `ps` to list all running processes.

~~~
ps
~~~
{. :bash}
This lists all processes currently running (including background processes that keep the instance running, but we don't need to care about that). If you've run the earlier command correctly you should be able to see the command you ran in the fourth column.
~~~
TO FILL
~~~
{: .output}


Flye is likely to take a *couple of hours* to finish assembling so you should go and do something else for a while. If you've run everything correctly you should be able to entirely log off from the instance (& even shut your computer down) with Flye still running.

In the meantime, if you wanted to read more about assembly and metagenomics there's a few papers and resources below with recommended reading.

> ## Recommended reading:
> While you're waiting for the assembly to finish you might want to read about the different programs devoted to [genome](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2874646/) and metagenome assembly, which can use different assembly strategies such as; Greedy extension, OLC and De Bruijn charts.
> MORE HERE
{: .callout}

## Determining if the assembly has finished

After leaving it at least a couple of hours (or longer!) we can check if Flye has finished using a combination of the two above commands.

First we can see if the command is still running using `ps`. If you don't see the `nohup` command running it likely means it has finished.
We can then check the output file from Flye.

~~~
cd analysis/assembly/
less flye.out
~~~
{: .bash}
If you navigate to the end of the file you should see something like:
~~~
TO FILL -- from INFO: >>>STAGE: finalize etc
~~~
{: .output}

## Assembly stats

If we `ls` in the assembly directory we can see the files that Flye has created.



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
