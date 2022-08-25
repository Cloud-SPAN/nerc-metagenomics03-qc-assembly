---
title: "Polishing with long reads"
teaching: 30
exercises: 10
questions:
- "Why do assemblies need to be polished?"
- "What is the difference between reads and contigs?"
- "How can we assemble a metagenome?"
objectives:
- "Understand what is an assembly."  
- "Run a metagenomics assembly workflow."
- "Use an enviroment in a bioinformatic pipeline."
keypoints:
- "Assembly groups reads into contigs."
- "De Brujin Graphs use Kmers to assembly cleaned reads"
- "MetaSPAdes is a metagenomes assembler."
- "Assemblers take FastQ files as input and produce a Fasta file as output."
---

Now we have generated a draft genome  
We need to polish it blah blah  
What is polishing?  
Why do we need to polish?  
How do we polish?  

## Polishing an assembly with long reads

We're first going to use the filtered raw long reads to polish the draft Flye assembly.
As with the assembly, we need to use polishing software that is especially written for long read raw reads.

[Medaka](https://github.com/nanoporetech/medaka) is a command line tool built by Oxford Nanopore Technologies which will polish an assembly by generating a consensus from raw Nanopore sequences using a recurrent neural network.

We will be using the command medaka_consensus, which is a pipeline that will first align the raw reads to the draft assembly, processes this alignment to generate a pileup which is presented to a recurrent neural network in order to produce a consensus sequence.

Medaka has been pre-installed on the instance so, first we can look at the help page for `medaka_consensus`.
Note: because this is a pipeline made up multiple steps we will be looking at the help documentation for `medaka_consensus` (no space!)

~~~
medaka_consensus -h
~~~
{: .bash}

> ## `medaka_consensus` Help
>
> ~~~
> medaka 1.7.0
>
> Assembly polishing via neural networks. Medaka is optimized
> to work with the Flye assembler.
>
> medaka_consensus [-h] -i <fastx> -d <fasta>
>
>     -h  show this help text.
>     -i  fastx input basecalls (required).
>     -d  fasta input assembly (required).
>     -o  output folder (default: medaka).
>     -g  don't fill gaps in consensus with draft sequence.
>     -r  use gap-filling character instead of draft sequence (default: None)
>     -m  medaka model, (default: r941_min_hac_g507).
>         Choices: r103_fast_g507 r103_hac_g507 r103_min_high_g345 r103_min_high_g360 r103_prom_high_g360 r103_sup_g507 r1041_e82_400bps_fast_g615 r1041_e82_400bps_hac_g615 r1041_e82_400bps_sup_g615 r104_e81_fast_g5015 r104_e81_hac_g5015 r104_e81_sup_g5015 r104_e81_sup_g610 r10_min_high_g303 r10_min_high_g340 r941_e81_fast_g514 r941_e81_hac_g514 r941_e81_sup_g514 r941_min_fast_g303 r941_min_fast_g507 r941_min_hac_g507 r941_min_high_g303 r941_min_high_g330 r941_min_high_g340_rle r941_min_high_g344 r941_min_high_g351 r941_min_high_g360 r941_min_sup_g507 r941_prom_fast_g303 r941_prom_fast_g507 r941_prom_hac_g507 r941_prom_high_g303 r941_prom_high_g330 r941_prom_high_g344 r941_prom_high_g360 r941_prom_high_g4011 r941_prom_sup_g507 r941_sup_plant_g610
>         Alternatively a .tar.gz/.hdf file from 'medaka train'.
>     -f  Force overwrite of outputs (default will reuse existing outputs).
>     -x  Force recreation of alignment index.
>     -t  number of threads with which to create features (default: 1).
>     -b  batchsize, controls memory use (default: 100).
> ~~~
> {: .output}
{: .solution}

* From this we can see that `-i` for the input basecalls (meaning Nanopore raw-reads) and `-d` for the assembly are required.  
* As Medaka uses recurrent neural networks we need to pick an appropriate model (`-m`) for the data we're using. From the [documentation](https://github.com/nanoporetech/medaka#models), Medaka models are named to indicate i) the pore type, ii) the sequencing device (MinION or PromethION), iii) the basecaller variant, and iv) the basecaller version, with the format: `{pore}_{device}_{caller variant}_{caller version}`. Medaka doesn't offer an exact model for our dataset while it is possible to train a model yourself we will not be doing that here and instead will use the closest available model. This is the model `r941_prom_fast_g303`, so we also need to add that to our command as that isn't medakas default.  
* Finally, to speed this step up we need to specify the number of threads with `-t`.

This gives us the command:
~~~
medaka_consensus -i ERR3152367_sub5_filtered.fastq -d assembly.fasta -m r941_prom_fast_g303 -t 4 &> medaka.out &
~~~
{: .bash}
Note, we have added `&> medaka.out &` to redirect the output and run the command in the background. Medaka shouldn't take as long as Flye did in the previous step (probably around 20 mins), but it's a good idea to run things in the background so that you can do other things while the program is running.

Similar to Flye, we can look in the output file (`medaka.out`) to check the progress of the command.
~~~
less medaka.out
~~~
{: .bash}
If the medaka command has been run correctly you will see something like this at the start of the output:
~~~
Checking program versions
This is medaka 1.7.0
Program    Version    Required   Pass     
bcftools   1.15.1     1.11       True     
bgzip      1.15.1     1.11       True     
minimap2   2.24       2.11       True     
samtools   1.15.1     1.11       True     
tabix      1.15.1     1.11       True     
Aligning basecalls to draft
Constructing minimap index.
[M::mm_idx_gen::0.515*0.99] collected minimizers
[M::mm_idx_gen::0.648*1.40] sorted minimizers
[M::main::0.877*1.29] loaded/built the index for 146 target sequence(s)
[M::mm_idx_stat] kmer size: 15; skip: 10; is_hpc: 0; #seq: 146
[M::mm_idx_stat::0.910*1.28] distinct minimizers: 2598367 (94.68% are singletons); average occurrences: 1.076; average spacing: 5.350; total length: 14953273
~~~
{: .output}
Medaka first looks for the other programs that it needs (known as dependencies) and their versions, which have also been pre-installed on the instance. Once it confirms they are preset it begins by aligning the raw reads (basecalls) to the assembly using minimap.
Once medaka has completed the end of the file will contain something like:
~~~
[20:51:16 - DataIndx] Loaded 1/1 (100.00%) sample files.
[20:51:16 - DataIndx] Loaded 1/1 (100.00%) sample files.
[20:51:16 - DataIndx] Loaded 1/1 (100.00%) sample files.
[20:51:16 - DataIndx] Loaded 1/1 (100.00%) sample files.
[20:51:16 - DataIndx] Loaded 1/1 (100.00%) sample files.
[20:51:16 - DataIndx] Loaded 1/1 (100.00%) sample files.
[20:51:16 - DataIndx] Loaded 1/1 (100.00%) sample files.
[20:51:17 - DataIndx] Loaded 1/1 (100.00%) sample files.
[20:51:17 - DataIndx] Loaded 1/1 (100.00%) sample files.
[20:51:17 - DataIndx] Loaded 1/1 (100.00%) sample files.
Polished assembly written to medaka/consensus.fasta, have a nice day.
~~~
{: .output}

Once medaka has completed we can navigate into the output directory and look at the files Medaka has generated.

~~~
cd medaka
ls
~~~
{: .bash}

~~~
calls_to_draft.bam  calls_to_draft.bam.bai  consensus.fasta  consensus.fasta.gaps_in_draft_coords.bed  consensus_probs.hdf
~~~
{: .output}

We can see that medaka has created multiple files, these are the following:

* `calls_to_draft.bam` - this is a BAM file contaning the alignment of the raw reads (basecalls) to the draft assembly
* `calls_to_draft.bam.bai` - this is an index file of the above BAM file
* `consensus.fasta` - this is the consensus sequence, or polished assembly in our case in FASTA format
* `consensus.fasta.gaps_in_draft_coords.bed` - this is a BED file containing information about the location of any gaps in the consensus sequence which can be used when visualising the assembly
* `consensus_probs.hdf` - this is a file that contains the output of the neural network calculations and is not an output for end-users, so we don't need to worry about this file

`consensus.fasta` is the file we're interested in, which in our case is the polished assembly.

> ## BAM and SAM Files
> A [SAM file](https://genome.sph.umich.edu/wiki/SAM), is a tab-delimited text file that contains information for each individual read and its alignment to the genome. While we do not have time to go into detail about the features of the SAM format, the paper by [Heng Li et al.](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full) provides a lot more detail on the specification.  
> The compressed binary version of SAM is called a BAM file. We use this version to reduce size and to allow for indexing, which enables efficient random access of the data contained within the file.  
> The file begins with a header, which is optional. The header is used to describe the source of data, reference sequence, method of alignment, etc., this will change depending on the aligner being used. Following the header is the alignment section. Each line that follows corresponds to alignment information for a single read. Each alignment line has 11 mandatory fields for essential mapping information and a variable number of other fields for aligner specific information. An example entry from a SAM file is displayed below with the different fields highlighted.  
> See [Genomics - Variant Calling](https://cloud-span.github.io/04genomics/01-variant_calling/index.html) for a deeper dive.  
{: .callout}

## Polishing with short reads
We will be using the program Pilon to further polish the draft assembly using the short reads. Something about pilon here...


Bioinformatics programs are not built equally. Some programs like Flye or Medaka will require very few input files as they will generate any that they need within the pipeline. Some programs however, require a lot of user input to generate the input files that are needed.  

Pilon is in the latter group of bioinformatics software, so we will need to do some pre-processing using other programs to create some of the inputs needed.

### Generating the Pilon input files

We will first use the program [BWA](https://github.com/lh3/bwa) to generate an alignment of the raw short reads against the draft genome.

We've previously covered aligning reads to a genome in [Genomics - Variant Calling](https://cloud-span.github.io/04genomics/01-variant_calling/index.html). We will be using very similar commands here.

We first need to index the polished assembly we got from Medaka. Indexing allows the aligner to quickly find potential alignment sites for query sequences in a genome, which saves time during alignment.
~~~
bwa index consensus.fasta
~~~
{: .bash}
This should only take a few seconds to complete so we don't need to run the job in the background.
Once the indexing is complete you should see an output like:
~~~
[bwa_index] Pack FASTA... 0.51 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 5.86 seconds elapse.
[bwa_index] Update BWT... 0.10 sec
[bwa_index] Pack forward-only FASTA... 0.10 sec
[bwa_index] Construct SA from BWT and Occ... 1.81 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index consensus.fasta
[main] Real time: 8.704 sec; CPU: 8.395 sec
~~~
{: .output}

We can then align the short reads to the draft assembly.

It is possible to chain together commands in unix using a process known as "piping". This allows the output from one command to be directly passed to another command for further processing. This is especially useful for situations where you may not need the intermediate file again. To do this we use the pipe `|` character.  

You can find the pipe (<kbd>|</kbd>) character on your keyboard, usually by typing <kbd>⇧ Shift</kbd> + <kbd>\</kbd>.

You can use multiple pipes in one command but data will only go from the left to the right.
I.e.
`command1 | command2 | command3 | .... |`

We will be using two pipes to join three different steps in order to align the raw reads to the draft assembly and then sort this alignment to generate a sorted BAM file.

First we will generate the alignment using BWA mem, then convert the alignment into BAM with `samtools view` and finally sort the alignment with `samtools sort`. We have run each of these commands separately in [Genomics - Variant Calling](https://cloud-span.github.io/04genomics/01-variant_calling/index.html), if you want to remind yourself of what they do in more detail.

~~~
bwa mem -t 4 consensus.fasta ../ERR3152367_sub5_filtered.fastq | samtools view - -Sb | samtools sort - -@4 -o test.bam
~~~
{: .bash}

~~~
~~~
{: .output}

We then need to run the following command to index the aligment. We haven't added this command to the above pipe
~~~
samtools index test.bam
~~~
{: .bash}

~~~
java -Xmx16G -jar $EBROOTPILON/pilon.jar --genome flye_sub5_med.fasta --unpaired test.bam --outdir test_pilon --changes --threads 4
~~~
{: .bash}
