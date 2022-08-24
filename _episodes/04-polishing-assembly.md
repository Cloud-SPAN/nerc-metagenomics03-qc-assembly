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

Medaka has been pre-installed on the instance so, first we can look at the help page for medaka_consensus.
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

From this we can see that `-i` for the input basecalls (meaning Nanopore raw-reads) and `-d` for the assembly are required.  
As Medaka uses recurrent neural networks we need to pick an appropriate model (`-m`) for the data we're using. From the [documentation](https://github.com/nanoporetech/medaka#models), Medaka models are named to indicate i) the pore type, ii) the sequencing device (MinION or PromethION), iii) the basecaller variant, and iv) the basecaller version, with the format: `{pore}_{device}_{caller variant}_{caller version}`. Medaka doesn't offer an exact model for our dataset while it is possible to train a model yourself we will not be doing that here and instead will use the closest available model. This is the model `r941_prom_fast_g303`, so we also need to add that to our command as that isn't medakas default.  
Finally, to speed this step up we need to specify the number of threads with `-t`.

This gives us the command:
~~~
medaka_consensus -i ERR3152367_sub5_filtered.fastq -d assembly.fasta -m r941_prom_fast_g303 -t 4 &> medaka.out &
~~~
{: .bash}
Note, we have added `&> medaka.out &` to redirect the output and run the command in the background. Medaka shouldn't take as long as Flye did in the previous step (probably around 20 mins), but it's a good idea to run things in the background so that you can do other things while the program is running.



## Polishing with short reads
Something about pilon here....

When doing bioinformatics you will come across software that requires a different amount of user input. Some software, such as Flye or Medaka, will accept your files as they come from different step with minimal preprocessing as they do most of the things required "under the hood". However some programs require you to generate particular files or process your files in a different way. Pilon which we will be using in this step requires us to create an indexed BAM file in order to polish the genome with short reads. You will also see BWA which we use to create this BAM file, also requires some pre-process of the medaka polished assembly.

BWA
~~~
bwa index consensus.fasta
~~~
{: .bash}

Introduce piping a command & BWA
~~~
bwa mem -t 4 consensus.fasta ../ERR3152367_sub5_filtered.fastq | samtools view - -Sb | samtools sort - -@4 -o test.bam
~~~
{: .bash}

~~~
samtools index test.bam
~~~
{: .bash}

~~~
java -Xmx16G -jar $EBROOTPILON/pilon.jar --genome flye_sub5_med.fasta --unpaired test.bam --outdir test_pilon --changes --threads 4
~~~
{: .bash}
