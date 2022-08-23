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
## Polishing an assembly with long reads
Something about medaka here...

We're first going to use the filtered raw long reads to polish the draft Flye assembly.
As with the assembly, we need to use polishing software that is especially written for long read raw reads.
[Medaka](https://github.com/nanoporetech/medaka) is a command line tool built by Oxford Nanopore Technologies to generate a consensus sequence from  




~~~
medaka_consensus -i ERR3152367_sub5_filtered.fastq -d assembly.fasta -o flye_sub5_medtest_annie_full -t 12
~~~
{: .bash}

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
