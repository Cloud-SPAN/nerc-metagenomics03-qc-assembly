---
title: "QC polished assembly"
teaching: 30
exercises: 10
questions:
- "Why would be want to QC an assembly?"
- "How can we perform QC on an assembly?"
- "What metrics can we compare between assemblies to understand the quality of an assembly?"
objectives:
- "Understand the terms N50, misassembly and largest contig."  
- "Understand what factors might effect the quality of an assembly."
- "Know how to use seqkit and QUAST to identify the quality of an assembly."
keypoints:
- "The N50 is the contig length of the 50 percentile. Which means that 50% of the contigs are at least this length in the assembly"
- "A misassembly is when a portion of the assembly is incorrectly put back together"
- "The largest contig is the longest contiguous piece in the assembly"
- "Seqkit can generate summary statistics that will tell us the N50, largest contig and the number of gaps"
- "QUAST can generate additional information in a report which can be used to identify misassemblies"
---

## Using seqkit to generate summary statistics of an assembly

When we QC'ed the nanopore reads by quality score, we used [Seqkit](https://bioinf.shenwei.me/seqkit/) and the command `seqkit seq`. This time we will be using the command `seqkit stats` instead. We also covered how to access the help documentation for seqkit in the assembly section of [this lesson](https://cloud-span.github.io/metagenomics01-qc-assembly/02-QC-quality-raw-reads/index.html).

~~~
seqkit stats --help
~~~
{: .bash}


~~~
simple statistics of FASTA/Q files

Tips:
  1. For lots of small files (especially on SDD), use big value of '-j' to
     parallelize counting.

Usage:
  seqkit stats [flags]

Aliases:
  stats, stat

Flags:
  -a, --all                  all statistics, including quartiles of seq length, sum_gap, N50
  -b, --basename             only output basename of files
  -E, --fq-encoding string   fastq quality encoding. available values: 'sanger', 'solexa', 'illumina-1.3+', 'illumina-1.5+', 'illumina-1.8+'. (default "sanger")
  -G, --gap-letters string   gap letters (default "- .")
  -h, --help                 help for stats
  -e, --skip-err             skip error, only show warning message
  -T, --tabular              output in machine-friendly tabular format

Global Flags:
      --alphabet-guess-seq-length int   length of sequence prefix of the first FASTA record based on which seqkit guesses the sequence type (0 for whole seq) (default 10000)
      --id-ncbi                         FASTA head is NCBI-style, e.g. >gi|110645304|ref|NC_002516.2| Pseud...
      --id-regexp string                regular expression for parsing ID (default "^(\\S+)\\s?")
  -w, --line-width int                  line width when outputing FASTA format (0 for no wrap) (default 60)
  -o, --out-file string                 out file ("-" for stdout, suffix .gz for gzipped out) (default "-")
      --quiet                           be quiet and do not show extra information
  -t, --seq-type string                 sequence type (dna|rna|protein|unlimit|auto) (for auto, it automatically detect by the first sequence) (default "auto")
  -j, --threads int                     number of CPUs. (default value: 1 for single-CPU PC, 2 for others) (default 2)
~~~
{: .output}
