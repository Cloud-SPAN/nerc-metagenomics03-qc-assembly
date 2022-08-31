---
title: "QC polished assembly"
teaching: 30
exercises: 10
questions:
- "Why would we QC an assembly?"
- "How can we perform QC on an assembly?"
- "What metrics can we compare between assemblies to understand the quality of an assembly?"
objectives:
- "Understand the terms N50, misassembly and largest contig."
- "Understand what factors might affect the quality of an assembly."
- "Use the help doucumentation
- "Apply seqkit to assess multiple assemblies"
- "Use metaQUAST to identify the quality of an assembly."
keypoints:
- "The N50 is the contig length of the 50 percentile. Which means that 50% of the contigs are at least this length in the assembly"
- "A misassembly is when a portion of the assembly is incorrectly put back together"
- "The largest contig is the longest contiguous piece in the assembly"
- "Seqkit can generate summary statistics that will tell us the N50, largest contig and the number of gaps"
- "metaQUAST can generate additional information in a report which can be used to identify misassemblies"
---

## Using seqkit to generate summary statistics of an assembly

After we finished the draft assembly we used `seqkit stats` to see some basic statistics about the assembly (see the episode on [Assembly](https://cloud-span.github.io/metagenomics01-qc-assembly/03-assembly/index.html)). We will be using it again here to get some more statistics for all three of the assemblies to compare the polishing process.

We can again review the help documentation for seqkit stats.
~~~
seqkit stats --help
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

### The N50 length

As we touched on in a previous lesson, the N50 length is a useful statistic when looking at sequences of varying length as it indicates that 50% of the total sequence is in segments that are that size or larger. See the webpage [What's N50?](https://www.molecularecologist.com/2017/03/29/whats-n50/) for a good explanation.

This is a useful statistic to describe an assembly as it indicates the average size of the contigs the assembly software has produced.

A higher N50 length means that more of the assembly is in longer fragments. That means the chunks of sequence produced by the assembler are, on average, larger.

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

Next, run the command on the original draft assembly (`assembly/assembly.fasta`) to calculate the N50 length and answer the questions below about the output.

> ## Hint: Seqkit stats N50 command
> ~~~
> seqkit stats -a assembly/assembly.fasta
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
>> file          format type num_seqs sum_len    min_len  avg_len    max_len     Q1     Q2       Q3    sum_gap  N50      Q20(%) Q30(%) GC(%)
>> assembly.fasta FASTA  DNA  146      14,953,273  3,164  102,419.7  6,068,630  7,364  13,415.5 35,259 0       2,976,503  0     0     52.48
>> ~~~
>> {: .output}
>> b) Comparing the header line from this command to the original command we can see we've now got statistics for Q1, Q2, Q3, sum_gap, N50, Q20(%) and Q30(%)  
>> c) The N50 length for this assembly is 2,976,503 bp, this tells us that 50% of the assembly is in fragments that are almost 3m bases long or longer!  
>> Bonus: `Q1`, `Q2`, `Q3` is the quartile range of sequence length, `sum_gap` is the total number of ambiguous bases (N's) in the sequence, N50 we have covered, Q20(%) is the percentage of bases with a PHRED score over 20, Q30(%) is the percentage of bases with a PHRED score over 30. GC(%) is the [guanine-cytosine content](https://en.wikipedia.org/wiki/GC-content) of the sequence.   
> {: .solution}
{: .challenge}


Instead of passing one FASTA file to `seqkit stats` at once we will can use all three FASTA files we have generated.

First we need to navigate into the analysis directory.
~~~
cd ~/analysis/
~~~
{: .bash}

Within the analysis directory, these three files are:
* Draft assembly generate by Flye in `assembly/assembly.fasta`
* Long-read polished assembly by Medaka in `medaka/consensus.fasta`
* Short-read polished assembly by Pilon in `pilon/pilon.fasta`

This makes our command:
~~~
seqkit stats -a assembly/assembly.fasta medaka/consensus.fasta pilon/pilon.fasta
~~~
{: .bash}


| file                              | format | type | num_seqs | sum_len  | min_len | avg_len  | max_len | Q1   | Q2      | Q3    | sum_gap | N50     | Q20(%) | Q30(%) | GC(%) |
|-----------------------------------|--------|------|----------|----------|---------|----------|---------|------|---------|-------|---------|---------|--------|--------|-------|
| assembly/assembly.fasta                    | FASTA  | DNA  | 146      | 14953273 | 3164    | 102419.7 | 6068630 | 7364 | 13415.5 | 35259 | 0       | 2976503 | 0      | 0      | 52.48 |
| medaka/consensus.fasta            | FASTA  | DNA  | 146      | 14973646 | 3142    | 102559.2 | 6074419 | 7299 | 13333   | 35173 | 0       | 2991855 | 0      | 0      | 52.4  |
| pilon/pilon.fasta | FASTA  | DNA  | 146      | 14970478 | 3142    | 102537.5 | 6073731 | 7299 | 13333   | 35169 | 0       | 2991264 | 0      | 0      | 52.4  |


> ## Exercise X: Comparing the Polishing
>
>
> > ## solution
> >
> >
> {: .solution}
{: .challenge}
