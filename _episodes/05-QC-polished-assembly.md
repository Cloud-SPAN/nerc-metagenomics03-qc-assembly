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
- "Use the help documentation to work out an appropriate flag for seqkit"
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

Next, run the command on the original draft assembly (`~/analysis/assembly/assembly.fasta`) to calculate the N50 length and answer the questions below about the output.

> ## Hint: Seqkit stats N50 command
> ~~~
> seqkit stats -a analysis/assembly/assembly.fasta
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


### Generating statistics for all three assemblies

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


> ## Exercise X: Comparing the Assemblies
> Using the seqkit output for all three assemblies, compare the statistics for each of the three assemblies. What has changed across the two rounds of polishing? (From assembly>medaka>pilon)
>
> > ## Solution
> > Between the original assembly and the medaka polished assembly:
> > - Total length, maximum length and average length have all increased as has the N50, the minimum length and GC content have decreased as has the quartile range of lengths.
> > Between the medaka polished assembly and the pilon polished assembly:
> > - The total length, average length, maximum length, Q3 and N50 have all decreased. However they are all still larger than the original, un-polished, assembly.
> {: .solution}
{: .challenge}

While we can compare the basic assembly statistics, these do not tell the full story as there will also have been changes to the overall sequences such as correcting individual base errors. 

## Using metaQUAST to further assess assembly Quality

We will use [MetaQUAST](http://quast.sourceforge.net/metaquast) to further evaluate our metagenomic assemblies. MetaQUAST is based on the QUAST genome quality tool but accounts for high species divesity and misassemblies.

As MetaQUAST assesses the quality of assemblies using alignments to close references we need to determine which references are appropriate for our data. MetaQUAST can automatically select reference genomes to align the assembly too, however it does not always pick the most appropriate references. As we know what organisms make up our metagenome we will be supplying a file containing the references we want to use instead. If you use MetaQUAST on your own data you could use the default references MetaQUAST selects or provide your own if you have an idea what organisms could be in your dataset.

### Making a file to list our reference Metagenomes

We need to generate a text file on the instance to pass to MetaQUAST. There are multiple ways of creating a text file on command line, we will be using the program Nano (no relation to Oxford Nanopore!) here.


> ## Text editors
>
> Text editors," like nano, "notepad" on Windows or "TextEdit" on Mac are used to edit any plain text files. Plain text files are those that contain only characters, not images or formatting.
>
> We are using nano because it is one of the least complex Unix text editors. However, many programmers use Emacs or Vim (both of which require more time to learn), or a graphical editor such as Gedit.
>
> No matter what editor you use, you will need to know where it searches for and saves files. If you start it from the shell, it will (probably) use your current working directory as its default location.
{: .callout}

To use nano we type the command followed by the name of the text file we want to generate.

~~~
nano reference_genomes.txt
~~~
{: .bash}

When you press enter your terminal should change. You should see a white bar at the top with `GNU nano 2.3.1` and some suggested commands at the bottom of the page.
There should also be a white box which indicates where your cursor is.

You should paste the following list of organism names into this file
~~~
Bacillus subtilis
Cryptococcus neoformans
Enterococcus faecalis
Escherichia coli str. K-12 substr. MG1655
Lactobacillus fermentum
Listeria monocytogenes EGD-e
Pseudomonas aeruginosa
Saccharomyces cerevisiae
Salmonella enterica subsp. enterica serovar Typhimurium str. LT2
Staphylococcus aureus
~~~
{: .bash}

**Note on pasting in GitBASH!**

Then press <kbd>Ctrl</kbd>+<kbd>O</kbd> to save the file. You will the be prompted with `File Name to Write: reference_genomes.txt` - as we named the file when we first used the command we don't need to change this name and can press enter to save the file.

> ## Control, Ctrl, or ^ Key
>
> The Control key is also called the "Ctrl" key. There are various ways
> in which using the Control key may be described. For example, you may
> see an instruction to press the <kbd>Ctrl</kbd> key and, while holding it down,
> press the <kbd>X</kbd> key, described as any of:
>
> * `Control-X`
> * `Control+X`
> * `Ctrl-X`
> * `Ctrl+X`
> * `^X`
> * `C-x`
>
> In `nano`, along the bottom of the screen you'll see `^G Get Help ^O WriteOut`.
> This means that you can use <kbd>Ctrl</kbd>-<kbd>G</kbd> to get help and <kbd>Ctrl</kbd>-<kbd>O</kbd> to save your
> file.
>
> If you are using a Mac, you might be more familiar with the `Command` key, which is labelled with a <kbd>⌘</kbd> .
> But you will often use the the `Ctrl` key when working in a Terminal.
{: .callout}

You should then be able to see this file when you `ls`
~~~
ls
~~~
{: .bash}
~~~
output of ls here
~~~
{: .output}

MetaQUAST command
~~~
metaquast.py --references-list reference_genomes.txt assembly/assembly.fasta medaka/consensus.fasta pilon/pilon.fasta
~~~
{: .bash}
