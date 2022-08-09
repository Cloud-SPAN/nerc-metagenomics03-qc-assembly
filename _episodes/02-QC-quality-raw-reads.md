---
title: "Assessing Read Quality, Trimming and Filtering"
teaching: 30
exercises: 20
questions:
- "How can I describe the quality of my data?"
objectives:
- "Explain how a FASTQ file encodes per-base quality scores."
- "Interpret a FastQC plot summarizing per-base quality across all reads."
- "Use `for` loops to automate operations on multiple files."
keypoints:
- "Quality encodings vary across sequencing platforms."
- "`for` loops let you perform the same set of operations on multiple files with a single command."
- "It is important to know the quality of our data to be able to make decisions in the subsequent steps."
title: "Trimming and Filtering"
teaching: 30
exercises: 25
questions:
- "How can we get rid of sequence data that doesn't meet our quality standards?"
objectives:
- "Clean FASTQ reads using Trimmomatic."
- "Select and set multiple options for command line bioinformatic tools."
- "Write `for` loops with two variables."
keypoints:
- "The options you set for the command-line tools you use are important!"
- "Data cleaning is essential at the beginning of metagenomics workflows."
- "Use Trimmomatic to get clean of reads without adapters or low quality bases."
- "Carefully fill the parameters and options required to call a function in the bash shell."
- "Automate repetitive workflows using for loops"
---


## Quality control

<img align="left" width="325" height="226" src="{{ page.root }}/fig/short_analysis_flowchart_crop1.png" alt="Analysis flow diagram that shows the steps: Sequence reads and Quality control." />

We will now assess the quality of the sequence reads contained in our FASTQ files.
We will be adapting the quality control workflow from [Cloud-SPAN Genomics](https://cloud-span.github.io/00genomics/) for the metagenomics dataset used in this course.  

You may want to revisit [Assessing Read Quality](https://cloud-span.github.io/03genomics/01-quality-control/index.html) or [Trimming and Filtering](https://cloud-span.github.io/03genomics/02-trimming/index.html) to remind yourself of key concepts.

We have two different types of sequencing data (short-read Illumina sequence and long-read Nanopore sequence) available for this metagenome and will be using them both in a hybrid approach to assemble and analyse this metagenome.

Because the two types of sequencing are different in length and quality, we need to use different programs for each of them that are built to handle the different strengths and weaknesses each technology provides.

<br clear="left"/>

> ## Reminder of the FASTQ format
> See [Genomics - Assessing Read Quality](https://cloud-span.github.io/03genomics/01-quality-control/index.html) for a more in-depth reminder about the FASTQ format.
>
> In the [FASTQ file format](https://en.wikipedia.org/wiki/FASTQ_format), each ‘read’ (i.e. sequence) is described in four lines of information:
> 1. The first line always starts with an '@' followed by the sequence identifier (also called the header) and may contain other information about the read such as the length.
> 2. The second line is the sequence of bases itself
> 3. The third line is a separator line which starts with a ‘+’ and may repeat the information from line 1
> 4. The fourth line is a string of characters representing the quality scores for each base
{: .callout}


### Illumina Quality control

We will first quality control the raw Illumina data.
We will be adapting the methods for short reads used in [Genomics - Assessing Read Quality](https://cloud-span.github.io/03genomics/01-quality-control/index.html) with our Illumina short read data.

> Bit about PHRED scores here


Rather than assessing every read in the raw data by hand we can use the program [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to visualise the quality of the whole sequencing file.

FastQC has been pre-installed on your isntance so we can pull up the help documentation for fastqc to remind ourselves of the paramaters avaiable

~~~
$ fastqc -h
~~~
{: .bash}

> ## FastQC seq help documentation
> ~~~
>
>             FastQC - A high throughput sequence QC analysis tool
>
> SYNOPSIS
>
>       fastqc seqfile1 seqfile2 .. seqfileN
>
>     fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
>            [-c contaminant file] seqfile1 .. seqfileN
>
> DESCRIPTION
>
>     FastQC reads a set of sequence files and produces from each one a quality
>     control report consisting of a number of different modules, each one of
>     which will help to identify a different potential type of problem in your
>     data.
>     
>     If no files to process are specified on the command line then the program
>     will start as an interactive graphical application.  If files are provided
>     on the command line then the program will run with no user interaction
>     required.  In this mode it is suitable for inclusion into a standardised
>     analysis pipeline.
>     
>     The options for the program as as follows:
>     
>     -h --help       Print this help file and exit
>     
>     -v --version    Print the version of the program and exit
>     
>     -o --outdir     Create all output files in the specified output directory.
>                     Please note that this directory must exist as the program
>                     will not create it.  If this option is not set then the
>                     output file for each sequence file is created in the same
>                     directory as the sequence file which was processed.
>                     
>     --casava        Files come from raw casava output. Files in the same sample
>                     group (differing only by the group number) will be analysed
>                     as a set rather than individually. Sequences with the filter
>                     flag set in the header will be excluded from the analysis.
>                     Files must have the same names given to them by casava
>                     (including being gzipped and ending with .gz) otherwise they
>                     won't be grouped together correctly.
>                     
>     --nano          Files come from nanopore sequences and are in fast5 format. In
>                     this mode you can pass in directories to process and the program
>                     will take in all fast5 files within those directories and produce
>                     a single output file from the sequences found in all files.                    
>                     
>     --nofilter      If running with --casava then don't remove read flagged by
>                     casava as poor quality when performing the QC analysis.
>                    
>     --extract       If set then the zipped output file will be uncompressed in
>                     the same directory after it has been created.  By default
>                     this option will be set if fastqc is run in non-interactive
>                     mode.
>                     
>     -j --java       Provides the full path to the java binary you want to use to
>                     launch fastqc. If not supplied then java is assumed to be in
>                     your path.
>                    
>     --noextract     Do not uncompress the output file after creating it.  You
>                     should set this option if you do not wish to uncompress
>                     the output when running in non-interactive mode.
>                     
>     --nogroup       Disable grouping of bases for reads >50bp. All reports will
>                     show data for every base in the read.  WARNING: Using this
>                     option will cause fastqc to crash and burn if you use it on
>                     really long reads, and your plots may end up a ridiculous size.
>                     You have been warned!
>                     
>     --min_length    Sets an artificial lower limit on the length of the sequence
>                     to be shown in the report.  As long as you set this to a value
>                     greater or equal to your longest read length then this will be
>                     the sequence length used to create your read groups.  This can
>                     be useful for making directly comaparable statistics from
>                     datasets with somewhat variable read lengths.
>                     
>     -f --format     Bypasses the normal sequence file format detection and
>                     forces the program to use the specified format.  Valid
>                     formats are bam,sam,bam_mapped,sam_mapped and fastq
>                     
>     -t --threads    Specifies the number of files which can be processed
>                     simultaneously.  Each thread will be allocated 250MB of
>                     memory so you shouldn't run more threads than your
>                     available memory will cope with, and not more than
>                     6 threads on a 32 bit machine
>                   
>     -c              Specifies a non-default file which contains the list of
>     --contaminants  contaminants to screen overrepresented sequences against.
>                     The file must contain sets of named contaminants in the
>                     form name[tab]sequence.  Lines prefixed with a hash will
>                     be ignored.
>
>     -a              Specifies a non-default file which contains the list of
>     --adapters      adapter sequences which will be explicity searched against
>                     the library. The file must contain sets of named adapters
>                     in the form name[tab]sequence.  Lines prefixed with a hash
>                     will be ignored.
>                     
>     -l              Specifies a non-default file which contains a set of criteria
>     --limits        which will be used to determine the warn/error limits for the
>                     various modules.  This file can also be used to selectively
>                     remove some modules from the output all together.  The format
>                     needs to mirror the default limits.txt file found in the
>                     Configuration folder.
>                     
>    -k --kmers       Specifies the length of Kmer to look for in the Kmer content
>                     module. Specified Kmer length must be between 2 and 10. Default
>                     length is 7 if not specified.
>                     
>    -q --quiet       Supress all progress messages on stdout and only report errors.
>    
>    -d --dir         Selects a directory to be used for temporary files written when
>                     generating report images. Defaults to system temp directory if
>                     not specified.
>                     
> BUGS
>
>     Any bugs in fastqc should be reported either to simon.andrews@babraham.ac.uk
>     or in www.bioinformatics.babraham.ac.uk/bugzilla/
> ~~~
> {: .output}
{: .solution}

FastQC can accept multiple file names as input, and on both zipped and unzipped files,
so we can use the `\*.fastq*`wildcard to run FastQC on all of the FASTQ files in this directory.

~~~
$ fastqc *.fastq*
~~~
{: .bash}

You will see an automatically updating output message telling you the
progress of the analysis. It will start like this:

~~~
Started analysis of JC1A_R1.fastq.gz                                               
Approx 5% complete for JC1A_R1.fastq.gz                                            
Approx 10% complete for JC1A_R1.fastq.gz                                           
Approx 15% complete for JC1A_R1.fastq.gz                                           
Approx 20% complete for JC1A_R1.fastq.gz                                           
Approx 25% complete for JC1A_R1.fastq.gz                                           
Approx 30% complete for JC1A_R1.fastq.gz                                          
Approx 35% complete for JC1A_R1.fastq.gz  
~~~
{: .output}

In total, it should take about five minutes for FastQC to run on all
four of our FASTQ files. When the analysis completes, your prompt
will return. So your screen will look something like this:

~~~
Approx 80% complete for JP4D_R2.fastq.gz
Approx 85% complete for JP4D_R2.fastq.gz
Approx 90% complete for JP4D_R2.fastq.gz
Approx 95% complete for JP4D_R2.fastq.gz
Analysis complete for JP4D_R2.fastq.gz
$
~~~
{: .output}

The FastQC program has created several new files within our
`data/untrimmed_fastq/` directory.

~~~
$ ls
~~~
{: .bash}

~~~
JC1A_R1_fastqc.html             JP4D_R1.fastq                   
JC1A_R1_fastqc.zip              JP4D_R1_fastqc.html
JC1A_R1.fastq.gz                JP4D_R1_fastqc.zip                   
JC1A_R2_fastqc.html             JP4D_R2_fastqc.html           
JC1A_R2_fastqc.zip              JP4D_R2_fastqc.zip
JC1A_R2.fastq.gz                JP4D_R2.fastq.gz       
~~~
{: .output}


For each input FASTQ file, FastQC has created a `.zip` file and a
`.html` file. The `.zip` file extension indicates that this is
actually a compressed set of multiple output files. We'll be working
with these output files soon. The `.html` file is a stable webpage
displaying the summary report for each of our samples.


### Nanopore quality control

Now we will be assessing the quality of the Nanopore raw reads which are in the file  `~/cs_course/data/nano_fastq/ERR3152367_sub5.fastq`.

As before, we can view the first complete read in one of the files from our dataset by using `head` to look at the first four lines.

~~~
 cd ~/cs_course/data/nano_fastq/
 head -n 4 ERR3152367_sub5.fastq
~~~
{: .bash}

~~~
@ERR3152367.34573250 d8c83b24-b46e-4f1a-836f-768f835acf68 length=320
GGTTGGTTATGTGCATGTTTTCAGTTACATATTGCATCTGTGGGAGCATATTCTTGTTTATGGGTTATGTGTTGGTGGTTGCATGTGGTGTGTTGTTGTGTTAACAAGTGTGGAACCTGTTCATTGGGTTATGAACAACGACACAAGTGTTGCGTGTTGAGCTAGTTAACGTGTGTGTTGTTATTCTTCTGAACCAGTTAACTTATTTGTTTTGTTGGGTGTGAAGCAGTGGGCGTGAAGGTGAGCGATGAAGCGGCGTTGTTCTGTTGCGTGTTTGATTGTGTTGTGTTGCGTGAAGAAGCGTCGTTGTTGGGTGGTTC
+
$$##$$###$#%###%##$%%$$###$#$$#$%##%&$$$$$$%#$$$$#$%#%$##$#$%#%$$#$$$%#$$#$%$$$$$#$%#$#$%$$$##$%%#&$#$#$$$$$%$$%$$%%$$#"$#$$$#&$$$$$#$$$$$######$#$#$$###$%###$$$$%$$&%$$$#$#$$%#%$##$##%#$&$$$$$#$$$%$$$##%#%$##$%%$$#$$$$%#%$###$$$####%$%%$$'$$%$$$$$%$#$$&$$%$#$##$%%$$%$$%%$%&'##$##%$#$$%$###$$$$$#$$$$#$&%##$$#$$%$$$%###
~~~
{: .output}


We can see that this read is longer than the Illumina reads we looked at earlier. The length of a raw read from Nanopore sequencing varies depends on the length of the length of the DNA strand being sequenced.


Line 4 shows us the quality score of this read.

~~~
$$##$$###$#%###%##$%%$$###$#$$#$%##%&$$$$$$%#$$$$#$%#%$##$#$%#%$$#$$$%#$$#$%$$$$$#$%#$#$%$$$##$%%#&$#$#$$$$$%$$%$$%%$$#"$#$$$#&$$$$$#$$$$$######$#$#$$###$%###$$$$%$$&%$$$#$#$$%#%$##$##%#$&$$$$$#$$$%$$$##%#%$##$%%$$#$$$$%#%$###$$$####%$%%$$'$$%$$$$$%$#$$&$$%$#$##$%%$$%$$%%$%&'##$##%$#$$%$###$$$$$#$$$$#$&%##$$#$$%$$$%###
~~~
{: .output}

Based on the PHRED quality scores, see [Genomics - Quality Control](https://cloud-span.github.io/03genomics/01-quality-control/index.html) for a reminder,
we can see that the quality score of the bases in this read are between 1-10.

> ## PHRED score reminder
>~~~
>Quality encoding: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ
>                   |         |         |         |         |
>Quality score:    01........11........21........31........41   
>~~~
> {: .output}                           
{: .callout}

Rather than looking at every sequence by hand we are going to use a program called [NanoPlot](https://github.com/wdecoster/NanoPlot), which is preinstalled on the instance, to create some plots for the whole sequencing file.

We first need to navigate to the `qc` directory we made earlier `cs_course/results/qc`.
~~~
cd ~/cs_course/results/qc/
~~~
{: .bash}

We are now going to run `NanoPlot` with the raw Nanopore sequencing file.
First we can look at the help documenation for NanoPlot to see what options are available.
~~~
NanoPlot --help
~~~
{: .bash}


> ## NanoPlot Help Documentation
> ~~~
> usage: NanoPlot [-h] [-v] [-t THREADS] [--verbose] [--store] [--raw] [--huge]m --downsample 10000
>                 [-o OUTDIR] [-p PREFIX] [--tsv_stats] [--maxlength N]
>                 [--minlength N] [--drop_outliers] [--downsample N]
>                 [--loglength] [--percentqual] [--alength] [--minqual N]
>                 [--runtime_until N] [--readtype {1D,2D,1D2}] [--barcoded]
>                 [--no_supplementary] [-c COLOR] [-cm COLORMAP]
>                 [-f {eps,jpeg,jpg,pdf,pgf,png,ps,raw,rgba,svg,svgz,tif,tiff}]
>                 [--plots [{kde,hex,dot,pauvre} [{kde,hex,dot,pauvre} ...]]]
>                 [--listcolors] [--listcolormaps] [--no-N50] [--N50]
>                 [--title TITLE] [--font_scale FONT_SCALE] [--dpi DPI]
>                 [--hide_stats]
>                 (--fastq file [file ...] | --fasta file [file ...] | --fastq_rich file [file ...] | --fastq_minimal file [file ...] | --summary file [file ...] | --bam file [file ...] | --ubam file [file ...] | --cram file [file ...] | --pickle pickle | --feather file [file ...])
>
> CREATES VARIOUS PLOTS FOR LONG READ SEQUENCING DATA.
>
> General options:
>   -h, --help            show the help and exit
>   -v, --version         Print version and exit.
>   -t, --threads THREADS
>                         Set the allowed number of threads to be used by the script
>   --verbose             Write log messages also to terminal.
>   --store               Store the extracted data in a pickle file for future plotting.
>   --raw                 Store the extracted data in tab separated file.
>   --huge                Input data is one very large file.
>   -o, --outdir OUTDIR   Specify directory in which output has to be created.
>   -p, --prefix PREFIX   Specify an optional prefix to be used for the output files.
>   --tsv_stats           Output the stats file as a properly formatted TSV.
>
> Options for filtering or transforming input prior to plotting:
>   --maxlength N         Hide reads longer than length specified.
>   --minlength N         Hide reads shorter than length specified.
>   --drop_outliers       Drop outlier reads with extreme long length.
>   --downsample N        Reduce dataset to N reads by random sampling.
>   --loglength           Additionally show logarithmic scaling of lengths in plots.
>   --percentqual         Use qualities as theoretical percent identities.
>   --alength             Use aligned read lengths rather than sequenced length (bam mode)
>   --minqual N           Drop reads with an average quality lower than specified.
>   --runtime_until N     Only take the N first hours of a run
>   --readtype {1D,2D,1D2}
>                         Which read type to extract information about from summary. Options are 1D, 2D,
>                         1D2
>   --barcoded            Use if you want to split the summary file by barcode
>   --no_supplementary    Use if you want to remove supplementary alignments
>
> Options for customizing the plots created:
>   -c, --color COLOR     Specify a valid matplotlib color for the plots
>   -cm, --colormap COLORMAP
>                         Specify a valid matplotlib colormap for the heatmap
>   -f, --format {eps,jpeg,jpg,pdf,pgf,png,ps,raw,rgba,svg,svgz,tif,tiff}
>                         Specify the output format of the plots.
>   --plots [{kde,hex,dot,pauvre} [{kde,hex,dot,pauvre} ...]]
>                         Specify which bivariate plots have to be made.
>   --listcolors          List the colors which are available for plotting and exit.
>   --listcolormaps       List the colors which are available for plotting and exit.
>   --no-N50              Hide the N50 mark in the read length histogram
>   --N50                 Show the N50 mark in the read length histogram
>   --title TITLE         Add a title to all plots, requires quoting if using spaces
>   --font_scale FONT_SCALE
>                         Scale the font of the plots by a factor
>   --dpi DPI             Set the dpi for saving images
>   --hide_stats          Not adding Pearson R stats in some bivariate plots
>
> Input data sources, one of these is required.:
>   --fastq file [file ...]
>                         Data is in one or more default fastq file(s).
>   --fasta file [file ...]
>                         Data is in one or more fasta file(s).
>   --fastq_rich file [file ...]
>                         Data is in one or more fastq file(s) generated by albacore, MinKNOW or guppy
>                         with additional information concerning channel and time.
>   --fastq_minimal file [file ...]
>                         Data is in one or more fastq file(s) generated by albacore, MinKNOW or guppy
>                         with additional information concerning channel and time. Is extracted swiftly
>                         without elaborate checks.
>   --summary file [file ...]
>                         Data is in one or more summary file(s) generated by albacore or guppy.
>   --bam file [file ...]
>                         Data is in one or more sorted bam file(s).
>   --ubam file [file ...]
>                         Data is in one or more unmapped bam file(s).
>   --cram file [file ...]
>                         Data is in one or more sorted cram file(s).
>   --pickle pickle       Data is a pickle file stored earlier.
>   --feather file [file ...]
>                         Data is in one or more feather file(s).
>
> EXAMPLES:
>     NanoPlot --summary sequencing_summary.txt --loglength -o summary-plots-log-transformed
>     NanoPlot -t 2 --fastq reads1.fastq.gz reads2.fastq.gz --maxlength 40000 --plots hex dot
>     NanoPlot --color yellow --bam alignment1.bam alignment2.bam alignment3.bam --downsample 10000
> ~~~
> {: .output}
{: .solution}


As our data is in FASTQ format we are going to use the `--fastq` flag to specify the file, we are also going to use `--outdir` to specify an output directory, we're also going to use the flag `--loglength` to produce the plots with a log scale and finally we're going to use `--threads` to run the program on more than one thread to speed it up.

For the `--fastq` flag: The raw Nanopore data is in the location `/cs_workshop/data/nano_fastq/ERR3152367_sub5.fastq`, so we are going to use the absolute path in the NanoPlot command.

For the `--outdir` flag: As we are already in our `qc` directory we are going to specify `nano_qc` so that NanoPlot will create a directory within this directory to put the files it generates. (Note: with NanoPlot you don't need to create this directory before running the command, however this depends on the program you are using.)

For the `--threads` flag: we are going to run this on 4 threads to allow NanoPlot to use more compute power to speed it up.

The `--loglength` flag doesn't require any additional information.

~~~
NanoPlot --fastq ~/cs_course/data/nano_fastq/ERR3152367_sub5.fastq --outdir nano_qc --threads 4 --loglength
~~~
{: .bash}

Now we have the command set up we can press enter and wait for NanoPlot to finish.

This will take a **_couple of minutes_**, you will know it is finished once your cursor has returned (i.e. you can type in the terminal again).  



Once NanoPlot has finished we can have a look at the output.
First we need to navigate into the directory NanoPlot created, then list the files.
~~~
cd nano_qc
ls
~~~
{: .bash}
~~~
Dynamic_Histogram_Read_length.html            NanoPlot-report.html
HistogramReadlength.png                       NanoPlot_20220804_1931.log
LengthvsQualityScatterPlot_dot.png            NanoStats.txt
LengthvsQualityScatterPlot_kde.png            Weighted_HistogramReadlength.png
LengthvsQualityScatterPlot_loglength_dot.png  Weighted_LogTransformed_HistogramReadlength.png
LengthvsQualityScatterPlot_loglength_kde.png  Yield_By_Length.png
LogTransformed_HistogramReadlength.png
~~~
{: .output}


We can see that NanoPlot has generated a lot of different files.

As most of these are image or HTML files we won't be able to view them using terminal - luckily the `NanoPlot-report.html` file contains all of the plots and information held in the other files so we only need to download that one onto our local computer. To do this we will use `scp` which we have used in previous modules (see [Genomics - Quality Control](https://cloud-span.github.io/03genomics/01-quality-control/index.html)).

In a new terminal window that's **_not_** logged into the instance, navigate to your Cloud-SPAN directory (that contains your pem file) using `cd`.
Once you're in the directory you want to download this file into we will use `scp` to download the file.

The command will look something like:
~~~
scp -i login-key-instanceNNN.pem csuser@instanceNNN.cloud-span.aws.york.ac.uk:~/cs_course/results/qc/nano_qc/NanoPlot-report.html .
~~~
{: .bash}
Remember to replace NNN with the instance number specific to you.
As the file is downloading you will see an output like:
~~~
TO FILL
~~~
{: .output}

Once the file has downloaded, using your file system (e.g. File explorer or Finder) you can find the file and double click it to open.
As this is a HTML file it should open up in your browser.  

If you had trouble downloading the file you can view it here [NanoPlot-report.html]({{ page.root }}/files/NanoPlot-report.html)

In the report we have Summary Statistics followed by plots showing the distribution of read lengths and also the read length against the quality of the reads.

Looking at the Summary Statistics table answer the following questions:

> ## Exercise X:
>
> 1. How many sequences are in this file?
> 2. How many bases are there in this entire file?
> 3. What is the length of the longest read in the file and its associated mean quality score?
>
>> ## Solution
>>   1. There are 692,758 sequences (also known as reads) in this file
>>   2. There are 3,082,258,211 bases (bp) in total in this FASTQ file
>>   3. The longest read in this file is 413,847 bp and it has a mean quality score of 3.7
> {: .solution}
{: .challenge}

> ## Quality Encodings Vary
>
> Although we've used a particular quality encoding system to demonstrate interpretation of
> read quality, different sequencing machines use different encoding systems. This means that,
> depending on which sequencer you use to generate your data, a `#` may not be an indicator of
> a poor quality base call.
>
> This mainly relates to older Solexa/Illumina data,
> but it's essential that you know which sequencing platform was
> used to generate your data, so that you can tell your quality control program which encoding
> to use. If you choose the wrong encoding, you run the risk of throwing away good reads or
> (even worse) not throwing away bad reads!
{: .callout}

> ## N50
> The N50 value is a useful statistic when looking at sequences of varying length as it indicates that 50% of the total sequence is in reads (i.e. chunks) that are that size or larger.
> For this FASTQ file 50% of the total bases are in reads that have a length of 5,373 bp or longer.
> See the webpage [What's N50?](https://www.molecularecologist.com/2017/03/29/whats-n50/) for a good explanation.
>We will be coming back to this statistic in more detail when we get to the assembly step.
{: .callout}

We can also look at some of the plots produced by NanoPlot.  
One useful plot is the plot titled
### Read lengths vs Average read quality plot using dots after log transformation of read lengths

<img align="left" width="816" height="785" src="{{ page.root }}/fig/02_lengthvsquality_log.png" alt="NanoPlot KDE plot with the title Read lengths vs Average read quality plot using dots after log transformation of read lengths">

This plot shows the read length of the sequences when compared to the average quality of the sequence.  
We can see that the majority of the sequences have a quality score of 4 and above, and many of those with an average quality score of 4 are shorter in length.
This means that for this dataset we can remove those with a lower quality score in order to improve the overall quality of the raw sequences before assembling the metagenome.

<br clear="left"/>

## Filtering Nanopore sequences by quality

We can use the program `seqkit seq` to create a new file containing only the sequences with an average quality above a certain value.

After returning to our home directory, we can view the `seqkit seq` help documentation with the following command:

~~~
cd ~/cs_course/
seqkit seq -h
~~~
{: .bash}

> ## Seqkit seq help documentation
> ~~~
> transform sequences (extract ID, filter by length, remove gaps...)
>
> Usage:
>   seqkit seq [flags]
>
> Flags:
>   -k, --color                     colorize sequences - to be piped into "less -R"
>   -p, --complement                complement sequence, flag '-v' is recommended to switch on
>       --dna2rna                   DNA to RNA
>   -G, --gap-letters string        gap letters (default "- \t.")
>   -h, --help                      help for seq
>   -l, --lower-case                print sequences in lower case
>   -M, --max-len int               only print sequences shorter than the maximum length (-1 for no limit) (default -1)
>   -R, --max-qual float            only print sequences with average quality less than this limit (-1 for no limit) (default -1)
>   -m, --min-len int               only print sequences longer than the minimum length (-1 for no limit) (default -1)
>   -Q, --min-qual float            only print sequences with average quality qreater or equal than this limit (-1 for no limit) (default -1)
>   -n, --name                      only print names
>   -i, --only-id                   print ID instead of full head
>   -q, --qual                      only print qualities
>   -b, --qual-ascii-base int       ASCII BASE, 33 for Phred+33 (default 33)
>   -g, --remove-gaps               remove gaps
>   -r, --reverse                   reverse sequence
>       --rna2dna                   RNA to DNA
>   -s, --seq                       only print sequences
>   -u, --upper-case                print sequences in upper case
>   -v, --validate-seq              validate bases according to the alphabet
>   -V, --validate-seq-length int   length of sequence to validate (0 for whole seq) (default 10000)
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


From this we can see that the flag `-Q` will `only print sequences with average quality qreater or equal than this limit (-1 for no limit) (default -1)`.

From the plot above we identified that many of the lower quality reads below 4 were shorter _more here_ so we should set the minimum limit to 4.

~~~
seqkit seq -Q 4 data/nano_fastq/ERR3152367_sub5.fastq > data/nano_fastq/ERR3152367_sub5_filtered.fastq
~~~
{: .bash}

We are using redirecting (`>`) to generate a new file `data/nano_fastq/ERR3152367_sub5_filtered.fastq` containing only the reads with an average quality of 7 or above.

We can now re-run NanoPlot on the filtered file to see how it has changed.

~~~
cd results/qc/ # move into the qc directory

NanoPlot --fastq ~/cs_course/data/nano_fastq/ERR3152367_sub5_filtered.fastq --outdir nano_qc_trim --threads 4 --loglength
~~~
{: .bash}

Once again, wait for the command to finish and then `scp` the `NanoPlot-report.html` to your local computer.

If you had trouble downloading the file you can view it here [NanoPlot-filtered-report.html]({{ page.root }}/files/NanoPlot-filtered-report.html)

<img align="left" width="816" height="785" src="{{ page.root }}/fig/02_lengthvsquality_filtered_log.png" alt="NanoPlot KDE plot of the filtered raw reads Read lengths vs Average read quality plot using dots after log transformation of read lengths">
<br clear="left"/>

**Compare the NanoPlot statistics of the Nanopore raw reads [before filtering]({{ page.root }}/files/NanoPlot-report.html) and [after filtering]({{ page.root }}/files/NanoPlot-filtered-report.html)  and answer the questions below.**

> ## Exercise X:
>
> 1. How many reads have been removed by filtering?
> 2. How many bases have been removed by filtering?
> 3. What is the length of the new longest read and its associated average quality score?
>
>> ## Solution
>>   1. Initially there were 692,758 reads in the filtered file there are 666,597 reads so 26,161 reads have been removed by the quality filtering
>>   2. Initially there were 3,082,258,211 bases and after filtering there are 	3,023,658,929 base which means filtering has removed 58,599,282 bases
>>   3. The longest read in the filtered file is 229,804bp and it has a mean quality score of 6.7
> {: .solution}
{: .challenge}
