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

## Bioinformatic workflows


<img align="right" width="325" height="506" src="{{ page.root }}/fig/short_analysis_flowchart.png" alt="Flow diagram that shows the steps: Sequence reads, Quality control, Assembly, Binning and Taxonomy" />

When working with high-throughput sequencing data, the raw reads you get off of the sequencer need to pass
through a number of  different tools in order to generate your final desired output.  

The use of this set of tools in a specified order is commonly referred to as a *workflow* or a *pipeline*.  

Here is an example of the workflow we will be using for our analysis with a brief
description of each step.  

1. Quality control - Assessing quality using FastQC and Trimming and/or filtering reads (if necessary)
2. Metagenome assembly
3. Binning
4. Taxonomic assignment

Workflows in bioinformatics often adopt a plug-and-play approach so the output of one tool can be easily used as input to another tool.
The use of standard data formats in bioinformatics (such as FASTA or FASTQ, which we will be using here) makes this possible.
The tools that are used to analyze data at different stages of the workflow are therefore built under the assumption that the data will be provided in a specific format.

<br clear="right"/>

## Quality control

<img align="left" width="325" height="226" src="{{ page.root }}/fig/short_analysis_flowchart_crop1.png" alt="Flow diagram that shows the steps: Sequence reads and Quality control." />


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


### Add Illumina FASTQC here

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
> usage: NanoPlot [-h] [-v] [-t THREADS] [--verbose] [--store] [--raw] [--huge]
> ~~~
> {: .output}
{: .solution}


As our data is in FASTQ format we are going to use the `--fastq` flag to specify the file, we are also going to use `--outdir` to specify an output directory and finally we're going to use `--threads` to run the program on more than one thread to speed it up.

For the `--fastq` flag: The raw Nanopore data is in the location `/cs_workshop/data/nano_fastq/ERR3152367_sub5.fastq`, so we are going to use the absolute path in the NanoPlot command.

For the `--outdir` flag: As we are already in our `qc` directory we are going to specify `nano_qc` so that NanoPlot will create a directory within this directory to put the files it generates. (Note: with NanoPlot you don't need to create this directory before running the command, however this depends on the program you are using.)

For the `--threads` flag: we are going to run this on 4 threads to allow NanoPlot to use more compute power to speed it up.

~~~
NanoPlot --fastq ~/cs_course/data/nano_fastq/ERR3152367_sub5.fastq --outdir nano_qc --threads 4
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
Dynamic_Histogram_Read_length.html      NanoPlot_20220804_1454.log
HistogramReadlength.png                 NanoStats.txt
LengthvsQualityScatterPlot_dot.png      Weighted_HistogramReadlength.png
LengthvsQualityScatterPlot_kde.png      Weighted_LogTransformed_HistogramReadlength.png
LogTransformed_HistogramReadlength.png  Yield_By_Length.png
NanoPlot-report.html
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

If you had trouble downloading the file you can view it here [NanoPlot-report.html]({{ page.root }}/data/NanoPlot-report.html)

In the report we have Summary Statistics for the raw sequences and some plots.

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
>>   2. The longest read in this file is 413,847 bp and it has a mean quality score of 3.7
> {: .solution}
{: .challenge}

The N50 is a useful statistic when looking at sequences of varying length as it indicates that 50% of the total sequence is in reads (i.e. chunks) that are that size or larger. So for this FASTQ file 50% of the total bases are in reads that have a length of 5,373 bp or longer.
See the webpage [What's N50?](https://www.molecularecologist.com/2017/03/29/whats-n50/) for a good explanation.

We will be coming back to this statistic in more detail when we get to the assembly step.

We can also look at some of the plots produced by NanoPlot.
One useful plot is "Read lengths vs Average read quality plot using a kernel density estimation"

<img align="left" width="510" height="490" src="{{ page.root }}/fig/02_lengthvsquality_all.png">

This is the plot that we generated using NanoPlot above, however because there's a few of long reads
<br clear="left"/>

<img align="right" width="500" height="500" src="{{ page.root }}/fig/02_lengthvsquality_small.png">

<br clear="right"/>


At this point, lets validate that all the relevant tools are installed. If you are using the AWS AMI then these _should_ be preinstalled.

~~~
$ fastqc -h
~~~
{: .bash}

~~~
            FastQC - A high throughput sequence QC analysis tool

SYNOPSIS

        fastqc seqfile1 seqfile2 .. seqfileN

    fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
           [-c contaminant file] seqfile1 .. seqfileN

DESCRIPTION

    FastQC reads a set of sequence files and produces from each one a quality
    control report consisting of a number of different modules, each one of
    which will help to identify a different potential type of problem in your
    data.

    If no files to process are specified on the command line then the program
    will start as an interactive graphical application.  If files are provided
    on the command line then the program will run with no user interaction
    required.  In this mode it is suitable for inclusion into a standardised
    analysis pipeline.

    The options for the program as as follows:

    -h --help       Print this help file and exit

    -v --version    Print the version of the program and exit

    -o --outdir     Create all output files in the specified output directory.                                                                                    
                    Please note that this directory must exist as the program
                    will not create it.  If this option is not set then the
                    output file for each sequence file is created in the same
                    directory as the sequence file which was processed.

    --casava        Files come from raw casava output. Files in the same sample
                    group (differing only by the group number) will be analysed
                    as a set rather than individually. Sequences with the filter
                    flag set in the header will be excluded from the analysis.
                    Files must have the same names given to them by casava
                    (including being gzipped and ending with .gz) otherwise they
                    won't be grouped together correctly.

    --nano          Files come from naopore sequences and are in fast5 format. In
                    this mode you can pass in directories to process and the program
                    will take in all fast5 files within those directories and produce
                    a single output file from the sequences found in all files.

    --nofilter      If running with --casava then don't remove read flagged by
                    casava as poor quality when performing the QC analysis.

    --extract       If set then the zipped output file will be uncompressed in
                    the same directory after it has been created.  By default
                    this option will be set if fastqc is run in non-interactive
                    mode.

    -j --java       Provides the full path to the java binary you want to use to
                    launch fastqc. If not supplied then java is assumed to be in
                    your path.

    --noextract     Do not uncompress the output file after creating it.  You
                    should set this option if you do not wish to uncompress
                    the output when running in non-interactive mode.

    --nogroup       Disable grouping of bases for reads >50bp. All reports will
                    show data for every base in the read.  WARNING: Using this
                    option will cause fastqc to crash and burn if you use it on
                    really long reads, and your plots may end up a ridiculous size.
                    You have been warned!

    -f --format     Bypasses the normal sequence file format detection and
                    forces the program to use the specified format.  Valid
                    formats are bam,sam,bam_mapped,sam_mapped and fastq

    -t --threads    Specifies the number of files which can be processed
                    simultaneously.  Each thread will be allocated 250MB of
                    memory so you shouldn't run more threads than your
                    available memory will cope with, and not more than
                    6 threads on a 32 bit machine

    -c              Specifies a non-default file which contains the list of
    --contaminants  contaminants to screen overrepresented sequences against.
                    The file must contain sets of named contaminants in the
                    form name[tab]sequence.  Lines prefixed with a hash will
                    be ignored.

    -a              Specifies a non-default file which contains the list of
    --adapters      adapter sequences which will be explicity searched against
                    the library. The file must contain sets of named adapters
                    in the form name[tab]sequence.  Lines prefixed with a hash
                    will be ignored.

    -l              Specifies a non-default file which contains a set of criteria
    --limits        which will be used to determine the warn/error limits for the
                    various modules.  This file can also be used to selectively
                    remove some modules from the output all together.  The format
                    needs to mirror the default limits.txt file found in the
                    Configuration folder.

   -k --kmers       Specifies the length of Kmer to look for in the Kmer content
                    module. Specified Kmer length must be between 2 and 10. Default
                    length is 7 if not specified.

   -q --quiet       Supress all progress messages on stdout and only report errors.

   -d --dir         Selects a directory to be used for temporary files written when
                    generating report images. Defaults to system temp directory if
                    not specified.

BUGS

    Any bugs in fastqc should be reported either to simon.andrews@babraham.ac.uk
    or in www.bioinformatics.babraham.ac.uk/bugzilla/
~~~
{: .output}

If FastQC is not installed then you would expect to see an error like
~~~
$ fastqc -h
~~~
{: .bash}

~~~
The program 'fastqc' is currently not installed. You can install it by typing:
sudo apt-get install fastqc
~~~
{: .error}

If this happens check with your instructor before trying to install it.

## Assessing quality using FastQC
In real life, you won't be assessing the quality of your reads by visually inspecting your
FASTQ files. Rather, you'll be using a software program to assess read quality and
filter out poor quality reads. We'll first use a program called [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to visualize the quality of our reads.
Later in our workflow, we'll use another program to filter out poor quality reads.

FastQC has a number of features which can give you a quick impression of any problems your
data may have, so you can take these issues into consideration before moving forward with your
analyses. Rather than looking at quality scores for each individual read, FastQC looks at
quality collectively across all reads within a sample. The image below shows one FastQC-generated plot that indicates a very high quality sample:

 <a href="{{ page.root }}/fig/03-02-03.png">
  <img src="{{ page.root }}/fig/03-02-03.png" alt="Graphic of boxplots, where are all of them are at the top of the y axis in the good range of scores." />
</a>

The x-axis displays the base position in the read, and the y-axis shows quality scores. In this
example, the sample contains reads that are 40 bp long. This is much shorter than the reads we
are working with in our workflow. For each position, there is a box-and-whisker plot showing
the distribution of quality scores for all reads at that position. The horizontal red line
indicates the median quality score and the yellow box shows the 1st to
3rd quartile range. This means that 50% of reads have a quality score that falls within the
range of the yellow box at that position. The whiskers show the absolute range, which covers
the lowest (0th quartile) to highest (4th quartile) values.

For each position in this sample, the quality values do not drop much lower than 32. This
is a high quality score. The plot background is also color-coded to identify good (green),
acceptable (yellow), and bad (red) quality scores.

Now let's take a look at a quality plot on the other end of the spectrum.

 <a href="{{ page.root }}/fig/03-02-04.png">
  <img src="{{ page.root }}/fig/03-02-04.png" alt="Graphic of boxplots, where the first ones are in the good range of scores of the y axis and extend to the acceptable and bad ranges of scores toward the right of the x axis" />
</a>

Here, we see positions within the read in which the boxes span a much wider range. Also, quality scores drop quite low into the "bad" range, particularly on the tail end of the reads. The FastQC tool produces several other diagnostic plots to assess sample quality, in addition to the one plotted above.

## Running FastQC  

We will now assess the quality of the reads that we downloaded. First, make sure you're still in the `untrimmed_fastq` directory

~~~
$ cd ~/dc_workshop/data/untrimmed_fastq/
~~~
{: .bash}

> ## Exercise 2: Looking at files metadata
>
> How would you see the size of the files in the `untrimmed_fastq\` directory?  
> (Hint: Look at the options for the `ls` command to see how to show file sizes.)  
> a) `ls -a`  
> b) `ls -S`  
> c) `ls -l`  
> d) `ls -lh`  
> e) `ls -ahls`  
>   
>> ## Solution
>>  
>> ~~~
>>   
>> a) No. The flag `-a` shows all of the contents, including hidden files and directories.  
>> b) No. The flag `-S` shows the content Sorted by size starting with the largest file.  
>> c) Yes. The flag `-l` shows the contents with metadata including file size.    
>> d) Yes. The flag `-lh` shows the content  with metadata in a human readable manner.  
>> e) Yes. The combination of all of the flags shows all of the contents with metadata including hidden files, sorted by size.  
>>
>> ~~~
>> {: .bash}
>>
>> ~~~
>> -rw-r--r-- 1 dcuser dcuser  24M Nov 26 21:34 JC1A_R1.fastq.gz                      
>> -rw-r--r-- 1 dcuser dcuser  24M Nov 26 21:34 JC1A_R2.fastq.gz                      
>> -rw-r--r-- 1 dcuser dcuser 616M Nov 26 21:34 JP4D_R1.fastq              
>> -rw-r--r-- 1 dcuser dcuser 203M Nov 26 21:35 JP4D_R2.fastq.gz   
>> ~~~
>> {: .output}
>>
>> There are four FASTQ files ranging from 24M (24MB) to 616M.
>>
> {: .solution}
{: .challenge}

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

We want to keep our data files and our results files separate, so we
will move these
output files into a new directory within our `results/` directory.

~~~
$ mkdir -p ~/dc_workshop/results/fastqc_untrimmed_reads
$ mv *.zip ~/dc_workshop/results/fastqc_untrimmed_reads/
$ mv *.html ~/dc_workshop/results/fastqc_untrimmed_reads/
~~~
{: .bash}

Now we can navigate into this results directory and do some closer
inspection of our output files.

~~~
$ cd ~/dc_workshop/results/fastqc_untrimmed_reads/
~~~
{: .bash}

## Viewing the FastQC results

If we were working on our local computers, we'd be able to look at
each of these HTML files by opening them in a web browser. However, these
files are currently sitting on our remote AWS instance, where our local
computer can't see them. Since we are only logging into the AWS instance
via the command line our remote computer it doesn't have any web browser
setup to display these files either. So, the easiest way to look at these webpage
summary reports is to transfer them to our local computers (i.e. your laptop).
To copy a file from a remote server to our own machines, we will use `scp`,
which we learned yesterday in the Introduction to the Command Line lesson.

First, open a new terminal in you local computer, we will make a new directory
on our computer to store the HTML files
we're transferring. Let's put it on our desktop for now. Open a new
tab in your terminal program (you can use the pull down menu at the
top of your screen or the Cmd+t keyboard shortcut) and type:

~~~
$ mkdir -p ~/Desktop/fastqc_html
~~~
{: .bash}

Now we can transfer our HTML files to our local computer using `scp`.

~~~
$ scp dcuser@ec2-34-238-162-94.compute-1.amazonaws.com:~/dc_workshop/results/fastqc_untrimmed_reads/*.html ~/Desktop/fastqc_html
~~~
{: .bash}

As a reminder, the first part
of the command `dcuser@ec2-34-238-162-94.compute-1.amazonaws.com` is
the address for your remote computer. Make sure you replace everything
after `dcuser@` with your instance number (the one you used to log in).

The second part starts with a `:` and then gives the absolute path
of the files you want to transfer from your remote computer. Don't
forget the `:`. We used a wildcard (`*.html`) to indicate that we want all of
the HTML files.

The third part of the command gives the absolute path of the location
you want to put the files in. This is on your local computer and is the
directory we just created `~/Desktop/fastqc_html`.

You should see a status output like this:

~~~
JC1A_R1_fastqc.html     100%  253KB 320.0KB/s   00:00     
JC1A_R2_fastqc.html     100%  262KB 390.1KB/s   00:00     
JP4D_R1_fastqc.html     100%  237KB 360.8KB/s   00:00     
JP4D_R2_fastqc.html     100%  244KB 385.2KB/s   00:00
~~~
{: .output}

Now we can go to our new directory and open the 4 HTML files.

Depending on your system,
you should be able to select and open them all at once via a right click menu
in your file browser.

> ## Exercise 3: Discuss the quality of sequencing files
>
> Discuss your results with a neighbor. Which sample(s) looks the best
> in terms of per base sequence quality? Which sample(s) look the
> worst?
>
>> ## Solution
>> All of the reads contain usable data, but the quality decreases toward
>> the end of the reads. File JC1A_R2_fastqc shows the lowest quality.
> {: .solution}
{: .challenge}

## Decoding the other FastQC outputs
We've now looked at quite a few "Per base sequence quality" FastQC graphs, but there are nine other graphs that we haven't talked about! Below we have provided a brief overview of interpretations for each of these plots. For more information, please see the FastQC documentation [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)

+ [**Per tile sequence quality**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/12%20Per%20Tile%20Sequence%20Quality.html): the machines that perform sequencing are divided into tiles. This plot displays patterns in base quality along these tiles. Consistently low scores are often found around the edges, but hot spots can also occur in the middle if an air bubble was introduced at some point during the run.
+ [**Per sequence quality scores**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/3%20Per%20Sequence%20Quality%20Scores.html): a density plot of quality for all reads at all positions. This plot shows what quality scores are most common.
+ [**Per base sequence content**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html): plots the proportion of each base position over all of the reads. Typically, we expect to see each base roughly 25% of the time at each position, but this often fails at the beginning or end of the read due to quality or adapter content.
+ [**Per sequence GC content**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html): a density plot of average GC content in each of the reads.  
+ [**Per base N content**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/6%20Per%20Base%20N%20Content.html): the percent of times that 'N' occurs at a position in all reads. If there is an increase at a particular position, this might indicate that something went wrong during sequencing.  
+ [**Sequence Length Distribution**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/7%20Sequence%20Length%20Distribution.html): the distribution of sequence lengths of all reads in the file. If the data is raw, there is often on sharp peak, however if the reads have been trimmed, there may be a distribution of shorter lengths.
+ [**Sequence Duplication Levels**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html): a distribution of duplicated sequences. In sequencing, we expect most reads to only occur once. If some sequences are occurring more than once, it might indicate enrichment bias (e.g. from PCR). If the samples are high coverage (or RNA-seq or amplicon), this might not be true.  
+ [**Overrepresented sequences**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/9%20Overrepresented%20Sequences.html): a list of sequences that occur more frequently than would be expected by chance.
+ [**Adapter Content**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/10%20Adapter%20Content.html): a graph indicating where adapater sequences occur in the reads.
+ [**K-mer Content**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/11%20Kmer%20Content.html): a graph showing any sequences which may show a positional bias within the reads.

## Working with the FastQC text output

Now that we've looked at our HTML reports to get a feel for the data,
let's look more closely at the other output files. Go back to the tab
in your terminal program that is connected to your AWS instance
(the tab label will start with `dcuser@ip`) and make sure you're in
our results subdirectory.   

~~~
$ cd ~/dc_workshop/results/fastqc_untrimmed_reads/
$ ls
~~~
{: .bash}

~~~
JC1A_R1_fastqc.html           JP4D_R1_fastqc.html                  
JC1A_R1_fastqc.zip            JP4D_R1_fastqc.zip                   
JC1A_R2_fastqc.html           JP4D_R2_fastqc.html                  
JC1A_R2_fastqc.zip            JP4D_R2_fastqc.zip
~~~
{: .output}

Our `.zip` files are compressed files. They each contain multiple
different types of output files for a single input FASTQ file. To
view the contents of a `.zip` file, we can use the program `unzip`
to decompress these files. Let's try doing them all at once using a
wildcard.

~~~
$ unzip *.zip
~~~
{: .bash}

~~~
Archive:  JC1A_R1_fastqc.zip                                                       
caution: filename not matched:  JC1A_R2_fastqc.zip                                 
caution: filename not matched:  JP4D_R1_fastqc.zip                      
caution: filename not matched:  JP4D_R2_fastqc.zip  
~~~
{: .output}

This didn't work. It identified the first file and then got a warning
message for each of the other `.zip` files. This is because `unzip`
expects to get only one zip file as input. We could go through and
unzip each file one at a time, but this is very time consuming and
error-prone. Someday you may have 500 files to unzip!

A more efficient way is to use a `for` loop like we learned in the Command Line lesson to iterate through all of
our `.zip` files. Let's see what that looks like and then we'll
discuss what we're doing with each line of our loop.

~~~
$ for filename in *.zip
> do
> unzip $filename
> done
~~~
{: .bash}

In this example, the input is the four filenames (one filename for each of our `.zip` files).
Each time the loop iterates, it will assign a file name to the variable `filename`
and run the `unzip` command.
The first time through the loop,
`$filename` is `JC1A_R1_fastqc.zip`.
The interpreter runs the command `unzip` on `JC1A_R1_fastqc.zip`.
For the second iteration, `$filename` becomes
`JC1A_R2_fastqc.zip`. This time, the shell runs `unzip` on `JC1A_R2_fastqc.zip`.
It then repeats this process for the other `.zip` files in our directory.


When we run our `for` loop, you will see output that starts like this:

~~~
Archive:  JC1A_R1_fastqc.zip                                            
creating: JC1A_R1_fastqc/                                            
creating: JC1A_R1_fastqc/Icons/                                      
creating: JC1A_R1_fastqc/Images/                                    
inflating: JC1A_R1_fastqc/Icons/fastqc_icon.png                      
inflating: JC1A_R1_fastqc/Icons/warning.png                          
inflating: JC1A_R1_fastqc/Icons/error.png                            
inflating: JC1A_R1_fastqc/Icons/tick.png                             
inflating: JC1A_R1_fastqc/summary.txt                                
inflating: JC1A_R1_fastqc/Images/per_base_quality.png                
inflating: JC1A_R1_fastqc/Images/per_tile_quality.png                
inflating: JC1A_R1_fastqc/Images/per_sequence_quality.png            
inflating: JC1A_R1_fastqc/Images/per_base_sequence_content.png       
inflating: JC1A_R1_fastqc/Images/per_sequence_gc_content.png         
inflating: JC1A_R1_fastqc/Images/per_base_n_content.png              
inflating: JC1A_R1_fastqc/Images/sequence_length_distribution.png
inflating: JC1A_R1_fastqc/Images/duplication_levels.png              
inflating: JC1A_R1_fastqc/Images/adapter_content.png                 
inflating: JC1A_R1_fastqc/fastqc_report.html                         
inflating: JC1A_R1_fastqc/fastqc_data.txt                            
inflating: JC1A_R1_fastqc/fastqc.fo  
~~~
{: .output}

The `unzip` program is decompressing the `.zip` files and creating
a new directory (with subdirectories) for each of our samples, to
store all of the different output that is produced by FastQC. There
are a lot of files here. The one we're going to focus on is the
`summary.txt` file.

If you list the files in our directory now you will see:
~~~
$ ls
~~~
{: .bash}

~~~
JC1A_R1_fastqc                  JP4D_R1_fastqc                                                     
JC1A_R1_fastqc.html             JP4D_R1_fastqc.html                                        
JC1A_R1_fastqc.zip              JP4D_R1_fastqc.zip                                                  
JC1A_R2_fastqc                  JP4D_R2_fastqc                                               
JC1A_R2_fastqc.html             JP4D_R2_fastqc.html                                             
JC1A_R2_fastqc.zip              JP4D_R2_fastqc.zip                                         
~~~
{: .output}

The `.html` files and the uncompressed `.zip` files are still present,
but now we also have a new directory for each of our samples. We can
see for sure that it's a directory if we use the `-F` flag for `ls`.

~~~
$ ls -F
~~~
{: .bash}

~~~
JC1A_R1_fastqc/                  JP4D_R1_fastqc/                                                     
JC1A_R1_fastqc.html             JP4D_R1_fastqc.html                                        
JC1A_R1_fastqc.zip              JP4D_R1_fastqc.zip                                                  
JC1A_R2_fastqc/                  JP4D_R2_fastqc/                                               
JC1A_R2_fastqc.html             JP4D_R2_fastqc.html                                             
JC1A_R2_fastqc.zip              JP4D_R2_fastqc.zip                                         
~~~
{: .output}

Let's see what files are present within one of these output directories.

~~~
$ ls -F JC1A_R1_fastqc/
~~~
{: .bash}

~~~
fastqc_data.txt  fastqc.fo  fastqc_report.html	Icons/	Images/  summary.txt
~~~
{: .output}

Use `less` to preview the `summary.txt` file for this sample.

~~~
$ less JC1A_R1_fastqc/summary.txt
~~~
{: .bash}

~~~
PASS    Basic Statistics        JC1A_R1.fastq.gz                     
FAIL    Per base sequence quality       JC1A_R1.fastq.gz             
PASS    Per tile sequence quality       JC1A_R1.fastq.gz             
PASS    Per sequence quality scores     JC1A_R1.fastq.gz             
WARN    Per base sequence content       JC1A_R1.fastq.gz             
FAIL    Per sequence GC content JC1A_R1.fastq.gz                     
PASS    Per base N content      JC1A_R1.fastq.gz                     
PASS    Sequence Length Distribution    JC1A_R1.fastq.gz             
FAIL    Sequence Duplication Levels     JC1A_R1.fastq.gz             
PASS    Overrepresented sequences       JC1A_R1.fastq.gz             
FAIL    Adapter Content JC1A_R1.fastq.gz  
~~~
{: .output}

The summary file gives us a list of tests that FastQC ran, and tells
us whether this sample passed, failed, or is borderline (`WARN`). Remember, to quit from `less` you must type `q`.


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

## Cleaning reads

In the last episode, we took a high-level look at the quality
of each of our samples using `FastQC`. We visualized per-base quality
graphs showing the distribution of the quality at each base across
all the reads from our sample. This information help us to determinate
the quality threshold we will accept and, thus we saw information about
which samples fail which quality checks. Some of our samples failed
quite a few quality metrics used by FastQC. HOwever, this does not mean,
that our samples should be thrown out! It's very common to have some
quality metrics fail, and this may or may not be a problem for your
downstream application. For our workflow, we will be removing some of
the low quality sequences to reduce our false-positive rate due to
sequencing error.

To accomplish this, we will use a program called
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). This
useful tool filters poor quality reads and trims poor quality bases
from the specified samples.

## Trimmomatic options

Trimmomatic has a variety of options to accomplish its task.
If we run the following command, we can see some of its options:

~~~
$ trimmomatic
~~~
{: .bash}

Which will give you the following output:
~~~
Usage:
       PE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...
   or:
       SE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] <inputFile> <outputFile> <trimmer1>...
   or:
       -version
~~~
{: .output}

This output shows us that we must first specify whether we have paired
end (`PE`) or single end (`SE`) reads. Next, we will specify which flags we
would like to run Trimmomatic with. For example, you can specify `threads`
to indicate the number of processors on your computer that you want Trimmomatic
to use. In most cases using multiple threads(processors) can help to run the
trimming faster. These flags are not necessary, but they can give you more control
over the command. The flags are followed by **positional arguments**, meaning
the order in which you specify them is important. In paired end mode, Trimmomatic
expects the two input files, and then the names of the output files. These files are
described below. While, in single end mode, Trimmomatic will expect one file
as input, after which you can enter the optional settings and lastly the
name of the output file.

| Option    | Meaning |
| ------- | ---------- |
|  \<inputFile1>  | Input forward reads to be trimmed. Typically the file name will contain an `_1` or `_R1` in the name.|
| \<inputFile2> | Input reverse reads to be trimmed. Typically the file name will contain an `_2` or `_R2` in the name.|
|  \<outputFile1P> | Output file that contains surviving pairs from the `_1` file. |
|  \<outputFile1U> | Output file that contains orphaned reads from the `_1` file. |
|  \<outputFile2P> | Output file that contains surviving pairs from the `_2` file.|
|  \<outputFile2U> | Output file that contains orphaned reads from the `_2` file.|

The last thing Trimmomatic expects to see is the trimming parameters:

| step   | meaning |
| ------- | ---------- |
| `ILLUMINACLIP` | Perform adapter removal. |
| `SLIDINGWINDOW` | Perform sliding window trimming, cutting once the average quality within the window falls below a threshold. |
| `LEADING`  | Cut bases off the start of a read, if below a threshold quality.  |
|  `TRAILING` |  Cut bases off the end of a read, if below a threshold quality. |
| `CROP`  |  Cut the read to a specified length. |
|  `HEADCROP` |  Cut the specified number of bases from the start of the read. |
| `MINLEN`  |  Drop an entire read if it is below a specified length. |
|  `TOPHRED33` | Convert quality scores to Phred-33.  |
|  `TOPHRED64` |  Convert quality scores to Phred-64. |

We will use only a few of these options and trimming steps in our
analysis. It is important to understand the steps you are using to
clean your data. For more information about the Trimmomatic arguments
and options, see [the Trimmomatic manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf).

<a href="{{ page.root }}/fig/03-03-01.png">
  <img src="{{ page.root }}/fig/03-03-01.png" alt="Diagram showing the parts of the sequence that are reviewed by each parameter, and the parts that are mantained or discarted at the end of the reviewing process" />
</a>

However, a complete command for Trimmomatic will look something like the command below. This command is an example and will not work, as we do not have the files it refers to:

~~~
$ trimmomatic PE -threads 4 SRR_1056_1.fastq SRR_1056_2.fastq\
              SRR_1056_1.trimmed.fastq SRR_1056_1un.trimmed.fastq\
              SRR_1056_2.trimmed.fastq SRR_1056_2un.trimmed.fastq\
              ILLUMINACLIP:SRR_adapters.fa SLIDINGWINDOW:4:20
~~~
{: .bash}

In this example, we've told Trimmomatic:

| code   | meaning |
| ------- | ---------- |
| `PE` | that it will be taking a paired end file as input |
| `-threads 4` | to use four computing threads to run (this will spead up our run) |
| `SRR_1056_1.fastq` | the first input file name. Forward |
| `SRR_1056_2.fastq` | the second input file name. Reverse |
| `SRR_1056_1.trimmed.fastq` | the output file for surviving pairs from the `_1` file |
| `SRR_1056_1un.trimmed.fastq` | the output file for orphaned reads from the `_1` file |
| `SRR_1056_2.trimmed.fastq` | the output file for surviving pairs from the `_2` file |
| `SRR_1056_2un.trimmed.fastq` | the output file for orphaned reads from the `_2` file |
| `ILLUMINACLIP:SRR_adapters.fa`| to clip the Illumina adapters from the input file using the adapter sequences listed in `SRR_adapters.fa` |
|`SLIDINGWINDOW:4:20` | to use a sliding window of size 4 that will remove bases if their phred score is below 20 |



> ## Multi-line commands
> Some of the commands we ran in this lesson are long! When typing a long
> command into your terminal, you can use the `\` character
> to separate code chunks onto separate lines. This can make your code more readable.
{: .callout}



## Running Trimmomatic

Now, we will run Trimmomatic on our data. Navigate to your
`untrimmed_fastq` data directory and verify that you are
located in the untrimmed_fastq directory:

~~~
$ cd ~/dc_workshop/data/untrimmed_fastq
$ pwd
~~~
{: .bash}

~~~
$ /home/dcuser/dc_workshop/data/untrimmed_fastq
~~~
{: .output}

You should have onlye four files in this directory. Those files corresponds
to the files of forward and reverse reads from samples JC1A and JP4D.
~~~
$ ls
~~~
{: .bash}

~~~
$ JC1A_R1.fastq.gz  JC1A_R2.fastq.gz  JP4D_R1.fastq  JP4D_R2.fastq.gz   
~~~
{: .output}

We are going to run Trimmomatic on one of our paired-end samples.
While using FastQC, we saw that Universal adapters were present
in our samples. The adapter sequences came with the installation of
Trimmomatic, so we will first copy these sequences into our current
directory.

~~~
$ cp ~/.miniconda3/pkgs/trimmomatic-0.38-0/share/trimmomatic-0.38-0/adapters/TruSeq3-PE.fa .
~~~
{: .bash}


We will also use a sliding window of size 4 that will remove bases if their
phred score is below 20 (like in our example above). We will also
discard any reads that do not have at least 25 bases remaining after
this trimming step. This command will take a few minutes to run.

Before, we unzipped one of our files to work with it. Let's compress the
file corresponding to the sample `JP4D` again before we run Trimmomatic.
~~~
gzip JP4D_R1.fastq
~~~
{: .bash}

~~~
$ trimmomatic PE JP4D_R1.fastq.gz JP4D_R2.fastq.gz \
JP4D_R1.trim.fastq.gz  JP4D_R1un.trim.fastq.gz \
JP4D_R2.trim.fastq.gz  JP4D_R2un.trim.fastq.gz \
SLIDINGWINDOW:4:20 MINLEN:35 ILLUMINACLIP:TruSeq3-PE.fa:2:40:15
~~~
{: .bash}


~~~
TrimmomaticPE: Started with arguments:
 JP4D_R1.fastq.gz JP4D_R2.fastq.gz JP4D_R1.trim.fastq.gz JP4D_R1un.trim.fastq.gz JP4D_R2.trim.fastq.gz JP4D_R2un.trim.fastq.gz SLIDINGWINDOW:4:20 MINLEN:35 ILLUMINACLIP:TruSeq3-PE.fa:2:40:15
Multiple cores found: Using 2 threads
Using PrefixPair: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
ILLUMINACLIP: Using 1 prefix pairs, 0 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Quality encoding detected as phred33
Input Read Pairs: 1123987 Both Surviving: 751427 (66.85%) Forward Only Surviving: 341434 (30.38%) Reverse Only Surviving: 11303 (1.01%) Dropped: 19823 (1.76%)
TrimmomaticPE: Completed successfully
~~~
{: .output}

> ## Exercise 1: What did Trimmomatic do?
>
> Use the output from your Trimmomatic command to answer the
> following questions.
>
> 1) What percent of reads did we discard from our sample?
> 2) What percent of reads did we keep both pairs?
>
>> ## Solution
>> 1) 1.76%    
>> 2) 66.85%
> {: .solution}
{: .challenge}

You may have noticed that Trimmomatic automatically detected the
quality encoding of our sample. It is always a good idea to
double-check this or to enter the quality encoding manually.

We can confirm that we have our output files:

~~~
$ ls JP4D*
~~~
{: .bash}

~~~
JP4D_R1.fastq.gz       JP4D_R1un.trim.fastq.gz	JP4D_R2.trim.fastq.gz
JP4D_R1.trim.fastq.gz  JP4D_R2.fastq.gz		JP4D_R2un.trim.fastq.gz
~~~
{: .output}

The output files are also FASTQ files. It should be smaller than our
input file, because we've removed reads. We can confirm this with this
command:

~~~
$ ls JP4D* -l -h
~~~
{: .bash}

~~~
-rw-r--r-- 1 dcuser dcuser 179M Nov 26 12:44 JP4D_R1.fastq.gz
-rw-rw-r-- 1 dcuser dcuser 107M Mar 11 23:05 JP4D_R1.trim.fastq.gz
-rw-rw-r-- 1 dcuser dcuser  43M Mar 11 23:05 JP4D_R1un.trim.fastq.gz
-rw-r--r-- 1 dcuser dcuser 203M Nov 26 12:51 JP4D_R2.fastq.gz
-rw-rw-r-- 1 dcuser dcuser 109M Mar 11 23:05 JP4D_R2.trim.fastq.gz
-rw-rw-r-- 1 dcuser dcuser 1.3M Mar 11 23:05 JP4D_R2un.trim.fastq.gz
~~~
{: .output}


We've just successfully run Trimmomatic on one of our FASTQ files!
However, there is some bad news. Trimmomatic can only operate on
one sample at a time and we have more than one sample. The good news
is that we can use a `for` loop to iterate through our sample files
quickly!

~~~
$ for infile in *_R1.fastq.gz
do
base=$(basename ${infile} _R1.fastq.gz)
trimmomatic PE ${infile} ${base}_R2.fastq.gz \
${base}_R1.trim.fastq.gz ${base}_R1un.trim.fastq.gz \
${base}_R2.trim.fastq.gz ${base}_R2un.trim.fastq.gz \
SLIDINGWINDOW:4:20 MINLEN:35 ILLUMINACLIP:TruSeq3-PE.fa:2:40:15
done
~~~
{: .bash}


Go ahead and run the `for` loop. It should take a few minutes for
Trimmomatic to run for each of our four input files. Once it's done,
take a look at your directory contents. You'll notice that even though we ran Trimmomatic on file `JP4D` before running the for loop, there is only one set of files for it. Because we matched the ending `_R1.fastq.gz`, we re-ran Trimmomatic on this file, overwriting our first results. That's ok, but it's good to be aware that it happened.

~~~
$ ls
~~~
{: .bash}

~~~
JC1A_R1.fastq.gz                     JP4D_R1.fastq.gz                                    
JC1A_R1.trim.fastq.gz                JP4D_R1.trim.fastq.gz                                
JC1A_R1un.trim.fastq.gz              JP4D_R1un.trim.fastq.gz                                
JC1A_R2.fastq.gz                     JP4D_R2.fastq.gz                                 
JC1A_R2.trim.fastq.gz                JP4D_R2.trim.fastq.gz                                 
JC1A_R2un.trim.fastq.gz              JP4D_R2un.trim.fastq.gz                                 
TruSeq3-PE.fa   
~~~
{: .output}

> ## Exercise 2: Adapter files
> We trimmed our FASTQ files with Nextera adapters,
> but there are other adapters that are commonly used.
> How can you check the full list of adapters that came with Trimmomatic installation?
>
> 1. `ls ~/.miniconda3/pkgs/trimmomatic-0.38-0/share/trimmomatic-0.38-0/adapters/`
> 2. `cp ~/.miniconda3/pkgs/trimmomatic-0.38-0/share/trimmomatic-0.38-0/adapters/`
> 3. `head ~/.miniconda3/pkgs/trimmomatic-0.38-0/share/trimmomatic-0.38-0/adapters/`
>
>> ## Solution
>> 1. Yes, the `ls` command will display the content of the adapters folder (came with
>> Trimmomatic installation).
>> 2. No, `cp` command is incomplete and it is used to copy adapters to your working
>> directory but will not show you the adapters list.
>> 3. No, `head` command is usually used to read files content.
> {: .solution}
>> ~~~
>> $ ls ~/.miniconda3/pkgs/trimmomatic-0.38-0/share/trimmomatic-0.38-0/adapters/
>> ~~~
>> {: .bash}
>>
>> ~~~
>> NexteraPE-PE.fa  TruSeq2-SE.fa    TruSeq3-PE.fa
>> TruSeq2-PE.fa    TruSeq3-PE-2.fa  TruSeq3-SE.fa
>> ~~~
>> {: .output}
>>
> {: .solution}
{: .challenge}

We've now completed the trimming and filtering steps of our quality
control process! Before we move on, let's move our trimmed FASTQ files
to a new subdirectory within our `data/` directory.

~~~
$ cd ~/dc_workshop/data/untrimmed_fastq
$ mkdir ../trimmed_fastq
$ mv *.trim* ../trimmed_fastq
$ cd ../trimmed_fastq
$ ls
~~~
{: .bash}

~~~
JC1A_R1.trim.fastq.gz    JP4D_R1.trim.fastq.gz                
JC1A_R1un.trim.fastq.gz  JP4D_R1un.trim.fastq.gz              
JC1A_R2.trim.fastq.gz    JP4D_R2.trim.fastq.gz                
JC1A_R2un.trim.fastq.gz  JP4D_R2un.trim.fastq.gz   
~~~
{: .output}

> ## Bonus Exercise (Advanced): Quality test after trimming
>
> Now that our samples have gone through quality control, they should perform
> better on the quality tests run by FastQC.
>
> Sort the following chunks of code and decide in which terminal
> (AWS or local) you should run them to be able to re-run
> FastQC on your trimmed FASTQ files and visualize the HTML files
> to see whether your per base sequence quality is higher after
> trimming.
> ~~~
> $ scp dcuser@ec2-34-203-203-131.compute-1.amazonaws.com:~/dc_workshop/data/trimmed_fastq/*.html ~/Desktop/fastqc_html/trimmed
> ~~~
> {: .bash}
>
> ~~~
> $ fastqc ~/dc_workshop/data/trimmed_fastq/*.fastq*
> ~~~
> {: .bash}
>
> ~~~
> $ mkdir ~/Desktop/fastqc_html/trimmed
> ~~~
> {: .bash}
>
>> ## Solution
>>
>> In your AWS terminal window do:
>>
>> ~~~
>> $ fastqc ~/dc_workshop/data/trimmed_fastq/*.fastq*
>> ~~~
>> {: .bash}
>>
>> In a terminal standing on your local computer do:
>>
>> ~~~
>> $ mkdir ~/Desktop/fastqc_html/trimmed
>> $ scp dcuser@ec2-34-203-203-131.compute-1.amazonaws.com:~/dc_workshop/data/trimmed_fastq/*.html ~/Desktop/fastqc_html/trimmed
>> ~~~
>> {: .bash}
>>
>> Then take a look at the html files in your browser.
>>
>> Remember to replace everything between the `@` and `:` in your scp
>> command with your AWS instance number.
>>
>> After trimming and filtering, our overall quality is much higher,
>> we have a distribution of sequence lengths, and more samples pass
>> adapter content. However, quality trimming is not perfect, and some
>> programs are better at removing some sequences than others. Trimmomatic
>> did pretty well though, and its performance is good enough for our workflow.
> {: .solution}
{: .challenge}
