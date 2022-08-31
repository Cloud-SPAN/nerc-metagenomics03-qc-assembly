---
title: "Introduction to Metagenomics"
teaching: 30
exercises: 5
questions:
- What is metagenomics?
- When should we use metagenomics?
- What does a metagenomics project look like?
objectives:
- Explain the difference between genomics, metagenomics and amplicon sequencing.
- Familiarise yourself with the metagenomics dataset used in this course.
keypoints:
- Genomics looks at the whole genome content of an organism
- Metagenomes contain multiple organisms within one sample unlike genomic samples.
- In metagenomes the organisms present are not usually present in the same abundance - except for mock communities.
- We can identify the organisms present in a sample using either amplicon sequencing or whole metagenome sequencing. Amplicon sequencing is cheaper and quicker, however it also limits the amount of downstream analysis that can be done with the data.
- Metagenomes can diff in their levels of complexity and this is determined by how many organisms are in the metagenome.
- Difference platforms allow us to perform different analyses. The suitability depends on the question you are asking.

---


## Background

Unlike genetics, which typically considers single or small sets of genes and their influence on an organism, genomics considers the whole genetic content of an organism. 

We study genomics through sequencing the DNA content. However, there are also lots of other "'omics" fields that are related, such as transcriptomics which covers RNA content and proteomics which focuses on protein content. Other common 'omics areas of study are outlined [here](https://en.wikipedia.org/wiki/Omics).

The key difference between metagenomics and genomics is that a genomic sample contains the genome of one organism, whereas a metagenomic sample contains the genomes of multiple coexisting microorganisms.

## Metagenomics
As mentioned above, a metagenome refers to a collections of genomic sequences from various (micro)organisms coexisting in a
given space. They are like snapshots that can give us information about the taxonomic and even metabolic or functional composition of the communities that we decide to study. Examples of well studied metagenomes include the human gut, the ocean floor and soil microbiomes. 

Because we are trying to analyse multiple instead of individual organisms, this introduces additional analysis complexity. This causes two additional analysis considerations. 
1. How we can assign genomic sequences to the multiple organisms that we know are present.
2. Identifying abundances of different organisms and using this abundance to separate the data into putative organisms, which we will cover in greater depth in the [assembly section of this course](https://cloud-span.github.io/metagenomics01-qc-assembly/03-assembly/index.html). This is because in most metagenomes the organisms in the community are not present in equal abundances.

Since metagenomes are mixtures of sequences that belong to different species, a metagenomic workflow is designed to answer two questions:
1. What species are represented in the sample?
2. What are they capable of doing?

To find which species are present in a niche, we have
to do a taxonomic assignation of the obtained sequences.
To find out their capabilities, we can
look at the genes directly encoded in the metagenome or the
genes associated with the species that we found. In order to
know which methodology we should use, it is important to
know what questions do we want to answer.

## Metagenomics sequencing approaches
There are two paths to obtain information from a complex sample:
1. **Whole genome metagenomics**  
2. **Amplicon/(16S) sequencing**.

Each is named after the sequencing methodology employed
and have particular use cases, with inherent advantages and disadvantages which are summarised below.

|        | Amplicon | Whole genome metagenomics |
|-------|--------------|
| Expense | Cheap | Expensive |
| Coverage depth | High | Lower - medium* |
| Taxonomy detection | Specific to amplicons used** | All in sample |
| Genome coverage | Only region amplified | All of genome |
| Turnaround time | Fast | Slower - more computational time for analysis needed |

__* This depends on the amount of sequencing done and the uniformity in abundance of the different organisms in the sample. There are often log fold differences in the abundances of different organisms in a metagenome. Given that only a small region is used to assign taxonomy, with amplicon-only sequencing there is lower resolution for closely related taxa, and so taxonomical assignment lower than genus is often not possible.__

__**Often amplicon sequencing is referred to as 16S, however this will only amplify and sequence bacteria. If you were amplifying funghi you would need to use ITS sequences, and for protozoa 18S sequences.__


### Whole metagenome sequencing

With **Whole genome Metagenomics**, we sequence random parts of the genomes present in a sample. We can search the origin of these pieces (_i.e._ their taxonomy) and try to find to what part of the genome they correspond to.

Depending on the complexity of the metagenome, and with enough sequences, it is even possible to obtain full individual genomes from a whole genome metagenome. These can identify not only the identities of the individuals in the sample, but also their functional abilities. 

For abundant organisms in your metagenome sample, there are likely to be enough data to generate reasonable genome coverage. However, this is not the case for low abundance organisms. Often deeper sequencing/ more total sequencing data is required to assemble the genomes of these organisms. 

Depending on the question your dataset is trying to answer and how many samples you will need to sequence, the cost of both preparing the samples and the computational effort required to analyse them can become prohibitively expensive quickly. This is especially the case when you are trying to include biological or technical replication in your experimental design. 

For more information on considerations for experimental design in sequencing see our other course [Statistical Design - Understanding Experimental design](https://cloud-span.github.io/experimental_design01-principles/02-design/index.html) or [Statistical Design - Statistical analysis](https://cloud-span.github.io/experimental_design01-principles/03-statistical-analysis/index.html).

### Amplicon sequencing

In comparison **Amplicon sequencing** tends to be cheaper,
which makes it more affordable to include additional replicates.  Only one small region is amplified through PCR, rather than the whole genome. This means that this region needs to be present in all the individuals in the community that you want to identify. 

The region varies between the organisms you are trying to profile. For bacteria, the 16S rRNA sequence is used to identify bacteria, whereas the ITs region is used for fungi, and the 18S rRNA sequence for protozoa. 

This means irrespective of depth, if an organism in the community does not contain the amplified region, it will not be present. This can be useful for removing host contamination, but this also limits our information about the organisms present to their relative abundance only.

Despite this, there are workflows such as [QIIME2](https://qiime2.org/), which are free and community led, which use database annotations of the reference versions of the organisms identified from the amplicon, to suggest what metabolic functions maybe present. The amplicon sequence is also limited because species may have genomic differences, but may be indistinguishable from the amplicon sequence alone. This means that amplicon sequencing can rarely resolve to less than a genus level.

<a href="{{ page.root }}/fig/analysis_flowchart_v3.png">
  <img src="{{ page.root }}/fig/analysis_flowchart_v3.png" width="325" height="880.5" alt="Flow chart that show the steps: Experimental design, Sampling, DNA extraction, Sequencing, Read quality, Assembly, Binning, Bin quality and Data analysis"  />
</a>


## Bioinformatic workflows


<img align="right" width="325" height="506" src="{{ page.root }}/fig/short_analysis_flowchart.png" alt="Flow diagram that shows the steps: Sequence reads, Quality control, Assembly, Binning and Taxonomy" />

When working with high-throughput sequencing data, the raw reads you get off of the sequencer need to pass
through a number of  different tools in order to generate your final desired output.  

The use of this set of tools in a specified order is commonly referred to as a *workflow* or a *pipeline*.  

Here is an example of the workflow we will be using for our analysis with a brief
description of each step.  

1. Sequence reads - obtaining raw reads from a sample via sequencing
2. Quality control - assessing quality using FastQC, and trimming and/or filtering reads (if necessary)
3. Metagenome assembly - piecing together genomes from reads using Flye, a long-read metagenome assembler
4. Binning - separating out genomes into 'bins' containing related contigs using Metabat2
5. Taxonomic assignment - assigning taxonomy to sequences/contigs using Krona and/or Pavion

Workflows in bioinformatics often adopt a plug-and-play approach so the output of one tool can be easily used as input to another tool.
The use of standard data formats in bioinformatics (such as FASTA or FASTQ, which we will be using here) makes this possible.
The tools that are used to analyze data at different stages of the workflow are therefore built under the assumption that the data will be provided in a specific format.

<br clear="right"/>


## Metagenome complexity

As the number of organisms increases in a community so does the complexity. See [Pimm , 1984](https://www.nature.com/articles/307321a0) for an explanation of community complexity.

A low-complexity microbial community is made up of fewer organisms and as a result the metagenomic data is usually easier to process and analyse than a high-complexity microbial community.

There's no official definition for what makes a community low or high complexity. But some examples of a low complexity microbial community include natural whey starter cultures, made up of around three different organisms, which are used in cheese production: see [Somerville _et al._, 2019](https://bmcmicrobiol.biomedcentral.com/articles/10.1186/s12866-019-1500-0).

> ## Defining our community
> Knowing now how we can define a metagenomic community and what organisms make up our dataset, how complex would you say this metagenome is? and why?
> > ## Solution
> > As we stated above there's no hard and fast rules for defining the complexity of a community.
> > However, given this community is made up of 10 organisms we could argue that it is a low complexity metagenome when compared to, for example, the natural soil metagenome which is potentially made up of thousands of different organisms.
> {: .solution}
{: .challenge}

## Mock metagenomes versus real samples

People often use publically available sequencing data to perform analyses or for testing an analysis pipeline on a smaller or similar dataset. 

However mock metagenomes are communities that are made by combining known species at known quantities to generate a mock sample. This mock community is often much simpler than a real metagenome, however it can be useful to know exactly what is in a metagenome, in order to test whether the result we get is sensible.

One reason for using these mock samples would be to test a workflow. The simplicity and quality of mock metagenome data means that performing analysis on these mock samples can be much easier computationally, and run more quickly.

However the main reason for using mock samples are to benchmark developing software tools. If you are aware of the level of complexity in a sample, and you know the content of the community, it is easier for you to identify false positives and false negatives. 

Some mock datasets, like the one we will be using for this course, have the organisms present with log fold differences in the abundance between the organisms. This is much more representative of a real community compared to a mock community, where all organisms are present at the same abundance. However it is useful to know that this is still an "idealised" dataset. In most cases the samples will not contain all known organisms, and often there will be degraded pieces of DNA in there.  

## Dataset used in this course

This course uses data from a mock metagenome community published from [Ultra-deep, long-read nanopore sequencing of mock microbial community standards](https://academic.oup.com/gigascience/article/8/5/giz043/5486468) which has long and short read sequencing data and has been used for benchmarking metagenome tools. 

This community contains eight bacteria (three gram positive and five gram negative) and two yeasts. These organisms are present in log abundances relative to each other. The organisms present and their known percentage of the total metagenome content is given in the table below.

This dataset contains whole metagenome sequencing. The short read data are generated on the Illumina platorm and the long reads are generated using Oxford Nanopore Technology's Nanopore platform. Other popular long read sequencing platforms exist, such as PacBio, though we will not be covering PacBio specific methodology here. Despite this, the same principle stages exist in the workflow, and often only different parameters are required to adapt analysis to different platforms.



|Species | Genomic DNA % | Genome coverage |
|-------|-----------|
| Listeria monocytogenes | 89.1 | 94.8 |
| Pseudomonas aeruginosa | 8.9 | 4.2 |
| Bacillus subtilis | 0.89 | 0.7 |
| Saccharomyces cerevisiae | 0.89 | 0.23 |
| Escherichia coli | 0.089 | 0.058 |
| Salmonella enterica | 0.089 | 0.059 |
| Lactobacillus fermentum | 0.0089 | 0.015 |
| Enterococcus faecalis | 0.00089 | 0.001 |
| Cryptococcus neoformans | 0.00089 | 0.00015 |
| Staphylococcus aureus | 0.000089 | 0.0001 |


## Differences between nanopore and illumina data

Elsewhere [(see here)](https://cloud-span.github.io/experimental_design01-principles/01-platform/index.html)we have covered information about how to use your research question to guide your preference of platform (Illumina/Nanopore)when designing an experiment. However typically the equivalent question for metagenomics is whether to sequence the whole genome, or whether amplicon sequencing will be preferable.

For non-metagenomic analyses, you can choose to do either a reference based or *de novo* approach. This will be dependent on whether there is a reasonable reference for your organism. 

However, for whole metagenome sequencing, reference genomes will exist that can be compared to organisms identified from your metagenome. However, this will be at the binning stage and a reference will not exist for your metagenome as a whole. Due to this, all of the assembly stages of the metagenome analysis pipeline are *de novo*. As a result, there is a bigger advantage to using long read sequencing over short read sequencing to assemble a metagenome if you were to choose only one method.

This will be covered in the [Genome Assembly](https://cloud-span.github.io/metagenomics01-qc-assembly/03-assembly/index.html) section of this course. There are pros and cons to each using both long and short reads, and so using them in combination is usually the best method. These pros and cons are irrespective of the application.

However, for metagenome analysis, if you were to use only short read sequencing for the assembly you would end up with a much more fragmented assembly to work with.

|        | Short reads | Long reads |
|-------|-----------|
| Technologies | Illumina | Nanopore and pacbio |
| Number of reads generated | 800 million paired end* | Depends on read length, and instrument, but usually only 10s of thousands** |
| Base accuracy | Very high | Low |
| Ease of assembly | Very difficult | Easier |
| Format output files | Fastq | Fastq, Fast5 |
| Read length | 150-300bp | Commonly 10-30kb*** |


__* As of July 2022, the NextSeq 550 high-output system runs were capable of generating upto [800 million paired-end reads](https://emea.illumina.com/systems/sequencing-platforms/nextseq/specifications.html) in one run__

__** There are different Nanopore instruments. The smaller instruments, like the minION, will generate far fewer reads. Larger instruments like the promethION will result in ~10-20k reads, but this will vary a lot between samples and their quality. They will never result in millions of reads like the Illumina platforms.__

__*** The read length will vary based on the DNA extraction method. If the DNA is very fragmented it will not result in very long reads. In many metagenomes bead beating is required to lyse cells, and so read length will still be longer than Illumina but shorter than non-metagenomic samples sequenced.__
