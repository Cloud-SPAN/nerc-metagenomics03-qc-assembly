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
-
---


## Background
## Metagenomics
Metagenomes are collections of genomic
sequences from various (micro)organisms that coexist in any
given space. They are like snapshots that can give us information
about the taxonomic and even metabolic or functional composition
of the communities that we decide to study. Thus, metagenomes
are usually employed to investigate the ecology of defining
characteristics of niches(*e.g.*, the human gut, or the ocean floor).

Since metagenomes are mixtures of sequences that belong to different species,
a metagenomic workflow is designed to answer two questions:
1. What species are represented in the sample?
2. What are they capable of doing?

To find which species are present in a niche, we have
to do a taxonomic assignation of the obtained sequences.
To find out their capabilities, we can
look at the genes directly encoded in the metagenome or the
genes associated with the species that we found. In order to
know which methodology we should use, it is important to
know what questions do we want to answer.

## Whole genome and amplicon
There are two paths to obtain information from a complex sample:
1. **Whole genome metagenomics**  
2. **Amplicon/(16S) sequencing**.

Each is named after the sequencing methodology employed
and have particular use cases, with inherent advantages and disadvantages.

|        | Amplicon | Whole genome metagenomics |
|-------|--------|-----------------|---------|
| Expense | Cheap | Expensive |
| Coverage depth | High | Lower - medium* |
| Taxonomy detection | Specific to amplicons used** | All in sample |
| Genome coverage | Only region amplified | All of genome |
| Turnaround time | Fast | Slower - more computational time for analysis needed |

*This depends on the amount of sequencing done and the uniformity in abundance of the different organisms in the sample. There are often log fold differences in the abundances of different organisms in a metagenome. Given that only a small region is used to assign taxonomy with amplicon sequencing only you have less resolution for closely related taxa, and so assignment to lower than genus is often not possible.

**Often amplicon sequencing is referred to as 16S, however this will amplify and sequence bacteria. If you were amplifying funghi you would need to use ITS sequences, and for protozoa 18S sequences.



With **Whole genome Metagenomics**, we sequence random parts (ideally all of them) of the
genomes present in a sample. We can search the origin of these
pieces (_i.e.,_ their taxonomy) and also try to find to what
part of the genome they correspond. Given enough pieces, it is even possible
to obtain full individual genomes from a Whole genome metagenome,
which could give us a bunch of information about the species
in our study. This, however, requires that we have a lot of genomic
sequences from one organism, and since the sequencing is done at random,
we usually have to sequence our community a lot (have a high sequencing depth)
to make sure that we obtain enough pieces of a given genome. This gets
exponentially harder when our species of interest is not very abundant.
It also requires that we have enough DNA to work with, which can be
difficult to obtain in certain cases. Finally, a lot of sequencing
means a lot of expenses, and because of this, making technical
and biological replicates can be prohibitively costly.   

On the contrary, **Amplicon sequencing** tends to be cheaper,
which makes it easier to duplicate and even triplicate
them without taking a big financial hit. This is because
amplicon sequencing is the collection of small genomic fragments
present in the community and amplified through PCR. If the
amplified region is present only once in every genome, ideally,
we wouldn't need to sequence the amplicon metagenome so thoroughly
because one sequence is all we need to get the information
about that genome, and by extension, about that species. On the other
hand, if a genome in the community lacks the region targeted by the
PCR primers, then no amount of sequencing can give us information
about that genome. This is why the most popular amplicon used for
this methodology are 16S amplicons for Bacteria since every known
bacteria have this particular region. Other regions can be chosen,
but they are used for very specific cases. However, even 16S amplicons
are limited to, well, the 16S region, so amplicon metagenomes cannot
directly tell us a lot about the metabolic functions found in each genome,
although educated guesses can be made by knowing which genes are
commonly found in every identified species. The information 16S gives us
is also limited because we there is a lot variation in the 16S variable
region compared to a genome, so we may be unable to identify several species
that belong to the same genus for instance.

<a href="{{ page.root }}/fig/analysis_flowchart_v3.png">
  <img src="{{ page.root }}/fig/analysis_flowchart_v3.png" alt="Flow chart that show the steps: Experimental design, Sampling, DNA extraction, Sequencing, Read quality, Assembly, Binning, Bin quality and Data analysis " />
</a>


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


## Metagenome complexity

As the number of organisms increases in a community so does the complexity. See [Pimm , 1984](https://www.nature.com/articles/307321a0) for an explanation of community complexity.
A low-complexity microbial community is made up of fewer organisms and as a result the metagenomic data is usually easier to process and analyse than a high-complexity microbial community.
There's no official definition for what makes a community low or highly complex. But some examples of a low complexity microbial community include natural whey starter cultures, made up of around three different organisms, which are used in cheese production see [Somerville _et al._, 2019](https://bmcmicrobiol.biomedcentral.com/articles/10.1186/s12866-019-1500-0).

> ## Defining our community
> Knowing now how we can define a metagenomic community and what organisms make up our dataset, how complex would you say this metagenome is? and why?
> > ## Solution
> > As we stated above there's no hard and fast rules for defining the complexity of a community.
> > However, given this community is made up of 10 organisms we could argue that it is a low complexity metagenome when compared to for example natural soil metagenome which are potentially made up of thousands of different organisms.
> {: .solution}
{: .challenge}

## Mock metagenomes aren't necessarily a good proxy for real world data

Something about this giving us the option to test the workflow in a short amount of time without needing too much compute
But this isn't representative of how long this workflow would take for real data.
