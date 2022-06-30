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

## Shotgun and amplicon
There are two paths to obtain information from a complex sample:
1. **Shotgun Metagenomics**  
2. **Amplicon/(16S) sequencing**.

Each is named after the sequencing methodology employed
and have particular use cases, with inherent advantages and disadvantages.

With **Shotgun Metagenomics**, we sequence random parts (ideally all of them) of the
genomes present in a sample. We can search the origin of these
pieces (_i.e.,_ their taxonomy) and also try to find to what
part of the genome they correspond. Given enough pieces, it is even possible
to obtain full individual genomes from a shotgun metagenome,
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

<a href="{{ page.root }}/fig/03-01-01.png">
  <img src="{{ page.root }}/fig/03-01-01.png" alt="Flow chart that show the steps: Experimental design, Sampling, DNA extraction, Sequencing, Read quality, Assembly, Binning, Bin quality and Data analysis " />
</a>
