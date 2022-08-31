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

When we QC'ed the nanopore reads by quality score, we used [Seqkit](https://bioinf.shenwei.me/seqkit/) and the command `seqkit seq`. This time we will be using the command `seqkit stats` instead.
