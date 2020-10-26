---
layout: lesson
root: .
---

A lot of genomics analysis is done using command-line tools for three reasons:   
1) you will often be working with a large number of files, and working through the command-line rather than 
through a graphical user interface (GUI) allows you to automate repetitive tasks,    
2) you will often need more compute power than is available on your personal computer, and 
connecting to and interacting with remote computers requires a command-line interface, and    
3) you will often need to customize your analyses, and command-line tools often enable more 
customization than the corresponding GUI tools (if in fact a GUI tool even exists).   

In a [previous lesson](https://kccg.github.io/shell-genomics/), you learned how to use the bash shell to interact with your computer through a command line interface. In this 
lesson, you will be applying this new knowledge to carry out a common genomics workflow - identifying variants among sequencing samples 
taken from multiple individuals within a population. We will be starting with a set of sequenced reads (`.fastq` files), performing
some quality control steps, aligning those reads to a reference genome, and ending by identifying and visualizing variations among these
samples. 

As you progress through this lesson, keep in mind that, even if you aren't going to be doing this same workflow in your research, 
you will be learning some very important lessons about using command-line bioinformatic tools. What you learn here will enable you to 
use a variety of bioinformatic tools with confidence and greatly enhance your research efficiency and productivity.

> ## Prerequisites
>
> This lesson assumes a working understanding of the bash shell. If you haven't already completed the [Shell Genomics](https://kccg.github.io/shell-genomics/) lesson, and aren't familiar with the bash shell, please review those materials
> before starting this lesson.
>
> This lesson also assumes some familiarity with biological concepts, including the structure of DNA, nucleotide abbreviations, and the 
> concept of genomic variation within a population. 
>
> This lesson is to be completed on the Garvan high performance cluster (Wolfpack). 
See the [previous lesson](https://kccg.github.io/shell-genomics/) for information on how to access the cluster, as well as information on how to access the practice data.
{: .prereq}
