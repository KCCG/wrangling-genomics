---
title: "Variant Calling Workflow"
teaching: 35
exercises: 25
questions:
- "How do I find sequence variants between my sample and a reference genome?"
objectives:
- "Understand the steps involved in variant calling."
- "Describe the types of data formats encountered during variant calling."
- "Use command line tools to perform variant calling."
keypoints:
- "Bioinformatic command line tools are collections of commands that can be used to carry out bioinformatic analyses."
- "To use most powerful bioinformatic tools, you'll need to use the command line."
- "There are many different file formats for storing genomics data. It's important to understand what type of information is contained in each file, and how it was derived."
---

We mentioned before that we are working with files from a long-term evolution study of an *E. coli* population (designated Ara-3). Now that we have looked at our data to make sure that it is high quality, and removed low-quality base calls, we can perform variant calling to see how the population changed over time. We care how this population changed relative to the original population, *E. coli* strain REL606. Therefore, we will align each of our samples to the *E. coli* REL606 reference genome, and see what differences exist in our reads versus the genome.

# Alignment to a reference genome

![workflow_align](../img/variant_calling_workflow_align.png)

We perform read alignment or mapping to determine where in the genome our reads originated from. There are a number of tools to
choose from and, while there is no gold standard, there are some tools that are better suited for particular NGS analyses. We will be
using the [Burrows Wheeler Aligner (BWA)](http://bio-bwa.sourceforge.net/), which is a software package for mapping low-divergent
sequences against a large reference genome. 

The alignment process consists of two steps:

1. Indexing the reference genome
2. Aligning the reads to the reference genome


# Setting up

First we download the reference genome for *E. coli* REL606. Although we could copy or move the file with `cp` or `mv`, most genomics workflows begin with a download step. 
In the interests of time, the file has already been downloaded and placed in one of your directories. 
The code below shows how this **would** have been acheived. 
Study this code and try to figure out where the reference genome data would be stored.
Also, quickly check the `man` page for the `curl` command to find out what the "-L" and "-o" options do (use <kbd>/</kbd> to search).

~~~
$ cd ~/course
$ mkdir -p data/ref_genome
$ curl -L -o data/ref_genome/ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
$ gunzip data/ref_genome/ecoli_rel606.fasta.gz
~~~
{: .bash}

Actually, a bacterial genome is not **that** big, at least compared to mamallian genomes, so downloading it would not have been **too** time consuming.
But don't download anything that you don't have to -- it's expensive and time-consuming.
In fact, there is a lot of big biological data already on /share/ClusterShare.
Navigate to /share/ClusterShare and poke around until you find it.
You'll notice that a lot of it is organised into directories based on user names.
That way you know who to ask if you have questions about a particular dataset.
If you do need to download a lot of data, especially reference data, and you think it might be useful to others, make a directory and store the data here.

I've placed the E.Coli reference genome data in my directory in biodata/contrib.
You'll need to refer to this data as part of your analysis.
But you won't need to actually modify this data.
So instead of copying the data to your project, let's create a `link`.
A `link` is like an "alias" on a Mac or a "shortcut" on Windows.
(If that doesn't mean anything to you then don't worry, it will make more sense after you see it.)

~~~
$ ln --symbolic /share/ClusterShare/biodata/contrib/johree/ecoli_ref_genome ~/course/data/ref_genome
$ cd ~/course/data
$ ls --classify
$ ls -l
$ cd ~/course/data/ref_genome
~~~
{: .bash}

In the first directory listing ("ls --classify") you should see an "@" symbol after "ref_genome".
This indicates that "ref_genome" is a link rather than a "real" directory.
In the second long form directory listing ("ls -l") you should see an arrow pointing from "ref_genome" to location where the real directory actually lives.

Note that we have used the "--symbolic" option with the `ln` command ("-s" for short).
This creates a "soft" link. 
If you delete the link then the original file or directory will remain untouched, whereas with a "hard" link (without the "--symbolic" option) both the source and the destination will be deleted.
Hopefully it is fairly obvious why that is generally not a good idea.

> ## Exercise 
> 
> We saved this file as `data/ref_genome/ecoli_rel606.fasta.gz` and then decompressed it. 
> What is the real name of the genome? 
> 
>> ## Solution
>> 
>> ~~~
>> $ head data/ref_genome/ecoli_rel606.fasta
>> ~~~
>> {: .bash}
>> 
>> The name of the sequence follows the `>` character. The name is `CP000819.1 Escherichia coli B str. REL606, complete genome`.
>> Keep this chromosome name (`CP000819.1`) in mind, as we will use it later in the lesson. 
> {: .solution}
{: .challenge}

We will also download a set of trimmed FASTQ files to work with. 
These are small subsets of our real trimmed data, and will enable us to run our variant calling workflow quite quickly. 
Again, I've saved you the trouble of actually downloading the files, but do take a second to see how it would be done.
Find the trimmed_fastq_small directory within your directory structure.

~~~
$ curl -L -o sub.tar.gz https://ndownloader.figshare.com/files/14418248
$ tar xvf sub.tar.gz
$ mv sub/ ~/course/data/trimmed_fastq_small
~~~
{: .bash}

> ## Tip
> In the example above, the `tar` command uses what is called "option stacking".
> At least, that's what I call it anyway...
> The options here are "-x", "-v" and "-f". 
> The `tar` command allows you to type "tar xvf" as a shorthand for "tar -x -v -f".
> Check out the `man` page for `tar` to see what these options do.
> As well as using <kbd>/</kbd> to search for "-x" etc, I also recommend using <kbd>g</kbd> in between searches to go back to the top of the page.
{: .callout}

You will also need to create directories for the results that will be generated as part of this workflow. We can do this in a single
line of code, because `mkdir` can accept multiple new directory
names as input.

~~~
$ mkdir -p results/sam results/bam results/bcf results/vcf
~~~
{: .bash}

Do you remember what the "-p" option does?
Look up the `man` page if you are not sure.
Use /-p to search the manual.


### Index the reference genome
Our first step is to index the reference genome for use by BWA. 
Indexing allows the aligner to quickly find potential alignment sites for query sequences in a genome, which saves time during alignment. 
Indexing the reference only has to be run once. 
The only reason you would want to create a new index is if you are working with a different reference genome or you are using a different tool for alignment.

Before we can use bwa we'll need to find and load the appropriate module.

> ## Exercise 
> Which bwa modules are available on the cluster?
> 
>> ## Solution
>> 
>> ~~~
>> $ modgrep bwa
>> ~~~
>> {: .bash}
>> aarsta/bwa/0.7.9a
>> aarsta/bwa-meth/0.09
>> aarsta/bwa-meth/0.10
>> aarsta/bwa-meth/git
>> aarsta/bwa-meth/git_aaron
>> evaben/bwa/gcc-7.3.0/0.7.15
>> gi/bwa/0.5.8c
>> gi/bwa/0.7.10
>> gi/bwa/0.7.12
>> gi/bwa/0.7.4
>> gi/bwa/0.7.5a
>> gi/bwa/0.7.6a
>> gi/bwa/0.7.8
>> gi/bwa/0.7.9a
>> kevyin/bwa/0.6.2
>> marcow/bwa/0.7.3a/gcc-4.4.6
>> marcow/bwa/gcc-4.4.6/0.7.3a
>> marsmi/bwa/0.7.17
>> pethum/bwa/gcc-4.4.6/0.7.12
>> pethum/bwakit/prebuilt/0.7.12
>> phuluu/bwa-meth/0.10
>> timpet/bwa-meth/git_aaron_mem_patch
>> ~~~
>> 
> {: .solution}
{: .challenge}

> ## Exercise 
> Which bwa module should we use?
>
> Check the [software requirements](https://datacarpentry.org/genomics-workshop/setup.html) for this workshop.
> 
>> ## Solution
>> 
>> marsmi/bwa/0.7.17
> {: .solution}
{: .challenge}

The command for indexing the reference genome, shown below, is quite simple.
But it is a good idea to consult reference manuals even when following a recipe.
That way you can understand what is happening, or notice other options that might be useful for what you are trying to do.
There might also be warnings, cautions or other "gotchas" that you need to avoid.

Have a quick look at the reference manual for `bwa`: http://bio-bwa.sourceforge.net/bwa.shtml
In particular, scan the section for "bwa index" in the [Commands and Options](http://bio-bwa.sourceforge.net/bwa.shtml#3) section.
What you will (hopefully) notice is that `bwa` has mutiple `subcommands` such as `index`, `mem` and `aln`.
Each one of these has its own set of options.
This kind of pattern is very common for bioinformatics tools.

On a small bacterial genome, indexing is a relatively quick process.
There is no problem just running the command interactively on a compute node requested with `qrsh`.
But we also want to practice making job scripts, and later we are going to wrap up our **entire** workflow into a job script.  

This is a good chance to clarify how `qrsh` and `qsub` work together. 
Request an interactive node with `qsub` before you run the rest of the commands in this chapter.
(No need to request extra memory.) 
But at the same time, create a new job script in `Atom` or `Sublime` (based on the job script template) and gradually keep track of the different steps as you go.
Your script will evolve as you work through this chapter.
You will load more modules, and run more commands.
The outputs from one command will be the input for another command.
You will make certain decisions about work directories and variables that you might end up reviewing.
It's OK to change your mind, but by having your script open in an editor it is easy to capture each additional step.
Once the script is "mature" you can run it using `qsub` to verify that you have correctly encapsulated the entire workflow.

Some things to think about as you work on your script.

> 1. Start with the job script template
> 2. Change the name of the job
> 3. Set the working directory 
> 4. Load all the necessary modules
> 5. Work out which parts to parameterize, so that your script can be easily adopted to other files, other samples, or other projects.

OK, back to working on the data...

~~~
$ qrsh
$ module load marsmi/bwa/0.7.17 
$ bwa index data/ref_genome/ecoli_rel606.fasta
$ exit
~~~
{: .bash}

The output should look something like this:
~~~
[bwa_index] Pack FASTA... 0.04 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 1.05 seconds elapse.
[bwa_index] Update BWT... 0.03 sec
[bwa_index] Pack forward-only FASTA... 0.02 sec
[bwa_index] Construct SA from BWT and Occ... 0.57 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index data/ref_genome/ecoli_rel606.fasta
[main] Real time: 1.765 sec; CPU: 1.715 sec
~~~
{: .output}

Transfer the relavent parts of the command above to the job script in your editor.
Include the module loading step. 
Parameterize some parts as appropriate.
For example, you could define a variable called REF_GENOME and then refer to this variable as part of the bwa command step.


### Align reads to reference genome

We will use the BWA-MEM algorithm, which is suited well to aligning accurate short-read transcriptomic Illumina data to genomic sequences. 
Alternatively, aligners such as minimap2 are well-suited for aligning noisy long-read data or short-read genomic Illumina data. 
The appropriate choice of aligner depending on the sequencing read types is crucial for down-stream high-quality genomic data analysis and some time should be spent choosing the best tool for the job.

An example of what a `bwa` command looks like is below. 
This command will not run, as we do not have the files `ref_genome.fa`, `input_file_R1.fastq`, or `input_file_R2.fastq`.

~~~
$ bwa mem ref_genome.fasta input_file_R1.fastq input_file_R2.fastq > output.sam
~~~
{: .bash}

Have a look at the [bwa options page](http://bio-bwa.sourceforge.net/bwa.shtml). 
While we are running bwa with the default parameters here, your use case might require a change of parameters.
*NOTE: Always read the manual page for any tool before using and make sure the options you use are appropriate for your data.*

We're going to start by aligning the reads from just one of the samples in our dataset (`SRR2584866`). 
Later, we'll be iterating this whole process on all of our sample files.

This is the commmand that we want to run.

~~~
$ qrsh
$ module load marsmi/bwa/0.7.17 
$ bwa mem data/ref_genome/ecoli_rel606.fasta \  
          data/trimmed_fastq_small/SRR2584866_1.trim.sub.fastq \  
          data/trimmed_fastq_small/SRR2584866_2.trim.sub.fastq \ 
          > results/sam/SRR2584866.aligned.sam
$ exit
~~~
{: .bash}

The output starts like this: 

~~~
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 77446 sequences (10000033 bp)...
[M::process] read 77296 sequences (10000182 bp)...
[M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR): (48, 36728, 21, 61)
[M::mem_pestat] analyzing insert size distribution for orientation FF...
[M::mem_pestat] (25, 50, 75) percentile: (420, 660, 1774)
[M::mem_pestat] low and high boundaries for computing mean and std.dev: (1, 4482)
[M::mem_pestat] mean and std.dev: (784.68, 700.87)
[M::mem_pestat] low and high boundaries for proper pairs: (1, 5836)
[M::mem_pestat] analyzing insert size distribution for orientation FR...
~~~
{: .output}

Once again, modify your job script to incude this command.
Parametise all of the input and output files by assigning them to variables.
Think carefully about what to do with the work directory.
The paths in the command above are `relative paths`. 
What directory are they relative to?
Rewrite the command to use variables rather than hard-coded file names.
This will help when we generalise our script to iterate over all files in a moment.


#### SAM/BAM format
The [SAM file](https://genome.sph.umich.edu/wiki/SAM),
is a tab-delimited text file that contains information for each individual read and its alignment to the genome. 
While we do not have time to go into detail about the features of the SAM format, the paper by 
[Heng Li et al.](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full) provides a lot more detail on the specification.

**The compressed binary version of SAM is called a BAM file.** 
We use this version to reduce size and to allow for *indexing*, which enables efficient random access of the data contained within the file.

The file begins with a **header**, which is optional. 
The header is used to describe the source of data, reference sequence, method of alignment, etc., this will change depending on the aligner being used. 
Following the header is the **alignment section**. 
Each line that follows corresponds to alignment information for a single read. 
Each alignment line has **11 mandatory fields** for essential mapping information and a variable number of other fields for aligner specific information. 
An example entry from a SAM file is displayed below with the different fields highlighted.

![sam_bam1](../img/sam_bam.png)


![sam_bam2](../img/sam_bam3.png)

> ## Exercise 
> Which samtools modules are available on the cluster?
> 
>> ## Solution
>> 
>> ~~~
>> $ modgrep samtools
>> aarsta/samtools/0.1.19
>> aledre/samtools/prebuilt/1.10
>> annsen/samtools/gcc-7.3.0/1.9
>> briglo/samtools/1.5
>> briglo/samtools/1.9
>> evaben/samtools/gcc-7.3.0/1.8
>> evaben7/samtools/1.9/gcc-8.2.0
>> gi/samtools/0.1.19
>> gi/samtools/1.0
>> gi/samtools/1.1
>> gi/samtools/1.2
>> marsmi/samtools/1.6
>> marsmi/samtools/1.7
>> vpethum/samtools/gcc-4.4.6/1.2
>> pethum/samtools/gcc-4.4.6/1.3
>> phuluu/samtools/1.4
>> vxiuque/samtools
>> ~~~
>> {: .bash}
> {: .solution}
{: .challenge}

> ## Exercise 
> Which samtools module should we use?
>
> Check the [software requirements](https://datacarpentry.org/genomics-workshop/setup.html) for this workshop.
> 
>> ## Solution
>> 
>> There are two modules for samtools 1.9 on the cluster.
>> Previous groups have successfully completed this task using briglo/samtools/1.9
> {: .solution}
{: .challenge}

We will convert the SAM file to BAM format using the `samtools` program with the `view` command and tell this command that the input is in SAM format (`-S`) and to output BAM format (`-b`).  
Run the command interactively.

~~~
$ qrsh
$ module load briglo/samtools/1.9
$ samtools view -S -b results/sam/SRR2584866.aligned.sam > results/bam/SRR2584866.aligned.bam
$ exit
~~~
{: .bash}

~~~
[samopen] SAM header is present: 1 sequences.
~~~
{: .output}

Edit your evolving job script to include this new step.
Load the extra module.
Note: you can do this in one step by tacking on another module at the end of the "module load" command.
You have a design decision when it comes to parameterisation.
You can either create two parameter variables (one for the input file and another for the output file) or you can create one parameter variable for the sample id.
Which do you think is a better option?
Is there a happy middle ground?


### Sort BAM file by coordinates

Next we sort the BAM file using the `sort` subcommand from `samtools`. 
`-o` tells the command where to write the output.

~~~
$ qrsh
$ module load briglo/samtools/1.9
$ samtools sort -o results/bam/SRR2584866.aligned.sorted.bam results/bam/SRR2584866.aligned.bam 
$ exit
~~~
{: .bash}

Our files are pretty small, so we won't see this output. If you run the workflow with larger files, you will see something like this:
~~~
[bam_sort_core] merging from 2 files...
~~~
{: .output}

Once again, add this new step to your ever evolving job script.

SAM/BAM files can be sorted in multiple ways, e.g. by location of alignment on the chromosome, by read name, etc. It is important to be aware that different alignment tools will output differently sorted SAM/BAM, and different downstream tools require differently sorted alignment files as input.

You can use samtools to learn more about this bam file as well.

~~~
$ qrsh
$ module load briglo/samtools/1.9
$ samtools flagstat results/bam/SRR2584866.aligned.sorted.bam
$ exit
~~~
{: .bash}

This will give you the following statistics about your sorted bam file:

~~~
351169 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
1169 + 0 supplementary
0 + 0 duplicates
351103 + 0 mapped (99.98% : N/A)
350000 + 0 paired in sequencing
175000 + 0 read1
175000 + 0 read2
346688 + 0 properly paired (99.05% : N/A)
349876 + 0 with itself and mate mapped
58 + 0 singletons (0.02% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
~~~
{: .output}

You could include this step in your job script if you want.
If you do then output will be saved in the job logs.
But it probably makes more sense to run this kind of command interactively, so you can leave it out of your job script if you want.

## Variant calling

A variant call is a conclusion that there is a nucleotide difference vs. some reference at a given position in an individual genome
or transcriptome, often referred to as a Single Nucleotide Polymorphism (SNP). The call is usually accompanied by an estimate of 
variant frequency and some measure of confidence. Similar to other steps in this workflow, there are a number of tools available for 
variant calling. In this workshop we will be using `bcftools`.

> ## Exercise 
> Which bcftools modules are available on the cluster?
> 
>> ## Solution
>> 
>> ~~~
>> modgrep bcftools
>> aarsta/bcftools/1.2
>> aarsta/bcftools/1.6
>> aledre/bcftools/prebuilt/1.10
>> annsen/bcftools/gcc-7.3.0/1.9
>> briglo/bcftools/1.9
>> evaben/bcftools/gcc-7.3.0/1.8
>> evaben7/bcftools/1.9/gcc-8.2.0
>> gi/bcftools/1.0
>> julyin/bcftools/1.3.1
>> pethum/bcftools/gcc-4.4.6/1.3
>> ~~~
>>
> {: .solution}
{: .challenge}

> ## Exercise 
> Which bcftools module should we use?
>
> Check the [software requirements](https://datacarpentry.org/genomics-workshop/setup.html) for this workshop.
> 
>> ## Solution
>> 
>> This time there are **three** modules for bcftools 1.9 on the cluster.
>> Can you spot them all?
>> We want evaben7/bcftools/1.9/gcc-8.2.0
> {: .solution}
{: .challenge}

There are a few things we need to do before actually calling the variants.

![workflow](../img/variant_calling_workflow.png)

### Step 1: Calculate the read coverage of positions in the genome

Do the first pass on variant calling by counting read coverage with 
[bcftools](https://samtools.github.io/bcftools/bcftools.html). 
We will use the command `mpileup`. 
The flag `-O b` tells bcftools to generate a bcf format output file, `-o` specifies where to write the output file, and `-f` flags the path to the reference genome:

~~~
$ qrsh
$ module load evaben7/bcftools/1.9/gcc-8.2.0
$ bcftools mpileup -O b -o results/bcf/SRR2584866_raw.bcf \
-f data/ref_genome/ecoli_rel606.fasta results/bam/SRR2584866.aligned.sorted.bam 
$ exit
~~~
{: .bash}

~~~
[mpileup] 1 samples in 1 input files
~~~
{: .output}

Once again, add this step to your growing `job script`.

We have now generated a file with coverage information for every base.

### Step 2: Detect the single nucleotide polymorphisms (SNPs)

Identify SNPs using bcftools `call`. 
We have to specify ploidy with the flag `--ploidy`, which is one for the haploid *E. coli*. 
`-m` allows for multiallelic and rare-variant calling, `-v` tells the program to output variant sites only (not every site in the genome), and `-o` specifies where to write the output file:

~~~
$ qrsh
$ module load evaben7/bcftools/1.9/gcc-8.2.0
$ bcftools call --ploidy 1 -m -v -o results/bcf/SRR2584866_variants.vcf results/bcf/SRR2584866_raw.bcf 
$ exit
~~~
{: .bash}

After running this command interactively, add it to your job script.

### Step 3: Filter and report the SNP variants in variant calling format (VCF)

Filter the SNPs for the final output in VCF format, using `vcfutils.pl`:

~~~
$ qrsh
$ module load evaben7/bcftools/1.9/gcc-8.2.0
$ vcfutils.pl varFilter results/bcf/SRR2584866_variants.vcf  > results/vcf/SRR2584866_final_variants.vcf
$ exit
~~~
{: .bash}

You know the drill by now.
Add this to your job script.


## Explore the VCF format:

This next step is exploratory.
You don't have to add it to your job script.

~~~
$ less -S results/vcf/SRR2584866_final_variants.vcf
~~~
{: .bash}

You will see the header (which describes the format), the time and date the file was
created, the version of bcftools that was used, the command line parameters used, and 
some additional information:

~~~
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##bcftoolsVersion=1.8+htslib-1.8
##bcftoolsCommand=mpileup -O b -o results/bcf/SRR2584866_raw.bcf -f data/ref_genome/ecoli_rel606.fasta results/bam/SRR2584866.aligned.sorted.bam
##reference=file://data/ref_genome/ecoli_rel606.fasta
##contig=<ID=CP000819.1,length=4629812>
##ALT=<ID=*,Description="Represents allele(s) other than observed.">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of reads supporting an indel">
##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of reads supporting an indel">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version=
##INFO=<ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">
##INFO=<ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">
##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
##INFO=<ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)">
##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=ICB,Number=1,Type=Float,Description="Inbreeding Coefficient Binomial test (bigger is better)">
##INFO=<ID=HOB,Number=1,Type=Float,Description="Bias in the number of HOMs number (smaller is better)">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">
##bcftools_callVersion=1.8+htslib-1.8
##bcftools_callCommand=call --ploidy 1 -m -v -o results/bcf/SRR2584866_variants.vcf results/bcf/SRR2584866_raw.bcf; Date=Tue Oct  9 18:48:10 2018
~~~
{: .output}

Followed by information on each of the variations observed: 

~~~
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  results/bam/SRR2584866.aligned.sorted.bam
CP000819.1      1521    .       C       T       207     .       DP=9;VDB=0.993024;SGB=-0.662043;MQSB=0.974597;MQ0F=0;AC=1;AN=1;DP4=0,0,4,5;MQ=60
CP000819.1      1612    .       A       G       225     .       DP=13;VDB=0.52194;SGB=-0.676189;MQSB=0.950952;MQ0F=0;AC=1;AN=1;DP4=0,0,6,5;MQ=60
CP000819.1      9092    .       A       G       225     .       DP=14;VDB=0.717543;SGB=-0.670168;MQSB=0.916482;MQ0F=0;AC=1;AN=1;DP4=0,0,7,3;MQ=60
CP000819.1      9972    .       T       G       214     .       DP=10;VDB=0.022095;SGB=-0.670168;MQSB=1;MQ0F=0;AC=1;AN=1;DP4=0,0,2,8;MQ=60      GT:PL
CP000819.1      10563   .       G       A       225     .       DP=11;VDB=0.958658;SGB=-0.670168;MQSB=0.952347;MQ0F=0;AC=1;AN=1;DP4=0,0,5,5;MQ=60
CP000819.1      22257   .       C       T       127     .       DP=5;VDB=0.0765947;SGB=-0.590765;MQSB=1;MQ0F=0;AC=1;AN=1;DP4=0,0,2,3;MQ=60      GT:PL
CP000819.1      38971   .       A       G       225     .       DP=14;VDB=0.872139;SGB=-0.680642;MQSB=1;MQ0F=0;AC=1;AN=1;DP4=0,0,4,8;MQ=60      GT:PL
CP000819.1      42306   .       A       G       225     .       DP=15;VDB=0.969686;SGB=-0.686358;MQSB=1;MQ0F=0;AC=1;AN=1;DP4=0,0,5,9;MQ=60      GT:PL
CP000819.1      45277   .       A       G       225     .       DP=15;VDB=0.470998;SGB=-0.680642;MQSB=0.95494;MQ0F=0;AC=1;AN=1;DP4=0,0,7,5;MQ=60
CP000819.1      56613   .       C       G       183     .       DP=12;VDB=0.879703;SGB=-0.676189;MQSB=1;MQ0F=0;AC=1;AN=1;DP4=0,0,8,3;MQ=60      GT:PL
CP000819.1      62118   .       A       G       225     .       DP=19;VDB=0.414981;SGB=-0.691153;MQSB=0.906029;MQ0F=0;AC=1;AN=1;DP4=0,0,8,10;MQ=59
CP000819.1      64042   .       G       A       225     .       DP=18;VDB=0.451328;SGB=-0.689466;MQSB=1;MQ0F=0;AC=1;AN=1;DP4=0,0,7,9;MQ=60      GT:PL
~~~
{: .output}

This is a lot of information, so let's take some time to make sure we understand our output.

The first few columns represent the information we have about a predicted variation. 

| column | info |
| ------- | ---------- |
| CHROM | contig location where the variation occurs | 
| POS | position within the contig where the variation occurs | 
| ID | a `.` until we add annotation information | 
| REF | reference genotype (forward strand) | 
| ALT | sample genotype (forward strand) | 
| QUAL | Phred-scaled probability that the observed variant exists at this site (higher is better) |
| FILTER | a `.` if no quality filters have been applied, PASS if a filter is passed, or the name of the filters this variant failed | 

In an ideal world, the information in the `QUAL` column would be all we needed to filter out bad variant calls.
However, in reality we need to filter on multiple other metrics. 

The last two columns contain the genotypes and can be tricky to decode. 

| column | info |
| ------- | ---------- |
| FORMAT | lists in order the metrics presented in the final column | 
| results | lists the values associated with those metrics in order | 

For our file, the metrics presented are GT:PL:GQ. 

| metric | definition | 
| ------- | ---------- |
| GT | the genotype of this sample which for a diploid genome is encoded with a 0 for the REF allele, 1 for the first ALT allele, 2 for the second and so on. So 0/0 means homozygous reference, 0/1 is heterozygous, and 1/1 is homozygous for the alternate allele. For a diploid organism, the GT field indicates the two alleles carried by the sample, encoded by a 0 for the REF allele, 1 for the first ALT allele, 2 for the second ALT allele, etc. |
| PL | the likelihoods of the given genotypes |
| GQ | the Phred-scaled confidence for the genotype | 
| AD, DP | the depth per allele by sample and coverage |

The Broad Institute's [VCF guide](https://www.broadinstitute.org/gatk/guide/article?id=1268) is an excellent place
to learn more about the VCF file format.

> ## Exercise
> 
> Use the `grep` and `wc` commands you've learned to assess how many variants are in the vcf file. 
>
>> ## Solution
>> 
>> ~~~
>> $ grep -v "#" results/vcf/SRR2584866_final_variants.vcf | wc -l
>> ~~~
>> {: .bash}
>> 
>> ~~~ 
>> 766
>> ~~~
>> {: .output}
>>
>> There are 766 variants in this file.
> {: .solution}
{: .challenge}

## Assess the alignment (visualization) - optional step

It is often instructive to look at your data in a genome browser. Visualization will allow you to get a "feel" for 
the data, as well as detecting abnormalities and problems. Also, exploring the data in such a way may give you 
ideas for further analyses.  As such, visualization tools are useful for exploratory analysis. In this lesson we 
will describe two different tools for visualization: a light-weight command-line based one and the Broad
Institute's Integrative Genomics Viewer (IGV) which requires
software installation and transfer of files.

In order for us to visualize the alignment files, we will need to index the BAM file using `samtools`:

~~~
$ qrsh
$ briglo/samtools/1.9
$ samtools index results/bam/SRR2584866.aligned.sorted.bam
$ # We'll exit in a moment.
~~~
{: .bash}

### Viewing with `tview`

[Samtools](http://www.htslib.org/) implements a very simple text alignment viewer based on the GNU
`ncurses` library, called `tview`. This alignment viewer works with short indels and shows [MAQ](http://maq.sourceforge.net/) consensus. 
It uses different colors to display mapping quality or base quality, subjected to users' choice. Samtools viewer is known to work with a 130 GB alignment swiftly. Due to its text interface, displaying alignments over network is also very fast.

In order to visualize our mapped reads, we use `tview`, giving it the sorted bam file and the reference file: 

~~~
$ # Still on interactive node with module loaded
$ samtools tview results/bam/SRR2584866.aligned.sorted.bam data/ref_genome/ecoli_rel606.fasta
~~~
{: .bash}

~~~
1         11        21        31        41        51        61        71        81        91        101       111       121
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATAC
..................................................................................................................................
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ..................N................. ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,........................
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ..................N................. ,,,,,,,,,,,,,,,,,,,,,,,,,,,.............................
...................................,g,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ....................................   ................
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,....................................   ....................................      ,,,,,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ....................................  ,,a,,,,,,,,,,,,,,,,,,,,,,,,,,,,,     .......
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, .............................  ,,,,,,,,,,,,,,,,,g,,,,,    ,,,,,,,,,,,,,,,,,,,,,,,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ...........................T.......   ,,,,,,,,,,,,,,,,,,,,,,,c,          ......
......................... ................................   ,g,,,,,,,,,,,,,,,,,,,      ...........................
,,,,,,,,,,,,,,,,,,,,, ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ,,,,,,,,,,,,,,,,,,,,,,,,,,,       ..........................
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   ................................T..  ..............................   ,,,,,,
...........................       ,,,,,,g,,,,,,,,,,,,,,,,,   ....................................         ,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,, ....................................  ...................................        ....
....................................  ........................  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,      ....
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
........................            .................................. .............................     ....
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,   ....................................        ..........................
...............................       ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ....................................
...................................  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ..................................
.................................... ,,,,,,,,,,,,,,,,,,a,,,,,,,,,,,,,,,,,        ,,,,,,,,,,,,,,,,,,,,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ............................ ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
~~~
{: .output}

The first line of output shows the genome coordinates in our reference genome. The second line shows the reference
genome sequence. The third line shows the consensus sequence determined from the sequence reads. A `.` indicates
a match to the reference sequence, so we can see that the consensus from our sample matches the reference in most
locations. That is good! If that wasn't the case, we should probably reconsider our choice of reference.

Below the horizontal line, we can see all of the reads in our sample aligned with the reference genome. Only 
positions where the called base differs from the reference are shown. You can use the arrow keys on your keyboard
to scroll or type `?` for a help menu. To navigate to a specific position, type `g`. A dialogue box will appear. In
this box, type the name of the "chromosome" followed by a colon and the position of the variant you would like to view
(e.g. for this sample, type `CP000819.1:50` to view the 50th base. Type `Ctrl^C` or `q` to exit `tview`. 

> ## Exercise 
> 
> Visualize the alignment of the reads for our `SRR2584866` sample. What variant is present at 
> position 4377265? What is the canonical nucleotide in that position? 
> 
>> ## Solution
>> 
>> ~~~
>> $ samtools tview ~/dc_workshop/results/bam/SRR2584866.aligned.sorted.bam ~/dc_workshop/data/ref_genome/ecoli_rel606.fasta
>> ~~~
>> {: .bash}
>> 
>> Then type `g`. In the dialogue box, type `CP000819.1:4377265`. 
>> `G` is the variant. `A` is canonical. This variant possibly changes the phenotype of this sample to hypermutable. It occurs
>> in the gene *mutL*, which controls DNA mismatch repair.
>> Exit the compute node when you are done.
> {: .solution}
{: .challenge}

### Viewing with IGV

[IGV](http://www.broadinstitute.org/igv/) is a stand-alone browser, which has the advantage of being installed locally and providing fast access. 
Web-based genome browsers, like [Ensembl](http://www.ensembl.org/index.html) or the [UCSC browser](https://genome.ucsc.edu/), are slower, but provide more functionality. 
They not only allow for more polished and flexible visualization, but also provide easy access to a wealth of annotations and external data sources. 
This makes it straightforward to relate your data with information about repeat regions, known genes, epigenetic features or areas of cross-species conservation, to name just a few.

In order to use IGV, we will need to transfer some files to our local machine. 
We know how to do this with `sshfs` or `scp`. 

Next, we need to open the IGV software. 
If you haven't done so already, you can download IGV from the [Broad Institute's software page](https://www.broadinstitute.org/software/igv/download), double-click the `.zip` file to unzip it, and then drag the program into your Applications folder (Mac). 

1. Open IGV.
2. Load our reference genome file (`ecoli_rel606.fasta`) into IGV using the **"Load Genomes from File..."** option under the **"Genomes"** pull-down menu.
3. Load our BAM file (`SRR2584866.aligned.sorted.bam`) using the **"Load from File..."** option under the **"File"** pull-down menu. 
4.  Do the same with our VCF file (`SRR2584866_final_variants.vcf`).

Your IGV browser should look like the screenshot below:

![IGV](../img/igv-screenshot.png)

There should be two tracks: one coresponding to our BAM file and the other for our VCF file. 

In the **VCF track**, each bar across the top of the plot shows the allele fraction for a single locus. The second bar shows
the genotypes for each locus in each *sample*. We only have one sample called here, so we only see a single line. Dark blue = 
heterozygous, Cyan = homozygous variant, Grey = reference.  Filtered entries are transparent.

Zoom in to inspect variants you see in your filtered VCF file to become more familiar with IGV. See how quality information 
corresponds to alignment information at those loci.
Use [this website](http://software.broadinstitute.org/software/igv/AlignmentData) and the links therein to understand how IGV colors the alignments.

Now that we've run through our workflow for a single sample, we want to repeat this workflow for our other five
samples. However, we don't want to type each of these individual steps again five more times. That would be very
time consuming and error-prone, and would become impossible as we gathered more and more samples. Luckily, we
already know the tools we need to use to automate this workflow and run it on as many files as we want using a
single line of code. Those tools are: wildcards, for loops, and bash scripts. We'll use all three in the next 
lesson. 

> ## Installing Software
> 
> It's worth noting that all of the software we are using for this workshop has been pre-installed on the cluster. 
> This saves us a lot of time - installing software can be a time-consuming and frustrating task.
> However, this does mean that you won't be able to start doing these analyses on your own computer. 
> You'll need to install the software first. 
> Look at the [setup instructions](http://www.datacarpentry.org/wrangling-genomics/setup.html) for more information on installing these software packages.
> Installing software on your laptop is actually _relatively_ straightforward because you usually have adminstrator privileges for your laptop.
> Installing on the cluster is harder, because you don't have administrator privileges and you often have to compile software from source.
> Most biofinformatics software will include an option for compiling from source, but it helps to be aware that this is a thing.
> Another challenge is that (as of August 2020) the operating system on the cluster is quite old, and so often the software that you want install will depend on a core operating system library that either doesn't exist on the cluster or which is out of date.
> There are ways around this but the cluster is about to get an overhaul, which should minimise this particular headache.
> For now, my key advise is to seek help if you experience difficulty software on the cluster.
> 
> The other option for installing software on the cluster is `conda`.
> This tool takes care of all the dependencies for you, and installs software in isolated `environments` to avoid conflicting requirements.
> There is a `miniconda` module available on the cluster so please don't go cluttering up your home directory by installing `conda` there.
> You can find a [brief guide](https://intranet.gimr.garvan.org.au/display/BINF/Using+Community+Conda+on+the+Cluster) to setting up `conda` in the [Community User Guides](https://intranet.gimr.garvan.org.au/display/BINF/Community+User+Guides) section of the Bioinformatics corner of Confluence.
> Note that you will need to run `conda init` the first time that you use `conda`.
> Then you will need to log out and back in again for the configuration changes to take effect.
> After that there is no need to load the `miniconda` module each time - you will be able to use and install `conda` software as needed.
> Just make a note of which environment you are in and PLEASE do not install anything in the `(base)` environment.
{: .callout}

> ## BWA Alignment options
> BWA consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence 
> reads up to 100bp, while the other two are for sequences ranging from 70bp to 1Mbp. BWA-MEM and BWA-SW share similar features such 
> as long-read support and split alignment, but BWA-MEM, which is the latest, is generally recommended for high-quality queries as it 
> is faster and more accurate. 
{: .callout}


