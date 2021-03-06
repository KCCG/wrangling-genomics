---
title: "Assessing Read Quality"
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
---

# Bioinformatic workflows

When working with high-throughput sequencing data, the raw reads you get off of the sequencer will need to pass
through a number of  different tools in order to generate your final desired output. The execution of this set of
tools in a specified order is commonly referred to as a *workflow* or a *pipeline*. 

An example of the workflow we will be using for our variant calling analysis is provided below with a brief
description of each step. 

![workflow](../img/variant_calling_workflow.png)


1. Quality control - Assessing quality using FastQC
2. Quality control - Trimming and/or filtering reads (if necessary)
3. Align reads to reference genome 
4. Perform post-alignment clean-up
5. Variant calling

These workflows in bioinformatics adopt a plug-and-play approach in that the output of one tool can be easily
used as input to another tool without any extensive configuration. Having standards for data formats is what 
makes this feasible. Standards ensure that data is stored in a way that is generally accepted and agreed upon 
within the community. The tools that are used to analyze data at different stages of the workflow are therefore 
built under the assumption that the data will be provided in a specific format.  

# Starting with Data

Often times, the first step in a bioinformatic workflow is getting the data you want to work with onto a computer where you can work with it. If you have outsourced sequencing of your data, the sequencing center will usually provide you with a link that you can use to download your data. Today we will be working with publicly available sequencing data.

We are studying a population of *Escherichia coli* (designated Ara-3), which were propagated for more than 50,000 generations in a glucose-limited minimal medium. We will be working with three samples from this experiment, one from 5,000 generations, one from 15,000 generations, and one from 50,000 generations. The population changed substantially during the course of the experiment, and we will be exploring how with our variant calling workflow. 

The data are paired-end, so we will download two files for each sample. We will use the [European Nucleotide Archive](https://www.ebi.ac.uk/ena) to get our data. The ENA "provides a comprehensive record of the world's nucleotide sequencing information, covering raw sequencing data, sequence assembly information and functional annotation." The ENA also provides sequencing data in the fastq format, an important format for sequencing reads that we will be learning about today. 

The commands below show how you would download the data from ENA. **You don't have to do this -- the files have already been downloaded for you.**
But take advantage of some of the time saved by **not** having to download the files to peruse this code and understand what you would normally have to do.

Here we are using the `-p` option for `mkdir`. This option allows `mkdir` to create the new directory, even if one of the parent directories doesn't already exist. It also supresses errors if the directory already exists, without overwriting that directory. 

It will take about 15 minutes to download the files.
**At least, it would if you did actually had to download them.**

~~~
mkdir -p ~/course/data/untrimmed_fastq/
cd ~/course/data/untrimmed_fastq

curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_2.fastq.gz 
~~~
{: .bash}


The data comes in a compressed format, which is why there is a `.gz` at the end of the file names. This makes it faster to transfer, and allows it to take up less space on our computer. Let's unzip one of the files so that we can look at the fastq format.

~~~
$ gunzip SRR2584863_1.fastq.gz 
~~~
{: .bash}

# Quality Control

We will now assess the quality of the sequence reads contained in our fastq files. 

![workflow_qc](../img/var_calling_workflow_qc.png)
## Details on the FASTQ format

Although it looks complicated (and it is), we can understand the
[fastq](https://en.wikipedia.org/wiki/FASTQ_format) format with a little decoding. Some rules about the format
include...

|Line|Description|
|----|-----------|
|1|Always begins with '@' and then information about the read|
|2|The actual DNA sequence|
|3|Always begins with a '+' and sometimes the same info in line 1|
|4|Has a string of characters which represent the quality scores; must have same number of characters as line 2|

We can view the first complete read in one of the files our dataset by using `head` to look at
the first four lines. 

~~~
$ head -n 4 SRR2584863_1.fastq 
~~~
{: .bash}

~~~
@SRR2584863.1 HWI-ST957:244:H73TDADXX:1:1101:4712:2181/1
TTCACATCCTGACCATTCAGTTGAGCAAAATAGTTCTTCAGTGCCTGTTTAACCGAGTCACGCAGGGGTTTTTGGGTTACCTGATCCTGAGAGTTAACGGTAGAAACGGTCAGTACGTCAGAATTTACGCGTTGTTCGAACATAGTTCTG
+
CCCFFFFFGHHHHJIJJJJIJJJIIJJJJIIIJJGFIIIJEDDFEGGJIFHHJIJJDECCGGEGIIJFHFFFACD:BBBDDACCCCAA@@CA@C>C3>@5(8&>C:9?8+89<4(:83825C(:A#########################
~~~
{: .output}

Line 4 shows the quality for each nucleotide in the read. Quality is interpreted as the 
probability of an incorrect base call (e.g. 1 in 10) or, equivalently, the base call 
accuracy (e.g. 90%). To make it possible to line up each individual nucleotide with its quality
score, the numerical score is converted into a code where each individual character 
represents the numerical quality score for an individual nucleotide. For example, in the line
above, the quality score line is: 

~~~
CCCFFFFFGHHHHJIJJJJIJJJIIJJJJIIIJJGFIIIJEDDFEGGJIFHHJIJJDECCGGEGIIJFHFFFACD:BBBDDACCCCAA@@CA@C>C3>@5(8&>C:9?8+89<4(:83825C(:A#########################
~~~
{: .output}

The numerical value assigned to each of these characters depends on the 
sequencing platform that generated the reads. The sequencing machine used to generate our data 
uses the standard Sanger quality PHRED score encoding, using Illumina version 1.8 onwards.
Each character is assigned a quality score between 0 and 41 as shown in 
the chart below.

~~~
Quality encoding: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ
                   |         |         |         |         |
Quality score:    01........11........21........31........41                                
~~~
{: .output}

Each quality score represents the probability that the corresponding nucleotide call is
incorrect. This quality score is logarithmically based, so a quality score of 10 reflects a
base call accuracy of 90%, but a quality score of 20 reflects a base call accuracy of 99%. 
These probability values are the results from the base calling algorithm and depend on how 
much signal was captured for the base incorporation. 

Looking back at our read: 

~~~
@SRR2584863.1 HWI-ST957:244:H73TDADXX:1:1101:4712:2181/1
TTCACATCCTGACCATTCAGTTGAGCAAAATAGTTCTTCAGTGCCTGTTTAACCGAGTCACGCAGGGGTTTTTGGGTTACCTGATCCTGAGAGTTAACGGTAGAAACGGTCAGTACGTCAGAATTTACGCGTTGTTCGAACATAGTTCTG
+
CCCFFFFFGHHHHJIJJJJIJJJIIJJJJIIIJJGFIIIJEDDFEGGJIFHHJIJJDECCGGEGIIJFHFFFACD:BBBDDACCCCAA@@CA@C>C3>@5(8&>C:9?8+89<4(:83825C(:A#########################
~~~
{: .output}

we can now see that there is a range of quality scores, but that the end of the sequence is
very poor (`#` = a quality score of 2). 

> ## Exercise
> 
> What is the last read in the `SRR2584863_1.fastq ` file? How confident
> are you in this read? 
> 
>> ## Solution
>> ~~~
>> $ tail -n 4 SRR2584863_1.fastq 
>> ~~~
>> {: .bash}
>> 
>> ~~~
>> @SRR2584863.1553259 HWI-ST957:245:H73R4ADXX:2:2216:21048:100894/1
>> CTGCAATACCACGCTGATCTTTCACATGATGTAAGAAAAGTGGGATCAGCAAACCGGGTGCTGCTGTGGCTAGTTGCAGCAAACCATGCAGTGAACCCGCCTGTGCTTCGCTATAGCCGTGACTGATGAGGATCGCCGGAAGCCAGCCAA
>> +
>> CCCFFFFFHHHHGJJJJJJJJJHGIJJJIJJJJIJJJJIIIIJJJJJJJJJJJJJIIJJJHHHHHFFFFFEEEEEDDDDDDDDDDDDDDDDDCDEDDBDBDDBDDDDDDDDDBDEEDDDD7@BDDDDDD>AA>?B?<@BDD@BDC?BDA?
>> ~~~
>> {: .output}
>> 
>> This read has more consistent quality at its end than the first 
>> read that we looked at, but still has a range of quality scores, 
>> most of them high. We will look at variations in position-based quality
>> in just a moment.
>> 
> {: .solution}
{: .challenge}

## Loading software on the cluster

There are a lot of software tools for doing bioinformatics.
Keeping all of these tools loaded in memory all of the time would waste a lot of resources for things that you might only need occassionally.
There is a LOT of software installed on the cluster, but to use it you first need to load it.

### Searching for software modules

To see a list of all the software modules installed on the cluster, type the following command:

~~~
$ module avail
~~~
{: .bash}

There might be a brief pause before you see any response. 
In fact, there could be quite a long pause. Be patient!
The first few lines of the output will look something like this:

~~~
--------------------------------- /share/ClusterShare/Modules/modulefiles/contrib ----------------------------------
aarsta/bcftools/1.2                                gi/trim_galore/0.3.7
aarsta/bcftools/1.6                                gi/trim_galore/0.4.0
aarsta/bedtools/2.20.1                             gi/trimmomatic/0.30
aarsta/bismark/0.13.0                              gi/trimmomatic/0.32
~~~
{: .output}

Actually, there are so many lines of output that it will probably scroll past you faster than you can keep up.
Before we learn how to interpret this output, let's use pipes and redirection to break it into bite sized pieces.
Type the following command (we'll worry about how it works in a moment):

~~~
module avail 2>&1 | less
~~~
{: .bash}

You should see the modules listed one per line, in a "pageable" format.
You can press the space bar to page down, or `b` to go back.
In fact, all of your favourite `less` shortcuts will work. 
For example, `/` to search, `n` to go to the next match and `q` to quit back to the command line.

Type "/fastqc" to see which versions of the `fastqc` tool are available. 
Then press "n" a few times to jump to the next match.
Press "q" to quit when you are done.

As we saw in the Shell Genomics workshop, anytime you have a command that produces way too much output for you to visually process, you can pipe the results through `less`.
But here we need to do some fancy pants redirection to get the pipe to work.
That because, in their wisdom, the authors of the `module` tool chose to print the results of "module avail" to the standard error output (STDERR) rather than to the standard output (STDOUT).
Now both error messages (STDERR) and regular output (STDOUT) both normally get displayed on the terminal, so the distinction is not very important most of the time.
But when you are piping and redirecting sometimes you need to be more specific.
Here `2>` means "redirect the standard error stream (STDERR)" and `&1` is a "file descriptor" for the standard output (STDOUT).
If all that gives you a headache, just remember "2>&1 |" as an idiom for "pipe the standard error stream to the next command". 

Or, better still, you can create an `alias`.
Actually, I saved you the trouble.

### Handy aliases

If you type "alias" at the command line, you should see an alias called "mods" which does just this.

~~~
$ alias
...
...
alias mods="module avail 2>&1 | less"
alias modgrep="module avail 2>&1 | grep"
...
~~~
{: .bash}

So from now on, just type "mods" to scroll through the available modules. 
There is another alias called `modgrep` that you can use when you are searching for a **particular** module.
Try it out to get a sense how it works.
Again, if you feel that "modgrep" is a dumb name then feel free to edit ~/.bash_aliases to change it to something more memorable.
You could even add `--ignore-case` after grep to make sure that "modgrep fastqc" matches "Fastqc" and "FastQC" as well as "fastqc".
(Most module creators use all lowercase, but there is no guarantee that this is consistent.)

### Finding the right module

Let's go back and have a look at the list of modules.
The basic pattern is "username/software/version".  
Sometimes you might see a group name instead of a username.
Occassionally you might see "username/software/compiler/version".
These module names all look like file paths.
There is a reason for this, which we'll get to in a moment.

Generally speaking, you will usually want the latest version of the software or at least the latest version that is already available on the cluster.
Sometimes you might want to use an older version in order to reproduce somebody's workflow, or follow a tutorial (like this one!)
At this point, you can basically ignore the compiler information. 
Just be sure to check all of the available versions, not just the first one that you see.
Bioinformatics tools evolve quite quickly, and there can be big differences even between minor version numbers.
In particular, make sure that you right versions of the tools for this workshop.


I'm not going to tell you exactly which modules to use, because figuring that out is half the trick.
Be warned, not every module is fully functional.
In some cases, they need a complementary module to also be loaded at the same time.
But the documentation on this is often sparse or non-existent.

In general, the "gi" group make pretty good modules.
You might also recognise the user names of bioinformaticians in your group.
In any case, if you are having trouble with a module don't be shy about reaching out to the person who made it.
That's why we include the usernames as the first part of the module name.

### Module location

Module names look like file paths because that's what they are -- paths relative to the base module directory.
And where is this base module directory of which I speak?
Let's find out.
~~~
$ echo $MODULEPATH
/share/ClusterShare/Modules/modulefiles/contrib:/usr/share/Modules/modulefiles:/etc/modulefiles
$ cd /share/ClusterShare/Modules/modulefiles/contrib
$ ls
~~~
{: .bash}

You should see a whole bunch of directories named after six-letter Garvan usernames, as well a couple of directories with shorter names for groups.
If you drill down a couple of more directories you will find the actual module definition files.

~~~
$ cd gi
$ ls
$ cd fastqc
$ ls
$ less 0.11.5
~~~
{: .bash}

Don't worry too much about the content of the module definition file.
We won't be making any modules as part of this course.
But if you ever do need to make a module, then you now know where to put the definition files.
And sometimes looking at a module file can give you a better idea of how it works.

Note that the module definition file is not the actual software.
That lives in another directory on the same volume (/share/ClusterShare) with a parallel directory structure.

~~~
$ cd /share/ClusterShare
$ ls
$ cd software
$ ls
$ cd contrib
$ ls
~~~
{: .bash}

Once again, you should see a lot of directories based on usernames.
You'll also see a bunch of other random stuff that probably shouldn't be there.
If you do need to install software on the cluster, please try to be tidy by keeping it all in the appropriate directory.

Let's drill down into where the software actually lives.

~~~
$ cd gi
$ ls
$ cd fastqc
$ ls
$ cd 0.11.5
$ ls
$ less README.txt
~~~
{: .bash}

In this directory you can find all of the files associated with the software that actually runs when you use a module.
Most of the time you can happily ignore all of these files.
But every now and then you might want to have a look at the "README" file or maybe check some of the configuration.
This is where to find it.

We won't be covering software installation in this course (unlesss we get some extra time at the end).
But note that you can install software here **without** making it into a module, if it is for your own personal use.
Of course, chances are that other people would benefit from your software too, so making a module is nice community-minded thing to do.
You don't necessarily have to make a module, but I **do** encourage you to install software on /share/ClusterShare rather than in your home directory.
Otherwise you risk blowing through your disk quota quite quickly.

### Compute nodes

Before we start loading software, let's pause to think about where we are and what we are doing.
Take a look at your prompt. 
Mine looks like this (except colourful).

~~~
johree@dice02:~:$
~~~

The "dice02" is there to remind me that I am currently logged in to the "login node" called "dice02".
Most of the time, there will be a bunch of other people also logged in to this computer.
If I start doing anything computationally intensive then the performance of the whole node will suffer, making life difficult for everyone else.
And bioinformatics software can be **very** computationally intensive.
So any time you want to do that kind of heavy lifting, make a habit of requesting a "compute node" first by using the `qrsh` command.
A "compute node" is your own private computer in the cloud, and you can request whatever resources you think you will need.
Just don't forget to release those resources when you are not using them by logging out from the compute node.

~~~
$ qrsh
$ pwd
$ ls
$ exit
~~~
{: .bash}

Notice how your prompt changes when you connect to a compute node.
There is no need to note down the name of the node, but looking at your prompt can be a helpful reminder of where you are.
Also note how all the files in your home directory and other volumes are still available to you from the compute node.
The CPU and memory have changed, but you can take all your files with you and you still have access to all the same data.
It's a bit like plugging an enormous USB memory stick into another computer.

### Loading modules

Once you have foundthe module you want, loading it is easy. Just type "module load".
For example,

~~~
$ qrsh
$ modgrep fastqc
danrod/fastqc/0.10.1
danrod/fastqc/0.11.4
gi/fastqc/0.10.1
gi/fastqc/0.11.1
gi/fastqc/0.11.2
gi/fastqc/0.11.3
gi/fastqc/0.11.5
leemar/fastqc/0.11.5
marcow/fastqc/prebuilt/0.10.1

$ module load gi/fastqc/0.11.5
$ module list
Currently Loaded Modulefiles:
  1) rocks-openmpi      2) gi/fastqc/0.11.5
  
$ module unload gi/fastqc/0.11.5
$ module list
Currently Loaded Modulefiles:
  1) rocks-openmpi
$ exit
~~~
{: .bash}

### Loading the fastqc module

We will use fastqc version 0.11.5 by the GI group (that's Derrick and Manuel, by the way).
Load the module and test that it is working by displaying the help menu.

~~~
$ qrsh
$ module load gi/fastqc/0.11.5
$ fastqc -h
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
{: .bash}

## Assessing Quality using FastQC
In real life, you won't be assessing the quality of your reads by visually inspecting your 
FASTQ files. Rather, you'll be using a software program to assess read quality and 
filter out poor quality reads. We'll first use a program called [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to visualize the quality of our reads. 
Later in our workflow, we'll use another program to filter out poor quality reads. 

FastQC has a number of features which can give you a quick impression of any problems your
data may have, so you can take these issues into consideration before moving forward with your
analyses. Rather than looking at quality scores for each individual read, FastQC looks at
quality collectively across all reads within a sample. The image below shows one FastQC-generated plot that indicates
a very high quality sample:

![good_quality](../img/good_quality1.8.png)

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

![bad_quality](../img/bad_quality1.8.png)

Here, we see positions within the read in which the boxes span a much wider range. Also, quality scores drop quite low into the "bad" range, particularly on the tail end of the reads. The FastQC tool produces several other diagnostic plots to assess sample quality, in addition to the one plotted above. 

## Running FastQC  

We will now assess the quality of the reads that we downloaded. First, make sure you're still in the `untrimmed_fastq` directory

~~~
$ cd ~/course/data/untrimmed_fastq/ 
~~~
{: .bash}

> ## Exercise
> 
>  How big are the files?
> (Hint: Look at the options for the `ls` command to see how to show
> file sizes.)
>
>> ## Solution
>>  
>> ~~~
>> $ ls -l -h
>> ~~~
>> {: .bash}
>> 
>> ~~~
>> -rw-rw-r-- 1 dcuser dcuser 545M Jul  6 20:27 SRR2584863_1.fastq
>> -rw-rw-r-- 1 dcuser dcuser 183M Jul  6 20:29 SRR2584863_2.fastq.gz
>> -rw-rw-r-- 1 dcuser dcuser 309M Jul  6 20:34 SRR2584866_1.fastq.gz
>> -rw-rw-r-- 1 dcuser dcuser 296M Jul  6 20:37 SRR2584866_2.fastq.gz
>> -rw-rw-r-- 1 dcuser dcuser 124M Jul  6 20:22 SRR2589044_1.fastq.gz
>> -rw-rw-r-- 1 dcuser dcuser 128M Jul  6 20:24 SRR2589044_2.fastq.gz
>> ~~~
>> {: .output}
>> 
>> There are six FASTQ files ranging from 124M (124MB) to 545M. 
>> 
> {: .solution}
{: .challenge}

FastQC can accept multiple file names as input, and on both zipped and unzipped files, so we can use the \*.fastq* wildcard to run FastQC on all of the FASTQ files in this directory.

~~~
$ fastqc *.fastq* 
~~~
{: .bash}

You will see an automatically updating output message telling you the 
progress of the analysis. It will start like this: 

~~~
Started analysis of SRR2584863_1.fastq
Approx 5% complete for SRR2584863_1.fastq
Approx 10% complete for SRR2584863_1.fastq
Approx 15% complete for SRR2584863_1.fastq
Approx 20% complete for SRR2584863_1.fastq
Approx 25% complete for SRR2584863_1.fastq
Approx 30% complete for SRR2584863_1.fastq
Approx 35% complete for SRR2584863_1.fastq
Approx 40% complete for SRR2584863_1.fastq
Approx 45% complete for SRR2584863_1.fastq
~~~
{: .output}

In total, it should take about five minutes for FastQC to run on all
six of our FASTQ files. When the analysis completes, your prompt
will return. So your screen will look something like this:

~~~
Approx 80% complete for SRR2589044_2.fastq.gz
Approx 85% complete for SRR2589044_2.fastq.gz
Approx 90% complete for SRR2589044_2.fastq.gz
Approx 95% complete for SRR2589044_2.fastq.gz
Analysis complete for SRR2589044_2.fastq.gz
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
SRR2584863_1.fastq        SRR2584866_1_fastqc.html  SRR2589044_1_fastqc.html
SRR2584863_1_fastqc.html  SRR2584866_1_fastqc.zip   SRR2589044_1_fastqc.zip
SRR2584863_1_fastqc.zip   SRR2584866_1.fastq.gz     SRR2589044_1.fastq.gz
SRR2584863_2_fastqc.html  SRR2584866_2_fastqc.html  SRR2589044_2_fastqc.html
SRR2584863_2_fastqc.zip   SRR2584866_2_fastqc.zip   SRR2589044_2_fastqc.zip
SRR2584863_2.fastq.gz     SRR2584866_2.fastq.gz     SRR2589044_2.fastq.gz
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
$ mkdir -p ~/course/results/fastqc_untrimmed_reads 
$ mv *.zip ~/course/results/fastqc_untrimmed_reads/ 
$ mv *.html ~/course/results/fastqc_untrimmed_reads/ 
~~~
{: .bash}

Now we can navigate into this results directory and do some closer
inspection of our output files.

~~~
$ cd ~/course/results/fastqc_untrimmed_reads/ 
~~~
{: .bash}

## Viewing the FastQC results

If we were working on our local computers, we'd be able to look at 
each of these HTML files by opening them in a web browser.

However, these files are currently sitting on the cluster.
And, since we are only logging into the cluster via the 
command line - it doesn't have any web browser setup to display 
these files either.

So the easiest way to look at these webpage summary reports will be 
to transfer them to our local computers (i.e. your laptop).

The simplest way to do this is using `sshfs` (or `sshfs-win`).
If you connect to the cluster via sshfs and then navigate to the "results" directory, all you have to do is double click on an HTML file.
It should open in a browser.

Sometimes you do want to copy the results to your laptop, and the instructions below remind you how to do that with the `scp` command.
But in the interests of time we'll just use sshfs.
**Skip down to the exercise below**

To transfer a file from a remote server to our own machines, we will
use `scp`, which we learned yesterday in the Shell Genomics lesson. 

First we
will make a new directory on our computer to store the HTML files
we're transferring. Let's put it on our desktop for now. Open a new
tab in your terminal program (you can use the pull down menu at the
top of your screen or the Cmd+t keyboard shortcut) and type: 

~~~
$ mkdir -p ~/Desktop/fastqc_html 
~~~
{: .bash}

Now we can transfer our HTML files to our local computer using `scp`.

~~~
$ scp username@dice01:/share/ScratchGeneral/username/course/results/fastqc_untrimmed_reads/*.html ~/Desktop/fastqc_html
~~~
{: .bash}

If you try this, then don't forget to swap "username" for your actual username.

The second part starts with a `:` and then gives the absolute path
of the files you want to transfer from your remote computer. Don't
forget the `:`. We used a wildcard (`*.html`) to indicate that we want all of
the HTML files. 

The third part of the command gives the absolute path of the location
you want to put the files. This is on your local computer and is the 
directory we just created `~/Desktop/fastqc_html`. 

You should see a status output like this:

~~~
SRR2584863_1_fastqc.html                      100%  249KB 152.3KB/s   00:01    
SRR2584863_2_fastqc.html                      100%  254KB 219.8KB/s   00:01    
SRR2584866_1_fastqc.html                      100%  254KB 271.8KB/s   00:00    
SRR2584866_2_fastqc.html                      100%  251KB 252.8KB/s   00:00    
SRR2589044_1_fastqc.html                      100%  249KB 370.1KB/s   00:00    
SRR2589044_2_fastqc.html                      100%  251KB 592.2KB/s   00:00  
~~~
{: .output}

Now we can go to our new directory and open the 6 HTML files. 

Depending on your system, 
you should be able to select and open them all at once via a right click menu
in your file browser.

> ## Exercise
> 
> Discuss your results with a neighbor. Which sample(s) looks the best
> in terms of per base sequence quality? Which sample(s) look the
> worst?
> 
>> ## Solution
>> All of the reads contain usable data, but the quality decreases toward
>> the end of the reads.
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
+ [**Sequence Duplication Levels**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html): A distribution of duplicated sequences. In sequencing, we expect most reads to only occur once. If some sequences are occurring more than once, it might indicate enrichment bias (e.g. from PCR). If the samples are high coverage (or RNA-seq or amplicon), this might not be true.  
+ [**Overrepresented sequences**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/9%20Overrepresented%20Sequences.html): A list of sequences that occur more frequently than would be expected by chance. 
+ [**Adapter Content**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/10%20Adapter%20Content.html): a graph indicating where adapater sequences occur in the reads.
+ [**K-mer Content**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/11%20Kmer%20Content.html): a graph showing any sequences which may show a positional bias within the reads.

## Working with the FastQC text output

Now that we've looked at our HTML reports to get a feel for the data,
let's look more closely at the other output files. 
Make sure you're in our results subdirectory.   

~~~
$ cd ~/course/results/fastqc_untrimmed_reads/ 
$ ls 
~~~
{: .bash}

~~~
SRR2584863_1_fastqc.html  SRR2584866_1_fastqc.html  SRR2589044_1_fastqc.html
SRR2584863_1_fastqc.zip   SRR2584866_1_fastqc.zip   SRR2589044_1_fastqc.zip
SRR2584863_2_fastqc.html  SRR2584866_2_fastqc.html  SRR2589044_2_fastqc.html
SRR2584863_2_fastqc.zip   SRR2584866_2_fastqc.zip   SRR2589044_2_fastqc.zip
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
Archive:  SRR2584863_1_fastqc.zip
caution: filename not matched:  SRR2584863_2_fastqc.zip
caution: filename not matched:  SRR2584866_1_fastqc.zip
caution: filename not matched:  SRR2584866_2_fastqc.zip
caution: filename not matched:  SRR2589044_1_fastqc.zip
caution: filename not matched:  SRR2589044_2_fastqc.zip
~~~
{: .output}

This didn't work. We unzipped the first file and then got a warning
message for each of the other `.zip` files. This is because `unzip` 
expects to get only one zip file as input. We could go through and 
unzip each file one at a time, but this is very time consuming and 
error-prone. Someday you may have 500 files to unzip!

A more efficient way is to use a `for` loop like we learned in the Shell Genomics lesson to iterate through all of
our `.zip` files. Let's see what that looks like and then we'll 
discuss what we're doing with each line of our loop.

~~~
$ for filename in *.zip
> do
> unzip $filename
> done
~~~
{: .bash}

In this example, the input is six filenames (one filename for each of our `.zip` files).
Each time the loop iterates, it will assign a file name to the variable `filename`
and run the `unzip` command.
The first time through the loop,
`$filename` is `SRR2584863_1_fastqc.zip`. 
The interpreter runs the command `unzip` on `SRR2584863_1_fastqc.zip`.
For the second iteration, `$filename` becomes 
`SRR2584863_2_fastqc.zip`. This time, the shell runs `unzip` on `SRR2584863_2_fastqc.zip`.
It then repeats this process for the four other `.zip` files in our directory.


When we run our `for` loop, you will see output that starts like this:

~~~
Archive:  SRR2589044_2_fastqc.zip
   creating: SRR2589044_2_fastqc/
   creating: SRR2589044_2_fastqc/Icons/
   creating: SRR2589044_2_fastqc/Images/
  inflating: SRR2589044_2_fastqc/Icons/fastqc_icon.png  
  inflating: SRR2589044_2_fastqc/Icons/warning.png  
  inflating: SRR2589044_2_fastqc/Icons/error.png  
  inflating: SRR2589044_2_fastqc/Icons/tick.png  
  inflating: SRR2589044_2_fastqc/summary.txt  
  inflating: SRR2589044_2_fastqc/Images/per_base_quality.png  
  inflating: SRR2589044_2_fastqc/Images/per_tile_quality.png  
  inflating: SRR2589044_2_fastqc/Images/per_sequence_quality.png  
  inflating: SRR2589044_2_fastqc/Images/per_base_sequence_content.png  
  inflating: SRR2589044_2_fastqc/Images/per_sequence_gc_content.png  
  inflating: SRR2589044_2_fastqc/Images/per_base_n_content.png  
  inflating: SRR2589044_2_fastqc/Images/sequence_length_distribution.png  
  inflating: SRR2589044_2_fastqc/Images/duplication_levels.png  
  inflating: SRR2589044_2_fastqc/Images/adapter_content.png  
  inflating: SRR2589044_2_fastqc/fastqc_report.html  
  inflating: SRR2589044_2_fastqc/fastqc_data.txt  
  inflating: SRR2589044_2_fastqc/fastqc.fo  
~~~
{: .output}

The `unzip` program is decompressing the `.zip` files and creating
a new directory (with subdirectories) for each of our samples, to 
store all of the different output that is produced by FastQC. There
are a lot of files here. The one we're going to focus on is the 
`summary.txt` file. 

If you list the files in our directory now you will see: 

~~~
SRR2584863_1_fastqc       SRR2584866_1_fastqc       SRR2589044_1_fastqc
SRR2584863_1_fastqc.html  SRR2584866_1_fastqc.html  SRR2589044_1_fastqc.html
SRR2584863_1_fastqc.zip   SRR2584866_1_fastqc.zip   SRR2589044_1_fastqc.zip
SRR2584863_2_fastqc       SRR2584866_2_fastqc       SRR2589044_2_fastqc
SRR2584863_2_fastqc.html  SRR2584866_2_fastqc.html  SRR2589044_2_fastqc.html
SRR2584863_2_fastqc.zip   SRR2584866_2_fastqc.zip   SRR2589044_2_fastqc.zip
~~~
{:. output}

The `.html` files and the uncompressed `.zip` files are still present,
but now we also have a new directory for each of our samples. We can 
see for sure that it's a directory if we use the `-F` flag for `ls`. 

~~~
$ ls -F 
~~~
{: .bash}

~~~
SRR2584863_1_fastqc/      SRR2584866_1_fastqc/      SRR2589044_1_fastqc/
SRR2584863_1_fastqc.html  SRR2584866_1_fastqc.html  SRR2589044_1_fastqc.html
SRR2584863_1_fastqc.zip   SRR2584866_1_fastqc.zip   SRR2589044_1_fastqc.zip
SRR2584863_2_fastqc/      SRR2584866_2_fastqc/      SRR2589044_2_fastqc/
SRR2584863_2_fastqc.html  SRR2584866_2_fastqc.html  SRR2589044_2_fastqc.html
SRR2584863_2_fastqc.zip   SRR2584866_2_fastqc.zip   SRR2589044_2_fastqc.zip
~~~
{: .output}

Let's see what files are present within one of these output directories.

~~~
$ ls -F SRR2584863_1_fastqc/ 
~~~
{: .bash}

~~~
fastqc_data.txt  fastqc.fo  fastqc_report.html	Icons/	Images/  summary.txt
~~~
{: .output}

Use `less` to preview the `summary.txt` file for this sample. 

~~~
$ less SRR2584863_1_fastqc/summary.txt 
~~~
{: .bash}

~~~
PASS    Basic Statistics        SRR2584863_1.fastq
PASS    Per base sequence quality       SRR2584863_1.fastq
PASS    Per tile sequence quality       SRR2584863_1.fastq
PASS    Per sequence quality scores     SRR2584863_1.fastq
WARN    Per base sequence content       SRR2584863_1.fastq
WARN    Per sequence GC content SRR2584863_1.fastq
PASS    Per base N content      SRR2584863_1.fastq
PASS    Sequence Length Distribution    SRR2584863_1.fastq
PASS    Sequence Duplication Levels     SRR2584863_1.fastq
PASS    Overrepresented sequences       SRR2584863_1.fastq
WARN    Adapter Content SRR2584863_1.fastq
~~~
{: .output}

The summary file gives us a list of tests that FastQC ran, and tells
us whether this sample passed, failed, or is borderline (`WARN`). Remember, to quit from `less` you must type `q`.

## Documenting Our Work

We can make a record of the results we obtained for all our samples
by concatenating all of our `summary.txt` files into a single file 
using the `cat` command. We'll call this `fastqc_summaries.txt` and move
it to `~/course/docs`.

~~~
$ cat */summary.txt > ~/course/docs/fastqc_summaries.txt 
~~~
{: .bash}

> ## Exercise
> 
> Which samples failed at least one of FastQC's quality tests? What
> test(s) did those samples fail?
>
>> ## Solution
>> 
>> We can get the list of all failed tests using `grep`. 
>> 
>> ~~~ 
>> $ cd ~/course/docs
>> $ grep FAIL fastqc_summaries.txt
>> ~~~
>> {: .bash}
>> 
>> ~~~
>> FAIL    Per base sequence quality       SRR2584863_2.fastq.gz
>> FAIL    Per tile sequence quality       SRR2584863_2.fastq.gz
>> FAIL    Per base sequence content       SRR2584863_2.fastq.gz
>> FAIL    Per base sequence quality       SRR2584866_1.fastq.gz
>> FAIL    Per base sequence content       SRR2584866_1.fastq.gz
>> FAIL    Adapter Content SRR2584866_1.fastq.gz
>> FAIL    Adapter Content SRR2584866_2.fastq.gz
>> FAIL    Adapter Content SRR2589044_1.fastq.gz
>> FAIL    Per base sequence quality       SRR2589044_2.fastq.gz
>> FAIL    Per tile sequence quality       SRR2589044_2.fastq.gz
>> FAIL    Per base sequence content       SRR2589044_2.fastq.gz
>> FAIL    Adapter Content SRR2589044_2.fastq.gz
>> ~~~
>> {: .output}
>> 
> {: .solution}
{: .challenge}


# Other notes  -- Optional 

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


> ## Same Symbols, Different Meanings
>
> Here we see `>` being used as a shell prompt, whereas `>` is also
> used to redirect output.
> Similarly, `$` is used as a shell prompt, but, as we saw earlier,
> it is also used to ask the shell to get the value of a variable.
>
> If the *shell* prints `>` or `$` then it expects you to type something,
> and the symbol is a prompt.
>
> If *you* type `>` or `$` yourself, it is an instruction from you that
> the shell should redirect output or get the value of a variable.
{: .callout}
