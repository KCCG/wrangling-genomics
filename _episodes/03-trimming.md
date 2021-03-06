---
title: "Trimming and Filtering"
teaching: 30
exercises: 25
questions:
- "How can I get rid of sequence data that doesn't meet my quality standards?"
objectives:
- "Clean FASTQ reads using Trimmomatic."
- "Select and set multiple options for command-line bioinformatic tools."
- "Write `for` loops with two variables."
keypoints:
- "The options you set for the command-line tools you use are important!"
- "Data cleaning is an essential step in a genomics workflow."
---

# Cleaning Reads

In the previous episode, we took a high-level look at the quality
of each of our samples using FastQC. We visualized per-base quality
graphs showing the distribution of read quality at each base across
all reads in a sample and extracted information about which samples
fail which quality checks. Some of our samples failed quite a few quality metrics used by FastQC. This doesn't mean,
though, that our samples should be thrown out! It's very common to have some quality metrics fail, and this may or may not be a problem for your downstream application. For our variant calling workflow, we will be removing some of the low quality sequences to reduce our false positive rate due to sequencing error.

We will use a program called
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) to
filter poor quality reads and trim poor quality bases from our samples.

First we will need to find the appropriate module and get it working.
Don't forget to get a compute node with `qrsh`.

> ## Exercise
>
> 1) Figure out which trimmomatic module to use.  
> 2) What happens when you try to load the module and type "trimmomatic" to test the installation?
>
>> ## Solution
>> 1) gi/trimmomatic/0.36
>> 2) Command not found
> {: .solution}
{: .challenge}

Trimmomatic is a software tool that has been written using the Java programming language.
As a general rule, Java programs are a bit of a headache.
Let's do some detective work by taking a look at the Quick Start guide at the Trimmomatic link above.

> ## Exercise
>
> What are the first three components of the command suggested in the Quick Start guide?
>
>> ## Solution
>> java -jar trimmomatic-xxx.jar
>> Here xxx is the version number.
> {: .solution}
{: .challenge}

To run trimmomatic, we are going to have to find the .jar file.

> ## Exercise
>
> Think back to the introductory overview of the cluster architecture.
> Which volume is used to store modules?
>
>> ## Solution
>> /share/ClusterShare
> {: .solution}
{: .challenge}

> ## Exercise
>
> Navigate (`cd`) to the directory in the solution to the previous exercise. 
> Follow your nose from there and try to locate the .jar file for trimmomatic. 
> Clue: think about the module name from the first exercise in this chapter.
>
>> ## Solution
>> /share/ClusterShare/software/contrib/gi/trimmomatic/0.36/trimmomatic.jar
> {: .solution}
{: .challenge}

Phew. We are finally ready to test the software.

> ## Exercise
>
> Run the following command to test the installation
> ~~~
> $ java -jar /share/ClusterShare/software/contrib/gi/trimmomatic/0.36/trimmomatic.jar
> ~~~
>
>> ## Solution
>> ~~~
>> Error occurred during initialization of VM
>> Could not allocate metaspace: 1073741824 bytes
>> ~~~
> {: .solution}
{: .challenge}

Grrrh. Did I mention that Java tools can be a headache?
This is because Java software runs inside a "virtual machine" (VM) -- a kind of computer inside a computer.
You don't need to understand all the technicalities, but you do need to allocate enough memory for the VM.
That's what the "metaspace" message is all about.

When you request a compute node with the `qrsh` command by default you get given a certain amount of memory.
By default (I think) you get two CPU cores and 4 GB of memory per core for a total of about 8 GB. 
For most applications, this is plenty.
But sometimes you need to request more.
After a bit of trial and error, it turns out that you need about 16 gigabytes of RAM.
You can request it like this (log out of the compute node first):

~~~
$ qrsh -l mem_requested=16000M
~~~
{: .bash}

Note that 16 GB is about 16000 MB.
Now you can finally run the basic trimmomatic command to verify that it is installed.

~~~
$ java -jar /share/ClusterShare/software/contrib/gi/trimmomatic/0.36/trimmomatic.jar
~~~
{: .bash}

~~~
Usage: 
       PE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...
   or: 
       SE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-quiet] <inputFile> <outputFile> <trimmer1>...
   or: 
       -version
~~~
{: .output}

A lot of commands output a "usage" summary like the one above when you type the command without any arguments.
It is also quite common for commands to display a usage summary when they can't match the arguments that you provided against the arguments that they were expecting.
If this happens, check the command arguments that you entered carefully against the command usage.

One final comment about Java tools is that it often helps to put bounds on the minimum and maximum amount of memory that can be give to the VM.
We do this with `-X` options, such as "-Xmx".
Google it if you want to know details, but at least be aware that this kind of tweaking is often part and parcel of getting Java tools to run.


## Trimmomatic Options

Trimmomatic has a variety of options to trim your reads, as summarised by the "Usage" above.

This output shows us that we must first specify whether we have paired end (`PE`) or single end (`SE`) reads.
Next, we specify what flag we would like to run. For example, you can specify `threads` to indicate the number of
processors on your computer that you want Trimmomatic to use. In most cases using multiple threads (processors) can help to run the trimming faster. These flags are not necessary, but they can give you more control over the command. The flags are followed by positional arguments, meaning the order in which you specify them is important. 
In paired end mode, Trimmomatic expects the two input files, and then the names of the output files. These files are described below. While, in single end mode, Trimmomatic will expect 1 file as input, after which you can enter the optional settings and lastly the name of the output file.

| option    | meaning |
| ------- | ---------- |
|  \<inputFile1>  | Input reads to be trimmed. Typically the file name will contain an `_1` or `_R1` in the name.|
| \<inputFile2> | Input reads to be trimmed. Typically the file name will contain an `_2` or `_R2` in the name.|
|  \<outputFile1P> | Output file that contains surviving pairs from the `_1` file. |
|  \<outputFile1U> | Output file that contains orphaned reads from the `_1` file. |
|  \<outputFile2P> | Output file that contains surviving pairs from the `_2` file.|
|  \<outputFile2U> | Output file that contains orphaned reads from the `_2` file.|

The last thing trimmomatic expects to see is the trimming parameters:

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

However, a complete command for Trimmomatic will look something like the command below. This command is an example and will not work, as we do not have the files it refers to:

~~~
$ java -Xmx4000M -jar /share/ClusterShare/software/contrib/gi/trimmomatic/0.36/trimmomatic.jar \
              PE -threads 4 SRR_1056_1.fastq SRR_1056_2.fastq  \
              SRR_1056_1.trimmed.fastq SRR_1056_1un.trimmed.fastq \
              SRR_1056_2.trimmed.fastq SRR_1056_2un.trimmed.fastq \
              ILLUMINACLIP:SRR_adapters.fa SLIDINGWINDOW:4:20
~~~
{: .bash}

In this example, we've told Trimmomatic:

| code   | meaning |
| ------- | ---------- |
| `PE` | that it will be taking a paired end file as input |
| `-threads 4` | to use four computing threads to run (this will spead up our run) |
| `SRR_1056_1.fastq` | the first input file name |
| `SRR_1056_2.fastq` | the second input file name |
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

Now we will run Trimmomatic on our data. To begin, navigate to your `untrimmed_fastq` data directory:

~~~
$ cd ~/course/data/untrimmed_fastq
~~~
{: .bash}

We are going to run Trimmomatic on one of our paired-end samples. 
While using FastQC we saw that Nextera adapters were present in our samples. 
The adapter sequences came with the installation of trimmomatic, and you can find them in your "untrimmed_fastq" directory as "NexteraPE-PE.fa".

We will also use a sliding window of size 4 that will remove bases if their
phred score is below 20 (like in our example above). We will also
discard any reads that do not have at least 25 bases remaining after
this trimming step. This command will take a few minutes to run.

~~~
$ java -Xmx4000M -jar /share/ClusterShare/software/contrib/gi/trimmomatic/0.36/trimmomatic.jar \
                PE SRR2589044_1.fastq.gz SRR2589044_2.fastq.gz \
                SRR2589044_1.trim.fastq.gz SRR2589044_1un.trim.fastq.gz \
                SRR2589044_2.trim.fastq.gz SRR2589044_2un.trim.fastq.gz \
                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 
~~~
{: .bash}


~~~
TrimmomaticPE: Started with arguments:
 SRR2589044_1.fastq.gz SRR2589044_2.fastq.gz SRR2589044_1.trim.fastq.gz SRR2589044_1un.trim.fastq.gz SRR2589044_2.trim.fastq.gz SRR2589044_2un.trim.fastq.gz SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
Multiple cores found: Using 2 threads
Using PrefixPair: 'AGATGTGTATAAGAGACAG' and 'AGATGTGTATAAGAGACAG'
Using Long Clipping Sequence: 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
Using Long Clipping Sequence: 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
Using Long Clipping Sequence: 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'
Using Long Clipping Sequence: 'CTGTCTCTTATACACATCTGACGCTGCCGACGA'
ILLUMINACLIP: Using 1 prefix pairs, 4 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Quality encoding detected as phred33
Input Read Pairs: 1107090 Both Surviving: 885220 (79.96%) Forward Only Surviving: 216472 (19.55%) Reverse Only Surviving: 2850 (0.26%) Dropped: 2548 (0.23%)
TrimmomaticPE: Completed successfully
~~~
{: .output}

> ## Exercise
>
> Use the output from your Trimmomatic command to answer the
> following questions.
>
> 1) What percent of reads did we discard from our sample?
> 2) What percent of reads did we keep both pairs?
>
>> ## Solution
>> 1) 0.23%
>> 2) 79.96%
> {: .solution}
{: .challenge}

You may have noticed that Trimmomatic automatically detected the
quality encoding of our sample. It is always a good idea to
double-check this or to enter the quality encoding manually.

We can confirm that we have our output files:

~~~
$ ls SRR2589044*
~~~
{: .bash}

~~~
SRR2589044_1.fastq.gz       SRR2589044_1un.trim.fastq.gz  SRR2589044_2.trim.fastq.gz
SRR2589044_1.trim.fastq.gz  SRR2589044_2.fastq.gz         SRR2589044_2un.trim.fastq.gz
~~~
{: .output}

The output files are also FASTQ files. It should be smaller than our
input file, because we've removed reads. We can confirm this:

~~~
$ ls SRR2589044* -l -h
~~~
{: .bash}

~~~
-rw-rw-r-- 1 user user 124M Jul  6 20:22 SRR2589044_1.fastq.gz
-rw-rw-r-- 1 user user  94M Jul  6 22:33 SRR2589044_1.trim.fastq.gz
-rw-rw-r-- 1 user user  18M Jul  6 22:33 SRR2589044_1un.trim.fastq.gz
-rw-rw-r-- 1 user user 128M Jul  6 20:24 SRR2589044_2.fastq.gz
-rw-rw-r-- 1 user user  91M Jul  6 22:33 SRR2589044_2.trim.fastq.gz
-rw-rw-r-- 1 user user 271K Jul  6 22:33 SRR2589044_2un.trim.fastq.gz
~~~
{: .output}

***(Please "exit" the compute node promptly as soon as the command finishes successfully.)***

We've just successfully run Trimmomatic on one of our FASTQ files!

## Submitting jobs to the cluster

We've finally got Trimmomatic working, but we had to request extra memory.
And if a class full of people all request exta memory at the same time, this will leave less resources for the people doing real computational analysis.
**(Polite reminder to exit the compute node if you haven't already.)**
Fortunately there is another, more efficient way to use the computing resources on the cluster.
This involves creating a `job script` and submitting a `job` to the cluster.

A `job script` is basically the same as a bash script, but with a few extra features.
Let's look at an example.

~~~
$ less --line-numbers ~/course/scripts/trim_one.sh
~~~
{: .bash}

The content of this job script file looks like this.

~~~
#$ -S /bin/bash
#$ -N trim_one 
#$ -wd /home/user/course/data/untrimmed_fastq
#$ -pe smp 4
#$ -l mem_requested=5G
#$ -M user@garvan.org.au
#$ -m bae

# === Modules ===
module load gi/trimmomatic/0.36
JAR_DIR=/share/ClusterShare/software/contrib/gi/trimmomatic/0.36/


# === Parameters ===
SAMPLE=SRR2589044                            # Sample number (as used in fastq files for that sample)
ADAPTER=ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 # Adapter sequences
# Log which sample we are processing
echo "Processing $SAMPLE with $ADAPTER using $NSLOTS cores"

# === Main script body ===
# Define input and output files
# Paired end input files 1 and 2
INPUT_1=${SAMPLE}_1.fastq.gz
INPUT_2=${SAMPLE}_2.fastq.gz
# Output files containing surviving reads (that were not removed)
SURVIVING_1=${SAMPLE}_1.trim.fastq.gz
SURVIVING_2=${SAMPLE}_2.trim.fastq.gz
# Output files containing orphan reads (that were trimmed out)
ORPHAN_1=${SAMPLE}_1un.trim.fastq.gz
ORPHAN_2=${SAMPLE}_2un.trim.fastq.gz
# Log input and output file definitions:
echo "Input files:" $INPUT_1 $INPUT_2
echo "Surviving trimmed reads:" $SURVIVING_1 $SURVIVING_2
echo "Orphaned reads:" $ORPHAN_1 $ORPHAN_2

# Run the trimmomatic command with the specified input files and adapter
# See trimmomatic manual for other options: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

java -Xmx4000M -jar $JAR_DIR/trimmomatic.jar PE -threads $NSLOTS \
      $INPUT_1 $INPUT_2 \
      $SURVIVING_1 $SURVIVING_2 \
      $ORPHAN_1 $ORPHAN_2 \
      SLIDINGWINDOW:4:20 MINLEN:25 $ADAPTER
~~~
{: .output}

This job script does exactly the same thing as the single command that we just ran a moment ago.
On the face of it, that might seem like a lot of extra work.
I'm going to try to convince you that it is worth it.

First, let's look at the macro structure of the job script.
At the top of the script, there are half a dozen lines that look like comments.
These lines starting with "#$" are actually instructions to the Sun Grid Engine about how to run the job.
We'll look at these more closely in a moment.

Right below the job instructions, we have a section called "=== Modules ===".
This is where we load the software needed for the job.
Usually this section will just consist of one or two "module load ..." commands.
Some jobs don't actually load any special software, in which case you can skip this section entirely.
But remember how much trouble we went to in order to get trimmomatic working?
Here we define a variable (JAR_DIR) to hold the path to the .jar file for trimmomatic. 
This extra line is not strictly necessary, but I hope you'll agree that they make the main body of the script clearer (see below).
This also makes the script easier to maintain if you want to upgrade to a newer version of trimmomatic in the future.

Next, we have a section called "=== Parameters ===".
This is the section that you would edit if you wanted to run this job again on a different sample, or in a different directory, or on samples that used different adapter sequences.
(Almost) everything else about the script can stay exactly the same, you just have to edit these two variables: SAMPLE and ADAPTER.
(You might also want to change the job name and the working directory.)
Different jobs will need different parameters, and sometimes you can skip this section entirely. 
But by putting all the parameters up near the top of the script we make it clear to the user (most likely your future self) what kind of information the script needs to do its job.
This makes it **much** easier to reuse and recycle your job scripts.

You can actually write your bash scripts and job scripts so that they process parameters (command line arguments) just like regular commands.
We won't go down that particular rabbit hole right now, but it is a neat thing to be able to do.
There is a concise explanation [here](https://tecadmin.net/tutorial/bash-scripting/bash-command-arguments/) or check out the "Shell" courses on DataCamp.

Note the `echo` statement at the end of the "Parameters" section.
Normally the `echo` command just echoes text to the terminal.
But in a job script these messages will be included in the job log.
Liberally including `echo` statements makes your jobs self-documenting, and also makes debugging **much** easier.

The next session is called "=== Main script body ===".
This is where the actual work gets done.
We start by creating variables to hold the names of all the different input and output files.
Again, this is not strictly necessary -- you could just hard-code the file names into the `trimmomatic` command.
But trimmomatic works with six different files (in paired-end mode anyway) and the names of these files are all related to the sample name in a systematic way.
By defining just one variable (SAMPLE) in the "Parameters" section, we can then calculate the names of all the files.
If you want to process a different sample, you only have to change one part of the script.
Not only is this less work, it also avoids copy-paste errors and ensures that output files are named consistently.
These output files will become the inputs for the next step of our workflow, so consistency is crucial -- especially when we try to automate the entire workflow at the end of this workshop.

A simpler job script might only deal with one or two files, say "INPUT" and "OUTPUT".
In this case, you might include these variables in the "Parameters" section.
The names and structure is flexible -- I'm just trying to give you an example of good style.

More `echo` messages to keep a record of exactly which files this job will process.
If something goes wrong you can compare these logs to the actual files in the directory.
Sometimes it is worth checking to make sure that the expected input files actually exist.
That way you can log a user-friendly messages if they are missing.

Finally, we run the same java command that we ran interactively a moment ago.
Hopefully the structure is this command is a lot clearer now that we have gone to the trouble of creating variables for all of the different files.

How do you run a `job script`? With the `qsub` command. Example below, but **don't try this just yet.**

~~~
$ qsub trim_one.sh
Your job 2608293 ("trim_one") has been submitted
~~~

The `qsub` command sends the `job script` (in this case "trim_one.sh") to the Sun Grid Engine (SGE).
The SGE then looks at all the `#$` instructions to figure out how much memory and other resources you are requesting (among other things).
If the requested resources are not available, then your job goes in a `queue` until other jobs finish running (or other people log out of their `qrsh` sessions!)
Then the SGE allocates the necessasry resources to your job, including a compute node, and runs the rest of the job script (everything below the `#$` lines) on the compute node.
As soon as the job finishes, the resources are released and made available for other jobs.

### Job instructions

The best explanation of job instructions that I have found is [Evan Ben's cheatsheet](https://intranet.gimr.garvan.org.au/pages/viewpage.action?pageId=74712562) on Confluence.
Here I will just explain the instructions in the job script above.
Open the script in nano, Atom or Sublime so that you can make a few changes as you read through.

| instruction | meaning | 
| ----------- | ------- |
| #$ -S /bin/bash | This tells the Sun Grid Engine to use bash for the shell. There are alternatives, but we won't worry about them here. |
| #$ -N trim_one | This gives the job a name, which can be anything you like and does **not** have to match the file name of the job script. |
| #$ -wd | This defines the work directory for the job. The SGE will `cd` to this directory before running any of the other commands. Use "-cwd" for "current working directory" when appropriate, but then you'll have to keep a record of what directory you were in. Annoyingly, you can't use bash shortcuts like `~` here or bash environment variables like `$HOME`. |
| #$ -pe smp 4 | This instruction requests four CPU cores. Use the $NSLOTS variable to tell the commands within your script to actually use all the cores. $NSLOTS is one of several variables created by SGE when you submit a job. |
| #$ -l mem_requested=5G | This requests extra memory (RAM). **Note** that the total amount of memory allocated will be this value multiplied by the number of cores, so do some calculations and don't go overboard. In this case we will get 4 x 5 = 20 GB. |
| #$ -M user@garvan.org.au | This specifies an email address for notifications. Change it to your address!!! |
| #$ -m bae | This specifies when notifications will be issued -- at the **b**eginning, when a job **a**borts or when it **e**nds. You decide how much you want to get spammed, but some kind of notification is helpful, especially for long-running jobs. |

If that seems like a lot to remember then don't worry, you don't have to.
You'll find a template file at `~/templates/job.sh`. 
You can start with this template file and just make a few minor changes, such as the job name.

### qrsh versus qsub

So, when do you use `qrsh` and when should you use `qsub` instead?
The line can be a bit blurry sometimes, but here are some rough rules of thumb.
* `qrsh` is great when you are doing something "one-off", such as installing some software or downloading a file.
* `qrsh` is also good for testing and development -- working out which module to use and how to use it, figuring out how much memory you need, or which options you need to use with a command.
* Once you know what module/commands/options/steps you need to do, then there is a lot to be gained from spelling it all out in a job script and submitting the job via `qsub`. You have a record of exactly what you did, and SGE provides logs of the output and any errors. 
* You can even use `qsub` jobs for relatively mundane tasks if you are confident that you know exactly how to do it. Examples include downloading a bunch of files, or moving files to or from another server.

**Note**: You can't run the `qsub` command from an interactive compute node that you requested with the `qrsh` command. 
Quit out of the interactive node, and then submit the job from the login node.

### Job logs

When you type a command directly into the command line the output usually appears in the terminal, unless it is redirected to a file.
Even when you run a `bash script`, the output from each step usually appears on the screen.
However with job scripts the output **never** appears on the screen -- it is always redirected to log files.

By the time you have read up to here, the job you submitted a moment ago should have finished executing.
If you remembered to change the email address, you will receive summary information via email.
But any messages produced by the job will be stored in one of four `job logs`.

> ## Exercise
>
> Navigate to the working directory for the job script.
>
>> ## Solution
>>
>> Look at the `#$ -wd` job instruction in the job script.
>> ~~~
>> cd /home/user/course/data/untrimmed_fastq
>> ls
>> ~~~
> {: .solution}
{: .challenge}

In the job working directory, you should see the following four files (along with the various .fastq files).
~~~
$ ls ~/course/untrimmed_fastq
trim_one.o1234567
trim_one.e1234567
trim_one.po1234567
trim_one.po1234567
~~~
{: .bash}

Here "trim_one" is of course the name of the job that we just submitted.
And the seven-digit number at the end is the job number.
The content of each of these files is summarized below:

| File | Content |
| ---- | ------- |
| *.o* | The `Standard output` stream from the job script, usually abbreviated as `STDOUT`. |
| *.e* | The `Standard error` stream from the job script, usually abbreviated as `STDERR`. |
| *.po* | Output associated with parallel processing. |
| *.pe* | Error messages associated with parallel processing. |

Note that "error" here doesn't necessarily mean that something has gone wrong.
Some programmers choose to send progress update messages to `STDERR` as well as actual error messages.
Note also that I've **never** seen anything interesting in either the `po` log or the `pe` log.
But I guess there is a first time for everything, so I make a point of checking anyway.

If you look inside the `o` file you should see the `echo` statements from the job script.
These messages can be very handy for debugging.
Always, always make sure that you check the `e` file as well.
Sometimes everything seems to have worked just fine, but then a message in the `e` file will alert you that all is not how it seems.

Occassioanlly, jobs will "hang".
You get a message that they started, but they never seem to finish.
You can use the following command to check on the status of your job.

~~~
$ qstat -u $USER
~~~
{: .bash}

In fact, this is such a handy little snippet of code that it is almost worth making an alias for it.
If you are cruising for time, add the following line to your ~/.bash_aliases file.

~~~
alias myjobs=qstat -u $USER
~~~
{: .bash}

There are other handy commands for checking up on jobs or deleting frozen jobs, etc.
Consult the [HPC guide](https://intranet.gimr.garvan.org.au/display/HPC/Wolfpack+Cluster+-+Getting+started) for details.
But let's get back to the bioinformatics.

## Trimming all the files

Trimmomatic can only operate on one sample at a time and we have more than one sample. 
The good news is that we can use a `for` loop to iterate through our sample files quickly! 

We unzipped one of our files before to work with it, let's compress it again before we run our for loop.

~~~
gzip SRR2584863_1.fastq 
~~~
{: .bash}

Now all of the samples are starting off the same.

### Making a job script

Open the job script we were looking at a moment ago. 
You can use `nano` if you like, but my strong recommendation is that you use `Atom` or `Sublime` via `sshfs`.
Change the email address to match **your email address**, and save the file so that you can use it as a template each time you want to make a new job script.
Then immediately save the file under a new name (such as "trim_all.sh") to prevent overwriting the template file.

We'll break down the process of creating the `job script` into small steps.
If you get stuck on a particular step, you can peek at the solution.
Or if you are running behind then feel free to skip to the end and just type out the solution.
You'll still learn something from seeing how it fits together, but I'd encourage you to try to figure it out step by step.

> ## Exercise
>
> In trim_one.sh we only processed one sample.  
> Now we want to process three samples.  
> What are the sample codes for the three samples?
>
>> ## Solution
>> ~~~
>> $ cd ~/course/data/untrimmed_fastq
>> $ ls
>> NexteraPE-PE.fa        SRR2584863_2.fastq.gz  SRR2584866_2.fastq.gz  SRR2589044_2.fastq.gz
>> SRR2584863_1.fastq.gz  SRR2584866_1.fastq.gz  SRR2589044_1.fastq.gz
>> ~~~
>> The sample codes are SRR2584863, SRR2584866 and SRR2589044
> {: .solution}
{: .challenge}

> ## Exercise
>
> One way to create a `for` loop is to just hard code the sample codes into the loop.  
> Write a `for` loop that loops over the three sample codes and simply `echo`s the codes to the terminal.  
> (No need to write a job script for this.
> A simple `for` loop typed directly into the terminal will do.)
>
>> ## Solution
>> You can simply write the sample codes, separated by spaces, as follows:
>> ~~~
>> $ for sample in SRR2584863 SRR2584866 SRR2589044; do echo $sample; done
>> ~~~
>> You can enter the `for` loop over multiple lines if you prefer.
> {: .solution}
{: .challenge}

> ## Exercise
>
> Hard-coding details like this in the main body of script is difficult to read, maintain and debug.
> We can parametize the list of sample codes as follows.
> (Still no need to put this in a job script, use the command line to experiment.)
> ~~~
> $ sample_list="SRR2584863 SRR2584866 SRR2589044"
> $ for sample in $sample_list; do echo $sample; done
> ~~~
> That's an improvement, but what if we didn't know the sample names or file names in advance?  
> 1) How can we use wildcards to get a list of the sample file names ending with "\_1.fastq.gz"?  
> 2) How can we extract the sample names from this list of file names? Hint: loop over the result from 1)
>
>> ## Solution
>> 1) Use `*` to match any number of characters. You can store the result in a variable. For example..
>> ~~~
>> $ file_list=*_1.fastq.gz
>> $ echo $file_list
>> SRR2584863_1.fastq.gz SRR2584866_1.fastq.gz SRR2589044_1.fastq.gz
>> ~~~
>> 
>> 2) Use `basename` to extract just the sample code
>> ~~~
>> $ file_list=*_1.fastq.gz
>> for file in $file_list; do 
>>     SAMPLE=$(basename $file _1.fastq.gz)
>>     echo $SAMPLE
>> done
>> ~~~
>> Here the `basename` command looks at each $file and strips off "_1.fastq.gz" from the end.
>> The `$( )` surrounding the call to the `basename` command means "Run the command inside the parentheses and return the result."
>> Here we assign that result to a variable called "SAMPLE".
> {: .solution}
{: .challenge}

> ## Exercise
>
> Now that we have figured out how to loop over all of our samples, let's make a job script.
> Merge the solution from the previous exercise with the content of "trim_one.sh".
> Think about which bits belong in the "Parameters" section, and which bits belong in the main script body.
> In particular, think carefully about what needs to go inside the `for` loop, and what can be left as a kind of prologue.
> 
> There is no shame in cheating here by peeking at the solution.
> If you do, try to understand how the solution works then replicate it yourself.
> If possible, wait until you have finished before checking your answer against the solution again.
> But don't feel too bad about having another peek if necessary.
>
>> ## Solution
>> Note that here I am defining the work directory as a variable (WORKDIR) rather than as a job instruction (#$ -wd).
>> That's because I want to define the work directory relative to my home directory, but I can't do that with the job instruction without hard coding the path.
>> I also wanted to show you an alternative way of handling the work directory issue.
>> Be aware that job logs will still be output to the directory defined by the job instruction ("#$ -wd" or "#$ -cwd").
>> ~~~
>> #$ -S /bin/bash
>> #$ -N trim_all
>> #$ -cwd
>> #$ -pe smp 4
>> #$ -l mem_requested=5G
>> #$ -M j.reeves@garvan.org.au
>> #$ -m bae
>> 
>> # === Modules ===
>> module load gi/trimmomatic/0.36
>> JAR_DIR=/share/ClusterShare/software/contrib/gi/trimmomatic/0.36/
>> 
>> # === Parameters ===
>> WORKDIR=$HOME/course/data/untrimmed_fastq/
>> ADAPTER=ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 # Adapter sequences
>> echo "Processing $WORKDIR with $ADAPTER using $NSLOTS cores"
>> 
>> #=== Main script body ===
>> # Move to work directory
>> cd $WORKDIR
>> echo "Trimming files in work directory:" $WORKDIR
>> # Get list of files to work with
>> FILE_LIST=*_1.fastq.gz
>> echo "Starting with the following file list:"  $FILE_LIST
>> 
>> # Loop over each sample
>> for file in $FILE_LIST; do
>>    SAMPLE=$(basename $file _1.fastq.gz)
>>    # Define input and output files
>>    # Paired end input files 1 and 2
>>    INPUT_1=${SAMPLE}_1.fastq.gz
>>    INPUT_2=${SAMPLE}_2.fastq.gz
>>    # Output files containing surviving reads (that were not removed)
>>    SURVIVING_1=${SAMPLE}_1.trim.fastq.gz
>>    SURVIVING_2=${SAMPLE}_2.trim.fastq.gz
>>    # Output files containing orphan reads (that were trimmed out)
>>    ORPHAN_1=${SAMPLE}_1un.trim.fastq.gz
>>    ORPHAN_2=${SAMPLE}_2un.trim.fastq.gz
>>    # Log which sample we are processing
>>    echo "================================"
>>    echo "Processing $SAMPLE with $ADAPTER"
>>    echo "Input files:" $INPUT_1 $INPUT_2
>>    echo "Surviving trimmed reads:" $SURVIVING_1 $SURVIVING_2
>>    echo "Orphaned reads:" $ORPHAN_1 $ORPHAN_2
>>    # Run the trimmomatic command with the specified input files and adapter
>>    # See trimmomatic manual for other options:
>>    # http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
>>    java -Xmx4000M -jar $JAR_DIR/trimmomatic.jar PE -threads $NSLOTS \
>>          $INPUT_1 $INPUT_2 \
>>          $SURVIVING_1 $SURVIVING_2 \
>>          $ORPHAN_1 $ORPHAN_2 \
>>          SLIDINGWINDOW:4:20 MINLEN:25 $ADAPTER
>> done
~~~
> {: .solution}
{: .challenge}

Go ahead and submit the job. 
It should take a few minutes for Trimmomatic to run for each of our six input files. 
Once it's done running, take a look at your directory contents. 
You'll notice that even though we ran Trimmomatic on file `SRR2589044` before running the for loop, there is only one set of files for it. 
Because we matched the ending `_1.fastq.gz`, we re-ran Trimmomatic on this file, overwriting our first results. 
That's ok, but it's good to be aware that it happened.

~~~
$ ls
~~~
{: .bash}

~~~
NexteraPE-PE.fa               SRR2584866_1.fastq.gz         SRR2589044_1.trim.fastq.gz
SRR2584863_1.fastq.gz         SRR2584866_1.trim.fastq.gz    SRR2589044_1un.trim.fastq.gz
SRR2584863_1.trim.fastq.gz    SRR2584866_1un.trim.fastq.gz  SRR2589044_2.fastq.gz
SRR2584863_1un.trim.fastq.gz  SRR2584866_2.fastq.gz         SRR2589044_2.trim.fastq.gz
SRR2584863_2.fastq.gz         SRR2584866_2.trim.fastq.gz    SRR2589044_2un.trim.fastq.gz
SRR2584863_2.trim.fastq.gz    SRR2584866_2un.trim.fastq.gz
SRR2584863_2un.trim.fastq.gz  SRR2589044_1.fastq.gz
~~~
{: .output}

We've now completed the trimming and filtering steps of our quality
control process! Before we move on, let's tidy up by moving our trimmed FASTQ files
to a new subdirectory within our `data/` directory.

~~~
$ cd ~/course/data/untrimmed_fastq
$ mkdir ../trimmed_fastq
$ mv *.trim* ../trimmed_fastq
$ cd ../trimmed_fastq
$ ls
~~~
{: .bash}

~~~
SRR2584863_1.trim.fastq.gz    SRR2584866_1.trim.fastq.gz    SRR2589044_1.trim.fastq.gz
SRR2584863_1un.trim.fastq.gz  SRR2584866_1un.trim.fastq.gz  SRR2589044_1un.trim.fastq.gz
SRR2584863_2.trim.fastq.gz    SRR2584866_2.trim.fastq.gz    SRR2589044_2.trim.fastq.gz
SRR2584863_2un.trim.fastq.gz  SRR2584866_2un.trim.fastq.gz  SRR2589044_2un.trim.fastq.gz
~~~
{: .output}

Actually, this kind of "tidying up" is quite common after running a job that produces a lot of output files.
Although the tidy up commands above are not particularly intensive and so don't really warrant a job script of their own, it is quite common practice to include this kind of tidy up step at the end of the job script.
If you like, copy and paste the `mkdir` and `mv` commands above at the bottom of your "trim_all.sh" job script (after the "done" statement).


Now that our samples have gone through quality control, they should perform better on the quality tests run by FastQC.
Later we will re-run FastQC on the **trimmed** FASTQ files and visualize the HTML files to see whether your per base sequence quality is higher after trimming.

