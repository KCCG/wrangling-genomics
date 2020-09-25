---
title: "Automating a Variant Calling Workflow"
teaching: 30
exercises: 15
questions:
- "How can I make my workflow more efficient and less error-prone?"
objectives:
- "Write a shell script with multiple variables."
- "Incorporate a `for` loop into a shell script."
keypoints:
- "We can combine multiple commands into a shell script to automate a workflow."
- "Use `echo` statements within your scripts to get an automated progress update."
---
The original Data Wrangling and Processing Genomics workshop finishes up by writing a script that encapsulates the entire workflow.
We've already been writing scripts in order to submit them as jobs.
But it's still worthwhile reviewing some of the fundamentals.

# What is a shell script?

You wrote a simple shell script in a [previous lesson](http://www.datacarpentry.org/shell-genomics/05-writing-scripts/) that we used to extract bad reads from our
FASTQ files and put them into a new file. 

Here's the script you wrote:

~~~
grep -B1 -A2 NNNNNNNNNN *.fastq > scripted_bad_reads.txt

echo "Script finished!"
~~~
{: .bash}

That script was only two lines long, but shell scripts can be much more complicated
than that and can be used to perform a large number of operations on one or many 
files. This saves you the effort of having to type each of those commands over for
each of your data files and makes your work less error-prone and more reproducible. 
For example, the variant calling workflow we just carried out had about eight steps
where we had to type a command into our terminal. 
Most of these commands were pretty long. 
If we wanted to do this for all six of our data files, that would be forty-eight steps. 
If we had 50 samples (a more realistic number), it would be 400 steps! 
You can see why we want to automate this.

We've also used `for` loops in previous lessons to iterate one or two commands over multiple input files. 
In these `for` loops, the filename was defined as a variable in the `for` statement, which enabled you to run the loop on multiple files. 
We will be using variable assignments like this in our new shell scripts.

Here's the `for` loop you wrote for unzipping `.zip` files: 

~~~
$ for filename in *.zip
> do
> unzip $filename
> done
~~~
{: .bash}

You also wrote a loop for running Trimmomatic on all of our `.fastq` sample files.
Here is the version from the original workship.
There are a few differences, but the structure is basically the same.
See if you can following the logic in the loop below.

~~~
$ for infile in *_1.fastq.gz
> do
>   base=$(basename ${infile} _1.fastq.gz)
>   trimmomatic PE ${infile} ${base}_2.fastq.gz \
>                ${base}_1.trim.fastq.gz ${base}_1un.trim.fastq.gz \
>                ${base}_2.trim.fastq.gz ${base}_2un.trim.fastq.gz \
>                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 
> done
~~~
{: .bash}

Notice that in this `for` loop, we used two variables, `infile`, which was defined in the `for` statement, and `base`, which was created from the filename during each iteration of the loop.

> ## Creating Variables
> Within the Bash shell you can create variables at any time (as we did
> above, and during the 'for' loop lesson). Assign any name and the
> value using the assignment operator: '='. You can check the current
> definition of your variable by typing into your script: echo $variable_name.
{: .callout}

In this lesson, we'll use two shell scripts to automate the variant calling analysis: one for FastQC analysis (including creating our summary file), and a second for the remaining variant calling. To write a script to run our FastQC analysis, we'll take each of the commands we entered to run FastQC and process the output files and put them into a single file with a `.sh` extension. The `.sh` is not essential, but serves as a reminder to ourselves and to the computer that this is a shell script.

# Analyzing Quality with FastQC

Let's write a `job script` for the quality control processing that we previously entered directly into the command line.
Open the job script template in Atom or Sublime (via sshfs) and immediately save it under a new name.

## === Job instructions ===

The first section of our job script consists of job instructions (staring wtih #$).
Review the instructions in the template, and change as appropriate.

> ## Exercise
> 1. What should we do about the number of cores (#$ -pe smp)?
> Look up the [fastqc documentation](https://manpages.debian.org/stretch/fastqc/fastqc.1.en.html) and see if this tool can use multiple threads/cores.
> 2. What about memory? Will you need to request extra memory?
>> ## Solution
>> 
>> 1. Yes, you can use the `--threads` option.
>> This option will process multiple files in parallel, using one file per thread.
>> So there is no point requesting more threads than the number of files to process.
>> 2. No. Fastqc needs 250 MB per thread (see the documentation) but SGE gives us 4 GB (=4000 MB) per core which is more than enough.
>> This command worked fine without extra memory when we ran it interactively before, so no need to request more memory as a job script.
>> In fact, you can delete the "#$ -l mem_requested" instruction if you want.
>> {: .bash}
> {: .solution}
{: .challenge}

## === Modules ===

Do you remember which module you used for FastQC?
Me neither. 
You can look back at chapter 2 or you can consult your `history`.

There are three good ways to quickly check your history.
1) Type <kbd>Control</kbd>-<kbd>R</kbd> and then start typing your search term (eg. "fastqc").
This will find the last time you typed that term.
As you type more letters, the search updates each time. 
Keep pressing <kbd>Control</kbd>-<kbd>R</kbd> repeatedly to scroll back through earlier matches.
2) Use the `hgrep` alias that I set up for you.
Type "alias" to see how it works under the hood.
Or just type "hgrep fastqc" to see a list of all the times you have ever used the search term "fastqc".
3) Pipe history through `less` by typing "history | less".
Then you can either page through your history with <kbd>space</kbd>, or jump to the most recent commands with <kbd>shift</kbd>-<kbd>g</kbd> and then scroll backwards with <kbd>b</kbd>. 
Alternatively, you can search your history by typing "/" followed by the search term (use "?" to search backwards) and tap <kbd>n</kbd> for the next match.
The advantage of this method is that you can see the context before and after the match as well, which is really helpful especially when you can't remember the name of a command but you **do** remember what you did just before or after that command.
This combination is so useful, you could even consider making an alias for it.
{: .callout}

Now that you have found the right module, enter a line in the script to load it.

## === Parameters ===

It might not be obvious what the parameters are this point.
But as you work your way through the script, look out for variables that belong here.

## === Main script body ===

There is a neat little hack that you can use in both job scripts and regular bash scripts.

~~~
set -e
~~~

The `set` command makes various "settings" to your bash environment.
The "-e" flag tells bash to terminate the whole script the moment any command in the script fails.
(More precisely, the moment that any command has a non-zero exit status. Google it if you are curious.)
This is so handy that you probably want to include it in pretty much all of your scripts.
You can put it right up near the top of the script, or at the beginning of the main script body.
Go ahead and update your template. I'll wait.

These next two lines should look something like this.
~~~
echo "Running FastQC in " $PWD
fastqc *.fastq*
~~~

The echo statement lets us know that the script managed to progress to this point.
By including a variable in the statement (in this case `PWD`) it also provides information that can be helpful for debugging.
If something goes wrong then knowing that you were in the wrong directory when the `fastqc` command executed might be just the clue you need.

> ## Exercise
> Which directory do you want to be in anyway?
>
>> ## Solution
>> The `.fastq` files that we want to analyse are in ~/course/data/untrimmed_fastq/
>> You might need to `cd` into that directory before you run the `fastqc` command.
>> This will depend on what you specified for the `#$ -wd` job instruction.
>> There is more than one way to do this -- just make sure that you think it through carefully.
>> You might want to define a variable to hold the location of the data files.
>> This might be a parameter.
>> (These are hints, not instructions. Find your own solution.)
>> {: .bash}
> {: .solution}
{: .challenge}

## == Tidying up ===

The `fastqc` command spits out a bunch of result files that you should tidy up neatly.

Our next line will create a new directory to hold our FastQC output files. 
You will need to use the `-p` option for `mkdir` again. 
It is a good idea to use this option in your shell scripts to avoid running into errors if you don't have the directory structure you think you do.

> ## Exercise
> 1. How can we create the following directory?
> ~~~
> ~/course/results/fastqc_untrimmed_reads
> ~~~
> Assume that you are not sure if ~/course/results/ exists.
> 
> 2. How can you parameterize this command?
>> ## Solution
>> 1. Use the `-p` option as follows.
>> ~~~
>> mkdir -p ~/course/results/fastqc_untrimmed_reads
>> ~~~
>> 2. There several ways to paramerize this.
>> a) Using an absolute path ...
>> ~~~
>> RESULT_DIR=~/course/results/fastqc_untrimmed_reads
>> mkdir -p $RESULT_DIR
>> ~~~
>> b) Or using a relative path ...
>> ~~~
>> RESULT_DIR=results/fastqc_untrimmed_reads
>> mkdir -p $RESULT_DIR
>> ~~~
>> This option assumes that your current directory is ~/course
>>
>> Whichever you decide on, include it in your job script rather than entering it in the terminal.
>> {: .bash}
> {: .solution}
{: .challenge}

As part of the tidy up, you need to move all of the files with a `.zip` or a `.html` extension to the directory we just created for storing our FastQC results. 
(Now you see why defining a variable was worth the effort.) 
Once again it is a good idea to include a status message.

~~~
echo "Saving FastQC results..."
mv *.zip $RESULT_DIR
mv *.html $RESULT_DIR
~~~
{: .output}

The next line moves us to the results directory where we've stored our output.
Variables may feel like extra work to set up in the first place, but often they save you a lot of time later.
More importantly, if you change something (such as where you want the results stored) you only have to change it in one place.

~~~
cd $RESULT_DIR
~~~
{: .output}

Hopefully by now you have realised that RESULTS_DIR is really a kind of parameter.
If necessary, move it to the "Parameters" block.

The next five lines should look very familiar. First we give ourselves a status message to tell us that we're unzipping our ZIP
files. Then we run our for loop to unzip all of the `.zip` files in this directory.

~~~
echo "Unzipping..."
for filename in *.zip
do
  unzip $filename
done
~~~
{: .output}

Next we concatenate all of our summary files into a single output file, with a status message to remind ourselves that this is 
what we're doing.

~~~
echo "Saving summary..."
cat */summary.txt > ~/course/docs/fastqc_summaries.txt
~~~
{: .output}

> ## Using `echo` statements
> 
> We've used `echo` statements to add progress statements to our script. 
> When you run a regular bash script, these statements will be printed to the terminal as the script runs, so that you can ee how far our script has progressed.
> If you include `echo` statements in a job script the messages won't be output to the screen but they will be save in the job log. 
> This can still be useful for debugging if things go wrong.
> The basic structure is as follows:
> 1. Write a comment explaining what the variable means eg: 
> ~~~
> # This variable holds the name of the input file
> ~~~
> 2. Define the variable, eg: 
> ~~~
> INPUT=/mydir/input.txt
> ~~~
> 3. Echo the value of the variable for debugging purposes, eg: 
> ~~~
> echo "Input file:" $INPUT
> ~~~
> 4. Use the variable, eg:
> ~~~
> mv $INPUT data/raw/
> ~~~
{: .callout}

Your full job script should now look something like this:

~~~
#!/bin/bash
#$ -S /bin/bash
#$ -N quality_control 
#$ -wd ~/course/
#$ -pe smp 4
#$ -l mem_requested=4G
#$ -M user@garvan.org.au
#$ -m bae

# === Modules ===
module load gi/fastqc/0.11.5

# === Parameters ===
DATA_DIR=~/course/data/untrimmed_fastq
RESULTS_DIR=~/course/results/fastqc_untrimmed_reads
DOC_DIR=~/course/docs

# === Main script body ===
echo "Running FastQC ..."
fastqc $DATA_DIR/*.fastq*

mkdir -p $RESULTS_DIR

# === Tidy up ===
echo "Saving FastQC results..."
mv *.zip $RESULTS_DIR
mv *.html $RESULTS_DIR

cd $RESULTS_DIR
echo "Unzipping..."
for filename in *.zip
do
  unzip $filename
done

echo "Saving summary..."
cat */summary.txt > $DOC_DIR/fastqc_summaries.txt
~~~
{: .output}

Save your file. 
We can now submit our job:

~~~
$ qsub read_qc.sh
~~~
{: .bash}

When you get a notification email telling you that your job is complete, check the output in the `<job-name>.o<job-number>` log.
It should look something like this.

~~~
Running FastQC ...
Started analysis of SRR2584866.fastq
Approx 5% complete for SRR2584866.fastq
Approx 10% complete for SRR2584866.fastq
Approx 15% complete for SRR2584866.fastq
Approx 20% complete for SRR2584866.fastq
Approx 25% complete for SRR2584866.fastq
. 
. 
. 
~~~
{: .output}


For each of your sample files, FastQC will ask if you want to replace the existing version with a new version. This is 
because we have already run FastQC on this samples files and generated all of the outputs. We are now doing this again using
our scripts. Go ahead and select `A` each time this message appears. It will appear once per sample file (six times total).

~~~
replace SRR2584866_fastqc/Icons/fastqc_icon.png? [y]es, [n]o, [A]ll, [N]one, [r]ename:
~~~
{: .output}


# Automating the Rest of our Variant Calling Workflow

=== Parking here for future editing ===

>> After trimming and filtering, our overall quality is much higher, 
>> we have a distribution of sequence lengths, and more samples pass 
>> adapter content. However, quality trimming is not perfect, and some
>> programs are better at removing some sequences than others. Because our
>> sequences still contain 3' adapters, it could be important to explore
>> other trimming tools like [cutadapt](http://cutadapt.readthedocs.io/en/stable/) to remove these, depending on your
>> downstream application. Trimmomatic did pretty well though, and its performance
>> is good enough for our workflow.
=== 

We can extend these principles to the entire variant calling workflow. To do this, we will take all of the individual commands that we wrote before, put them into a single file, add variables so that the script knows to iterate through our input files and write to the appropriate output files. This is very similar to what we did with our `read_qc.sh` script, but will be a bit more complex.

Download the script from [here](https://raw.githubusercontent.com/datacarpentry/wrangling-genomics/gh-pages/files/run_variant_calling.sh). Download to `~/course/scripts`.

~~~
curl -O https://raw.githubusercontent.com/datacarpentry/wrangling-genomics/gh-pages/files/run_variant_calling.sh
~~~
{: .bash}

Our variant calling workflow has the following steps:

1. Index the reference genome for use by bwa and samtools.
2. Align reads to reference genome.
3. Convert the format of the alignment to sorted BAM, with some intermediate steps.
4. Calculate the read coverage of positions in the genome.
5. Detect the single nucleotide polymorphisms (SNPs).
6. Filter and report the SNP variants in VCF (variant calling format).

Let's go through this script together:

~~~
$ cd ~/course/scripts
$ less run_variant_calling.sh
~~~
{: .bash}

The script should look like this:

~~~
set -e
cd ~/dc_workshop/results

genome=~/course/data/ref_genome/ecoli_rel606.fasta

bwa index $genome

mkdir -p sam bam bcf vcf

for fq1 in ~/course/data/trimmed_fastq_small/*_1.trim.sub.fastq
    do
    echo "working with file $fq1"

    base=$(basename $fq1 _1.trim.sub.fastq)
    echo "base name is $base"

    fq1=~/dc_workshop/data/trimmed_fastq_small/${base}_1.trim.sub.fastq
    fq2=~/dc_workshop/data/trimmed_fastq_small/${base}_2.trim.sub.fastq
    sam=~/dc_workshop/results/sam/${base}.aligned.sam
    bam=~/dc_workshop/results/bam/${base}.aligned.bam
    sorted_bam=~/dc_workshop/results/bam/${base}.aligned.sorted.bam
    raw_bcf=~/dc_workshop/results/bcf/${base}_raw.bcf
    variants=~/dc_workshop/results/bcf/${base}_variants.vcf
    final_variants=~/dc_workshop/results/vcf/${base}_final_variants.vcf 

    bwa mem $genome $fq1 $fq2 > $sam
    samtools view -S -b $sam > $bam
    samtools sort -o $sorted_bam $bam
    samtools index $sorted_bam
    bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam
    bcftools call --ploidy 1 -m -v -o $variants $raw_bcf 
    vcfutils.pl varFilter $variants > $final_variants
   
    done
~~~
{: .output}

Now, we'll go through each line in the script before running it.

First, notice that we change our working directory so that we can create new results subdirectories
in the right location. 

~~~
cd ~/course/results
~~~
{: .output}

Next we tell our script where to find the reference genome by assigning the `genome` variable to 
the path to our reference genome: 

~~~
genome=~/course/data/ref_genome/ecoli_rel606.fasta
~~~
{: .output}

Next we index our reference genome for BWA: 

~~~
bwa index $genome
~~~
{: .output}

And create the directory structure to store our results in: 

~~~
mkdir -p sam bam bcf vcf
~~~
{: .output}

Then, we use a loop to run the variant calling workflow on each of our FASTQ files. The full list of commands
within the loop will be executed once for each of the FASTQ files in the 
`data/trimmed_fastq_small/` directory. 
We will include a few `echo` statements to give us status updates on our progress.

The first thing we do is assign the name of the FASTQ file we're currently working with to a variable called `fq1` and
tell the script to `echo` the filename back to us so we can check which file we're on.

~~~
for fq1 in ~/dc_workshop/data/trimmed_fastq_small/*_1.trim.sub.fastq
    do
    echo "working with file $fq1"
~~~
{: .bash}

We then extract the base name of the file (excluding the path and `.fastq` extension) and assign it
to a new variable called `base`. 
~~~
    base=$(basename $fq1 _1.trim.sub.fastq)
    echo "base name is $base"
~~~
{: .bash}

We can use the `base` variable to access both the `base_1.fastq` and `base_2.fastq` input files, and create variables to store the names of our output files. This makes the script easier to read because we don't need to type out the full name of each of the files: instead, we use the `base` variable, but add a different extension (e.g. `.sam`, `.bam`) for each file produced by our workflow.


~~~
    #input fastq files
    fq1=~/course/data/trimmed_fastq_small/${base}_1.trim.sub.fastq
    fq2=~/course/data/trimmed_fastq_small/${base}_2.trim.sub.fastq
    
    # output files
    sam=~/dc_workshop/results/sam/${base}.aligned.sam
    bam=~/dc_workshop/results/bam/${base}.aligned.bam
    sorted_bam=~/course/results/bam/${base}.aligned.sorted.bam
    raw_bcf=~/course/results/bcf/${base}_raw.bcf
    variants=~/course/results/bcf/${base}_variants.vcf
    final_variants=~/course/results/vcf/${base}_final_variants.vcf     
~~~
{: .bash}


And finally, the actual workflow steps:

1) align the reads to the reference genome and output a `.sam` file:

~~~
    bwa mem $genome $fq1 $fq2 > $sam
~~~
{: .output}

2) convert the SAM file to BAM format:

~~~
    samtools view -S -b $sam > $bam
~~~
{: .output}

3) sort the BAM file:

~~~
    samtools sort -o $sorted_bam $bam 
~~~
{: .output}

4) index the BAM file for display purposes:

~~~
    samtools index $sorted_bam
~~~
{: .output}

5) calculate the read coverage of positions in the genome:

~~~
    bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam 
~~~
{: .output}

6) call SNPs with bcftools:

~~~
    bcftools call --ploidy 1 -m -v -o $variants $raw_bcf 
~~~
{: .output}

7) filter and report the SNP variants in variant calling format (VCF):

~~~
    vcfutils.pl varFilter $variants  > $final_variants
~~~
{: .output}



> ## Exercise
> It's a good idea to add comments to your code so that you (or a collaborator) can make sense of what you did later. 
> Look through your existing script. Discuss with a neighbor where you should add comments. Add comments (anything following
> a `#` character will be interpreted as a comment, bash will not try to run these comments as code). 
{: .challenge}


Now we can run our script:

~~~
$ bash run_variant_calling.sh
~~~
{: .bash}


> ## Exercise
>
> The samples we just performed variant calling on are part of the long-term evolution experiment introduced at the 
> beginning of our variant calling workflow. From the metadata table, we know that SRR2589044 was from generation 5000,
> SRR2584863 was from generation 15000, and SRR2584866 was from generation 50000. How did the number of mutations per sample change
> over time? Examine the metadata table. What is one reason the number of mutations may have changed the way they did?
> 
> Hint: You can find a copy of the output files for the subsampled trimmed FASTQ file variant calling in the 
> `~/.solutions/wrangling-solutions/variant_calling_auto/` directory.
> 
>> ## Solution
>> 
>> ~~~
>> $ for infile in ~/course/results/vcf/*_final_variants.vcf
>> > do
>> >     echo ${infile}
>> >     grep -v "#" ${infile} | wc -l
>> > done
>> ~~~
>> {: .bash}
>> 
>> For SRR2589044 from generation 5000 there were 10 mutations, for SRR2584863 from generation 15000 there were 25 mutations, 
>> and SRR2584866 from generation 766 mutations. In the last generation, a hypermutable phenotype had evolved, causing this
>> strain to have more mutations. 
> {: .solution}
{: .challenge}


> ## Bonus Exercise
> 
> If you have time after completing the previous exercise, use `run_variant_calling.sh` to run the variant calling pipeline 
> on the full-sized trimmed FASTQ files. You should have a copy of these already in `~/dc_workshop/data/trimmed_fastq`, but if 
> you don't, there is a copy in `~/.solutions/wrangling-solutions/trimmed_fastq`. Does the number of variants change per sample?
{: .challenge} 



