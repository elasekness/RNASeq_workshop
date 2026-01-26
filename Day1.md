## Tutorial 1

**Objective:** Successfully log on to Google Cloud Platform (GCP) Console and become familiar with working in a command-line (Linux) environment. Use various tools to download and manipulate data files.

<br>

## Opening a terminal window

We are working on remote servers or virtual machines (VMs) hosted by Google.
To open a terminal window where we will work, we need to log on to the GCP console with our Wadsworth credentials,
start the VM associated with our Project and `ssh` to those servers.
`ssh` is a secure shell protocol for safely connecting to remote services.
To `ssh` you need an account on the server with a username and password as well as the IP address of the server.

- Navigate to the [GCP console](https://console.cloud.google.com) and log on with your Wadsworth credentials
- Once connected to the console, click the Compute Engine link and start your VM by expanding the three dot icon to the right
- Click the SSH button once the VM has started. This will open a terminal window on your VM.

<br>

To perform any operations in a Linux environment, we need to tell the computer what to do by typing specific commands into our terminal window.
In general, the syntax will be: `Command File`,
where `Command` is the function or operation you want to perform on the` File` or `Directory`
that you specify. Sometimes specifying a command is all we need.

<br>

Examine the contents of your directory with:
 
	ls

 > **`ls`** = list

<br>


## Directory structure and navigation

Directories are folders with specific locations on the server.
To navigate from one directory to another, we must specify the location of the directory or its path.
Directories have a forward slash `/` after their names.  Thus to get to a subdirectory within a directory
you would specify the path by stringing together the directory names, separated by `/`. We'll practice below

<br>

First, determine the name and location of your home directory.

	echo $HOME

> **`echo`** = repeat, **`$HOME`** = variable name for your home directory and its location on the server. <br>
> Try the **`echo`** command with other variable names, such as **`$SHELL`** or **`$PATH`**.

<br>

You can also print your working directory (where you are).
 
	pwd

> **`pwd`** = print working directory

<br>

Make a directory within your home directory called `fastq`, where we will store your fastq files.

	mkdir fastq

> **`mkdir`** = make directory

<br>

Now try making two directories at the same time, named `reference` and `rnaseq_workshop`.

<br>

Rename a file or directory. You can rename a file or directory with the `mv` (move) command.

	mv reference reference_db

> **`mv`** = move. <br>
> Renaming files with `mv` will overwrite any existing file.  

<br>

You can also mv a file to a different directory or move one directory into another.
Let's make `reference_db` and `fastq` subdirectories of `rnaseq_workshop`.  

	mv reference_db fastq rnaseq_workshop

> Notice that you can move multiple files or directories at the same time, with the last directory being the final location of those previously listed.

<br>


Change (move) to a different directory.

	cd rnaseq_workshop/fastq

> **`cd`** = change directory.
> Because we are in our home directories and fastq is in a subdirectory of our home directory, we need to specify the path to get there. <br>
> You can think of changing directories as physically moving from one directory to another, which means your point of reference has changed. This will become more evident in the upcoming examples. <br>
> Use **`../`** or **`..`** to move up one directory or **`../../`** to move up two directories (back to your home directory). <br>
> What does **`cd`** alone do? <br>

<br>

Now try copying files from one directory to another.  Here we will copy the fastq files required for our future RNASeq analysis.
We will copy R1 and R2 paired-end read files from two conditions (a control and an experimental treatment), each having two replicates, from `/home/rnaseq_workshop/fastq` to your `fastq` directory. 
To copy them, we can specify an absolute path, which is the location
of these files with respect to the root directory (i.e. going through the entire filesystem to get to your file) or a relative path, which is
where these files are located with respect to your current working directory (i.e. where you are when you enter `pwd`). If it's not, make sure `fastq` is your current working dirctory.

<br>

Let's try using absolute and relative paths.

	cp /home/rnaseq_workshop/fastq/wt-1_R1.fastq.gz .

> **`cp`** = copy. <br>
> Here we used the absolute path to copy the fastq file `wt-1_R1.fastq.gz` to `.`, which represents our current location.
> Since we cd'ed into `fastq` in our previous command, this should be our current location.
> It's always good to check that your commands worked. `ls` your current directory to see that the fastq file copied.

<br>

Now use a relative path to copy the fastq file to where you are. This will simply overwrite the exisiting file.

	cp ../../../../rnaseq_workshop/fastq/wt-1_R1.fastq.gz .

> Notice that we had to move up four directories to get to `rnaseq_workshop`. <br>
> Also notice that we must always specify an end location for our copied files, which is our current location: `.`. 
> If we specified `new-fastq_R1.fastq.gz` instead of `.`, the contents of `wt-1_R1.fastq.gz` would be written to this new file in our current location.

<br>

You could continue to copy each file individually, or we could copy all our fastq files simultaneously by employing metacharacters.

	cp /home/rnaseq_workshop/fastq/*gz .

> The **star/asterisk** has special meaning.  It represents any character any number of times
> Given the syntax above, I am copying any file with a `gz` ending from the fastq subdirectory to my current directory.
> Other metacharacters are described at the end of this tutorial.

<br>

You can also copy the entire rnaseq_workshop directory and all of its contents to your home directory. 
**Don't do this.**

	cp -r /home/rnaseq_workshop/ .
	
> Notice that we need the `-r` recursive argument to specify that we are copying a directory and all the content within that directory, including subdirectories.

## Permissions


Directories and files have specific permissions associated with them or things that you are allowed to do to them.
There are three permissions: 1) the ability to read (r) 2) the ability to write (w) and 3) the ability to execute (a script or program; x).
There are three sets of permissions representing 1) the User (you), 2) the group (multiple users),
and 3) Others (Everyone else in the world).

<br>

Change directories with the `cd rna_workshop` command (if you are not already in this directory) directory and examine the permissions of your 'fastq' directory.

	ls -lh fastq

> **`ls -l`** = long list, providing you information on when the **`fastq`** directory was created and its associated permissions. <br>
> the **`-h`** specifies that the output be in human-readable format.

<br>

Change the permissions associated with your `fastq` directory with the 'chmod' command (change mode).

	chmod 775 fastq

> Permissions are represented by three-digit (for user, group, and other) octal numbers.
> Here we are allowing the user and group universal permissions (7 = read, write, and execute) and all others
> the ability to read and write only (5).
> For more information on permissions, see: [Linux permissions](https://www.guru99.com/file-permissions.html#linux_file_ownership)

<br>


## Manipulating files (making, (re)moving, editing, and viewing)

There are many ways to make and view files in a Linux OS (operating system).
We can redirect output from a command that prints its output to your screen (STDOUT) to a file instead,
we can generate files on the fly by opening them in a text editor, or we can copy an existing file.
Similarly, we can view and edit files in a text editor or we can print their contents to the screen with various command-line options.
As a general rule, it's always good to examine some of the content of your file to ensure you've generated the results you want in the
format you want or that you are using the correct file and format for downstream applications.

<br>

Redirect STDOUT to a file.

	ls /usr/bin > programs.txt

> The path **`/usr/bin`** specifies the location where various SHELL commands are found. When you type a command, **`/usr/bin`** is one of the locations your computer searches to find and execute the command. <br>
> Was **`/usr/bin`** part of your **`$PATH`**? <br>
> Here we are redirecting the STDOUT from the **`ls`** command to a file named **`programs.txt`**. The **`>`** sign is responsible for the redirection.

<br>

Scroll through the contents of your file.

	more programs.txt

> Scroll through line-by-line with the enter key.  Scroll through page-by-page with the space key. <br>
> Do you notice that the file contains some of the commands you have just used? <br>
> Exit with **`control-c`**.

<br>

Display the first ten lines of your file.

	head programs.txt

> **`head`** displays the first ten lines by default but you can specify the number of lines with a flag. <br>
> For example, **`head -200 programs.txt`** will display the first two hundred lines of your file. <br>
> In general, most commands have additional arguments (or flags) associated with them. <br>
> You can see the different usage statements by typing the command with a **`-help`** or **`--help`** option depending on the program.

<br>

Display the last ten lines of your file.

	tail programs.txt

<br>

Print the entire contents of your file to your screen.

	cat programs.txt

> **`cat`** = concatenate.  The **`cat`** command can also join the contents of multiple files together.

<br>


View and edit the contents of a file with a text editor.

	nano programs.txt

> Nano, emacs, vim, and vi are all text editors. <br>
> You can make an empty file on the fly by typing **`nano`** or **`nano filename`**.  This will open a blank text editor screen. <br>
> Save your changes with **`control-o`**. <br>
> To exit the text editor, use **`control-x`**. <br>
> More information on Nano commands can be found here: 
> [Nano](https://www.howtogeek.com/howto/42980/the-beginners-guide-to-nano-the-linux-command-line-text-editor/)

<br>

Remove a file. Let's remove files we don't need.

	rm programs.txt

> **`rm`** = remove.  Remember, a removed file cannot be restored. <br>
> Can you remove a file from a different directory without having to change directories? <br>
> How would you remove a directory?

<br>


## Notes on working in a Linux environment

* Use **autocomplete** for speedier typing and to avoid typos.  Autocomplete will fill in the unique part of a command or file. For example, if I had only one file in my directory that began with a “b,” I could type `b` and then press the tab key to autocomplete the name of the directory.
* Hit the up-arrow key to recall a command you entered previously.
* Everything in Linux is case sensitive.
* Avoid creating file names/directories with spaces. Use underscores, dashes, or periods instead to separate multi-part names.
* Spaces need to be escaped (more on that later).  For example, if you tried to make a directory called “my directory”, you would get two directories, `my` and `directory`.
* Can’t find a command? Try **`which`** to see if the command is in your path – whether the command is in a location that the computer searches for executing commands.

<br>

## Databases and obtaining sequences


There are several sequence databases – NCBI, JGI, EMBL - that you might encounter.
We will focus on NCBI.  Our goals are to download reference genome and reference coding sequence files to which we will map our reads. We can accomplish this several ways.  Here will will download sequences directly to our VMs with the E-Utilities and datasets functions installed on your VM.


Navigate to NCBI’s homepage: [https://www.ncbi.nlm.nih.gov/](https://www.ncbi.nlm.nih.gov/)


> Notice that there are options to submit sequences, download sequences, and even analyze data. <br>
> There is also an option to submit sequences to the BLAST server. <br>
> BLAST is an alignment tool to look for sequences that are similar (a proxy for homology) to your queries. <br>


Choose `Assembly` under the top left pull-down menu (set to `All Databases` by default). We can use BOOLEAN search terms (such as AND, OR, and NOT) to make our search more specific.
For now, simply type `pseudomonas aeruginosa UCBPP-PA14`  How many entries are returned?
If we didn't know the genome or assembly accession, we could use the filter drop-down menu to narrow our search. From the `filters` menu, specify that you want genomes with annotation and those 
that are completely assembled.  This now returns one assembly. Clicking the assembly link associated with this strain brings you to a page with additional information, including various ways to download the genome and associated data.  Notice that there are GenBank and RefSeq versions of this assembly even though it is the same genome.

We can download the genome in fasta file format to our computers using the `download` option at the top left of the page.
However, we would then need to transfer this file from our computers to a GCP storage bucket and from the storage bucket to our VM. 
An easier method to obtain these sequences would be to use one of NCBI's tools for interacting with their databases.

<br>

Return to your VM terminal, navigate to your `reference_db` directory, and type:

	efetch --help

> This brings up a long menu of options for the efetch tool, which can be used to download a variety of data in different formats from NCBI. <br>
> **`esearch`** and **`efetch`** are part of Entrez Direct suite of [Eutilities]([https://www.ncbi.nlm.nih.gov/books/NBK179288/] that allow you to search NCBI databases from your terminal window. <br>
> Relevant arguments are typically the database **`-db`** we want to search, the format **`-format`** of the data, and the **`-id`** of our query, which is the accession of a sequence of interest.

<br>

We can use the genome accession associated with the assembly to download the Pa14 genome.  You can obtain this accession by going back to the assembly page for Pa14 and clicking the link "view GenBank sequences." This takes you to the genome page where the accession listed is as `CP000438`. 

	efetch -db nuccore -id CP000438 -format fasta > pa14.fasta
	
> This fetches the nucleotide sequence for the Pa14 genome from the core nucleotide database in fasta format.
> Remember that STDOUT is output from a command that gets printed to your screen while the **`>`** symbol redirects this output to a file (`pa14.fasta`).

<br>

We need to specify a different output format to download the coding sequences associated with this file.

	efetch -db nuccore -id CP000438 -format fasta_cds_na > pa14_cds.fasta
	
> If we didn't have a genome/nucleotide accession, we could query the nucleotide core database with the **`esearch`** command first and pass that information to **`efetch`** using a pipe.

<br>

An easier and more efficient method to download a genome and its coding sequences is to use NCBI's [Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/getting_started/) command-line tool. A help menu will appear if you type **`datasets`** without any arguments.  
Typing **`datasets download`** will give you additional information regarding the download function. Note the option to download a genome by its assembly accession or taxonomy.  


	datasets download genome accession GCA_000014625.1 --include genome,cds
	
> This command will download the GenBank version of the Pa14 assembly as well as the coding sequences in fasta format.
> Data are saved to a zipped folder entitled `ncbi_dataset` by default.  To inflate its contents, simply type **`unzip ncbi_dataset`**.
> If you **`ls`** the directory structure, you'll see that your files are here: `ncbi_dataset/data/GCA_000014625.1/`.

<br>

Move the genome and coding sequence files to `reference_db`.

	mv ncbi_dataset/data/GCA_000014625.1/*fna .
	
> This command assumes you are in your `reference_db` directory.
> The genome and coding sequence fasta files end with `.fna`, making it easy to move both.

As you can see, there are usually multiple ways to solve a problem in bioinformatics.

<br>


## Manipulating data: Parsing files, modifying content, and piping


Often, we want to extract information from and/or alter the contents of a file. We might also want to determine some basic features of the file. For example, we might want to know how many sequences are in our fasta file without having to open and scroll through a file. You will probably use the following commands most frequently to parse and manipulate files – grep, sed, cut, paste, sort and uniq. These commands can perform simple routines such as search and replace but combined with regular expressions, these tools are incredibly powerful and efficient.

Piping is specified by **`|`** and simply pipes the STDOUT from one command to another so that you can string multiple operations together on one line. Sometimes the most challenging bioinformatics operations are wrangling your data into the proper format.

<br>


Count how many sequences are in the `cds_from_genomic.fna` file.

	grep -c ">" cds_from_genomic.fna

> **`grep`** = global regular expression print.  Grep searches a file line-by-line for patterns that you specify and returns every line containing that pattern. <br>
> The **`-c`** option counts the number of lines that contain the search pattern instead of returning the lines. <br>
> Try **`grep`** without the **`-c`** argument to see the difference. <br>
> Combined with metacharacters, **`grep`** is a powerful way to search your document for complicated patterns.

<br>

Grab the first 5 header lines from your fasta file with `grep` and a pipe.

	grep  ">" cds_from_genomic.fna | head -5

> Here we are using a pipe, designated by **`|`** to capture the output of grep and pass it to another command (**`head`**).
> Piping is a really useful skill to learn for parsing and processing data more efficiently. <br>
> Note that you can string many pipes together, if necessary. <br>
> As is the case for most operations conducted in Linux, there are multiple ways to do things. <br>
> Use the manual page for grep to find an alternative way to obtain the first five header lines (**`man grep`**).

<br>

Count the number of sequences in the fasta file using a pipe to **`wc`** instead.

	grep ">" cds_from_genomic.fna | wc -l

 > Here we are passing the output of **`grep`** to the word count command, **`wc`**.  The **`-l`** argument specifies that we want to count lines. <br>
 
<br>

You'll notice that the definition lines or sequence names are quite long and include the chromosome accession, the coding sequence accession, the locus tag, the protein description, and more.
For our purposes, we want to define our genes by the locus_tag.  Without the use of regular expressions and the **`sed`** command, modifying this file would be incredibly difficult and tedious.
Use **`sed`** to rename the definition lines of your fasta file so that only the locus tags remain. But first, **`grep`** the definition lines so that you can view your changes more easily.


	grep ">" cds_from_genomic.fna | sed "s/>.*locus_tag=\(PA14_.....\)\].*/>\1/"

> **`sed`** = stream editor.  **`sed`** is essentially a search and replace function. <br>
> Like **`grep`**, we can search for complicated patterns when we use this command with regular expressions. Unlike **`grep`**, we can replace these complicated patterns with another. <br>
> The syntax for the search and replace command is **`'s/search pattern/replacement pattern/'`** where the 's' stands for substitute. <br>
> Although the definition lines are unique, the formatting is consistent, allowing us to search for a general pattern and replace it with a new one. <br>
> In other words, we can use a regular expression to replace all of the patterns displayed in the sequence names, without having to search for each pattern individually, by employing our special metacharacters. <br>
> Like **`grep`**, **`sed`** will search for your pattern line by line. It will make a replacement once (unless you specify otherwise, see the manual page). <br>
> Here we start our search by specifying that we are looking for lines with a `>` followed by any character `.` any number of times `*`.  This is followed by the `locus_tag=` expression. <br>
> To search for a literal period or star, we would have to exit out of or escape the metacharacter with a backslash: **`\.`** and **`\*`**.
> The **`\`** in front of the parentheses give the parentheses a special meaning.  It is a way of grouping information and saving it for recall later. <br>
> With our special parenthetical notation, we are saving the locus tag information, which we notate as **`PA14_.....`**, where the periods represent any character once for the five digits associated with the tag. <br>  
> We could also specify the digits of the locus tag like this: **`[0-9]*`** where brackets specify a range of numbers or characters and the asterisk specifies any number of times. <br>
> The final part of the expression: `\].*` encompasses all information after the locus tag in the definition line. Anything not saved in the parentheses will be replaced by our replacement expression. <br>
> The replacement expression adds the **`>`** back to the definition line and recalls the information grouped in the parenthesis with **`\1`**, where the backslash essentially makes the 1 a special character recalling the stored information.<br>

If we didn't **`grep`** the definitions lines first,  **`sed`** would print your entire document to STDOUT, making it hard to find the changes, which are only in the definition lines.

* How would we check that changes have been made to all of our definition lines and that we aren't missing any anomalies?

Once you know your **`sed`** command has worked to reformat all definition lines, you can actually implement the change.

	sed "s/>.*locus_tag=\(PA14_.....\)\].*/>\1/" cds_from_genomic.fna > pa14_cds.fna
	
> We've written the modified content to a new file but you can overwrite the content of the original file in place: <br>
> **`sed -i "s/>.*locus_tag=\(PA14_.....\)\].*/>\1/" cds_from_genomic.fna > pa14_cds.fna`** <br>
> However, writing in place should only be performed if you are sure the modified content is correct.

<br>


We now have our fastq files and our reference genome and coding sequence files.  We can proceed to read mapping and quantitation.

## More piping

Let's try some more complicated parsing of our data using various Bash commands and pipes. 

## Bash for loops


Bash for loops are basically little shell scripts that can be excecuted from the command line (BASH is the command-line language we are using). Like all loops, they allow you to automate iterative processes. For example, instead of opening 200 hundred fasta files and manually changing the definition lines in each, I can run a for loop that will open each fasta file and make the changes that I specify.

<br>

The basic syntax is:

	for FILE in *common_file_ending; do command $FILE; done

> The interpretation of this code is: <br>
> For every file that ends in some common ending (such as .txt or .gz), perform (do) some command on that file until there are no more files on which to operate, whereby “done” will exit us from the loop. <br>
> The $ in front of FILE indicates that $FILE is a variable, a placeholder which is referring to each file that enters the loop, just as x is a variable that represents 2 in the equation x = 2. <br>
> The **`for`**, **`in`**, **`do`**, and **`done`** are required parts of the for-loop syntax.

<br>

We will try some examples in class.

<br>


## Regular expressions (regex) and special characters (metacharacters)

Regular expressions are search terms that incorporate special characters to make searches more powerful (both broader and more specific).
Metacharacters have a special meaning and include:

- `*` 	Star is a greedy metacharacter meaning match anything any number of times
- `[]` Brackets are often used to specify a range of numbers or letters to include in a search pattern
- `.` 	A period represents any character once
- `?` 	Match one character
- `$`		End of Line
- `^`		Beginning of line

Usually, we escape special characters with a backslash to interpret them literally.

<br>

