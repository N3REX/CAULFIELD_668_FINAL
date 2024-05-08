# CAULFIELD_668_FINAL

My final project is the usage of Anvi’o in order to visualize and try to track changes in select contigs across metagenome time series. These metagenomic paired end reads were taken from a Crohn's disease patient before, during and on the downtrend of a flareup. This is known to us by the level of calprotectin in the fecal samples, which is indicative of how much neutrophils are in the gut at the time (i.e. the immune response). Our time points in order are represented by the file names **18894** (before), **19027** (during), and **19029** (after).

The contigs we're interested in are viruses pulled from anaerobic bateria that Cole was able to incubate. Frankly since things are still being figured out outside of this class final project, I and the other people in my lab aren't completely sure of the big picture yet, including these lil' guys.

For this pipeline/workflow you'll need two things:

Your contig(s) of interest: in this case ```sussytigs.fasta``` 

Your reads: ```X_R1.fastq``` and ```X_R2.fastq``` 



## On Anvi'o and Docker
You'll be working inside the directory where you park your datafiles (make sure you have space since these take up ~150GB alone, and some later generated files are also huge). You do not need to use Docker like I am in this walkthrough/pseudo-tutorial , but **if you're using an Apple computer with an M1/M2 chip, you will not be able to download all the needed dependencies**. You can download/install Anvi'o [here](https://anvio.org/install/)- the official spelling/use of anvi'o is with a lowercase "a" which throws me off, so I'll be ignoring it. 
The program is well made and has a [pretty helpful github page where people Q&A for issues they find](https://github.com/merenlab/anvio/issues/). My only issue with the program is that it occasionally uses late 2000s/early 2010s internet lingo in it's readouts which can be very cringy and infuriating when you've been trying to troubleshoot a command that isn't working and you're working off of  ```':3 Oopsies, I has an error```.

Docker is pretty neat since it's essentially a VM which doesn't emulate hardware, and for a lot of programs like this and QIIME, you can just download a container that has everything you need.

### Docker Installation
You'll first need to [install docker](https://www.docker.com/get-started/), which is free (as long as you aren't making like 20 million USD a year or something). Once you have it installed, and have it running, go to your command line. This first command will pull the container onto your local machine, and will eat up 15GB of space, which is chump-change for what we're doing.
```
docker pull meren/anvio:7
```
Then you can run the container- this script will run it in the current working directory/
```
docker run --rm -it -v `pwd`:`pwd` -w `pwd` -p 8080:8080 meren/anvio:7
```
If you need to restart/reopen the docker container in your terminal, you can navigate to the docker desktop app and use either of the following commands (the top one you'd use the container's ID code, and the bottom the weird name that gets generated for the container.
```
docker restart c19224c72a3a 
docker exec -ti tender_northcutt sh
```

## Contigs

First we need to do a little prepwork with our contigs of interest. This first command will essentially clean up our fasta files into something we can work with (unique identifiers and the removal of any weird things that my interfere). Do note that complex headers ***are*** something that I ran into issues with, something I'll touch on again.
```
anvi-script-reformat-fasta sussytigs.fasta -o contig-fixed.fa -l 0 
```
This is an optional step I kept in so I could easily grab a ```contigs.fa``` and know I had the file which is transformed into ```contigs.db```.
```
mv contig-fixed.fa contigs.fa
```
Here's where we make our ```contigs.db``` file which we'll be using later to make our anvi-profile!
```
anvi-gen-contigs-database -f contigs.fa -o contigs.db -n 'Suspect contigs database'
```


## Prepping the metagenome

Now we move onto the bulk of the pipeline- these steps will take a while if you're not on a nice setup. The first step is making an index using our contigs to allign our reads to. This ```contofint.fa``` file is essentially ```contigs.fa``` but but with the headers cut down; for example, ```>phi119.2contig2 length=13894nt depth=0.29x``` was shifted to just ```>phi119.2contig2```. This was done to satisfy some header issue that was popping up with Bowtie2.

```
bowtie2-build contofint.fa contofint
```
Here we map our forward and reverse reads to the index, producing a .sam file. This process took the longest for me to run. 
```
bowtie2 -x contofint -1 18894_nh_1.fastq -2 18894_nh_r2.fastq -S 18894.sam
```
Conversion of the sam to bam, then sorting/indexing the bam. Though the sorting step here can take a while, after this we are good to move onto the Anvi'o profile!
```
samtools view -b -o 18894.bam 18894.sam 
samtools sort -o 18894_sort.bam 18894.bam 
samtools index 18894_sort.bam
```


## The Whole Enchilada
Now onto the main course, Anvi'o itself. We'll use that ```_.bam``` file from earlier, and the ```contigs.db``` from the start to create a profile. The ```--cluster-contigs``` flag is "optional" in that you can either specify it here as I'm doing or use another flag on the actual ```anvi-interactive``` to force it. It's required to make the phylogenetic tree bit of your visualization so if you don't do either of these then you'll get an error.
```
anvi-profile -i 18894_sort.bam -c contigs.db -o 18894_sussytig --cluster-contigs
```
The next command will reuse the contigs database file you make along with the ```PROFILE.db``` file *within* the profile output from the previous line. If you don't specify this file and instead just call the profile folder (in the case below just ```19029_sussytig```), you will get an another error.
```
anvi-interactive -p 18894_sussytig/PROFILE.db -c contigs.db 
```



### Repeating for the other two time points

During the flare-up. We get to reuse ```contigs.db``` but still need to make new sam/bam and profiles.
```
bowtie2 -x contofint -1 19027_nh.1.fastq -2 19027_nh.2.fastq -S 19027.sam 
samtools view -b -o 19027.bam 19027.sam 
samtools sort -o 19027_sort.bam 19027.bam 
samtools index 19027_sort.bam
anvi-profile -i 19027_sort.bam -c contigs.db -o 19027_sussytig --cluster-contigs
anvi-interactive -p 19027_sussytig/PROFILE.db -c contigs.db 
```
After the flare-up. Same deal as the 19027.
```
bowtie2 -x contofint -1 19029_NH_R1.fastq -2 19029_NH_R2.fastq -S 19029.sam
samtools view -b -o 19029.bam 19029.sam 
samtools sort -o 19029_sort.bam 19029.bam 
samtools index 19029_sort.bam
anvi-profile -i 19029_sort.bam -c contigs.db -o 19029_sussytig --cluster-contigs
anvi-interactive -p 19029_sussytig/PROFILE.db -c contigs.db
```
## So what does it mean?
When you open up the server link provided to you, you should see an interface with a multi-layered ring. At this point its completely up to personal preference how you tweak the parameters to change the visualization. On the right side of the page, having the mouse tab open will let you see various values when you hover over them on the figure.  You **can** export your figure as an .svg, or save it's state so you can come back to it later, but given that the bin state saving doesn't seem to work consistently, I wouldn't reccomend relying on this feature. If you import more data you can show stuff in that quarter of the circle that is missing, but I opted to increase the radial coverage to fully display all the data labels. To share details with my labmates more easially I used some third party illustration program to copy paste the mouse output next to each contig, as you can see here.
![Modified screenshot of 18894's anvi'o with the mouse tab for each contig displayed.](https://github.com/N3REX/CAULFIELD_668_FINAL/blob/main/anvint_18894_v1.png)
![Modified screenshot of 19027's anvi'o with the mouse tab for each contig displayed.](https://github.com/N3REX/CAULFIELD_668_FINAL/blob/main/anvint_19027_v1.png)
![Modified screenshot of 19029's anvi'o with the mouse tab for each contig displayed.](https://github.com/N3REX/CAULFIELD_668_FINAL/blob/main/anvint_19029_v1.png)

### Things of interest
In the first time point we see Phi119.2 Contig 2 being the most detected contig with a mere 23%- this number would jump up along with the other contigs during the flare up, with Phi29E2.1 Contig 3 having as high as 48%. These numbers are pretty great for viral genes, but the real stand-out result we see is in 19029 *after the flareup*. Phi29E1.1 Contig 2, which up until after the flare-up had a detection of zero, explodes to a detection of **98%**. Its something incredibly interesting and is something I plan on looking into further to explore why we see this really outlying result. So far I've done a little conserved domain searching and found that this contig has an unusual N<sub>2</sub>0 reductase, but I'm still looking into it.
