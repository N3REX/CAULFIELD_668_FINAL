# CAULFIELD_668_FINAL

My final project is the usage of Anvi’o in order to visualize and try to track changes in select contigs across metagenome time series. These metagenomic paired end reads were taken from a Crohn's disease patient before, during and on the downtrend of a disease flareup. This is known via the level of calprotectin in the fecal samples, which is indicative of how much neutrophils are in the gut at the time (i.e. the immune response). Our time points in order are represented by the file names **18894** (before), **19027** (during), and **19029** (after).

The contigs we're interested in are viruses pulled from anaerobic bateria that Cole was able to incubate. Frankly since things are still being figured out outside of this class final project, I and the other people in my lab aren't completely sure of the big picture yet, including these lil' guys.

For this pipeline/workflow you'll need two things:

Your contig(s) of interest: in this case ```sussytigs.fasta```

Your reads: ```X_R1.fastq``` and ```X_R2.fastq``` 


### Docker Installation/opening
You'll be working inside the directory where you park your datafiles (make sure you have space since these take up ~150GB alone, and some later generated files are also huge). You do not need to use Docker like I am in this walkthrough/pseudo-tutorial , but **if you're using an Apple computer with an M1/M2 chip, you will not be able to download all the needed dependencies**. You can download/install Anvi'o [here](https://anvio.org/install/)- the official spelling/use of anvi'o is with a lowercase "a" which throws me off, so I'll be ignoring it. The program is pretty well made and has a [pretty helpful github page where people Q&A for issues they get](https://github.com/merenlab/anvio/issues/). My only issue with the program is that it occasionally uses late 2000s/early 2010s internet lingo in it's readouts which can be very cringe and infuriating when you've been trying to troubleshoot a command that isn't working and you're working off of  ```':3 Oopsies, I has an error```.

Docker is pretty neat since it's essentially a VM which doesn't emulate hardware, and for a lot of programs like this and QIIME, you can just download a container that has everything you need.
```
docker pull meren/anvio:7

docker run --rm -it -v `pwd`:`pwd` -w `pwd` -p 8080:8080 meren/anvio:7

docker restart c19224c72b4b 
docker exec -ti tender_northcutt sh
```

### Contigs
```
anvi-script-reformat-fasta sussytigs.fasta -o contig-fixed.fa -l 0 

mv contig-fixed.fa contigs.fa

anvi-gen-contigs-database -f contigs.fa -o contigs.db -n 'Suspect contigs database'
```
*Had to reduce contig.fasta file headers to single word

### Prepping the metagenome
```
bowtie2-build contofint.fa contofint

bowtie2 -x contofint -1 18894_nh_1.fastq -2 18894_nh_r2.fastq -S 18894.sam 

samtools view -b -o 18894.bam 18894.sam 
samtools sort -o 18894_sort.bam 18894.bam 
samtools index 18894_sort.bam

```
### The Whole Enchilada
Now onto the main course, Anvi'o itself. We'll use that ```_.bam``` file from earlier, and the ```contigs.db``` from the start to create a profile. The ```--cluster-contigs``` flag is "optional" in that you can either specify it here as I'm doing or use another flag on the actual ```anvi-interactive``` to force it. It's required to make the phylogenetic tree bit of your visualization so if you don't do either of these then you'll get an error.
```

anvi-profile -i 18894_sort.bam -c contigs.db -o 18894_sussytig --cluster-contigs
```
The next command will reuse the contigs database file you make along with the ```PROFILE.db``` file *within* the profile output from the previous line. If you don't specify this file and instead just call the profile folder (in the case below just ```19029_sussytig```), you will get an another error.
```
anvi-interactive -p 18894_sussytig/PROFILE.db -c contigs.db 
```


### Repeating for the other two time points

```
bowtie2 -x contofint -1 19027_nh.1.fastq -2 19027_nh.2.fastq -S 19027.sam 
samtools view -b -o 19027.bam 19027.sam 
samtools sort -o 19027_sort.bam 19027.bam 
samtools index 19027_sort.bam
anvi-profile -i 19027_sort.bam -c contigs.db -o 19027_sussytig --cluster-contigs
anvi-interactive -p 19027_sussytig/PROFILE.db -c contigs.db 
```

```
bowtie2 -x contofint -1 19029_NH_R1.fastq -2 19029_NH_R2.fastq -S 19029.sam
samtools view -b -o 19029.bam 19029.sam 
samtools sort -o 19029_sort.bam 19029.bam 
samtools index 19029_sort.bam
anvi-profile -i 19029_sort.bam -c contigs.db -o 19029_sussytig --cluster-contigs
anvi-interactive -p 19029_sussytig/PROFILE.db -c contigs.db

```

