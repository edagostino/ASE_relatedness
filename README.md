# ASE_relatedness

I used the following to generate allele-specific expression data at variant coding sites in the *Drosophila melanogaster* genome. Note that this shows the steps for the head expression data; the body expression data was done in parallel, identically. Also, note that some of these were submitted as SLURM jobs, some were submitted via SLURM wrappers (left in the code here), and some were just run on the login node.

##### Copying over the files

First, I copied over the directory Luisa made to a location of my own:

```bash
cp -r /Genomics/ayroleslab2/lamaya/bigProject/GeneExpression /Genomics/ayroleslab2/emmanuel/relatedness_ase
```

I also copied over the genome to a new directory I added to my copy of`GeneExpression`:

```bash
cp /Genomics/grid/users/lamaya/genomes/dmel_genome/dmel-all-chromosome-r6.14.fa /Genomics/ayroleslab2/emmanuel/relatedness_ase/GeneExpression/genome
```

##### Subsetting and indexing the `.bam`s

I then took the list of the subset of "good" sample bams, `SamplesHeadCTRL.txt`, and the list of file paths, `allRNAbamfiles_CTRLhead_bigproject_noReseq_Nov1521`, and made a new list of file paths, `good_head_bams_list.txt`, to just the good samples. The `grep -e` option `grep`s files by either of the two ways they could be named in the `.bam`s.

```bash
#We just want the "plate_well" which is just the second column of SamplesHeadCTRL.txt
awk -F'[ ]' '{print $2}' < SamplesHeadCTRL.txt > plate_well_SamplesHeadCTRL.txt

#Read in this new list using the _ as the deliminter
while IFS='_' read -r plate well

do
#Grep ones that match it either way with -e
#Note that we need to put the slashes in front, and either . or _ (depending on the pattern) after
grep -e "/${plate}h_RNA_${well}\." -e "/${plate}-RNA-${well}_" < allRNAbamfiles_CTRLhead_bigproject_noReseq_Nov1521 > good_head_bam_temp.txt
#Append each one to a master list
cat good_head_bam_temp.txt >> good_head_bams_list.txt
#rm temp
rm good_head_bam_temp.txt
done < plate_well_SamplesHeadCTRL.txt```
```

Next, I copied the "good" `.bam`s to a new directory I made, `head_bams` , so as not to mess up Luisa's files.

```bash
for file_with_path in `cat good_head_bams_list.txt`
do
cp $file_with_path ./head_bams
done
```

I also generated a text file, `head_bam_list_for_bcftools`, that's just a list of the`.bam`s  preceded by "`./head_bams/`" (i.e., the paths to the `.bam`s from my main `GeneExpression` working directory.) I indexed these bams, too:

``` bash
for bam in `cd head_bams; ls *bam`
do
sbatch --wrap "
cd head_bams
module load samtools
samtools index $bam"
done
```

##### Using `bcftools mpileup` to count allele-specific expression

From that directory, I used `bcftools mpileup` to get counts of allele-specific reads at each site (the `-a AD` flag does this perfectly). To parallelize this, I used a wrapper to `mpileup` separately for each chomosome, and get six outfiles. Also, note, at Julien's advice, I'm discarding anything that maps to a scaffold that's not one of the six Muller elements, and I kept all of Luisa's quality filters (except the SNP one, since I'm not interested in a predefined list of specific sites).

```bash
for chr in X 2L 2R 3L 3R 4
do
sbatch --time=7-0 --wrap "
module load bcftools/1.9
bcftools mpileup -Oz -B -q 60 -d 1000 -r ${chr} -a AD -f /Genomics/ayroleslab2/emmanuel/relatedness_ase/GeneExpression/genome/dmel-all-chromosome-r6.14.fa -b head_bam_list_for_bcftools -o ${chr}.headbam.vcf.gz"
done
```

I also indexed the chromosome VCFs I generated.

```bash
for vcf in X.headbam.vcf.gz 2L.headbam.vcf.gz 2R.headbam.vcf.gz 3L.headbam.vcf.gz 3R.headbam.vcf.gz 4.headbam.vcf.gz
do
sbatch --wrap "
module load bcftools
bcftools index ${vcf}"
done
```

I concatenated these into one VCF and indexed it.

```bash
module load bcftools
bcftools concat -Oz -o head_all.vcf.gz X.headbam.vcf.gz 2L.headbam.vcf.gz 2R.headbam.vcf.gz 3L.headbam.vcf.gz 3R.headbam.vcf.gz 4.headbam.vcf.gz
bcftools index head_all.vcf.gz
```

##### Subsetting to just variants in coding regions

We're only really interested in coding regions, since we're working with expression `.bam`s, and since we only care about allele-specific expression of genes. I copied the `.gff` file from FlyBase (making sure to use the same release as the reference genome Luisa provided). https://ftp.flybase.org/genomes/dmel/dmel_r6.14_FB2017_01/gff/dmel-all-r6.14.gff.gz. I also made a "`gff`" directory in my main `GeneExpression ` directory, and put it there.

In that sub-directory, I used `awk` to make a version of the `.gff` that 1) only contains annotations of the chromosomes (as opposed to the other scaffolds), 2) only contains the coding sequences, and 3) only contains the chromosome, start position, and end position of each coding sequence (so, like a `.bed` file, but not zero-indexed).

```bash
zcat dmel-all-r6.14.gff.gz | awk '{ if ((($1 == "X") || ($1 == "2L") || ($1 == "2R") || ($1 == "3L") || ($1 == "3R") || ($1 == "4")) && ($3 == "CDS")) { print $1 "\t" $4 "\t" $5} }' > dmel-CDS-chr-regions.txt
```

Now that I have this regions file, I used `bcftools view` to subset my VCF from the last section to just the regions it contains.

```bash
module load bcftools
bcftools view -R ./gff/dmel-CDS-chr-regions.txt -Oz -o head_CDS.vcf.gz head_all.vcf.gz
```

