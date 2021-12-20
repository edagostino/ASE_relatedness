ASE relatedness project

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
done < plate_well_SamplesHeadCTRL.txt
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

I am also interested in sites that are biallelic (since invariant sites aren't informative, and multiallelic sites are most likely due to genotyping or barcoding errors) and that are SNPs. Again, using `bcftools view` (and, in retrospect I should have combined this with the above step):

```bash
module load bcftools
bcftools view -m3 -M3 -e 'INFO/INDEL=1' -Oz -o head_CDS_biallelic_SNPs.vcf.gz head_CDS.vcf.gz
```

Note that `-m3 -M3` is different from the standard `-m2 -M2` used for finding biallelic sites from the output of `bcftools call` . This is because, without callling variants, bcftools doesn't assign an alternate allele definitively, instead listing alternate alleles as (e.g.) "A, <\*>". <\*> represents an alternate allele not known (see section 5.5 of the VCF 4.3 specification). Therefore, we want sites with three alleles: the reference allele, the "unofficial" alternate allele, and the <\*> (hence `-m3 -M3`). The `-e 'INFO/INDEL=1'` flag drops indels as well.

##### Editing VCF to make ASE table

We now have the VCF that we'll conduct our downstream analyses on, but we want to createa a .txt file that's clearer and easier to manipulate. To do this, we first create a header with `bcftools query -l`, by printing five headers we want as well as the list of samples. We use `sed` to delete the file paths for clarity, and then use `tr` to tab-delimit our list. 

```bash
module load bcftools
printf "CHROM\nPOS\nREF\nALT\nDP\n$(bcftools query -l head_CDS_biallelic_SNPs.vcf.gz | sed s'|\./head_bams/||')" | tr "\n" "\t" > final_head_header.txt
```

Next, we create the body of our text file with`bcftools query -f`, taking the columns we mentioned in our previous command, and printing allele depths at each sample. Note allele depths are represented as, say, `5,1,0` (5 reads mapping to the reference allele, one mapping to the primary inferred alternate allele, and zero mapping to <*>). While the first two are variable, the last number is always zero, since if there were a mapping to an alternate base besides the reference or primary alternate, it wouldn't be a biallelic SNP, and we're only using those.

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP[\t%SAMPLE=%AD]\n' head_CDS_biallelic_SNPs.vcf.gz | sed 's|\./head_bams/[^=]*=||g' | sed 's|,0\t|\t|g' | sed 's|,0\n|\n|g' > final_head_rows.txt
```

For clarity, we use three `sed`s . The first one deletes the sample names at each site, since we have those from the header we made above, and the second and third one delete the always-zero mappings to <*>. We do this twice because sometimes that zero is followed by a tab (when it's within a line) and sometimes it's followed by a newline character (at the end of a line).

Last, we just concatenate our header to this. Yay!

```bash
cat final_head_header.txt final_head_rows.txt > final_head_all.txt
```

##### Redoing the filtering

After conversations with Julien and Luisa, I redid some filtering. First, I redownloaded Luisa's gene_pos_dmel.txt. To manipulate it, I want to delete the first row, rename the chromosomes (i.e., "2L" instead of "chr2L"), and only use regions on the main chromosomes.

```bash
cat gene_pos_dmel.txt | tail -n +2 | sed 's/chr//' | awk '{ if (($2 == "X") || ($2 == "2L") || ($2 == "2R") || ($2 == "3L") || ($2 == "3R") || ($2 == "4")) {print $2 "\t" $3 "\t" $4}}' > formatted_gene_pos_dmel.txt
```

I used a series of bcftools view commands to subset the VCF. First, I used the newly-formatted table to filter to only the gene regions. Since we're working with expression data, everything should map, so this is more of a sanity check, which it passed: the VCF before this filtering step was 33GB, and after, it was 32G: I'm guessing just a few .bams may have mapped to the wrong place. I will redo this to annotate with FlyBase gene names instead of just filter by their regions, which will be more informative, too.

```bash
bcftools view -R ./gene_positions/genic_regions_for_view.txt -Oz -o head_genic_regions.vcf.gz head_all.vcf.gz
```

Next, I filtered to only use SNPs. As before, note that `-m3 -M3` is different from the standard `-m2 -M2` used for finding biallelic sites from the output of `bcftools call` . This is because, without callling variants, bcftools doesn't assign an alternate allele definitively, instead listing alternate alleles as (e.g.) "A, <\*>". <\*> represents an alternate allele not known (see section 5.5 of the VCF 4.3 specification). Therefore, we want sites with three alleles: the reference allele, the "unofficial" alternate allele, and the <\*> (hence `-m3 -M3`). The `-e 'INFO/INDEL=1'` flag drops indels as well.

I should have combined these bcftools views to be more efficient, but I did them separately because I was working through what I needed at the time.

```bash
bcftools view -m3 -M3 -e 'INFO/INDEL=1' -Oz -o head_genic_regions_biallelic_SNPs.vcf.gz head_genic_regions.vcf.gz
```

Finally, I used bcftools view to do some filtering. `N_PASS` is quite cool, and relatively new to bcftools. Basically, it lets you specify the number of samples needed to pass an expression check. I told it that the expression check I needed those samples to pass was `FORMAT/AD[:0]` (the number of alleles mapping to the reference) plus`FORMAT/AD[:1]` (the number of alleles mapping to the first alternate, which is really the only alternate for SNPs) had to be more than 9, and I told it that I need more than 49 samples to pass that expression check for a site to include that site. I also needed the total depth at that site, across samples, to be between 5k and 100k.

```bash
bcftools view -i'N_PASS(FORMAT/AD[:0]+FORMAT/AD[:1]>9)>49 && INFO/DP>5000 && INFO/DP<100000' -Oz -o final_head_all_filt.vcf.gz head_genic_regions_biallelic_SNPs.vcf.gz
```

