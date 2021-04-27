# README.md

In this repository you will find scripts & workflows for the benchmarking of genomic data, particularly variant annotation files (vcfs).
Although all steps can be used generate figures and analyze all genotype information, they primarily serve to interpret results of the ExomePLUS project 
at the McGill Genome Centre. 

With manuscript publication set sometime this year (2021), no figures are available for reference just yet. I will update with the DOI upon release.

The code shown below was originally part of an automated #SLURM array job on ComputeCanada's Beluga server. As such, it relies on some configured modules.
Since they are all open source, they can also be installed locally and the commands can be run with any Bash-like shell.

# Analysis #1: [Metrics by Sample](metrics_by_sample)

- Apply PASS and FAIL filters based on hard-filters recommendations from GATK:

```
module load bcftools
for filepath in <input_dir>/*.vcf.gz
do
filename=$(basename $filepath)
bcftools filter -s FAIL -i '(TYPE="snp" & INFO/QD>=2 & INFO/FS<=60 & INFO/MQ>=40 & INFO/SOR<=3 & \
(INFO/MQRankSum="." | INFO/MQRankSum>=-12.4) & (INFO/ReadPosRankSum="." | INFO/ReadPosRankSum>=-8.0)) | \
(TYPE~"indel" & INFO/QD>=2 & INFO/FS<=200 & INFO/SOR<=10 & (INFO/ReadPosRankSum="." | \
INFO/ReadPosRankSum>=-20))' $filepath -Oz -o <output_dir>/${filename}
done
```
- Run [hap.py](https://github.com/Illumina/hap.py) on vcfs to generate metrics
- Format metrics into desired .csv input with [this](metrics_by_sample/metric_tables.R)
- Use [metric_figure.R](metrics_by_sample/metric_figure.R) for data visualization (creates ggplots for PASS & FAIL)
- Note: If only PASS variants are desired:
```
for filepath in <input_dir>/*.vcf.gz
do
filename=$(basename $filepath)
bcftools view -i "%FILTER!='FAIL'" $filepath -o <output_dir>/$filename
done
```

# Analysis #2: [Metrics by Depth Bins](metrics_by_depth_bins)

**Coverage depth must already be annotated with "DP_BIN"**
•	Aggregate with happy output (TP/FP/FN):
```
module load bcftools
module load tabix
for filepath in <input_dir>*.vcf.gz
do
filename=$(basename $filepath)
tabix $filepath
bcftools annotate -a $filepath -c "DP_BIN" <input_dir>${filename}| gzip -c > <output_dir>/${filename}
done
```
•	Take only PASS:
```
for filepath in <input_dir>/*.vcf.gz
do
filename=$(basename $filepath)
bcftools view -i "%FILTER!='FAIL'" $filepath -o <output_dir>/$filename
done
```
•	Make count tables with [depth_bins_count_tables.R](metrics_by_depth_bins/depth_bins_count_tables.R)
```
module load nixpkgs/16.09  
module load gcc/8.3.0
module load r/4.0.0

for file in <input_dir>/*
do
Rscript depth_bins_count_tables.R $file <output_dir>
done
```

# Analysis #3: [Analyses by Gene](by_gene)
- Extract coding (CDS) and UTR regions based on [GTF annotion](https://www.gencodegenes.org/human/release_34lift37.html) 
of 'basic genes' in the form of BED files
-	Annotate GC content in all vcfs
- Create bam files for each vcf
- Intersect BED with VCFs with findOverlap.R **Note** you will need to do this by chunks, perhaps parallelizing by chromosome
You can then recombine with combine_bychrom.R
- Generating Results
- Outputs: boxplots by average, median, stats_summary, and coverage~gc content 
`Rscript plot_by_gene <INPUT_DIR> <OUTPUT_DIR>`
