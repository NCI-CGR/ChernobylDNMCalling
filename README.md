# ChernobylDNMCalling

### I. Description
This code is used for calling de novo mutations in familial trios. It accepts BAM files for parents and offspring as input and outputs putative de novo mutation calls. The pipeline has been adapted from DNM calling methods developed by [Wong et al](https://www.nature.com/articles/ncomms10486) and the [Epi4K consortium](https://www.nature.com/articles/nature12439).

![image](https://user-images.githubusercontent.com/2903359/110140571-39a1a180-7da2-11eb-9afe-de1d97099f01.png)

### II. Dependencies
* GATK (HaplotypeCaller, PhaseByTransmission, ReadBackedPhasing)
* bcftools
* bedtools
* R (packages: sqldf, stringr, data.table)
* Python
* tabix

Optional (step 4-6):
* [SnpEff](https://pcingola.github.io/SnpEff/)
* liftOver

### III. Reference data
* Reference assembly
* dbSNP VCF
* interval_list for excluding any regions from variant calling
* [RepeatMasker repeating elements](http://genome.ucsc.edu/cgi-bin/hgTrackUi?g=rmsk)

Optional (step 4-6):
* [Known fragile sites in human chromosomes (bed file)](https://webs.iiitd.edu.in/raghava/humcfs/download.html)
* liftOver chain files
* SnpEff database (uses local version after downloading for first use)

### IV. Input
#### params.config
The file needs to be modified to specify working directories, location of reference data and tools, and location of all trio BAM files (assumed to be in a single directory).

#### combined_manifest.txt
A tab-delimited file in the following format needs to be provided. It lists Sample_ID (corresponding BAM files should be named <Sample_ID>.bam and should contain the same sample ID), Subject_ID (identifier for trio members, should be unique within each trio), Family (unique identifier for the trio), and Member for identifying each member within a given trio (can be fa, mo, c1..n for father, mother, and child 1..n)

```
Sample_ID	Subject_ID  Family	Member
SC400001	t0005c2	    t0005	c1
SC400002	t0005c1	    t0005	c2
SC400003	t0005fa	    t0005	fa
SC400004	t0005mo	    t0005	mo
SC400008	t0008c1	    t0008	c1
SC400009	t0008fa	    t0008	fa
SC400006	t0008mo	    t0008	mo
```

### V. Workflow and output
The code is designed to be run on SGE. However, it can be modified to run in other environments.

**Note:** Please source `params.config` before running the following steps.

**Step 1:** `create_manifest.py` uses `combined_manifest.txt` above to create a folder for each trio (separate folders for each offspring) and deposits a manifest.txt file into each trio folder. The input BAM and output data directories for the trio are inherited from `params.config`.

```
outId	t0008c1
child	SC400008
mom	    SC400006
dad	    SC400009
bamDir	/projects/dnm/bams
dataDir	/projects/dnm/trios/t0008c1
```

**Step 2:** `submit_hc_qsub.sh` calls `runHc_b38_qsub.sh` for each trio using parameters defined in the trio's manifest.txt and submits a job to the cluster for variant calling using HaplotypeCaller. GATK Queue was used for scattering and gathering HC jobs at the time this code was developed. However, this variant calling step can modified to use current GATK best practices.

**Step 3:**  `submit_filt_qsub.sh` calls `runPtFilt.sh` for each trio to create a ped file from manifest.txt, run PhaseByTransmission on trio VCF, filter MIEs using `runFilt.r` to extract putative DNMs with different filtering strategies, annotated with repeat information, and create a BED file for viewing the variants more conveniently in IGV. ReadBackedPhasing is then run on PhaseByTransmission output VCF using only the intervals as defined by putative DNM sites. Parent of origin is then determined with `extractparentoforigin.pl` (credit Laurent Francioli) using RbP VCF and DNM sites as the input.

**Step 4** (Optional): `submit_hg19.sh` lifts over b38 DNM coordinates to hg19 for comparison against data from alignment to hg19.

**Step 5** (Optional): `submit_snpeff.sh` annotates putative DNM sites (with the most relaxed filters) using SnpEff and source data from RefSeq and Ensembl.

**Step 6** (Optional): `runTriNuc.sh` extracts tri-nucleotide sequence around DNM sites using output from Step 5 and also annotates known fragile sites (from HumCFS). 

Output files for a single trio (using manifest.txt and BAM files as input):
```
├── hc.log
├── hg19
│   ├── hg19_input.txt
│   ├── hg19_out.txt
│   ├── hg19_unmapped_centromere.anno.txt
│   ├── hg19_unmapped.txt
│   ├── liftOver_cgemsIII_7ac6_ef760.bed
│   ├── liftOver_cgemsIII_7ac6_ef760.bedmapped
│   └── liftOver_cgemsIII_7ac6_ef760.bedunmapped
├── manifest.txt
├── RBP
│   ├── t0008c1.gqdp.filt.pt.vcf
│   ├── t0008c1.gqdp.filt.pt.vcf.idx
│   ├── t0008c1.mdnm.relaxed.dnm.buffer.list
│   ├── t0008c1.mdnm.relaxed.dnm.list
│   ├── t0008c1.pt.rbp.po.out
│   ├── t0008c1.pt.rbp.vcf
│   └── t0008c1.pt.rbp.vcf.idx
├── runHc_b38_qsub.sh.e9200091
├── runHc_b38_qsub.sh.o9200091
├── snpEff
│   ├── epi.bed
│   ├── epi.humcfs.txt
│   ├── epi.list
│   ├── epi.snpeff.ensembl.txt
│   ├── epi.snpeff.refseq.txt
│   ├── epi.trinuc.txt
│   ├── epi.vcf
│   ├── epi.vcf.idx
│   ├── snpEff_genes.txt
│   └── snpEff_summary.html
├── t0008c1.mdnm.relaxed.bed
├── t0008c1.mie.pt.filt.common.txt
├── t0008c1.mie.pt.filt.epi.anno.txt
├── t0008c1.mie.pt.filt.epi.txt
├── t0008c1.mie.pt.filt.homalt.txt
├── t0008c1.mie.pt.filt.mdnm.relaxed.txt
├── t0008c1.mie.pt.filt.mdnm.txt
├── t0008c1.mie.pt.txt
├── t0008c1.ped
├── t0008c1.pt.vcf
├── t0008c1.pt.vcf.idx
├── t0008c1.raw_variants.vcf
├── t0008c1.raw_variants.vcf.idx
├── t0008c1.raw_variants.vcf.out
└── variantCaller.jobReport.txt
```
