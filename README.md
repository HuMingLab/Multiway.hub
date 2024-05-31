# Multi-way Hub Identification

**Last updated: May. 31, 2024**

## Prerequisites

### Packages
* **R** (tested on 4.3.1)
    * **Pairtools** (tested on 1.0.2)
* (other packages required as dependencies)

### Data

* **Sam/Bam files for each single-cell in a directory**

## Installation

Download the github repository. Install all the requirements.

## Test Run

1. Make sure the example data (LC716 dataset) is downloaded.

2. Change reference genome dependency files path location in the Rscript "Multiway.Contacts.Step1.R" , "Hub.Step2.R" . We have provided files for mm10 i.e mm10.blacklist.bed , mm10.chrom.sizes , F_GC_M_MboI_10Kb_el.mm10.txt , 10KbBin_mm10_blacklist.txt.

Change the path of above mentioned files in following variables in R script "Multiway.Contacts.Step1.R" :

```
mm10.chrom.size = "./mm10.chrom.sizes"   ######## Path to mm10 chrom sizes file ########
mm10.blacklist.region = "./mm10.blacklist.bed" ####### Path to mm10 Black list regions #########
```
Similarly change the path of above mentioned files in following variables in R script "Hub.Step2.R" :

```
mm10.blacklist = '10KbBin_mm10_blacklist.txt'     ##### Path to mm10 Blacklist region
annotation = 'F_GC_M_MboI_10Kb_el.mm10.txt'       ##### Path to mm10 10Kb MboI annotation file
```

3. Set following two path for dataset in script "Multiway.Contacts.Step1.R" :

```
dir1 = "~/LC716/bam.sc/"      #### Path to bam/sam files ######
dir2 = "~/LC716/chimeric.sc.pairs.mapq40/"     #### Path to save all intermediate files ######****
```

Set following path for dataset in script "Hub.Step2.R". We have provided example for "cellnames" , "color.loc" in example.data directory. :
```
dir<-'~/LC716/chimeric.sc.pairs.mapq40/'   #### Path to directory where all intermediate result files were saved in step1 (dir2 in step1)

cellnames = "./example.data/LC716_MAPQ40_cell_names_Shreya_light.txt"     ###### Path to example data cellnames file
color.loc = './example.data/color_20_YangOrder.txt'                       ###### Path to Color assigned to each celltype 

metadata = './052124_LC716_MAPQ40_metadata_MultiContact.txt'      ###### Filename to save multiway Contacts metadata file
final.result = './hub_20celltypes_MAPQ40.txt'                     ###### Filename to save final hub file
```

4. Run command below :

```bash
Rscript Multiway.Contacts.Step1.R
```
Then run second step when first step is complete :

```bash
Rscript Hub.Step2.R
```

### 2. File format requirements

#### 2.1 Sam file

```
> head CTTCTAAAGAACTAAC.OGC.int.sam
@SQ	SN:chr1	LN:195471971
@SQ	SN:chr2	LN:182113224
@SQ	SN:chr3	LN:160039680
@SQ	SN:chr4	LN:156508116
@SQ	SN:chr5	LN:151834684
@SQ	SN:chr6	LN:149736546
@SQ	SN:chr7	LN:145441459
@SQ	SN:chr8	LN:129401213
@SQ	SN:chr9	LN:124595110
@SQ	SN:chr10	LN:130694993
@SQ	SN:chr11	LN:122082543
@SQ	SN:chr12	LN:120129022
@SQ	SN:chr13	LN:120421639
@SQ	SN:chr14	LN:124902244
@SQ	SN:chr15	LN:104043685
@SQ	SN:chr16	LN:98207768
@SQ	SN:chr17	LN:94987271
@SQ	SN:chr18	LN:90702639
@SQ	SN:chr19	LN:61431566
@SQ	SN:chrX	LN:171031299
@SQ	SN:chrY	LN:91744698
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem -SP5M -T0 -t16 /projects/ps-renlab/share/bwa_indices/mm10.fa /home/y2xie/scratch/62.schiC_MouseCortex_NovaSeq_230723/02.trimmed/.///LC716_R1_BC_cov_val_1.fq.gz /home/y2xie/scratch/62.schiC_MouseCortex_NovaSeq_230723/02.trimmed/.///LC716_R3_BC_cov_val_2.fq.gz
CTTCTAAAGAACTAAC:A01535:363:HFY77DSX7:4:1101:28004:1110	65	chr11	63157895	60	151M	=	63144151	-13745	TGTCTGTACAGAGGTTTGGGGCTTCACCAGCATATCTGCACCAGGTTGCTTCATAGAAGAGTTAGGCAGGCTCATCTTGCAGCCTGGGTCCTCACTGCCTGTCTTGAATGACCTAGGAAGAAGTGTTATGCTGAAAGAATGAGAAATGAGG	:F:::FFFFFFFFF:FF::FFFFFFFFFFFFFFFF,FFFF:,:FFF:F:FFFFFFF,FFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFF:FFF:FFFFFFFFFFFFFFFF:FF:F,FFFFFFFFF,F:F:FFF:FFFFFFFFF	NM:i:0	MD:Z:151	MC:Z:84M66S	AS:i:151	XS:i:0
CTTCTAAAGAACTAAC:A01535:363:HFY77DSX7:4:1101:28004:1110	129	chr11	63144151	60	84M66S	=	63157895	13745	CTATTGTGCTCACTAGCTTCTTATGGTGCTAGTGCAAAAGTGGGCACACTCTAGCCTTAGGGGCCAAAGGAAAGACACCACTTGTTCTAGGCGGAGGTGGGCATTTACAACACTAGGCAGGGACGATGTGTCTGAGCCTCCCAGCTCAGA	F,,,:FFFFFF:FF:F,,,F,:,,FFFFF:,FFF,FF,,FFFFF:FFF:FF:,F::,F,FF:F:F:,,FF,,:F,F,F,,,::FFFFF,FF:FF,FFFFFFF,,,,F,:,FF::,FFF:FFF:,FF,F,FFFFFFFFF:FFF,FF,F:FF	NM:i:1	MD:Z:46A37	MC:Z:151M	AS:i:79	XS:i:19	SA:Z:chr11,63158153,-,65M85S,60,1;
```

#### 2.2 Metadata file with columns "Cellname","Celltype"
```
> head LC716_MAPQ40_cell_names_Shreya_light.txt
AAACGAACAATCAGGG ITL6GL
AAACGAACATCCGTGG D12MSN
AAACGAAGTTAGGCTT ITL23GL
AAACGAATCTCAGATG CTGL
AAACTCGAGGACTAGC MGL
AAACTCGAGTCTGCTA CTGL
AAACTCGAGTGAATAC ITL45GL
AAACTCGCACTTTGGA MGL
AAACTCGCATAAAGTG ITL5GL
AAACTCGCATCCATAG ASC
```

#### 2.2 Color assignment file with columns "Celltype","Color"
```
Subclass	Subclass.color
ASC	#A13F9E
OGC	#985695
OPC	#B9A0D6
MGL	#DF8DB1
VLMC	#A78BDF
ITL23GL	#E59D9F
ITL45GL	#DDA58F
ITL5GL	#E1A462
ITL6GL	#9E4954
```

## Contact Us

For any question, contact Ming Hu (hum@ccf.org), Shreya Mishra (mishras10@ccf.org), or submit an issue on GitHub.
