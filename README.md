# SCIP: Suite for CNV Interpretation and Prioritization

### Where are all the files required by SCIP?
**hg19**
- Annotation files: https://drive.google.com/file/d/1rWwkJ-eFDL1xT0DzraNi9oQ-NtC62TCe/view?usp=sharing
- SCIP Filtration & Prioritization Module scripts: under `filtration_prioritization_hg19` and `pipeline_config/pipeline_config.txt`
- Alignment for a reference sample: ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR194/ERR194147 (`NA12878_S1.bam`)
- SCIP Visualization Module scripts: under `visualization_hg19`

**hg38**
- Annotation files: https://drive.google.com/file/d/1N4dt_UZ3CMw7CtkmU-KQDZe2tFmnt9kd/view?usp=sharing
- SCIP Filtration & Prioritization Module scripts: under `filtration_prioritization_hg38` and `pipeline_config/pipeline_config_hg38.txt`
- Alignment for a reference sample: ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334 (`NA12878.final.cram`, conversion to BAM may be required)
- SCIP Visualization Module scripts: under `visualization_hg38`

### Step-by-step Instructions for Using SCIP

**Preparations (Filtration & Prioritization Modules)**

1.	Create an empty directory on a UNIX-based operating system (e.g., CentOS), hereafter referred to as `./`, then create `./d1emp_server`, `./d1stat`, `./user_data`, and `./app_temp_file`.
2.	Download the required annotation files at the following links, then decompress to `./` using `tar -xvzf [package name]`. 
    - https://drive.google.com/file/d/1rWwkJ-eFDL1xT0DzraNi9oQ-NtC62TCe/view?usp=sharing (hg19)
    - https://drive.google.com/file/d/1N4dt_UZ3CMw7CtkmU-KQDZe2tFmnt9kd/view?usp=sharing (hg38)
3.	Obtain the OMIM `genemap2.txt` file from https://www.omim.org/downloads, registration required. Place the file in the `./hg19_files` or `./hg38_files` directory.
4.	Download the SCIP backend scripts from GitHub, place them at `./`. For each genome build, there are 13 scripts: 4 for the filtration module and 9 for the prioritization module.
5.	For each sample, prepare a tab-delimited file named `[name].unfiltered_CNV.txt` that stores all CNVs. Place the file in `./user_data`.
    - a. This file contains seven columns. The first four columns are chromosome, start position, end position, and type (DEL/DUP), respectively. The fifth and sixth columns are not used by SCIP and may contain free text. The last column is sample ID. 
    - b. To allow SCIP recognize members of the same family, we recommend naming samples using the following format: `ABC-123-001`. `ABC-123` is the family ID that is identical for all family members. The last three digits denote family relationship, with 001, 002, 003, 004+ denotes proband, mother, father and additional family members, respectively.
    - c. For example:
    > 1 620001 635000 DUP . . SAM-001-001

**Filtration & Prioritization Modules Configuration File**

6.	Download the SCIP `pipeline_config.txt` (hg19) or `pipeline_config_hg38.txt` file from GitHub, place it at `./`. Modify various entries in this file according to steps 7–16. 
7.	Specify `SAMPLE_ID` as the path of the file containing conversion information between sample ID and alignment file name. 
    - a. Two tab-separated columns. Each row is a sample.
    - b. The second column is the sample ID; the first column is the name prefix of the alignment file for this sample.
    - c. For example, if alignment file for sample SAM-001-001 is alignment001.bam, specify the following line in this file:
    > alignment001 SAM-001-001
    - d. Example file available at `./hg[19/38]_files/demo/sample_id.txt`.
8.	Specify `CURRENT_PATH` to the absolute path of `./`.
9.	Specify `ALIGNMENT_PATH` to the path containing alignment files. Each sample must have its own subdirectory. For example, alignment file for sample001 should be stored at: `ALIGNMENT_PATH/sample001/sample001.bam`. Corresponding index files (bai or crai) must also be available.
10.	Specify `OMIM` as the path of the OMIM `genemap2.txt` file. 
11.	Specify `expression_file` as the path of the file that contains expression level per gene in the tissue-of-interest. This information is optional. If a user does not wish to provide this file, specify the path to an empty file. 
    - a. Two tab-separated columns. Each row is a gene. The first column is HGNC gene symbol, and the second column is expression level. For example:
    > A1BG 7.627
    - b. Example file available at `./hg[19/38]_files/demo/gene_exp_demo.txt`.
12.	Specify `GO_terms` as the path of the file that contains GO terms information per gene. This information is optional. If a user does not wish to provide this file, specify the path to an empty file.
    - a. Two tab-separated columns. Each row is a gene. The first column is HGNC gene symbol; the second column is GO-terms-of-interest for this gene, separated by a vertical bar. For example:
    > DGKA Signaling|CalciumIon
    - b. Example file available at `./hg[19/38]_files/demo/go_terms_genes_demo.txt`.
13.	Specify `REF_BAM` as the path of the whole-genome alignment file of a reference sample. 
    - a. `NA12878_S1.bam`: ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR194/ERR194147 (for hg19).
    - b. `NA12878.final.cram`: ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334 (for hg38). Based on our experience, it might be necessary to convert the CRAM to BAM.
    - c. Corresponding index files (bai or crai) must also be available.
14.	Specify `cohort_CNV` as the path of the file containing CNVs detected in internal cohorts. If one would like to study transmission pattern, CNVs from family members must be included in this file. This information is optional. If a user does not wish to provide this file, specify the path to an empty file. 
    - a. Six tab-separated columns. Each row is a CNV.
    - b. The first four columns are chromosome, start position, end position, and type (DEL/DUP), respectively. The fifth column is sample ID, and the last column is the algorithm used to identify this variant. For example:
    > 1 10001 50001 DUP SAM-001-001 ERDS
    - c. See Step 5(b) for recommendations regarding sample ID format.
15.	Specify `gene_interest` as the path of the file containing a list of candidate genes. This information is optional. If a user does not wish to provide this file, specify the path to an empty file.
    - a. One column, HGNC gene symbols. 
    - b. Example file available at `./hg[19/38]_files/demo/candidate_gene_list_demo.txt`.
16.	Specify `search_terms` as the path of the file containing search terms to be listed in the Genes table in Section 6 of the SCIP Visualization Module.
    - a. Two tab-separated columns. The first column is the search term (for performing the Google search), based on your disease(s)-of-interest. The second column is abbreviation (for display in the Genes table). If abbreviation is not needed, the two columns can be identical. For example: 
    > developmental delay     DD
    - b. Example file available at `./hg[19/38]_files/demo/search_terms_demo.txt`.

*The above set-up steps are required only once. Periodic updates of some annotation file (e.g., OMIM, GenCC, ClinGen) may be recommended. See Table S4 of the SCIP manuscript for details.*

**Running the SCIP Variant Filtration Module**

17.	Run the four SCIP_filt scripts sequentially. Denote the `[name]` specified in step 5 with the `-n` flag. See the example below.
    - a. The script will update its progress by printing `SCIP Filtration Module script 01/02/03 processing hg19/hg38 chr[1-22,X]`. No errors are anticipated. 
    - b. The `[name].filtered.txt` in the `./user_data` directory contains CNVs that passed the SCIP Filtration Module. It was reformatted into `[name].hg19/38. variant_list.txt` for the Prioritization Module. 

```
perl SCIP_filt_01_hg38.pl -n MSG-4093
perl SCIP_filt_02_hg38.pl -n MSG-4093
perl SCIP_filt_03_hg38.pl -n MSG-4093
perl SCIP_filt_04_hg38.pl -n MSG-4093
```

**Running the SCIP Prioritization Module**

18.	Run `perl SCIP_pri_01_hg19/38.pl`. Use the `-n` flag to specify `[name].[hg19/hg38]`, also always specify `-u 1`. The `-u` flag is reserved for parallelization (currently disabled). See line 11 of the script for information. 
    - a. Note that `samtools` and `R` must be installed.
    - b. `samtools` may occasionally report warnings, e.g., `the index file is older than the data file`, `protocol not supported`, and/or `failed to open reference` (especially for CRAM files). These warnings are likely harmless and may be ignored.
    - c. The script will update its progress by printing each CNV analyzed and date/time. It will print whether it is generating new or reusing SAM and depth files (see step 20). No errors are anticipated except mentioned in (b) above. 
19.	Run `perl SCIP_pri_02_hg19/38.pl`. Use the `-n` flag to specify `[name].[hg19/hg38]`.
    - a. The SCIP Prioritization Module generates files to be used by the Visualization Module, stored in the `./app_temp_file` directory. Priority scores are stored in the `[name].[hg19/hg38].pipeline_summary.txt` file in the `./user_data` directory.
20.	Optional clean-up. SCIP stores SAM and depth information extracted from alignment files in the `./d1emp_server` directory. These files are no longer required after step 18 and may be removed. 
    - a. If the CNVs are re-analyzed in the future (e.g., using up-to-date annotation files), keeping these files allow SCIP to use them instead of querying the alignment BAM/CRAM files again. Therefore, we recommend keeping them (unless disk space is an issue).

**Preparations (Visualization Module)**

21.	The SCIP Visualization Module can be run on Windows or macOS computers. Create an empty directory on the computer (hereafter referred to as the working directory). 
22.	Download the `SCIP_interface.R` (hg19) or `SCIP_interface_hg38.R` script and the `interface_config.txt` files from GitHub, place them in the working directory. 
23.	Install `R`, `RStudio`, and R packages `shiny`, `DT` and `plotrix` (and any dependencies). 
24.	The `./` directory of the Filtration & Prioritization Modules needs to be accessible to the computer running the Visualization Module. In our institution, this was done by mounting the computer server running the Filtration & Prioritization Modules as a network volume using the SMB protocol on the PC running the Visualization Module.
    - a. When the network volume option is not available, the user can download all files in the `./app_temp_file` and `./user_data` directories to a local directory. Note that the relative location of `app_temp_file` and `user_data` must be maintained (i.e., they must sub-directories within the same directory).
25.	In the `interface_config.txt` file:
    - a. Modify the `TEMP_DIR` entry to the path of the app_temp_file.
    - b. Modify the `ROOT_DIR` entry to the path of the working directory (see step 22).
    - c. Modify the `USER` entry with your name/email. This is optional and is only used to track interpretations across multiple users. Otherwise, enter any string.

*The above set-up steps are required only once.*

**Running the SCIP Visualization Module**

26.	Modify the `LIST_NAME` entry in the interface_config.txt to `[name].[hg19/hg38]` (see steps 5, 18).
27.	Open `SCIP_interface.R` (or `SCIP_interface_hg38.R`) in `RStudio`. Click “Run App”, then at the top of the pop-up window, click “Open in Browser”. We recommend using Google Chrome.
    - a. RStudio sometimes will produce a warning, which can be discarded unless it is a fatal error (e.g., the Visualization Module crashes).
28.	CNVs are now ready for manual review using the SCIP Visualization Interface. 
