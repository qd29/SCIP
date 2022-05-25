# SCIP Usage Instructions
**This corresponds to Supplementary Text 1 of the SCIP paper.**

**When using SCIP for the first time, please follow the <a href='https://github.com/qd29/SCIP/blob/main/setup_instructions.md'>SCIP Setup Instructions</a> for initial setup steps.** This instruction assumes that all setup steps have been successfully completed. We use `./` to denote the directory created for SCIP in the <a href='https://github.com/qd29/SCIP/blob/main/setup_instructions.md'>SCIP Setup Instructions</a>.

**Prerequisites:**<br>
**Filtration & Prioritization Modules:** (a) A UNIX-based operating system (tested with CentOS 7). (b) The following software in the `$PATH` environment variable (version tested in parenthesis): `Perl` (v5.16), `R` (v3.5.1), `samtools` (v1.10), `bedtools` (v2.26), and `tabix` (v0.2.5).<br>
**Visualization Module:** Any operation system supported by RStudio and has a web browser (tested with Windows 10 and macOS Big Sur).


### Prepare Input Files (Variant Filtration & Prioritization Modules)
1. For each sample, prepare a tab-delimited file named `[name].unfiltered_CNV.txt` that includes all CNVs. Place the file in `./user_data`.

    - **a.** This file has 7 columns. The first four columns are chromosome, start position, end position, and type (DEL/DUP), respectively. The fifth and sixth columns are not used by SCIP and may contain free text. The last column is sample ID. Example format:<br>
      >1 620001 635000 DUP . . SAM-001-001

    - **b.** To allow SCIP to recognize members of the same family, we recommend naming samples using the following format: ABC-123-001. ABC-123 is the family ID that is identical for all family members. The last three digits denote family relationship, with 001, 002, 003, 004+ denoting proband, mother, father, and additional family members, respectively.

2. In the Filtration & Prioritization Modules Configuration File (`pipeline_config.txt` (hg19) or `pipeline_config_hg38.txt`), specify `SAMPLE_ID` as the path to the file containing conversion information between sample ID and alignment file name. If this file already exists, you may append information about the current sample to this file. 

    - **a.** Two tab-separated columns. Each row is a sample. The second column is the sample ID; the first column is the name prefix of the alignment file for this sample.

    - **b.** For example, if alignment file for sample `SAM-001-001` is `alignment001.bam`, specify the following line in this file:<br>
      >alignment001 SAM-001-001

    - **c.** Example file available at `./hg[19/38]_files/demo/sample_id.txt`.

3. In the Filtration & Prioritization Modules Configuration File (`pipeline_config.txt` (hg19) or `pipeline_config_hg38.txt`), specify `cohort_CNV` as the path to the file containing CNVs detected in internal cohorts. If this file already exists, you may append information about the current family to this file. **This file is optional, however, if one would like to study transmission pattern, CNVs from family members must be included in this file.** 

    - **a.** Optional information. If a user does not wish to provide this file, specify the path to an empty file. 

    - **b.** Six tab-separated columns. Each row is a CNV. The first four columns are chromosome, start position, end position, and type (DEL/DUP), respectively. The fifth column is sample ID, and the last column is the algorithm used to identify this variant. For example:<br>
      >1 10001 50001 DUP SAM-001-001 ERDS

### Run the SCIP Variant Filtration & Prioritization Modules
4. Run the following command. Denote the `[name]` specified in step 1 with the -n flag. For example (for the script name, change `hg19` to `hg38` as appropriate): `perl SCIP_backend_hg19.pl -n SAM-001-001`

    - **a.** Expected outputs on screen. Filtration Module: the following information will be printed on screen - `SCIP Filtration Module script 01/02/03 processing hg19/hg38 chr[1-22,X]`. Prioritization Module: the current date/time, name of the CNV being analyzed, and whether it generates new / reuses SAM and depth files will be printed. 

    - **b.** `samtools` and `tabix` may occasionally report warnings, e.g., the index file is older than the data file, protocol not supported, and/or failed to open reference (especially for CRAM files). No other errors/warnings are anticipated. 

5. **(Optional, Advanced Users Only)** Clean-up. SCIP stores SAM and depth information extracted from alignment files in the `./d1temp_server` directory. These files are no longer required after step 4 and may be removed. 

    - **a.** If the CNVs are re-analyzed in the future (e.g., using up-to-date annotation files), keeping these files allow SCIP to use them instead of querying the alignment BAM/CRAM files again. Therefore, we recommend keeping them (unless disk space poses an issue).

### Use the SCIP Visualization Module
6. Modify the `LIST_NAME` entry in the `interface_config.txt` to `[name].[hg19/hg38]` (see step 1).

8. Open `SCIP_interface.R` (or `SCIP_interface_hg38.R`) in RStudio. Click “Run App”, then at the top of the pop-up window, click “Open in Browser”. We recommend using Google Chrome.

    - **a.** RStudio sometimes produces a warning, which can be discarded unless it is a fatal error (e.g., the Visualization Module crashes).

**CNVs are now ready for manual review using the SCIP Visualization Interface.**
