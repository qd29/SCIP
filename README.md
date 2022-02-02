# SCIP: Suite for CNV Interpretation and Prioritization

### Step-by-step Instructions

**Preparations (Filtration & Prioritization Modules)**
1.	Create an empty directory on a UNIX-based operating system (e.g., CentOS), hereafter referred to as `./`, then create `./d1emp_server`, `./d1stat`, `./user_data`, and `./app_ temp_file`.
2.	Download the required annotation files at the following links, then decompress to `./` using `tar -xvzf [package name]`. 
https://drive.google.com/file/d/1rWwkJ-eFDL1xT0DzraNi9oQ-NtC62TCe/view?usp=sharing (hg19)
https://drive.google.com/file/d/1N4dt_UZ3CMw7CtkmU-KQDZe2tFmnt9kd/view?usp=sharing (hg38)
3.	Obtain the OMIM `genemap2.txt` file from https://www.omim.org/downloads, registration required. Place the file in the `./hg19_files` or `./hg38_files` directory.
4.	Download the SCIP backend scripts from GitHub, place them at `./`. For each genome build, there are 13 scripts: 4 for the filtration module and 9 for the prioritization module.
5.	For each sample, prepare a tab-delimited file named `[name].unfiltered_CNV.txt` that stores all CNVs. Place the file in `./user_data`.
    - a. This file contains seven columns. The first four columns are chromosome, start position, end position, and type (DEL/DUP), respectively. The fifth and sixth columns are not used by SCIP and may contain free text. The last column is sample ID. 
    - b. To allow SCIP recognize members of the same family, we recommend naming samples using the following format: `ABC-123-001`. `ABC-123` is the family ID that is identical for all family members. The last three digits denote family relationship, with 001, 002, 003, 004+ denotes proband, mother, father and additional family members, respectively.
    - c. For example:
    > 1 620001 635000 DUP . . SAM-001-001
