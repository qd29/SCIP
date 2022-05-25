# SCIP Setup Instructions
**This corresponds to Supplementary Text 2 of the SCIP paper.**

### Filtration & Prioritization Modules
**Prerequisites:** (a) A UNIX-based operating system (tested with CentOS 7). (b) The following software in the `$PATH` environment variable (version tested in parenthesis): `Perl` (v5.16), `R` (v3.5.1), `samtools` (v1.10), `bedtools` (v2.26), and `tabix` (v0.2.5).

1. Download the required annotation file at the following links. Place the file in a folder you intend to install the SCIP backend modules.

    - https://drive.google.com/file/d/1rWwkJ-eFDL1xT0DzraNi9oQ-NtC62TCe/view (hg19)
    - https://drive.google.com/file/d/1N4dt_UZ3CMw7CtkmU-KQDZe2tFmnt9kd/view (hg38)

2. At the command line, copy the following code block (depending on the genome build) and press the enter/return key to execute. The SCIP Filtration & Prioritization Modules will be installed in a new folder called `SCIP_backend`. For simplicity, this folder will be hereafter referred to as `./`.

    For hg19:
    ```
    mkdir -p ./SCIP_backend
    mv ./SCIP_hg19_files.tar.gz ./SCIP_backend
    cd ./SCIP_backend
    tar -xvzf ./SCIP_hg19_files.tar.gz
    mkdir -p d1temp_server d1stat user_data app_temp_file
    git clone https://github.com/qd29/SCIP.git
    mv ./SCIP/filtration_prioritization_hg19/* ./
    rm -rf ./SCIP
    ```
    
    For hg38:
    ```
    mkdir -p ./SCIP_backend
    mv ./SCIP_hg38_files.tar.gz ./SCIP_backend
    cd ./SCIP_backend
    tar -xvzf ./SCIP_hg38_files.tar.gz
    mkdir -p d1temp_server d1stat user_data app_temp_file
    git clone https://github.com/qd29/SCIP.git
    mv ./SCIP/filtration_prioritization_hg38/* ./
    rm -rf ./SCIP
    ```
    
3. Obtain the OMIM `genemap2.txt` file from https://www.omim.org/downloads, registration required. Place the file in the `hg19_files` or `hg38_files` folder under `./`.

### Filtration & Prioritization Modules Configuration File
**Configuration File Name:** `pipeline_config.txt` (hg19) or `pipeline_config_hg38.txt`

4. Specify `ALIGNMENT_PATH` to the path containing alignment files. Each sample must have its own subdirectory. For example, alignment file for sample001 should be stored at:

    >ALIGNMENT_PATH/sample001/sample001.bam

Corresponding index files (`bai` or `crai`) must also be available.

5. Specify `REF_BAM` as the path to the whole-genome alignment file of a reference sample.

    - NA12878_S1.bam: ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR194/ERR194147 (hg19)
    - NA12878.final.cram: ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334 (hg38)
    - **a.** Corresponding index files (`bai` or `crai`) must also be available. For the hg19 BAM file, EBI did not provide an index, thus you need to generate the index with the `samtools index` command.
    - **b.** For CRAM files, you may need to convert it to BAM and/or set the `$REF_PATH` and `$REF_CACHE` environment variables, so that reads can be printed using the `samtools view [file]` command.

**Setup of the SCIP Filtration & Prioritization Modules is now complete.** These steps are required only once. To run the Filtration & Prioritization Modules, go to step 1 of the <a href='https://github.com/qd29/SCIP/blob/main/usage_instructions.md'>SCIP Usage Instructions</a>. 

**Advanced Users:** (1) Consider periodically update the annotation files (e.g., OMIM, GenCC, ClinGen). (2) The `expression_file`, `GO_terms`, `gene_interest`, and `search_terms` entries in the configuration file can be customized by the user. See Table S4 of the SCIP paper for formatting details. 

### Visualization Module

6. The SCIP Visualization Module can be run on Windows or macOS computers. Create an empty directory on the computer (hereafter referred to as the working directory). 

7. Download the `SCIP_interface.R` (hg19) or `SCIP_interface_hg38.R` script and the `interface_config.txt` files from GitHub, place them in the working directory. 

8. Install R, RStudio, and R packages `shiny`, `DT` and `plotrix` (and any dependencies).
 
9. The `./` directory of the Filtration & Prioritization Modules needs to be accessible to the computer running the Visualization Module. At our institution, this was done by mounting the computer server running the Filtration & Prioritization Modules as a network volume on the PC running the Visualization Module.

    - **a.** When the network volume option is not available, the user can download all files in the `./app_temp_file` and `./user_data` directories to a local directory. Note that the relative location of `app_temp_file` and `user_data` must be maintained (i.e., they must sub-directories within the same directory).

10.	In the `interface_config.txt` file:

    - **a.** Modify the `TEMP_FILE_DIR` entry to the path to the `app_temp_file`.

    - **b.** Modify the `ROOT_DIR` entry to the path to the working directory (see step 6).

    - **c.** Modify the `USER` entry with your identifier (e.g., name/email). This is optional and is only used to track interpretations across multiple users.

**Setup of the SCIP Visualization Module is now complete.** These steps are required only once. **Exception:** if you use the approach described in step 9a, you will need to re-download the `./app_temp_file` and `./user_data` directories every time you have new samples. To run the Visualization Module, go to step 6 of the <a href='https://github.com/qd29/SCIP/blob/main/usage_instructions.md'>SCIP Usage Instructions</a>. 
