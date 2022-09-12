# Instructions for Updating SCIP Annotation Files

The following annotation files used by the SCIP Prioritization Module (Table S3 of the SCIP manuscript) need to be updated regularly.

Place them in the `hg19_files` (for hg19) or `hg38_files` (for hg38) folder, then modify the SCIP pipeline configuration file `pipeline_config.txt` (for hg19) or `pipeline_config_hg38.txt` (for hg38) to specify the paths to the new files.

### 1. ClinGen Dosage Curation Files (`ClinGen_dosage_region` and `ClinGen_dosage_gene`)

The <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8993917/'>current guidelines</a> recommend updating these files at least quarterly. Available at https://ftp.clinicalgenome.org/.

For hg19:</br>
https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh37.tsv </br>
https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh37.tsv

For hg38:</br>
https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv</br>
https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv

### 2. ClinVar CNV File (`ClinVar_CNV`)

The <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8993917/'>current guidelines</a> recommend updating this file at least quarterly. 

Path forthcoming.

### 3. OMIM and GenCC Gene-disease Association Files (`OMIM` and `GenCC`)

The <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8993917/'>current guidelines</a> recommend updating this file at least quarterly.

The updated `genemap2.txt` file must be directly obtained from OMIM.</br>
Path forthcoming for the updated GenCC file.

### 4. Non-coding Pathogenic Regions (`noncoding`)

Update this file if new pathogenic non-coding regions are discovered.</br>
There is no update yet to this file.
