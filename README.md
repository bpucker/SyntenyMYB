# SyntenyMYB
synteny analysis of MYBs

## 1) Prepare data sets ##

### a) Retrieve GFF3 file and matching CDS FASTA file ###
Extraction of CDS data set based on GFF3 and genomic FASTA file via get_pep_from_gff.py.
Adjust transcript IDs if data sets are retrieved from the NCBI.

### b) Retrieve annotated MYBs ###
If MYB IDs are matching the IDs in the CDS file, no additional processing is required. A text file should contain these IDs (one per line).

A FASTA file of the CDS or peptide sequences can be used to identify the corresponding IDs in the CDS file via myb_mapper.py.




## 2) Run synteny identificaiton via JCVI ##
All file locations can be specified via config file.


## 3) Downstream analyses ##
