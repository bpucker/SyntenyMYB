# SyntenyMYB
synteny analysis of MYBs

## 1) Prepare data sets ##

### a) Retrieve GFF3 file and matching CDS FASTA file ###
It is important that the sequence IDs in the FASTA file are matching the transcript/mRNA IDs in the GFF3 file. Extraction of CDS data set based on GFF3 and genomic FASTA file via get_pep_from_gff.py (translation step can be skipped) is possible to ensure this.

If datasets were retrieved from the NCBI, it might be necessary to modify the mRNA IDs. These are "rna0", "rna1", "rna2", .... which can cause issues in some of the JCVI analysis steps. A script is available to adjust transcript IDs in the CDS FASTA file and in the GFF3 file.


```
python3 myb_mod_ncbi_input.py
--gff <GFF3_FILENAME>
--cds <CDS_FILENAME>
--name <NAME>
--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
```

`--gff` specifies a GFF3 file. It is important that the IDs of mRNAs in this file are matching the IDs of the coding sequences in the FASTA file.

`--cds` specifies a FASTA file containing coding sequences. It is important that the IDs in this file are matching the IDs of GFF3 file.

`--name` specifies the prefix of all output filenames.

`--out` specifies the output folder. This folder will be created if it does not exist already.


### b) Retrieve annotated MYBs ###
Option 1: If MYB IDs are matching the IDs in the CDS file, no additional processing is required. A text file should contain these IDs (one per line).

Option 2: A FASTA file containing the MYB CDS or peptide sequences can be used to identify the corresponding IDs in the CDS file via myb_mapper.py.

```
python3 myb_mapper.py
--myb <MYB_SEQUENCE_FILENAME>
--cds <CDS_FILENAME>
--name <NAME>
--out <FULL_PATH_TO_OUTPUT_DIRECTORY>

optional:
--mode <BLAST_TYPE>[blastn]
```

`--myb` specifies a FASTA file containing all MYB sequences. These sequences can be CDS or peptide, but it is important to select the appropriate mode for the mapping to the CDS IDs.

`--cds` specifies a FASTA file containing coding sequences.

`--name` specifies the prefix of all output filenames.

`--out` specifies the output folder. This folder will be created if it does not exist already.

`--mode` specifies the BLAST method to use for the mapping. Options are blastn (for CDS) and tblastn (for peptides). The default is blastn.




## 2) Run synteny identification via JCVI ##
All file locations can be specified via config file.

myb_syn1.py


## 3) Downstream analyses ##
