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

```
python3 myb_syn1.py
option1:
--config
--out

option2:
--gff <COMMA_SEPARATED_LIST_OF_GFF_FILES>
--cds <COMMA_SEPARATED_LIST_OF_CDS_FILES>
--feature <COMMA_SEPARATED_LIST_OF_FEATURE_NAMES>(mRNA|transcript)
--ID <COMMA_SEPARATED_LIST_OF_ID_NAMES>[ID]
--specs <COMMA_SEPARATED_LIST_OF_SPECIES_NAMES>
--MYBs <COMMA_SEPARATED_LIST_OF_ID_FILES>
--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
```

`--config` specifies a config file containing paths to all files and also information about all other options. Information of each species is given in one row.

`--out` specifies the output folder. This folder will be created if it does not exist already.

`--gff` specifies a GFF3 file.

`--cds` specifies a FASTA file containing the coding sequences of a species. This file needs to match the GFF3 file for the JCVI analysis of synteny.

`--feature` specifies the feature tag that is used in the GFF3 file to indicate a mRNA entry. This is usually "mRNA" or "transcript".

`--ID` specifies the tag that is used to indicate the ID field in the last column of a GFF3 file.

`--specs` specifies a the names of all the species analysed in one batch.

`--MYBs` specifies the MYB ID containing text files. These IDs need to match the IDs in the CDS file.



## 3) Downstream analyses ##
