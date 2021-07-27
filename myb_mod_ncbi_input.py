### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python3 myb_mod_ncbi_input.py
					--gff <GFF_FILE>
					--cds <CDS_FILE>
					--out <OUTPUT_FOLDER>
					--name <SPEC_NAME>
					"""

import os, sys, re, subprocess

# --- end of imports --- #


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def main( arguments ):
	"""! @brief run everything """
	
	gff_file = arguments[ arguments.index('--gff')+1 ]
	cds_file = arguments[ arguments.index('--cds')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	name = arguments[ arguments.index('--name')+1 ]
	
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	cds_seqs = load_sequences( cds_file )
	
	cds_output_file = output_folder + name  + ".cds.fasta"
	gff_output_file = output_folder + name + ".gff"
	
	with open( cds_output_file, "w" ) as out:
		for seq in cds_seqs.keys():
			out.write( '>' + name + "_" + seq + "\n" + cds_seqs[ seq ] + "\n" )
	
	with open( gff_output_file, "w" ) as out:
		with open( gff_file, "r" ) as f:
			line = f.readline()
			while line:
				matches = re.findall( "rna\d+;", line )
				if len( matches ) > 0:
					out.write( line.replace( matches[0], name + "_" + matches[0] ) )
				else:
					out.write( line )
				line = f.readline()


if '--gff' in sys.argv and '--cds' in sys.argv and '--out' in sys.argv and '--name' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
