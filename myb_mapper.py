### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python3 myb_mapper.py
					--myb <MYB_CDS_FASTA_FILE>
					--cds <CDS_FILE>
					--out <OUTPUT_FOLDER>
					--name <SPEC_NAME>
					
					optional:
					--mode <BLAST_TYPE>(blastn|tblastn)[blastn]
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


def load_best_blast_hit( blast_result_file ):
	"""! @brief load best blast hit per query """
	
	best_hits = {}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				data = best_hits[ parts[0] ]
				if float( parts[-1] ) > data['score']:
					del best_hits[ parts[0] ]
					best_hits.update( { parts[0]: { 'score': float( parts[-1] ), 'subject': parts[1] } } )
			except:
				best_hits.update( { parts[0]: { 'score': float( parts[-1] ), 'subject': parts[1] } } )
			line = f.readline()
	return best_hits


def main( arguments ):
	"""! @brief run everything """
	
	myb_seq_file = arguments[ arguments.index('--myb')+1 ]
	cds_file = arguments[ arguments.index('--cds')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	name = arguments[ arguments.index('--name')+1 ]
	
	if '--mode' in arguments:
		mode = arguments[ arguments.index('--mode')+1 ]
		if mode not in [ "blastn", "tblastn" ]:
			mode = "blastn"
	else:
		mode = "blastn"
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	query_seqs = load_sequences( myb_seq_file )
	print ( "number of query sequences: " + str( len( query_seqs.keys() ) ) )
	
	blast_result_file = output_folder + "blast_results.txt"
	if not os.path.isfile( blast_result_file ):
		if mode == "blastn":
			p = subprocess.Popen( args="blastn -query " + myb_seq_file + " -subject " + cds_file + " -out " + blast_result_file + " -outfmt 6 -evalue 0.00001 -word_size 30", shell=True )
		else:
			p = subprocess.Popen( args="tblastn -query " + myb_seq_file + " -subject " + cds_file + " -out " + blast_result_file + " -outfmt 6 -evalue 0.00001", shell=True )
		p.communicate()
	
	blast_results = load_best_blast_hit( blast_result_file )
	print ( "number of BLAST hits: " + str( len( blast_results.keys() ) ) )
	
	mapping_table = output_folder + name + ".txt"
	with open( mapping_table, "w" ) as out:
		for each in blast_results.keys():
			out.write( blast_results[ each ]['subject'] + "\t" + each + "\n" )


if '--myb' in sys.argv and '--cds' in sys.argv and '--out' in sys.argv and '--name' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
