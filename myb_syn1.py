### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python3 myb_syn1.py
					option1:
					--config <CONFIG_FILE>
					--out <OUTPUT_FOLDER>
					
					option2:
					--gff <COMMA_SEPARATED_LIST_OF_GFF_FILES>
					--cds <COMMA_SEPARATED_LIST_OF_CDS_FILES>
					--feature <COMMA_SEPARATED_LIST_OF_FEATURE_NAMES>(mRNA|transcript)
					--ID <COMMA_SEPARATED_LIST_OF_ID_NAMES>(ID)
					--specs <COMMA_SEPARATED_LIST_OF_SPECIES_NAMES>
					--MYBS <COMMA_SEPARATED_LIST_OF_ID_FILES>
					--out <OUTPUT_FOLDER>
					"""

import os, sys, re, subprocess
import numpy as np

# --- end of imports --- #

def prepare_customized_blocks_file( genes, input_file, output_file ):
	"""! @brief select blocks of interest """
	
	data = []
	with open( input_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				#ID = re.findall( "[a-z]{4}", parts[0] )[0]
				data.append( { 'ID': parts[0], 'line': line } )
			except IndexError:
				print ( line )
			line = f.readline()
	
	# --- get min and max index --- #
	min_index = len( data )
	max_index = 0
	for gene in genes:
		for idx, entry in enumerate( data ):
			if gene == entry['ID']:
				if idx < min_index:
					min_index = idx + 0
				if idx > max_index:
					max_index = idx + 0
	
	if min_index < max_index:
		with open( output_file, "w" ) as out:
			lines = []
			for each in data[ min_index: max_index+1 ]:
				lines.append( each['line'] )
			out.write( "".join( lines ) )
	else:
		print ( "ERROR (BP): min index >= max index (candidate genes missing?)" )
	

def merge_block_files( block_file_per_spec, merged_block_file ):
	"""! @brief merge single block files """
	
	# --- read all data --- #
	mapping_table = {}
	sample_names = []
	genes = []
	for idx, filename in enumerate( block_file_per_spec ):
		sample_names.append( filename )
		tmp = {}
		with open( filename, "r" ) as f:
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				if idx == 0:
					genes.append( parts[0] )
				tmp.update( { parts[0]: parts[1:] } )
				line = f.readline()
		mapping_table.update( { filename: tmp } )
	# --- generate output file --- #
	with open( merged_block_file, "w" ) as out:
		for gene in genes:
			new_line = [ gene ]
			for sample in sample_names:
				new_line.append( "\t".join( mapping_table[ sample ][ gene ] ) )
			out.write( "\t".join( new_line ) + "\n" )


def modify_lifted_anchor_file( lifted_anchor_file ):
	"""! @brief modify lifted anchor file to remove entries of empty contigs """
	
	with open( lifted_anchor_file, "r" ) as f:
		content = f.read()
	content = content.replace( "###\n###\n", "###\n" ).replace( "###\n###\n", "###\n" ).replace( "###\n###\n", "###\n" ).replace( "###\n###\n", "###\n" ).replace( "###\n###\n", "###\n" ).replace( "###\n###\n", "###\n" ).strip()
	
	with open( lifted_anchor_file, "w" ) as out:
		out.write( content )


def load_config( config_file ):
	"""! @brief load config file """
	
	gff_files = []
	cds_files = []
	feature_types = []
	ID_types = []
	spec_names = []
	MYB_files = []
	
	with open( config_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != "#":
				parts = line.strip().split('\t')
				spec_names.append( parts[0] )
				ID_types.append( parts[1] )
				feature_types.append( parts[2] )
				cds_files.append( parts[3] )
				gff_files.append( parts[4] )
				MYB_files.append( parts[5] )
				
			line = f.readline()
	
	return gff_files, cds_files, feature_types, ID_types, spec_names, MYB_files


def load_MYB_IDs( spec_names, MYB_files ):
	"""! @brief load MYB IDs from given files """
	
	if len( spec_names ) != len( MYB_files ):
		sys.exit( "ERROR: number of species names does not match number of provided MYB ID files" )
	
	MYB_IDs = {}
	for idx, name in enumerate( spec_names ):
		IDs = []
		with open( MYB_files[ idx ], "r" ) as f:
			line = f.readline()
			while line:
				if "\t" in line:
					IDs.append( line.strip().split('\t')[0] )
				else:
					IDs.append( line.strip() )
				line = f.readline()
		MYB_IDs.update( { name: IDs } )
	return MYB_IDs


def calculate_conservation_per_gene( column1, column2, stretch=5 ):
	"""! @brief calculate conservation or each gene based on identified syntenic genes """
	
	conservation_per_gene = {}
	for idx, gene in enumerate( column1 ):
		if column2[idx] == ".":
			conservation_per_gene.update( { gene: 0 } )	#set to 0 if corresponding genes is missing
		else:
			x = max( [ idx-stretch, 0 ] )
			y = min( [ idx+stretch, len( column1 )-1 ] )
			score = int( 100.0*( (y-x) - column2[ x: y ].count('.') ) / (y-x) )
			conservation_per_gene.update( { gene: score } )
	return conservation_per_gene


def load_mapping_table( my_block_file ):
	"""! @brief load block file into mapping table """
	
	mapping_table = {}
	column1 = []
	column2 = []
	with open( my_block_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			mapping_table.update( { parts[0]: parts[1] } )
			column1.append( parts[0] )
			column2.append( parts[1] )
			line = f.readline()
	conservation_per_gene = calculate_conservation_per_gene( column1, column2 )
	return mapping_table, conservation_per_gene


def main( arguments ):
	"""! @brief run everything """
	
	if '--config' in arguments:
		config_file = arguments[ arguments.index('--config')+1 ]
		gff_files, cds_files, feature_types, ID_types, spec_names, MYB_files = load_config( config_file )
	else:
		gff_files = arguments[ arguments.index('--gff')+1 ].split(',')
		cds_files = arguments[ arguments.index('--cds')+1 ].split(',')
		feature_types = arguments[ arguments.index('--feature')+1 ].split(',')	#tag in third column
		ID_types = arguments[ arguments.index('--ID')+1 ].split(',')	#tag in last field
		spec_names = arguments[ arguments.index('--specs')+1 ].split(',') #species names appear as labels in figure
		MYB_files = arguments[ arguments.index('--MYBs')+1 ].split(',')

	cscore_cutoff = 0.1	#needs to be a value between 0 and 1
	number_of_matching_regions = 1	#get one corresponding region for each region in reference (needs to be increased for polyploids)

	output_folder = arguments[ arguments.index('--out')+1 ]

	# --- generate output folder and set it as working directory (important, because relative paths are used!) --- #
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	os.chdir( output_folder )


	# --- production of bed files --- #
	bed_files = []
	for idx, filename in enumerate( gff_files ):
		bed_file = spec_names[ idx ] + ".bed"
		if not os.path.isfile( bed_file ):
			p = subprocess.Popen( args="python3 -m jcvi.formats.gff bed --type " + feature_types[ idx ] + " --key=" + ID_types[ idx ] + " " + filename + " -o " + bed_file, shell=True )
			p.communicate()
		bed_files.append( bed_file )


	# --- generate clean CDS files --- #
	clean_cds_files = []
	for idx, filename in enumerate( cds_files ):
		cds_file = spec_names[ idx ] + ".cds"
		if not os.path.isfile( cds_file ):
			p = subprocess.Popen( args="python3 -m jcvi.formats.fasta format " + filename + " " + cds_file, shell=True )
			p.communicate()
		clean_cds_files.append( cds_file )


	# --- generate one block file per species --- #
	block_file_per_spec = {}	#could be useful to perform comparison across multiple species (might not be relevant)
	for idx1, spec1 in enumerate( spec_names ):
		for idx2, spec2 in enumerate( spec_names ):
			if idx2 != idx1:	#calculating all connections only in one direction to avoid redundancy
				# --- identification of putative orthologs --- #
				anchor_file = ".".join( [ spec_names[idx1], spec2 ] ) + ".anchors "
				lifted_anchor_file = ".".join( [ spec_names[idx1], spec2 ] ) + ".lifted.anchors"
				if not os.path.isfile( lifted_anchor_file ):
					p = subprocess.Popen( args="python3 -m jcvi.compara.catalog ortholog " + " ".join( [ spec_names[idx1], spec2 ] ) + " --cscore=" + str( cscore_cutoff ) + " --no_strip_names", shell=True )
					p.communicate()

				# --- generate more succint anchors file --- #
				modify_lifted_anchor_file( lifted_anchor_file )
				block_file = ".".join( [ spec_names[idx1], spec2 ] ) + ".i" + str( number_of_matching_regions ) + ".blocks"
				if not os.path.isfile( block_file ):
					p = subprocess.Popen( args="python3 -m jcvi.compara.synteny mcscan " + bed_files[idx1] + " " + lifted_anchor_file + " --iter=" + str( number_of_matching_regions ) + " -o " + block_file, shell=True )
					p.communicate()
				
				try:	#this block might not be relevant
					block_file_per_spec[ spec_names[idx1] ].append( block_file )
				except KeyError:
					block_file_per_spec.update( { spec_names[idx1]: [ block_file ] } )
	
	
	# --- get MYB IDs per species --- #
	MYB_IDs = load_MYB_IDs( spec_names, MYB_files )
	
	# --- iterate over species and check conservation of MYBs --- #
	detailed_output_file = output_folder + "detailed_output.txt"
	datamatrix = {}
	with open( detailed_output_file, "w" ) as details:
		details.write( "\t".join( [ "Species1", "Species2", "MYB1", "MYB2", "Conservation" ] ) + "\n" )
		for idx1, spec1 in enumerate( spec_names ):
			datamatrix.update( { spec1: {} } )
			for idx2, spec2 in enumerate( spec_names ):
				if idx1 != idx2:
					block_file= spec1 + "." + spec2 + ".i1.blocks"
					if os.path.isfile( block_file ):
						pairs, conservation_per_gene = load_mapping_table( block_file )
					else:
						#sys.exit( "ERROR: block file missing for species combi: " + spec1 + "\t" + spec2 )
						print( "ERROR: block file missing for species combi: " + spec1 + "\t" + spec2 )
					MYBs = MYB_IDs[ spec1 ]
					for MYB in MYBs:
						try:
							details.write( "\t".join( list( map( str, [ spec1, spec2, MYB, pairs[ MYB ], conservation_per_gene[ MYB ] ] ) ) ) + "\n" )
							try:
								datamatrix[ spec1 ][ MYB ].append( conservation_per_gene[ MYB ] )
							except KeyError:
								datamatrix[ spec1 ].update( { MYB : [ conservation_per_gene[ MYB ] ] } )
						except KeyError:
							print ( MYB )
	summary_output_file = output_folder + "summary_output.txt"
	with open( summary_output_file, "w" ) as out:
		out.write( "\t".join( [ "Species", "MYB", "AverageConservation", "IndividualConservationValues" ] ) + "\n" )
		for spec in list( sorted( datamatrix.keys() ) ):
			for gene in list( sorted( datamatrix[ spec ].keys() ) ):
				out.write( "\t".join( list( map( str, [ spec, gene, np.mean( datamatrix[ spec ][ gene ] ), ",".join( list( map( str, datamatrix[ spec ][ gene ] ) ) ) ] ) ) ) + "\n" )
				

	# # --- merge block files of all species --- #
	# merged_block_file = spec_names[0] + ".blocks"
	# #os.popen( "python -m jcvi.formats.base join " + " ".join( block_file_per_spec ) + " --noheader | cut -f1,2,4,6 > " + merged_block_file )
	# merge_block_files( block_file_per_spec, merged_block_file )	#MY SOLUTION, but might enfores plot with wrong locus


if '--gff' in sys.argv and '--cds' in sys.argv and '--feature' in sys.argv and '--ID' in sys.argv and '--specs' in sys.argv and '--chromosomes' in sys.argv:
	main( sys.argv )
elif '--config' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
