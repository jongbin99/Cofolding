import sys, os



def get_fingerprints(smiles_file, outname):
	pwd = os.getcwd() + '/'
	if os.path.isfile(pwd+outname+'.fp'):
		print('The fingerprint for these molecules already exists!')

	else:
		os.system('python ~jklyu/zzz.github/ChemInfTools/utils/teb_chemaxon_cheminf_tools/generate_chemaxon_fingerprints.py ' + smiles_file + ' ' + outname)
        	os.system('~jklyu/zzz.github/ChemInfTools/utils/convert_fp_2_fp_in_16unit/convert_fp_2_fp_in_uint16 ' + outname + '.fp ' + smiles_file+ ' ' + outname)

def main():
	smiles_file_hits = sys.argv[1]
	outname_hits = sys.argv[2]
	max_tc = sys.argv[3]
	max_clusters = sys.argv[4]

	get_fingerprints(smiles_file_hits, outname_hits)
	print('this one worked')
	os.system('nohup ~jklyu/zzz.github/ChemInfTools/utils/best_first_clustering_uint16/best_first_clustering_uint16 '
+ outname_hits + '_uint16.fp ' + outname_hits + '_uint16.count ' + smiles_file_hits + ' ' + max_tc + ' ' + max_clusters + ' > bfc_log &')




main()
