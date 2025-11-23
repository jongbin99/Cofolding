import sys, os



def get_fingerprints(smiles_file, outname):
	pwd = os.getcwd() + '/'
	if os.path.isfile(pwd+outname+'.fp'):
		print('The fingerprint for these molecules already exists!')

	else:
		os.system('/nfs/home/ttummino/zzz.scripts/freechem/freechem-19.15.r4/bin/generatemd c ' + smiles_file + ' -k ECFP -2 > ' + outname + '.fp')
		os.system("sed -i -e 's/^/fingerprint = /'  " + outname + '.fp' )
		os.system('~jklyu/zzz.github/ChemInfTools/utils/convert_fp_2_fp_in_16unit/convert_fp_2_fp_in_uint16 ' + outname + '.fp ' + smiles_file+ ' ' + outname)

def main():
	smiles_file_hits = sys.argv[1]
	outname_hits = sys.argv[2]

	smiles_file_ligands = sys.argv[3]
	outname_ligands = sys.argv[4]

	get_fingerprints(smiles_file_hits, outname_hits)
	get_fingerprints(smiles_file_ligands, outname_ligands)

	os.system('nohup ~jklyu/zzz.github/ChemInfTools/utils/cal_Tc_matrix_uint16/cal_Tc_matrix_uint16 ' 
+ outname_hits + '_uint16.fp ' + smiles_file_hits + '  ' + outname_hits + '_uint16.count '
+ outname_ligands + '_uint16.fp ' + smiles_file_ligands + '  ' + outname_ligands + '_uint16.count '  + 'max_tc > log ')




main()
