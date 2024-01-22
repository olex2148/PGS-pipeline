#!/usr/bin/env python

"""Workflow to parse and compute polygenic scores

gwf workflow to parse all GWAS sumstats in data folder and generate PRS from them 
using LDpred2-auto and lassosum as implemented in LDpred2

Example:
	gwf run
	
Todo:
	* Add arguments to gwf command to specify folder with input data and a 
	  shrinkage coefficient if adjustment is relevant

@Author: Ole S Hansen
@Date: 

"""

from gwf import Workflow, AnonymousTarget
import glob
from datetime import date
from code.modpath import modpath

gwf = Workflow()

### Defining gwf templates
def parse_sumstats(inputfile):
	'''
	Template for running the r script "prepare_and_parse_sumstats.R" which parses GWAS summary statistics
	'''

	# Using modpath() to create name of output file from inputfile - keeps basename,
	# but gets another path and another suffix
	outputfile = modpath(inputfile, parent=('steps/parsed_sumstats'), suffix=('_parsed.rds'))

	# Defining inputs, outputs and ressources
	inputs = [inputfile]
	outputs = [outputfile]
	options = {
		'memory': '10g',
		'walltime': '00:30:00',
		'cores': '2'
	}
	
	# Command to be run in the terminal
	spec = f'''
	
	cd ~/NCRR-PRS/faststorage/osh/PGS/pgs_workflow/

	Rscript code/parser.R {inputfile} {outputfile}
	
	'''

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def compute_pgs(inputfile):
	'''
	Template for running the r script "pgs_model.R" which computes PGS models using parsed sumstats
	'''
	today = date.today()

	base = modpath(inputfile, parent=(''), suffix=('_parsed.rds', ''))      # Getting the base name from the inputfile 
	base_path = f'results/{base}/{base}'                                    # New path with sumstat-specific folder (and filename without suffix)

	model_out = f'{base_path}_raw_models.rds'                               # Models, parameters, and scores are put in a folder specific to the sumstats
	scores_out = f'{base_path}_scores.rds'
	parameters_out = f'{base_path}_auto_paramters.rds'
	foelgefil = f'results/følgefiler/{today}/{base}_følgefil.xlsx'          # Følgefil is put in a separate folder specific to the batch run (date in folder name)

	inputs = [inputfile]
	outputs = [model_out, scores_out, parameters_out, foelgefil]
	options = {
		'memory': '60g',
		'walltime': '12:00:00',
		'cores': '23'
	}

	spec = f'''

	mkdir -p results/{base}
	mkdir -p results/følgefiler/{today}

	Rscript code/pgs_model.R {inputfile} {base_path}
	
	'''
	
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

### Defining targets

# Input
sumstats_path = input("Path to summary statistics: ")
sumstats = glob.glob(sumstats_path)

# Submitting jobs for every sumstat file in data folder (parsing and computing PGSs)
for file in sumstats:
	base_name = modpath(file, parent='', suffix='')

	a = gwf.target_from_template(
				name=f'parse_{base_name}',
				template=parse_sumstats(
						inputfile=file
				)
	)

	b =  gwf.target_from_template(
				name=f'PGS_model_{base_name}',
				template=compute_pgs(
						inputfile=a.outputs[0]
				)
	)
