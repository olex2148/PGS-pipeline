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
from code.aux.modpath import modpath

gwf = Workflow()

### Defining gwf templates
def munge_sumstats(inputfile):
	'''
	Template for running the r script "prepare_and_parse_sumstats.R" which parses GWAS summary statistics
	'''

	# Using modpath() to create name of output file from inputfile - keeps basename,
	# but gets another path and another suffix
	munged_sumstats = modpath(inputfile, parent=('steps/munged_sumstats'), suffix=('_munged.rds'))
	
	base_name = modpath(inputfile, parent=(''), suffix=(''))      # Getting the base name from the inputfile 

	# Defining inputs, outputs and ressources
	inputs = [inputfile]
	outputs = [munged_sumstats]
	options = {
		'memory': '50g',
		'walltime': '01:00:00',
		'cores': '23'
	}
	
	# Command to be run in the terminal
	spec = f'''
	
	cd ~/NCRR-PRS/faststorage/osh/PGS/pgs_workflow/
	
	mkdir -p results/{base_name}

	Rscript code/munge_sumstats.R {inputfile} {base_name}
	
	'''

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def compute_pgs(inputfile):
	'''
	Template for running the r script "pgs_model.R" which computes PGS models using parsed sumstats
	'''
	today = date.today()

	base = modpath(inputfile, parent=(''), suffix=('_munged.rds', ''))      # Getting the base name from the inputfile 
	base_path = f'results/{base}/{base}'                                    # New path with sumstat-specific folder (and filename without suffix)

	model_out = f'{base_path}_raw_models.rds'                               # Models, parameters, and scores are put in a folder specific to the sumstats
	scores_out = f'{base_path}_scores.rds'
	parameters_out = f'{base_path}_auto_parameters.rds'
	# foelgefil = f'results/foelgefiler/{today}/{base}_foelgefil.xlsx'          # Foelgefil is put in a separate folder specific to the batch run (date in folder name)

	inputs = [inputfile]
	outputs = [model_out, scores_out, parameters_out]
	options = {
		'memory': '60g',
		'walltime': '12:00:00',
		'cores': '23'
	}

	spec = f'''

	mkdir -p results/foelgefiler/{today}

	Rscript code/pgs_model.R {inputfile} {base_path}
	
	'''
	
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

### Defining targets

# Input
sumstats = glob.glob("data/ipsych_test/*")

# Submitting jobs for every sumstat file in data folder (parsing and computing PGSs)
for file in sumstats:
	base_name = modpath(file, parent='', suffix='')

	a = gwf.target_from_template(
				name=f'parse_{base_name}',
				template=munge_sumstats(
						inputfile=file
				)
	)

	b =  gwf.target_from_template(
				name=f'ldpred2_{base_name}',
				template=compute_pgs(
						inputfile=a.outputs[0]
				)
	)
