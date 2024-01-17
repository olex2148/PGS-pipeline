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
import os, re, glob
from datetime import date
gwf = Workflow()

### Helper functions
def modpath(p, parent=None, base=None, suffix=None):
	'''
	A function to split path strings into path, basename and suffix
	for easy manipulation of input and output names in the workflow.
	'''
	par, name = os.path.split(p)
	name_no_suffix, suf = os.path.splitext(name)
	if type(suffix) is str:
		suf = suffix
	if parent is not None:
		par = parent
	if base is not None:
		name_no_suffix = base

	new_path = os.path.join(par, name_no_suffix + suf)
	if type(suffix) is tuple:
		assert len(suffix) == 2
		new_path, nsubs = re.subn(r'{}$'.format(suffix[0]), suffix[1], new_path)
		assert nsubs == 1, nsubs
	return new_path

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
sumstats = glob.glob('data/*')

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