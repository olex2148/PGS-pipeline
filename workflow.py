#!/usr/bin/env python

"""Workflow to parse and compute polygenic scores

gwf workflow to parse all GWAS sumstats in data folder and generate PRS from them 
using LDpred2-auto and lassosum as implemented in LDpred2

@Author: Ole S Hansen
@Date: 

"""

from gwf import Workflow, AnonymousTarget
from datetime import date
from code.aux.modpath import modpath
from code.aux.input_paths import work_dir, sumstat_folder

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
	working_dir = work_dir
	options = {
		'memory': '20g',
		'walltime': '00:10:00',
		'cores': '23'
	}
	
	# Command to be run in the terminal
	spec = f'''
	
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
	
	inputs = [inputfile]
	outputs = [f'{base_path}_raw_models.rds', f'{base_path}_scores.rds', f'{base_path}_auto_parameters.rds']
	working_dir = work_dir
	options = {
		'memory': '30g',
		'walltime': '02:00:00',
		'cores': '23'
	}

	spec = f'''

	mkdir -p results/foelgefiler/{today}

	Rscript code/pgs_model.R {inputfile} {base_path}
	
	'''
	
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

### Defining targets

# Naming functions
def get_munge_name(idx, target):
  filename = modpath(target.inputs[0], parent='', suffix='')
  return f'munge_{filename}'

def get_pgs_name(idx, target):
  filename = modpath(target.inputs[0], parent='', suffix=('_munged.rds', ''))
  return f'ldpred2_{filename}'

# Input
sumstats = gwf.glob(sumstat_folder)

# Mapping over the inputs
parse_sumstats = gwf.map(munge_sumstats, sumstats, name=get_munge_name)
compute_scores = gwf.map(compute_pgs, parse_sumstats.outputs, name=get_pgs_name)

