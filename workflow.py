#!/usr/bin/env python

"""Workflow to parse and compute polygenic scores

gwf workflow to parse all GWAS sumstats in data folder and generate PRS from them 
using LDpred2-auto and lassosum as implemented in LDpred2

@Author: Ole S Hansen
@Date: 

"""
import os, json
from gwf import Workflow, AnonymousTarget
from datetime import date
from code.aux.modpath import modpath

gwf = Workflow()
paths = json.load(open('data/paths.json'))

### Defining gwf templates
def munge_sumstats(inputfile):
	'''
	Template for running the r script "prepare_and_parse_sumstats.R" which parses GWAS summary statistics
	'''

	# Name of folder to be created in steps and results
	folder_name = os.path.split(inputfile)[0].split("/")[-1]
	
	# Using modpath() to create name of output file from inputfile - keeps basename,
	# but gets another path and another suffix
	munged_sumstats = modpath(inputfile, parent=(f'steps/munged_sumstats/{folder_name}'), suffix=('_munged.rds'))
	foelgefil = modpath(inputfile, parent=(f'results/foelgefiler/{folder_name}'), suffix=('_foelgefil.xlsx'))
	
	# Defining inputs, outputs and ressources
	inputs = [inputfile]
	outputs = [munged_sumstats, foelgefil]
	working_dir = paths['work_dir']
	options = {
		'memory': '20g',
		'walltime': '00:10:00',
		'cores': '2'
	}
	
	# Command to be run in the terminal
	spec = f'''
	
	mkdir -p results/foelgefiler/{folder_name}
 	mkdir -p steps/munged_sumstats/{folder_name}

	Rscript code/munge_sumstats.R {inputfile} {munged_sumstats} {foelgefil}
	
	'''

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def compute_pgs(inputfile, foelgefil):
	'''
	Template for running the r script "pgs_model.R" which computes PGS models using parsed sumstats
	'''
	# Name of folder to be created in steps and results
	folder_name = os.path.split(inputfile)[0].split("/")[-1]

	base_name = modpath(inputfile, parent=(''), suffix=('_munged.rds', ''))             # Getting the base name from the inputfile 
	output_path = f'results/{folder_name}/{base_name}/{base_name}'                      # New path with sumstat-specific folder (and filename without suffix)

	working_dir = paths['work_dir']
	inputs = [inputfile]
	outputs = [
		f'{output_path}_raw_models.rds', 
		f'{output_path}_scores.rds', 
		f'{output_path}_auto_parameters.rds'
	]
	options = {
		'memory': '10g',
		'walltime': '01:00:00',
		'cores': '23'
	}

	spec = f'''

 	mkdir -p results/{folder_name}/{base_name}

	Rscript code/pgs_model.R {inputfile} {output_path} {foelgefil}
	
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
sumstats = gwf.glob(paths['sumstat_folder'])

# Mapping over the inputs
parse_sumstats = gwf.map(munge_sumstats, sumstats, name=get_munge_name)
compute_scores = gwf.map(compute_pgs, parse_sumstats.outputs, name=get_pgs_name)

