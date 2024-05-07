#!/usr/bin/env python

"""Workflow to parse and compute polygenic scores

gwf workflow to parse all GWAS sumstats in data folder and generate PRS from them 
using LDpred2-auto and lassosum as implemented in LDpred2

@Author: Ole S Hansen
@Date: 

"""
import os, json
from gwf import Workflow, AnonymousTarget
from gwf.workflow import collect
from datetime import date
from code.functions.modpath import modpath

gwf = Workflow()
paths = json.load(open('data/paths.json'))

### Defining gwf templates
def munge_sumstats(inputfile):
	'''
	Template for running the r script "munge_sumstats.R" which munges GWAS summary statistics and performs quality control
	'''

	# Name of folder to be created in steps and results
	folder_name = os.path.split(inputfile)[0].split("/")[-1]
	
	# Using modpath() to create name of output file from inputfile - keeps basename,
	# but gets another path and another suffix
	base_name = modpath(inputfile, parent='', suffix='')
	
	res_folder = f'results/{folder_name}/{base_name}'
	munged_folder = f'steps/munged_sumstats/{folder_name}'
	
	munged_sumstats = f'{munged_folder}/{base_name}_munged.rds'
	
	# Defining inputs, outputs and ressources
	inputs = {'sumstats': inputfile}
	outputs = [munged_sumstats]
	working_dir = paths['work_dir']
	options = {
		'memory': '60g',
		'walltime': '00:30:00',
		'cores': '2',
		'account': 'NCRR-PRS'
	}
	
	# Command to be run in the terminal
	spec = f'''
	
 	mkdir -p {res_folder}
 	mkdir -p {munged_folder}

	Rscript code/munge_sumstats.R {inputfile} {munged_sumstats} {res_folder}
	
	'''

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def compute_pgs(inputfile):
	'''
	Template for running the r script "pgs_model.R" which computes PGS models using munged sumstats
	'''
	# Name of folder to be created in steps and results
	folder_name = os.path.split(inputfile)[0].split("/")[-1]

	base_name = modpath(inputfile, parent=(''), suffix=('_munged.rds', ''))             # Getting the base name from the inputfile 
	output_path = f'results/{folder_name}/{base_name}/{base_name}'                      # New path with sumstat-specific folder (and filename without suffix)
	model_info = f'{output_path}_model_info.csv'

	working_dir = paths['work_dir']
	inputs = {'munged_sumstats': inputfile}
	outputs = {
		'models': f'{output_path}_raw_models.rds', 
		'scores': f'{output_path}_scores.rds', 
		'auto_parameters': f'{output_path}_auto_parameters.rds',
		'model_info': model_info
	}
	options = {
		'memory': '40g',
		'walltime': '12:00:00',
		'cores': '23',
		'account': 'NCRR-PRS'
	}

	spec = f'''

	Rscript code/pgs_model.R {inputfile} {output_path} {model_info}
	
	'''
	
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def concat_files(paths):
  folder_name = os.path.split(paths[0])[0].split('/')[-2]
  output_path = f'results/{folder_name}/{folder_name}_model_info.csv'
  list_str = ' '.join(paths)
  command_string = "awk 'FNR == 1 && NR != 1 {next} {print}'"
  
  inputs = {'paths': paths}
  outputs = {'concatenated_model_info': output_path}
  options = {
    'memory': '2g',
    'walltime': '00:05:00'
  }
  
  spec = f'''
  
  {command_string} {list_str} > {output_path}
  
  '''
  
  return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

### Defining targets

# Naming functions
def get_munge_name(idx, target):
  filename = modpath(target.inputs['sumstats'], parent='', suffix='')
  return f'munge_{filename}'

def get_pgs_name(idx, target):
  filename = modpath(target.inputs['munged_sumstats'], parent='', suffix=('_munged.rds', ''))
  return f'ldpred2_{filename}'

# Input
sumstats = gwf.glob(paths['sumstat_folder'])

# Mapping over the inputs
munge_sumstats = gwf.map(munge_sumstats, sumstats, name=get_munge_name)
compute_pgs = gwf.map(compute_pgs, munge_sumstats.outputs, name=get_pgs_name)

gwf.target_from_template(
  name='concatenate_model_info',
  template=concat_files(
    paths=collect(compute_pgs.outputs, ['model_info'])['model_infos']
  )
)
