#!/usr/bin/env python

"""Workflow to parse and compute polygenic scores

gwf workflow to parse all GWAS sumstats in data folder and generate PRS from them 
using LDpred2-auto and lassosum as implemented in LDpred2

AU sumary version

@Author: Ole S Hansen, Ainhoa Sanchez
@Date: 19/11/26

"""
import os, json
from gwf import Workflow, AnonymousTarget
from gwf.workflow import collect
from datetime import date
from code.functions.modpath import modpath
from pathlib import Path


gwf = Workflow()
paths = json.load(open('data/paths.json'))
run_model_qc = paths.get("run_model_qc", True)

###Code to use the Sumstats Library. This can be deactivated by applying "run_model_qc"=false (Noa)

def read_documentation(inputfile):
  #Read the file
  from openpyxl import load_workbook
  workbook = load_workbook(library_file, read_only=True, data_only=True, )
  
  try:
    #Check to get the PGS_code that uses ldpred2 for that Trait_code
    rows = workbook["04_PGS"].iter_rows(values_only=True)
    headers = next(rows)
    
    trait_idx = headers.index("Trait_code")
    pgs_idx = headers.index("PGS_code")
    method_idx = headers.index("PGS_method")
    
    trait_to_pgs = {}
    for row in rows:
      method = row[method_idx]
      if str(method).strip().casefold() != "ldpred2-auto": continue #We will only pick the pgs code if the method is ldpred2
      
      trait = row[trait_idx]
      pgs = row[pgs_idx]
      if trait is None or str(trait).strip() in ("", "NA"): continue #in case of errors putting the name
      if pgs is None or str(pgs).strip() in ("", "NA"): continue
      
      trait = str(trait).strip()
      pgs = str(pgs).strip()
      if trait in trait_to_pgs and trait_to_pgs[trait] != pgs: raise ValueError( f"Multiple LDpred-auto PGS_code values for {trait}")
      trait_to_pgs[trait] = pgs
          
    #File_name to Trait_code to PGS_code
    rows = workbook["02_Sumstats"].iter_rows(values_only=True)
    headers = next(rows)
    
    filename_idx = headers.index("File_name")
    trait_idx = headers.index("Trait_code")
    n_eff_idx = headers.index("N_eff")
    n_sumstats_idx = headers.index("N_sumstats")
    
    pgs_codes = {}
    seen_pgs_codes = set()
    sample_sizes = {}
  
    for row in rows:
      filename = row[filename_idx]
      trait = row[trait_idx]
    
      if filename is None or str(filename).strip() in ("", "NA"): #in case of errors putting the name
          continue
      filename = str(filename).strip()
      sample_sizes[filename] = (row[n_eff_idx], row[n_sumstats_idx])
      
      if trait is None or str(trait).strip() in ("", "NA"): 
          raise ValueError(f"Missing Trait_code for {filename}")
      trait = str(trait).strip()
      if trait not in trait_to_pgs: continue
      
      pgs = trait_to_pgs[trait]
      if filename in pgs_codes:
          raise ValueError(f"Duplicate File_name: {filename}")
      if pgs in seen_pgs_codes:
          raise ValueError(f"Duplicate PGS_code: {pgs}")
      pgs_codes[filename] = pgs
      seen_pgs_codes.add(pgs)
  finally:
    workbook.close()
    
  return pgs_codes, sample_sizes
    
if run_model_qc:
  au_sumstats_library = Path(paths["au_sumstats_library"])
  library_dir = au_sumstats_library / "supporting_info_sumstats"
  library_name = "Sumstats_PGS_iPSYCH_documentation.xlsx"
  library_file = max(library_dir.glob(f"????_??_??_{library_name}"), default=library_dir / library_name) #Takes the last library date, if not the original one.
  pgs_codes, sample_sizes = read_documentation(library_file)



### Defining gwf templates
def munge_sumstats(inputfile):
	'''
	Template for running the r script "munge_sumstats.R" which munges GWAS summary statistics and performs quality control
	'''

	# Name of folder to be created in steps and results
	folder_name = os.path.split(inputfile)[0].split("/")[-1]
	documentation = paths.get("run_model_qc", True)
	documented_n_eff = "NA"
	documented_n = "NA"
	
	#Change the name of the file for the PGS code. Noa
	if documentation:
	    base_name = pgs_codes[os.path.basename(inputfile)]
	    
	    #Use the n_eff or n_sumstats from the documentation if n_eff is not in the sumstats	    
	    n_eff_value, n_value = sample_sizes[os.path.basename(inputfile)]
	    if n_eff_value is not None: documented_n_eff = str(n_eff_value).strip() or "NA"
	    if n_value is not None: documented_n = str(n_value).strip() or "NA"
	    
	else:
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
		'memory': '100g',
		'walltime': '01:00:00',
		'cores': '2'
	}
	
	# Command to be run in the terminal
	spec = f'''
	
 	mkdir -p {res_folder}
 	mkdir -p {munged_folder}

	Rscript code/munge_sumstats.R {inputfile} {munged_sumstats} {res_folder} {documented_n_eff} {documented_n} {str(documentation).lower()}
	
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
		'cores': '23'
	}

	spec = f'''
	

	Rscript code/pgs_model.R {inputfile} {output_path} {model_info}
	
	'''
	
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

#Noa's code to use the QC into de model_info files and add the values to "AU_sumstats_library/supporting_info_sumstats/Sumstats_PGS_iPSYCH_documentation.xlsx". 
#It can be activated or deactivated in paths.json: "run_model_qc"


def check_model_qc(inputfile):
	'''
	Template for running the r script "check_parameters.R" which adds the QC values to "model_info.csv"
	'''
	folder_name = os.path.split(inputfile)[0].split("/")[-1]
	output_path = modpath(inputfile, suffix =("_model_info.csv", "_model_info_qc.csv"),) 
	model_info = f'{output_path}_model_info.csv'

	inputs = {'model_info': inputfile}
	outputs = {'model_info_qc': output_path}
	options = {
    'memory': '2g',
    'walltime': '00:05:00'
    }   
    
	spec = f'''

	Rscript code/check_parameters.R {inputfile} {output_path}
	
	'''
	
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def update_qc_library(inputfile):
	'''
	Template for running the r script "check_parameters.R" which adds the QC values to "model_info.csv and make a copy of the results and cleaned(munged) summary statistics into the AU_Sumstats_library folder
	'''
	pgs_code = modpath(qc_file, parent="", suffix=("_model_info_qc.csv", ""))
	group = Path(qc_file).parent.parent.name
	
	source_folder = Path(qc_file).parent
	scores_file = str(Path(qc_file).parent / f"{pgs_code}_scores.rds")
	munged_file = (f"steps/munged_sumstats/{group}/{pgs_code}_munged.rds")
	scores_dir = au_sumstats_library / "scores"
	cleaned_dir = au_sumstats_library / "cleaned_sumstats"
	scores_copy = str(scores_dir/ pgs_code / f"{pgs_code}_scores.rds")
	munged_copy = str(cleaned_dir / f"{pgs_code}_munged.rds")
  
	run_date = date.today().strftime("%Y_%m_%d")
	#the library id locked untill is rewrite it
	marker = modpath(inputfile, suffix=("_model_info_qc.csv", f"_{run_date}_library_updated.txt"))
	lock_file = str(library_dir / "qc_library.lock")

	inputs = {'model_info_qc': inputfile,
	        'scores': scores_file,
            'munged_sumstats': munged_file}
	outputs = {'updated_library': marker, 
	        'library_scores': scores_copy,
	        'library_munged_sumstats': munged_copy}
	options = {
    'memory': '2g',
    'walltime': '00:10:00'
    }   
    
	spec = f'''
	set -e
	cp -r {source_folder} {scores_dir}/
	cp {munged_file} {munged_copy}
	
	flock {lock_file}	python code/functions/update_qc_library.py {inputfile} {library_dir} {library_name} {run_date}
	touch {marker}
	'''
	
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)




def concat_files(paths, run_model_qc=False):
  folder_name = os.path.split(paths[0])[0].split('/')[-2]
  output_path = f'results/{folder_name}/{folder_name}_model_info.csv'
  output_path_qc = f'results/{folder_name}/{folder_name}_model_info_qc.csv' #Noa,QC file
  warn_fail_path = f'results/{folder_name}/{folder_name}_model_info_qc_warn_fail.csv' #Noa, QC file with FAIL/WARNING
  list_str = ' '.join(paths)
  command_string = "awk 'FNR == 1 && NR != 1 {next} {print}'"
  
  inputs = {'paths': paths}
  options = {
    'memory': '2g',
    'walltime': '00:05:00'
  }
  
  #Add if/else 
  if run_model_qc:
    output_path_qc = f'results/{folder_name}/{folder_name}_model_info_qc.csv' #Noa,QC file
    warn_fail_path = f'results/{folder_name}/{folder_name}_model_info_qc_warn_fail.csv' #Noa, QC file with FAIL/WARNING
    outputs = {'concatenated_model_info': output_path_qc,
               'warn_fail_model_info': warn_fail_path}
    #Look into the last collumn and filter only WARN and FAIL
    filter_command = """awk -F '\\t' 'NR == 1 || $NF ~ /^(FAIL|WARN)/'"""
  
  
    spec = f'''
      
    {command_string} {list_str} > "{output_path_qc}"
    {filter_command} "{output_path_qc}" > {warn_fail_path}
      
    '''
  
  else:  
    outputs = {'concatenated_model_info': output_path}
  
  
    spec = f'''
  
    {command_string} {list_str} > {output_path}
  
    '''
  
  return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

## End of addition




### Defining targets

# Naming functions
def get_munge_name(idx, target):
  if paths.get("run_model_qc", True):
    filename = pgs_codes[os.path.basename(target.inputs["sumstats"])]
  else:
    filename = modpath(target.inputs['sumstats'], parent='', suffix='')
  return f'munge_{filename}'

def get_pgs_name(idx, target):
  filename = modpath(target.inputs['munged_sumstats'], parent='', suffix=('_munged.rds', ''))
  return f'pgs_{filename}'
  
def get_qc_name(idx, target):
  filename = modpath(target.inputs['model_info'], parent='', suffix=('_model_info.csv', ''))
  return f'model_qc_{filename}'

# Input
sumstats = gwf.glob(paths['sumstat_folder'])

# Mapping over the inputs
munge_sumstats = gwf.map(munge_sumstats, sumstats, name=get_munge_name)
compute_pgs = gwf.map(compute_pgs, munge_sumstats.outputs, name=get_pgs_name)

#Noa code to use the QC
model_info_files = collect( compute_pgs.outputs, ["model_info"], )["model_infos"]
files_to_concat = model_info_files
if run_model_qc:
  qc_targets = gwf.map(check_model_qc, model_info_files, name=get_qc_name,)
  qc_files = collect(qc_targets.outputs, ['model_info_qc'])['model_info_qcs']
  for qc_file in qc_files: 
    pgs_code = modpath(qc_file, parent="", suffix=("_model_info_qc.csv", "")) 
    gwf.target_from_template(name=f"update_library_{pgs_code}", template=update_qc_library(inputfile=qc_file))
  files_to_concat = qc_files

### end of addition

    
concat_target = gwf.target_from_template(
  name='concatenate_model_info',
  template=concat_files(
    paths=files_to_concat, run_model_qc=run_model_qc, #Ole's: paths=collect(compute_pgs.outputs, ['model_info'])['model_infos']
  )
)
