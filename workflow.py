## GWF workflow


'''
Author: Ole S Hansen
gwf workflow to parse all GWAS sumstats in data folder and generate PRS from them using LDpred2-auto
'''

from gwf import Workflow, AnonymousTarget
import os, re, glob
gwf = Workflow()

### Helper functions
def modpath(p, parent=None, base=None, suffix=None):
  '''
  A function to split path strings into path, basename and suffix
  for easy manipulation of input and output names in the workflow.

  Rename these..
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
Rscript code/prepare_and_parse_sumstats.R {inputfile} {outputfile}
'''

return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)










