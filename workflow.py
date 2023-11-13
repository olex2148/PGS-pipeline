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
