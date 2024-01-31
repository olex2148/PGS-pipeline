#!/usr/bin/env python

"""Modpath function

A function to split path strings into path, basename and suffix for easy manipulation of input and output names in the workflow.

"""

import os, re

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
