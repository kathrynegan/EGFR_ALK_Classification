# -*- coding: utf-8 -*-

"""
author@esilgard

Copyright (c) 2015-2017 Fred Hutchinson Cancer Research Center

Licensed under the Apache License, Version 2.0:
http://www.apache.org/licenses/LICENSE-2.0

Refactored by: kathrynegan

run.py was adapted from the following module(s) in esilgard's master branch:
svm_pipeline.py
final_output.py
"""
import os
import sys
import csv
from utils.gentest_classifier import GenTestClassifier
from utils.vectorizer import Vectorizer


TEXT = 'full_path_text'  # name of pathology report field
ACC = 'accession_number_hosp'  # name of accession number field
PAT = 'patient_id'  # name of patient ID field
TUMOR = 'tumor_record'
REC = 'source_id'  # name of record ID field
TOTAL = 20000  # total number of records (estimated)
MARKERS = ['EGFR', 'ALK']  # markers to process

# pt_file = 'random_50_patients'
# rd_file = 'random_200_records'
# skip_file = 'skip_patients'
try:
	pt_subset = set(open(pt_file + '.txt').read().strip().split())
except NameError:
	pt_subset = set()
try:
	rd_subset = set(open(rd_file + '.txt').read().strip().split())
except NameError:
	rd_subset = set()
try:
	skip_set = set(open(skip_file + '.txt').read().strip().split())
except NameError:
	skip_set = set()


def run_pipeline():
	"""
	Pipeline for classification of EGFR and ALK test use, result, and method.
	Writes record level and patient level results to separate files.
	"""
	dirs = get_dirs()
	cases = process_records(dirs)
	process_patients(cases, dirs['case level'])


def get_dirs():
	""" Create or verify directories for models, input, and output.
	Assign paths for files.
	Returns:
		dict (str:str) : type of file mapped to file path
	"""
	dirs = {}
	home = os.path.dirname(os.path.normpath(os.path.realpath(__file__)))
	# get input file
	try:
		dirs['input'] = sys.argv[1]
	except IndexError:
		raise IndexError('Provide path to input file as first argument.')
	# establish output dir and files
	output_dir = os.path.join(home, 'output')
	try:
		os.mkdir(output_dir)
	except OSError:
		pass
	flag = '_' + rd_file if rd_subset else ''  # compile flag to mark output type
	flag += '_' + pt_file if pt_subset else ''
	flag += '_' + skip_file if skip_set else ''
	flag = 'all' if not flag else flag
	dirs['record level'] = os.path.join(
		output_dir, 'record_level_output{}.txt'.format(flag))
	dirs['case level'] = os.path.join(
		output_dir, 'case_level_output{}.txt'.format(flag))
	# establish model files
	model_dir = os.path.join(home, 'models')
	if not os.path.exists(model_dir):
		sys.stderr.write('Model directory not found.\nExiting...\n')
		sys.exit()
	dirs['model'] = model_dir
	return dirs


def process_records(dirs):
	""" Processes pathology reports, writes record-level results to file and
	iteratively determines patient level genetic testing status.
	Args:
		dirs (dict str:str) : type of file mapped to file path
	Returns:
		dict (str:str:(str, str)) :
			patient ID mapped to gen marker and status with deciding report ID
	"""
	classifier = GenTestClassifier(dirs['model'])
	vectorizer = Vectorizer()
	cases = {}
	with open(dirs['input'], 'r') as fin, open(dirs['record level'], 'w') as fout:
		reader = csv.reader(fin, delimiter='\t')
		headers = reader.next()
		row_length = len(headers)
		# write headers
		fout.write('\t'.join(
			headers[:get(headers, TEXT)] + headers[get(headers, TEXT) + 1:]))
		for marker in MARKERS:
			for cat in 'Reported', 'Result', 'Method':
				fout.write('\t{} {}'.format(marker, cat))
		fout.write('\n')
		sys.stderr.write('Log based on {} total records\n'.format(TOTAL))
		mark = 10
		num_processed = 0
		for row in reader:
			if len(row) != row_length:
				message = 'Differing row lengths detected. ' +\
					'Please check input data. [row index={}]\n'.format(reader.line_num)
				raise IOError(message)

			cases = process_row(fout, headers, row, cases, classifier, vectorizer)

			num_processed += 1
			percentage = num_processed * 100.0 / TOTAL
			if percentage > mark:
				sys.stderr.write('{}% of records processed...\n'.format(mark))
				mark += 10
	sys.stderr.write('100% of records processed\n')
	sys.stderr.write(
		'Record level results written to:\n{}\n'.format(dirs['record level']))
	return cases


def process_row(fout, headers, row, cases, classifier, vectorizer):
	""" Processes a single row in input file. Classifies report and
	resolves output, writing output to record-level file. Updates patient
	level result status.
	Args:
		fout (file) : open filestream to output file
		headers (list of str) : list of headers
		row (list of str) : row as list
		cases (dict str:str:(str, str)) :
			patient ID mapped to gen marker and status with deciding report ID
		classifier (Classifier) : classifier object
		vectorizer (Vectorizer) : vectorizing object
	Returns:
		cases (dict str:str:(str, str)) :
			patient ID mapped to gen marker and status with deciding report ID,
			updated according to result for this row
	"""
	accession = row[get(headers, ACC)]
	tumor = row[get(headers, TUMOR)]
	record = row[get(headers, REC)]
	patient = row[get(headers, PAT)]
	text = row.pop(get(headers, TEXT))
	case = '{}_{}'.format(patient, tumor)
	if pt_subset and patient not in pt_subset:
		return cases
	if rd_subset and record not in rd_subset:
		return cases
	if skip_set and case in skip_set:
		return cases
	fout.write('\t'.join(row))
	for marker in MARKERS:
		cases.setdefault(case, {}).setdefault(marker, ('Unknown', 'N/A'))
		vector = vectorizer.make_vector(text, accession, marker)
		reported, result, method = classifier.classify(vector)
		fout.write('\t' + '\t'.join([reported, result, method]))
		# take only first positive
		if cases[case][marker] == 'Positive':
			continue
		# take any ALK result or an EGFR result by mutational analysis
		if result == 'Positive':
			if marker == 'ALK' or method == 'Mutational Analysis':
				cases[case][marker] = (result, record)
		if cases[case][marker] == 'Negative':
			continue
		if result == 'Negative':
			if marker == 'ALK' or method == 'Mutational Analysis':
				cases[case][marker] = (result, record)
	fout.write('\n')
	return cases


def process_patients(cases, file):
	""" Write patient-level results to file.
	Args:
		cases (dict str:str:(str, str)) :
			patient ID and tumor ID mapped to gen marker
			and status with deciding report ID
		file (str) : path to case level output file
	"""
	with open(file, 'w') as f:
		f.write('patient_id')
		for marker in MARKERS:
			f.write('\t{} Result\t{} Record ID'.format(marker, marker))
		f.write('\n')
		for case in sorted(cases):
			f.write(case)
			for marker in MARKERS:
				f.write('\t{}\t{}'.format(*cases[case][marker]))
			f.write('\n')
	sys.stderr.write(
		'Case level results written to:\n{}\n'.format(file))


def get(headers, field):
	""" Returns index of the given data field name according to its
	positing in headers. If data field cannot be found, exits
	with error message.
	Args:
		headers (list of str) : list of headers
		field (str) : data field name to find in headers
	Returns:
		int : index of given data field
	"""
	try:
		return headers.index(field)
	except ValueError:
		sys.stderr.write(
			'ERROR: Data field {} not found in input.\n'.format(field))
		sys.exit()


if __name__ == "__main__":
	run_pipeline()
