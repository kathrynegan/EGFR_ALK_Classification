# -*- coding: utf-8 -*-

"""
author@esilgard

Copyright (c) 2015-2017 Fred Hutchinson Cancer Research Center

Licensed under the Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0

Refactored by: kathrynegan

gentest_classifier.py was adapted from the following module(s) in esilgard's master branch:
make_arrays.py
decoder.py
"""
import os
import numpy as np
from scipy.sparse import dok_matrix
from sklearn.externals import joblib


class GenTestClassifier:

	def __init__(self, model_dir):
		""" Initializes GenTestClassifier instance. """
		self.algorithms = {}
		for algorithm in os.listdir(model_dir):
			model = Model()
			with open(os.path.join(model_dir, algorithm, 'features.txt'), 'r') as f:
				for line in f.readlines():
					feature, index = line.split()
					index = int(index)
					model.mapping[feature] = index
					model.num_features = max(model.num_features, index + 1)
			# joblib is json for large, sparse numpy arrays
			model.model = joblib.load(
				os.path.join(model_dir, algorithm, 'model.pkl'))
			self.algorithms[algorithm] = model

	def classify(self, vector):
		""" Returns the label for each of results reported,
		results positive, and method of test for the given vector.
		Args:
			vector (list of str) : list of features as strings in vector
		Returns:
			reported (str) :
				'Results Reported' if results of a gen test were stated in this path report
				'Not Reported' otherwise
			result (str) :
				'Positive' if the results were a positive test
				'Negative' if the results were a negative test
				'N/A' if no results were reported or if the results were ambiguous
			method (str) : the test method that was used to obtain results (e.g. FISH)
		"""
		reported, result, method = 'Not Reported', 'N/A', 'N/A'
		if 'NO_KEYWORD_IN_TEXT' in vector:
			return reported, result, method
		if self._classify('svm_reported', vector) != 'Reported':
			return reported, result, method
		reported = 'Results Reported'
		result = self._classify('positive', vector)
		method = self._classify('method', vector)
		return reported, result, method

	def _classify(self, algorithm, vector):
		""" Returns the label for this vector given the algorithm.
		Args:
			algorithm (str) : name of algorithm
			vector (list of str) : list of features as strings in vector
		Returns:
			str : instance label
		"""
		# instantiate empty numinstances x numfeatures matrix
		model = self.algorithms[algorithm]
		matrix = dok_matrix((1, model.num_features), dtype=np.float64)
		# populate matrix with binary values representing features in instance
		for feature in vector:
			try:
				matrix[0, model.mapping[feature]] = 1
			except KeyError:
				continue
		matrix.tocsc()
		output = model.model.predict(matrix)[0]
		output = self._translate_output(algorithm == 'method', output)
		return output

	def _translate_output(self, is_method, output):
		""" Converts the output of the SVM to a meaningful label.
		Args:
			is_method (bool) :
				True if the SVM output is for a test method, False otherwise
			output (int) : output of SVM
		Returns:
			str : SVM output mapped to label
		"""
		method_mapping = {
			0: 'Mutational Analysis',
			1: 'IHC',
			2: 'FISH',
			3: 'OTHER',
			4: 'NONE'}
		result_mapping = {
			0: 'Unknown',
			1: 'Positive',
			2: 'Negative',
			3: 'Insufficient',
			4: 'Reported',
			5: 'Not Reported'}
		if is_method:
			return method_mapping[output]
		return result_mapping[output]


class Model:
	""" Tiny class to conveniently store a model and its associated data. """

	def __init__(self):
		self.model = None
		self.mapping = {}
		self.num_features = 0

	@property
	def model(self):
		return self.model

	@model.setter
	def model(self, model):
		self.model = model

	@property
	def mapping(self):
		return self.mapping

	@mapping.setter
	def mapping(self, mapping):
		self.mapping = mapping

	@property
	def num_features(self):
		return self.num_features

	@num_features.setter
	def num_features(self, num_features):
		self.num_features = num_features
