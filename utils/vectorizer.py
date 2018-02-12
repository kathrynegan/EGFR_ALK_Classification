"""
CREATED BY: Kathryn Egan

"""
import re
import os
import json


class Vectorizer:
	""" Stores data and behavior for cleaning and processing a
	pathology report as part of  EGFR/ALK classification. """

	def __init__(self):
		""" Initializies Vectorizor instance by compiling regexes. """
		self.test_patterns = self._compile_patterns(
			'condensed_patterns.json', True, r'[\W\^]', r'[\W$]')
		self.other_patterns = self._compile_patterns(
			'other_kw_patterns.json', False, r'[\W\^]', r'[\W$]')
		self.section_patterns = self._compile_patterns(
			'section_patterns.json', True, r'^', r'$')
		self.substitutions = self._compile_substitutions()
		self.cytology_pattern = re.compile(
			r'(cytoprep)|(cytolog)', flags=re.IGNORECASE)
		self.insufficient_pattern = re.compile(
			r'insufficient (tumor|sample)?|technical diff}(tumor|sample) insufficient',
			flags=re.IGNORECASE)
		self.accession_pattern = re.compile(
			r'[\W]([\(]?[A-Z]{1,2}[\- ]?[\d]{2,4}[\- ]{1,3}[\d]{2,8}[\)]?)[\W]')
		self.stop_list = re.compile(
			r'[\s\^](TO|THE|FOR|A|AN|AS|THIS|THAT|THESE|' +
			r'THEY|IN|OF|ON|OR|BY)( THE|A|AN)?[\s\$]')

	def _compile_patterns(self, file_name, uppercase, pre, post):
		"""
		Helper method to compile patterns for each regex dictionary with
		appropriate "cushions" before and after primary capture groups.
		Allows for optomizing pattern matching while allowing for slightly
		better readibility in json docs. Only uppercase all patterns
		(based on boolean flag) as long as no regex character classes
		are used in pattern e.g. we don't want [\w] to turn into [\W].
		Args:
			file_name (str) : name of file to read in
			uppercase (bool) : add uppercased pattern
			pre (str) : regex to capture before pattern
			post (str) : regex to capture after pattern
		Returns:
			dict str:regex: key mapped to regex
		"""
		home = os.path.dirname(os.path.normpath(os.path.realpath(__file__)))
		file = os.path.join(home, 'patterns', file_name)
		pattern_dict = json.load(open(file, 'r'))  # test patterns
		compiled = []
		for subin, pattern_list in pattern_dict.items():
			subin = subin.replace('<newline>', '\n')
			for pattern in pattern_list:
				# make sure match pattern is isolated from alphanumeric characters
				subout1 = r'{}({}){}'.format(pre, pattern, post)
				subout1 = re.compile(subout1, re.MULTILINE)
				compiled.append((subout1, subin))
				if uppercase:
					subout2 = r'{}({}){}'.format(pre, pattern.upper(), post)
					subout2 = re.compile(subout2, re.MULTILINE)
					compiled.append((subout2, subin))
		return compiled

	def _compile_substitutions(self):
		"""
		Map string to replace to replacement string, ordered by hierarchy.
		- condense coordinated test instance
		- (eg test 1 and test 2 are pending)
		- we dont want to break on test2
		- strip out long trails of  ____ -
		- there ARE single underscores IN the standardized features
		Returns:
			dict (index: str: str) : dictionary of index mapped to
				replacement string mapped to string to replace
		"""
		replacements = [
			'OTHER_TEST', 'PUBLICATION', 'TEST_INSTANCE', 'IHC', 'PATHOLOGIST',
			'BLOCK_ACC', 'SPECIFIC_MUT', 'MUT_ANALYSIS', 'FISH', 'AUTHOR']
		substitutions = {}
		for index, replacement in enumerate(replacements):
			subin = r' {} '.format(replacement)
			subout = re.compile(
				r'{}[,.\(\):\-;andor{} ]{{1,}}{}'.format(*[replacement] * 3))
			substitutions[index] = (subout, subin)
		substitutions[index + 1] = (re.compile(
			r'TEST_INSTANCE[,.\(\):\-;andor ]{1,}OTHER_TEST'), 'TEST_INSTANCE')
		substitutions[index + 2] = (re.compile(
			r'[0-9]{2}[\-\\\/][0-9]{2}[\-\\\/][0-9]{2,5}'), ' DATE ')
		substitutions[index + 3] = (re.compile(
			r'($|\s)[A-H][\)]?[.]'), ' SPECIMEN_LABEL ')
		substitutions[index + 4] = (re.compile(
			r'[\"\(\\\)\-\/\']'), ' ')
		substitutions[index + 5] = (re.compile(
			r'[.,;:\?]'), ' PUNCTUATION ')
		substitutions[index + 6] = (re.compile(
			r'[\[\]]'), ' ')
		return substitutions

	def make_vector(self, text, accession, marker):
		""" Initial/main method for vector creation. Requires a file that
		(at minimum) contains the unique id of the report (instance), the
		pathology report accession number and the raw text of the pathology
		report. Assumed to be a tab delimited file.
		"""
		self.accession = re.sub(r'[\- ]', '', accession)
		self.text = self._get_text(text)
		self.vector = []
		self.marker = marker
		self._cytology_report()
		self._positive_test()
		self._test_instance()
		self._standardize()
		self._other_accession()
		self._insufficient()
		self._substitute()
		self._stop_list()
		self._ngrams()
		self._test_mentions()
		return self.vector

	def _get_text(self, text):
		""" Returns ascii-only version of text. Subs non-ascii characters
		with white space and truncates multiple sequential space characters
		to one space character.
		Args:
			text (str) : text to remove non-ascii characters from
		Returns:
			str : processed text
		"""
		ascii_only = []
		for char in text:
			if ord(char) >= 128:
				char = ' '
			ascii_only.append(char)
		ascii_only = ''.join(ascii_only)
		ascii_only = re.sub(r' +', ' ', ascii_only)
		text = re.sub(r' +', ' ', ascii_only)
		text = text.replace('<newline>', '\n')
		return text

	def _cytology_report(self):
		""" Adds whether report is cytology related.
		Cytology reports unlikely to have reliable tests. """
		result = 1 if self.cytology_pattern.search(self.text) else 0
		self.vector.append('CYTO_RELATED_REPORT')

	def _positive_test(self):
		""" Adds existence of one or more positive genetic tests to vector. """
		for pattern, test in self.test_patterns:
			if test != self.marker:
				continue
			# the '+' get lost in the [\W] buffer around tests - pull them out
			# not pulling out '-', since it's ambiguous; minus or just a dash?
			# strips off the "cushion" from the end of the test instance pattern
			trunc_pattern = re.compile(pattern.pattern[:-5] + r'[\s]*[\+]')
			if trunc_pattern.search(self.text):
				self.vector.append('post_window=POSITIVE')
				self.vector.append('post_window=TEST_INSTANCE_POSITIVE')
				break

	def _test_instance(self):
		""" Replaces instances of a genetic test with generic placeholder. """
		for pattern, test in self.test_patterns:
			if test != self.marker:
				continue
			self.text = pattern.sub(' TEST_INSTANCE ', self.text)

	def _standardize(self):
		""" Standardizes patterns in report. """
		# kinda gross hack - this runs through twice to catch overlapping patterns
		# (because of [\W] buffer in pattern match)
		for i in range(2):
			# replace all other tests in the text
			for subout, subin in self.test_patterns:
				self.text = subout.sub(' OTHER_TEST ', self.text)
			for subout, subin in self.section_patterns:
				self.text = subout.sub(' {} '.format(subin), self.text)
			for subout, subin in self.other_patterns:
				self.text = subout.sub(' {} '.format(subin), self.text)

	def _other_accession(self):
		""" Adds whether other/previous pathology reports are mentioned to vector,
		subs accession numbers out of report text. """
		# find all accession numbers
		for match in self.accession_pattern.finditer(self.text):
			accession = re.sub(r'[()\- ]', '', match.group(1))
			# string to sub if accession number matches the one for this report
			if accession == self.accession:
				subin = ' THIS_ACC_NUM '
			else:
				subin = ' OTHER_ACC_NUM '
			# escape parentheses
			subout = match.group(1)
			subout = re.sub(r'\(', '\\(', subout)
			subout = re.sub(r'\)', '\\)', subout)
			self.text = re.sub(subout, subin, self.text)
		self.vector.append('OTHER_ACC_NUM_IN_TEXT')

	def _insufficient(self):
		""" Adds presence of insufficient samples to
		vector EVEN when test name isn't mentioned. """
		result = 1 if self.insufficient_pattern.search(self.text) else 0
		self.vector.append('INSUFFICIENT')

	def _substitute(self):
		""" Makes iterative substitutions in text. """
		for index in sorted(self.substitutions):
			subout, subin = self.substitutions[index]
			self.text = subout.sub(subin, self.text)

	def _stop_list(self):
		""" Removes stop list items from report. Run twice to catch
		downstream patterns. """
		self.text = self.text.upper()
		self.text = re.sub(self.stop_list, ' ', self.text)
		self.text = re.sub(self.stop_list, ' ', self.text)

	def _ngrams(self):
		""" Adds ngrams to vector. """
		tokens = self.text.strip().split()
		for index, token in enumerate(tokens):
			if token != 'TEST_INSTANCE':
				continue
			self.vector.append(self.marker)
			self._add_section(tokens, index)
			start, end = self._get_window(tokens, index)
			if index > 0 and index > start:
				self.vector.append('immediately_pre_window=' + tokens[index - 1])
			for i in reversed(range(start, index)):
				self.vector.append('pre_window=' + tokens[i])
				for j in 1, 2, 3:
					# if statement contains index per original
					if index - j > start:
						self.vector.append('pre_window={}_{}'.format(tokens[i - j], tokens[i]))
			if index < len(tokens) - 1 and index < end - 1:
				self.vector.append('immediately_post_window=' + tokens[index + 1])
			for i in range(index + 1, min(len(tokens), end)):
				self.vector.append('post_window=' + tokens[i])
				for j in 1, 2, 3:
					if i < end - j:
						self.vector.append('post_window={}_{}'.format(tokens[i], tokens[i + j]))

	def _get_window(self, tokens, index):
		""" Returns start and end indexes for window around given index.
		Args:
			tokens (list of str) : tokenized pathology report
			index (int) : current index
		Returns:
			start (int) : start index for this window
			end (int) : end index for this window
		"""
		window_start = max(index - 10, 0)
		window_end = min(len(tokens), index + 10)
		window = tokens[window_start:window_end]
		# break window size for new sections or other tests, etc
		breaks = [
			i for i, t in enumerate(window)
			if t in ('_SECTION_', 'PUNCTUATION', 'SPECIMEN_LABEL', 'OTHER_TEST')]
		breaks.append(0)
		start = window_start + max([i + 1 for i in breaks if i < len(window) / 2])
		breaks[-1] = len(window)
		# start is applied to end per original
		end = window_start + min([i for i in breaks if i > len(window) / 2])
		return start, end

	def _add_section(self, tokens, index):
		""" Adds section name to vector.
		Args:
			tokens (list of str) : tokenized pathology report
			index (int) : current index
		"""
		# look backwards from test instance to find name of section
		prev = ''
		for curr in reversed(tokens[:index]):
			if curr == '_SECTION_':
				self.vector.append('SECTION=' + prev)
				break
			prev = curr

	def _test_mentions(self):
		""" Adds number of test mentions to vector. """
		test_mentions = self.text.count('TEST_INSTANCE')
		self.vector.append('COUNT_TEST_INSTANCE')
		if not test_mentions:
			self.vector.append('NO_KEYWORD_IN_TEXT')
