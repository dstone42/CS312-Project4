#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random
import sys
sys.setrecursionlimit(10000)

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		pass

	def printMat(mat):
		# Print function to better visualize matrices
		for i in range(len(mat)):
			print(mat[i])

	def NucScore(self, aa1, aa2):
		# Return -5 for an indel
		if aa1 == "-" or aa2 == "-":
			return INDEL
		elif aa1 == aa2:
			return MATCH
		else:
			return SUB

	def score(self, seq1, seq2):
		# mat[seq1][seq2]
		# Initialize matrices
		scoreMat = [0] * (len(seq1) + 1)
		for i in range(len(seq1) + 1):
			scoreMat[i] = [0] * (len(seq2) + 1)
		backMat = [""] * (len(seq1) + 1)
		for i in range(len(seq1) + 1):
			backMat[i] = [""] * (len(seq2) + 1)
		# Initialize first row and column for indels
		for i in range(len(scoreMat)):
			scoreMat[i][0] = i * 5
		for i in range(len(scoreMat[0])):
			scoreMat[0][i] = i * 5
		# Build backtrack matrix
		for i in range(1, len(seq1) + 1):
			for j in range(1, len(seq2) + 1):
				down = scoreMat[i - 1][j] + self.NucScore(seq1[i - 1], "-")
				right = scoreMat[i][j - 1] + self.NucScore("-", seq2[j - 1])
				diagonal = scoreMat[i - 1][j - 1] + self.NucScore(seq1[i - 1], seq2[j - 1])
				scoreMat[i][j] = min([down, right, diagonal])
				if scoreMat[i][j] == down:
					backMat[i][j] = "down"
				elif scoreMat[i][j] == right:
					backMat[i][j] = "right"
				elif scoreMat[i][j] == diagonal:
					backMat[i][j] = "diagonal"
		return backMat, scoreMat[len(seq1)][len(seq2)]

	def out(self, backMat, v1, v2, i, j):
		if i == 0 or j == 0:
			if i > 0:
				o = self.out(backMat, v1, v2, i - 1, j)
				o[0] += v1[i - 1]
				o[1] += "-"
				return o
			elif j > 0:
				o = self.out(backMat, v1, v2, i, j - 1)
				o[0] += "-"
				o[1] += v2[j - 1]
				return o
			else:
				return ["", ""]
		elif backMat[i][j] == "right":
			o = self.out(backMat, v1, v2, i, j - 1)
			o[0] += "-"
			o[1] += v2[j - 1]
			return o
		elif backMat[i][j] == "down":
			o = self.out(backMat, v1, v2, i - 1, j)
			o[0] += v1[i - 1]
			o[1] += "-"
			return o
		else:
			o = self.out(backMat, v1, v2, i - 1, j - 1)
			o[0] += v1[i - 1]
			o[1] += v2[j - 1]
			return o

		# string1 = ''
		# string2 = ''
		# while i > 0 and j > 0:
		# 	if i == 0 or j == 0:
		# 		if i > 0:
		# 			string1 += v1[i - 1]
		# 			string2 += '-'
		# 			i -= 1
		# 		elif j > 0:
		# 			string1 += '-'
		# 			string2 += v2[j - 1]
		# 			j -= 1
		# 	elif backMat[i][j] == "right":
		# 		string1 = '-'
		# 		string2 = v2[j - 1]
		# 		j -= 1
		# 	elif backMat[i][j] == "down":
		# 		string1 += v1[i - 1]
		# 		string2 += '-'
		# 		i -= 1
		# 	else:
		# 		string1 += v1[i - 1]
		# 		string2 += v2[j - 1]
		# 		i -= 1
		# 		j -= 1
		# return [string1[::-1], string2[::-1]]

	def stringIndex(self, j, k):
		return j + (k - MAXINDELS)

	def bandedScore(self, seq1, seq2):
		# mat[seq1][seq2]
		# Initialize matrices
		seq1 = '-' + seq1
		seq2 = '-' + seq2
		scoreMat = [0] * (len(seq2))
		for i in range(len(seq2)):
			scoreMat[i] = [0] * ((MAXINDELS * 2) + 1)
		backMat = [""] * (len(seq2))
		for i in range((len(seq2))):
			backMat[i] = [""] * ((MAXINDELS * 2) + 1)
		# Initialize first row and column for indels
		for i in range(len(scoreMat[0])):
			if self.stringIndex(0, i) < 0:
				scoreMat[0][i] = math.inf
			else:
				scoreMat[0][i] = self.stringIndex(0, i) * 5
		# Build backtrack matrix
		for i in range(1, len(seq2)):
			for j in range(0, (MAXINDELS * 2) + 1):
				if not self.stringIndex(i, j) < 0 and not self.stringIndex(i, j) > len(seq1) - 1:
					if j + 1 == len(scoreMat[0]):
						down = math.inf
					else:
						down = scoreMat[i - 1][j + 1] + self.NucScore(seq1[self.stringIndex(i, j)], "-")
					if j == 0:
						right = math.inf
					else:
						right = scoreMat[i][j - 1] + self.NucScore("-", seq2[i])
					diagonal = scoreMat[i - 1][j] + self.NucScore(seq1[self.stringIndex(i, j)], seq2[i])
					scoreMat[i][j] = min([down, right, diagonal])
					if scoreMat[i][j] == down:
						backMat[i][j] = "down"
					elif scoreMat[i][j] == right:
						backMat[i][j] = "right"
					elif scoreMat[i][j] == diagonal:
						backMat[i][j] = "diagonal"
				else:
					scoreMat[i][j] = math.inf
		for i in reversed(range(len(scoreMat[len(seq2) - 1]))):
			if not scoreMat[len(seq2) - 1][i] == math.inf:
				break
		return backMat, scoreMat[len(seq2) - 1][i], i

	def bandedOut(self, backMat, v1, v2, i, j):
		if i == 0 or j == 0:
			if i > 0:
				o = self.bandedOut(backMat, v1, v2, i - 1, j)
				o[0] += v1[self.stringIndex(i, j) - 1]
				o[1] += "-"
				return o
			elif j > 0:
				o = self.bandedOut(backMat, v1, v2, i, j - 1)
				o[0] += "-"
				o[1] += v2[j - 1]
				return o
			else:
				return ["", ""]
		elif backMat[i][j] == "right":
			o = self.bandedOut(backMat, v1, v2, i, j - 1)
			o[0] += "-"
			o[1] += v2[j - 1]
			return o
		elif backMat[i][j] == "down":
			o = self.bandedOut(backMat, v1, v2, i - 1, j + 1)
			o[0] += v1[self.stringIndex(i, j) - 1]
			o[1] += "-"
			return o
		else:
			o = self.bandedOut(backMat, v1, v2, i - 1, j)
			o[0] += v1[self.stringIndex(i, j) - 1]
			o[1] += v2[j - 1]
			return o

# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you 
# how many base pairs to use in computing the alignment

	def align( self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length

		subseq1 = seq1[:align_length]
		subseq2 = seq2[:align_length]

		if banded:
			if len(subseq1) < len(subseq2):
				switch = subseq1
				subseq1 = subseq2
				subseq2 = switch
			# Score the sequences and return the back matrix
			backMat, score, start = self.bandedScore(subseq1, subseq2)
			# Go through the back matrix and construct the sequence
			seqs = self.bandedOut(backMat, subseq1, subseq2, len(subseq2), start)
		else:
			# Score the sequences and return the back matrix
			backMat, score = self.score(subseq1, subseq2)
			# Go through the back matrix and construct the sequence
			seqs = self.out(backMat, subseq1, subseq2, len(subseq1), len(subseq2))

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
# 		score = random.random()*100;
# 		alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
# 			len(seq1), align_length, ',BANDED' if banded else '')
# 		alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
# 			len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################					
		
		return {'align_cost':score, 'seqi_first100':seqs[0], 'seqj_first100':seqs[1]}


