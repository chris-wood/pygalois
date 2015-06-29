import pickle
import sys
from FFArithmetic import *

class FFNormalBasisArithmetic(FFArithmetic):
	''' Class that implements normal basis arithmetic for finite fields.
	'''
	def __init__(self, field, baseField):
		self.field = field
		self.base = baseField
