import numpy as np

class Cell:

	def __init__(self):
		self.particles = np.array([])

	def add(self, p):
		self.particles = np.append(self.particles, p)