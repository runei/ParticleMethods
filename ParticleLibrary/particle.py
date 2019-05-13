import numpy as np

class ParticleProperties:

	def __init__(self):
		pass

	def setZeros(self):
		pass

	def __add__(self, properties):
		return 0

class Particle:

	def __init__(self, position, properties):
		self.position = np.array(position)
		self.position_change = np.array([])
		self.properties = properties
		self.properties_change = ParticleProperties()

	def interact(self, q):
		#return [position_change, properties_change]
		return [0, ParticleProperties()]

	def evolve(self):
		pass