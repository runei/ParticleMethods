import numpy as np
import particle
import cell

class CellList:

	def __init__(self, particles, rc):
		# self.dimensions = np.array(dimensions) # sizes of each dimension, Lx, Ly...
		self.particles = particles
		self.cells = dict()
		self.rc = rc
		for p in self.particles:
			cell_index = self.getCellIndex(p)
			if not cell_index in self.cells:
				self.cells[cell_index] = cell.Cell()
			self.cells[cell_index].add(p)

	def getCellIndex(self, p):
		return tuple(np.floor(p.position / self.rc))

	def getAdjacentCell(self, cell_index):
		if not isinstance(cell_index, tuple):
			print("Cell index should be tuple\n")
			return None
		array_cell_index = np.asanyarray(cell_index)
		result = np.array([])
		if len(array_cell_index) == 1:
			result = np.array([array_cell_index[0] - 1, array_cell_index[0] + 1])
			result = result[result >= 0]
		else:
			dimensions = tuple([3] * len(array_cell_index))
			result = np.indices(dimensions).T.reshape(-1, len(array_cell_index)) - 1
			result = result[~np.all(result == 0, axis = 1)]
			result = result + array_cell_index
			result = result[~np.any(result < 0, axis=1)]
		result = map(tuple, result)
		result = [i for i in result if i in self.cells]
		return result

	def shortRangeInteractions(self):
		for p in self.particles:
			p.position_change[p.position_change != 0] = 0
			p.properties_change.setZeros()
			cell_index = self.getCellIndex(p)
			for q in self.cells[cell_index].particles:
				[pos_change, prop_change] = p.interact(q)
				p.position_change += pos_change
				p.properties_change += prop_change
			for adj_cell_index in self.getAdjacentCell(cell_index):
				for q in self.cells[adj_cell_index].particles:
					[pos_change, prop_change] = p.interact(q)
					p.position_change += pos_change
					p.properties_change += prop_change
					# [p.position_change, p.properties_change] += p.interact(q)
		for p in self.particles:
			p.evolve()
