import numpy as np
import particle
import cell

class VerletList:

	def __init__(self, cell_list, skin):
		self.cell_list = cell_list
		self.skin = skin
		self.verlet = dict()
		for p in self.cell_list.particles:
			self.verlet[p] = cell.Cell()
			cell_index = self.cell_list.getCellIndex(p)
			for q in self.cell_list.cells[cell_index].particles:
				#np.linalg.norm calculate distance
				if np.linalg.norm(p.position - q.position) <= (self.cell_list.rc + self.skin):
					self.verlet[p].add(q)
			for adj_cell_index in self.cell_list.getAdjacentCell(cell_index):
				for q in self.cell_list.cells[adj_cell_index].particles:
					if np.linalg.norm(p.position - q.position) <= (self.cell_list.rc + self.skin):
						self.verlet[p].add(q)

	# def getCellIndex(self, p):
		# return tuple(np.floor(p.position / (self.cell_list.rc + self.skin)))