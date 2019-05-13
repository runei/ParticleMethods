import numpy as np
import cell_list
import particle
import random
import matplotlib.pyplot as plt
import verlet_list

def printCellList(cell_list):
	random_cell = random.sample(cell_list.cells, 1)[0]
	adj_cells = cell_list.getAdjacentCell(random_cell)
	x = []
	y = []
	colors = []
	# print(cell_list.cells.items())
	for key, cell in cell_list.cells.items():
		if key == random_cell:
			color = 'green'
		elif key in adj_cells:
			color = 'yellow'
		else:
			color = 'blue'
		for p in cell.particles:
			x.append(p.position[0])
			y.append(p.position[1])
			colors.append(color)

	plt.scatter(x, y, c=colors)
	# plt.show()

def printVerletList(verlet_list):
	random_cell = random.sample(verlet_list.cell_list.cells, 1)[0]
	random_particle = np.random.choice(verlet_list.cell_list.cells[random_cell].particles)
	# print(random_particle)
	adj_cells = cell_list.getAdjacentCell(random_cell)
	x = []
	y = []
	colors = []
	for p in verlet_list.cell_list.particles:
		if p == random_particle:
			color = 'green'
		elif p in verlet_list.verlet[random_particle].particles:
			color = 'red'
		else:
			color = 'blue'
		x.append(p.position[0])
		y.append(p.position[1])
		colors.append(color)

	plt.scatter(x, y, c=colors)
	plt.show()

if __name__ == '__main__':
	n = 10000
	random_pos = np.random.rand(n, 2)
	particles = [particle.Particle(pos, particle.ParticleProperties()) for pos in random_pos]
	rc = 0.05
	cell_list = cell_list.CellList(particles, rc)
	printCellList(cell_list)


	verlet_list = verlet_list.VerletList(cell_list, 0.05)
	printVerletList(verlet_list)