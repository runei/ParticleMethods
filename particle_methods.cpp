#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include <numeric>
#include <algorithm>
#include <random>
#include <fstream>
#include <string>

//====================================================================================

namespace GeneralFunctions
{
	double calculateForce(const double, const double, const double);
	double calculateDistance(const std::vector<double> &, const std::vector<double> &);
	double calculateAcceleration(const double, const double);
	std::vector<double> calculateDimensionsDistance(const std::vector<double> &, const std::vector<double> &);
}

double GeneralFunctions::calculateAcceleration(const double m, const double f)
{
	return f / m;
}

double GeneralFunctions::calculateForce(const double m1, const double m2, const double r)
{
	static const double gravit_const = 0.002;//6.674 * std::pow(10, -11);

	return gravit_const * ((m1 * m2) / std::pow(r, 2));
}

double GeneralFunctions::calculateDistance(const std::vector<double> & pos1, const std::vector<double> & pos2)
{
	double result = 0.0;

	for (size_t i = 0; i < pos1.size(); ++i)
	{
		result += std::pow(pos1[i] - pos2[i], 2);
	}

	return std::sqrt(result);
}

std::vector<double> GeneralFunctions::calculateDimensionsDistance(const std::vector<double> &v1, const std::vector<double> &v2)
{
	//TO-DO this looks wrong, look for stl solution
	std::vector<double> result(v1.size());

	for (size_t i = 0; i < result.size(); ++i)
	{
		result[i] = v1[i] - v2[i];
	}

	return result;
}

//====================================================================================

class ParticleProperties
{
private:
	double m_mass;
public:
	//required
	ParticleProperties(double m) : m_mass(m) {}
	ParticleProperties() {}
	double mass() const { return m_mass; };
	void mass(const double m) { m_mass = m; };
	void setAllValues(const double val) { m_mass = val; };
	void add(const ParticleProperties & prop) { m_mass += prop.mass(); };
};

//====================================================================================

class Particle
{
protected:
	std::vector<double> m_position; // vector size is the # of dimensions os the particle
	std::vector<double> m_position_change;

	ParticleProperties m_properties;
	ParticleProperties m_properties_change;

public:
	std::vector<double> getPosition() const { return m_position; };
	ParticleProperties getProperties() const { return m_properties; };
	void setPositionChange(const double);
	void addPositionChange(const std::vector<double> &);
	ParticleProperties getPropertiesChange() const { return m_properties_change; };
	void setPropertiesChange(const double val) { m_properties_change.setAllValues(val); };
	void addPropertiesChange(const ParticleProperties & prop) { m_properties_change.add(prop); };

	Particle(const std::vector<double> & position, const ParticleProperties & properties) : m_position(position), m_position_change(position.size()), m_properties(properties), m_properties_change(0.0) {};

	std::tuple<std::vector<double>, ParticleProperties> interact(const Particle &);
	void evolve();
};

void Particle::setPositionChange(const double val) {
	std::fill(m_position_change.begin(), m_position_change.end(), val);
}

void Particle::addPositionChange(const std::vector<double> & pos)
{
	// assert(m_position_change.size() == pos.size());

	std::vector<double> result;
	result.reserve(m_position_change.size());

	std::transform(m_position_change.begin(), m_position_change.end(), pos.begin(), std::back_inserter(result), std::plus<double>());

	m_position_change = result;
}

std::tuple<std::vector<double>, ParticleProperties> Particle::interact(const Particle & p)
{
	std::vector<double> position_change(getPosition().size(), 0);
	ParticleProperties properties_change(0.0);

	if (getProperties().mass() > 0 && p.getProperties().mass() && p.getProperties().mass() > getProperties().mass())
	{
		const double dist = GeneralFunctions::calculateDistance(getPosition(), p.getPosition());
		double force = GeneralFunctions::calculateForce(p.getProperties().mass(), getProperties().mass(), dist);

		if (force > 1.0)
		{
			force = 1.0;
			properties_change.mass(- getProperties().mass());
		}

		for (size_t i = 0; i < position_change.size(); ++i)
		{
			// std::cout << getPosition()[i] << " " << force << " ";// exit(0);
			const double diff = p.getPosition()[i] - getPosition()[i];
			position_change[i] += force * diff;
			/*if (getPosition()[i] > p.getPosition()[i])
			{
				position_change[i] += -force * getPosition()[i];
			}
			else
			{
				position_change[i] += force * getPosition()[i];
			}
			std::cout << " " << position_change[i] << "\n";// exit(0);
			// position_change[i] = (temp[i] / sum) * acc;
			// std::cout << position_change[i] << " ";*/
		}

		/*const double acc = GeneralFunctions::calculateAcceleration(getProperties().mass(), force);


		std::vector<double> temp(GeneralFunctions::calculateDimensionsDistance(p.getPosition(), getPosition()));

		if (dist > 0.0)
		{
			// std::cout << dist << "\n";
			if (acc > dist)
			{
				// position_change = temp;
				properties_change.mass(- getProperties().mass());
				// properties_change.mass(p.getProperties().mass());
			}
			else
			{
				const double sum = std::accumulate(temp.begin(), temp.end(), 0.0);

				for (size_t i = 0; i < position_change.size(); ++i)
				{
					position_change[i] = (temp[i] / sum) * acc;
					std::cout << position_change[i] << " ";
				}
				std::cout << "\n";
			}
		}*/

		// std::cout << force / getProperties().mass << "\n";
	}

	return std::make_tuple(position_change, properties_change);
}

void Particle::evolve()
{
	for (size_t i = 0; i < m_position.size(); ++i)
	{
		// std::cout << m_position[i] << " " ;
		m_position[i] += m_position_change[i];
		// std::cout << m_position[i] << " " << m_position_change[i] << "\n";
	}
	// exit(0);
	m_properties.add(m_properties_change);
}

//====================================================================================

namespace ListUtilities
{
	std::vector<Particle> createRandomParticles(const double, const int, const int);
	void print(const std::vector<Particle> &);
}

std::vector<Particle> ListUtilities::createRandomParticles(const double max_dim_size, const int n_particles, const int n_dim)
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);

// std::cout << max_dim_size << "\n";
	std::vector<Particle> particles;

    for (int i = 0; i < n_particles; ++i)
    {
    	std::vector<double> position(n_dim);
    	for (size_t i = 0; i < position.size(); ++i)
    	{
    		position[i] = dis(gen);
			// std::cout << position[i] << "\n";
    	}
    	particles.push_back(Particle(position, ParticleProperties(dis(gen)*5)));
    }

	return particles;
}

void ListUtilities::print(const std::vector<Particle> & particles)
{
	static int file_number = 0;
	std::string name = "output//res" + std::to_string(file_number) + ".txt";
	std::ofstream res(name, std::ofstream::out);
	for (auto & p : particles)
	{
		for (auto & pos : p.getPosition())
		{
			res << pos << " ";
		}
		res << p.getProperties().mass() << "\n";
	}
	res.close();
	++file_number;
}

//====================================================================================

class GenericList
{
private:
	std::vector<Particle> m_particles;
	std::vector<double> m_dimensions;
public:
	void nextTimeInterval();
	void run();
	GenericList(const std::vector<Particle> & particles, const std::vector<double> & dimensions) : m_particles(particles), m_dimensions(dimensions) {};
	GenericList(const std::vector<Particle> & particles) : m_particles(particles) {};
	void destroyParticles();
};

void GenericList::nextTimeInterval()
{
	for (Particle & p : m_particles)
	{
		p.setPositionChange(0.0);
		p.setPropertiesChange(0.0);
		for (Particle & q : m_particles)
		{
			auto [pos_change_temp, prop_change_temp] = p.interact(q);
			p.addPropertiesChange(prop_change_temp);
			p.addPositionChange(pos_change_temp);
		}
	}
	for (Particle & p : m_particles)
	{
		p.evolve();
		// std::cout << p.getProperties().mass() << "\n";
	}
}

void GenericList::destroyParticles()
{
	size_t i = 0;
	while (i < m_particles.size())
	{
		if (m_particles[i].getProperties().mass() == 0.0)
		{
			m_particles.erase(m_particles.begin() + i);
		}
		else
		{
			++i;
		}
	}
	// m_particles.erase(
		// std::remove_if(m_particles.begin(), m_particles.end(), [](Particle p) {
			// return p.getProperties().mass() == 0.0;
		// })
	// );
}

void GenericList::run()
{
	ListUtilities::print(m_particles);
	for (int i = 0; i < 20; ++i)
	{
		nextTimeInterval();
		destroyParticles();
		ListUtilities::print(m_particles);
	}
}

//====================================================================================

class Cell
{
private:
	std::vector<Particle> m_particles;
public:
	Cell();
	void add(const Particle &);
	std::vector<Particle> getParticles() const { return m_particles; };
};

void Cell::add(const Particle & p)
{
	m_particles.push_back(p);
}


//====================================================================================

class CellList
{
private:
	std::vector<double> m_dimensions; // sizes of each dimension, Lx, Ly...
	std::vector<Cell> m_cells;
	std::vector<Particle> m_particles;
	const double rc;

	int getCellIndex(const Particle &);
	std::vector<int> adjacentCellIndex(const int);
public:
	CellList(const std::vector<double> &, const std::vector<Particle> &, const double);
	void nextTimeInterval();
};

CellList::CellList(const std::vector<double> & dimensions, const std::vector<Particle> & particles, const double _rc) :
		m_dimensions(dimensions),
		m_particles(particles),
		m_cells(std::accumulate(m_dimensions.begin(), m_dimensions.end(), 1, std::multiplies<int>())),
		rc(_rc)
{
	for (Particle & p : m_particles)
	{
		const int cell_index = getCellIndex(p);
		m_cells[cell_index].add(p);
	}
}

int CellList::getCellIndex(const Particle & p)
{
	int dim_mult = m_cells.size();
	int cell_index = 0;
	for (size_t i = 0; i < m_dimensions.size(); ++i)
	{
		dim_mult /= m_dimensions[i];
		cell_index += dim_mult * std::floor(p.getPosition()[i] / rc);
	}
	return cell_index;
};

std::vector<int> CellList::adjacentCellIndex(const int cell_index)
{
	const int n_adjacent_cells = std::pow(3, m_dimensions.size()) - 1;
	std::vector<int> result(n_adjacent_cells);

	

	return result;
}

void CellList::nextTimeInterval()
{
	for (Particle & p : m_particles)
	{
		p.setPositionChange(0.0);
		p.setPropertiesChange(0.0);
		const int cell_index = getCellIndex(p);
		for (Particle & q : m_cells[cell_index].getParticles())
		{
			auto [pos_change_temp, prop_change_temp] = p.interact(q);
			p.addPropertiesChange(prop_change_temp);
			p.addPositionChange(pos_change_temp);
		}
	}
}

//====================================================================================

int main()
{
	auto particles = ListUtilities::createRandomParticles(1.0, 5, 2);
	GenericList list(particles);
	list.run();
	// ListUtilities::print(particles);
	// GenericList g_list(particles);

	return 0;
}