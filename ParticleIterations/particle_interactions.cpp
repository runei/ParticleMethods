#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <random>
#include <tuple>
#include <algorithm>
#include <functional>

//====================================================================================

template<class T>
struct InteractResult
{
	std::vector<double> position_change;
	T properties_change;
};

//====================================================================================

template<class T, class U>
class Particle
{
protected:
	std::vector<double> m_position; // vector size is the # of dimensions os the particle
	std::vector<double> m_position_change;

	U m_properties;
	U m_properties_change;

public:
	std::vector<double> getPosition() const { return m_position; };
	U getProperties() const { return m_properties; };
	void setPositionChange(const double);

	Particle(const std::vector<double> & position, const U & properties) : m_position(position), m_position_change(position.size()), m_properties(properties) {};

	virtual std::tuple<std::vector<double>, U> interact(const T &) = 0;
	virtual void evolve() = 0;
};

template<class T, class U>
void Particle<T, U>::setPositionChange(const double val) {
	std::fill(m_position_change.begin(), m_position_change.end(), val);
}

//====================================================================================

struct StarProperties
{
	double mass;

	StarProperties(double m) : mass(m) {}
	StarProperties() {}
};

//====================================================================================

class Star: public Particle<Star, StarProperties>
{
private:
	double calculateForce(const double, const double, const double);
	double calculateDistance(const std::vector<double> &, const std::vector<double> &);
	double calculateAcceleration(const double, const double);
	std::vector<double> calculateDimensionsDistance(const std::vector<double> &, const std::vector<double> &);
public:
	Star(const std::vector<double> & _position, const StarProperties & _properties) : Particle(_position, _properties) {};

	std::tuple<std::vector<double>, StarProperties> interact(const Star &) override;
	void evolve() override;
};

std::vector<double> Star::calculateDimensionsDistance(const std::vector<double> &v1, const std::vector<double> &v2)
{
	//TO-DO this looks wrong, look for stl solution
	std::vector<double> result(v1.size());

	for (size_t i = 0; i < result.size(); ++i)
	{
		result[i] = v1[i] - v2[i];
	}

	return result;
}

double Star::calculateAcceleration(const double m, const double f)
{
	return f / m;
}

double Star::calculateForce(const double m1, const double m2, const double r)
{
	static const double gravit_const = 0.005;//6.674 * std::pow(10, -11);

	return gravit_const * ((m1 * m2) / std::pow(r, 2));
}

double Star::calculateDistance(const std::vector<double> & pos1, const std::vector<double> & pos2)
{
	double result = 0.0;

	for (size_t i = 0; i < pos1.size(); ++i)
	{
		result += std::pow(pos1[i] - pos2[i], 2);
	}

	return std::sqrt(result);
}

std::tuple<std::vector<double>, StarProperties> Star::interact(const Star & s)
{
	std::vector<double> position_change(getPosition().size(), 0);
	StarProperties properties_change(0.0);

	if (s.getProperties().mass > getProperties().mass)
	{
		const double dist = calculateDistance(getPosition(), s.getPosition());
		const double force = calculateForce(s.getProperties().mass, getProperties().mass, dist);
		const double acc = calculateAcceleration(getProperties().mass, force);

		std::vector<double> temp(calculateDimensionsDistance(s.getPosition(), getPosition()));

		if (acc > dist)
		{
			position_change = temp;
			properties_change.mass = s.getProperties().mass;
		}
		else
		{
			const double sum = std::accumulate(temp.begin(), temp.end(), 0.0);

			for (size_t i = 0; i < position_change.size(); ++i)
			{
				position_change[i] = (temp[i] / sum) * acc;
				std::cout << position_change[i] << "\n";
			}
		}

		// std::cout << force / getProperties().mass << "\n";
	}

	return std::make_tuple(position_change, properties_change);
}

void Star::evolve()
{

	m_properties.mass += m_properties_change.mass;
}

//====================================================================================

template<class T, class U>
class GenericList
{
private:
	std::vector<double> m_dimensions;
	std::vector<Particle<T, U> &> particles;
public:
	void nextTimeStep();
};

template<class T, class U>
void GenericList<T, U>::nextTimeStep()
{
	for (auto p : particles)
	{
		
	}
}

//====================================================================================

//====================================================================================
/*
class Cell
{
	std::vector<Particle> particles;
public:
	Cell();
	void add(const Particle &);
};

void Cell::add(const Particle & p)
{
	this->particles.push_back(p);
}

//====================================================================================

class CellList
{
	std::vector<int> dimensions; // sizes of each dimension, Lx, Ly...
	std::vector<Cell> cells;
	const int rc;

	int getCellIndex(Particle &);
public:
	CellList(const std::vector<int> &, const std::vector<Particle> &, const int);
};

CellList::CellList(const std::vector<int> & _dimension_sizes, const std::vector<Particle> & particles, const int _rc) :
		dimensions(_dimension_sizes),
		cells(std::accumulate(dimensions.begin(), dimensions.end(), 0, std::multiplies<int>())),
		rc(_rc)
{
	for (Particle p : particles)
	{
		const int cell_index = getCellIndex(p);
		cells[cell_index].add(p);
	}
}

int CellList::getCellIndex(Particle & p)
{
	int dim_mult = cells.size();
	int cell_index = 0;
	for (size_t i = 0; i < dimensions.size(); ++i)
	{
		dim_mult /= dimensions[i];
		cell_index += dim_mult * std::floor(p.getPosition()[i] / rc);
	}
	return cell_index;
}

//====================================================================================


*/
int main()
{

	Star star1({0.1}, StarProperties(1));
	Star star2({0.7}, StarProperties(2));

	star1.interact(star2);

	// Universe env(5, 2, 10);
	// std::cout << "asuidha\n";
	// env.print2D();

	return 0;
}