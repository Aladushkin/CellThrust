#ifndef __kerA__Adhesion__
#define __kerA__Adhesion__

#include <thrust\host_vector.h>
//функтор для поиска количества адгезионных связей
struct FindAdhesionContacts
{
	double3 ligand;

	FindAdhesionContacts(double3 _ligand) : ligand(_ligand){}

	__host__ __device__ std::pair<int, double3> operator() (const double3* a, const std::pair<int,double3> b) const {
		double3 nullLigand;  nullLigand.x = 0;  nullLigand.y = 0;  nullLigand.z = 0;
		std::pair<int, double3> resultPair; resultPair.first = 0; resultPair.second= nullLigand;

		if (a != NULL)
		{
			int n_bonds;
			double distance = sqrt((ligand.x - a->x)*(ligand.x - a->x) + (ligand.y - a->y)*(ligand.y - a->y) + (ligand.z - a->z)*(ligand.z - a->z));

			double k_f = 100 * exp(-0.002*(distance - 0.2)*(distance - 0.2) / 2 / (1.38*pow(10, -8)) / 310); // k0_f*exp(-sigma_ts*(l - l0)*(l - l0) / 2 / K_B / T);
			double k_r = 5 * exp(-0.001*(distance - 0.2)*(distance - 0.2) / 2 / (1.38*pow(10, -8)) / 310); // k0_r*exp((sigma_b-sigma_ts)*(l - l0)*(l - l0) / 2 / K_B / T)
			double P_f = 1 - exp(-k_f); //-k_f
			double P_r = 1 - exp(-k_r);

			if (b.first == 0) // связи отсутствуют
			{
				if (distance <= 0.5) // если лиганд достаточно близко
				{
					// создание новой связи
					double N_ran1 = 0.01 * (rand() % 101);
					if (P_f > N_ran1)
					{
						resultPair.first = 1;
						resultPair.second.x = ligand.x; resultPair.second.y = ligand.y; resultPair.second.z = ligand.z;
						return resultPair;
					}
				}
				return resultPair; // нет лигандов поблизости
			}
			else // связи уже существуют
			{
				n_bonds = b.first;
				if (ligand.x == b.second.x && ligand.y == b.second.y && ligand.z == b.second.z) // среди всех лигандов находим тот самый единственный соединеннный с рецептором
				{
					if (distance <= 0.5) // если лиганд достаточно близко
					{
						double P_f = 1 - exp(-k_f);
						double P_r = 1 - exp(-k_r);

						// разрыв существующих связей
						double N_ran2 = 0.01 * (rand() % 101);
						if (P_r > N_ran2)
						{
							n_bonds = n_bonds - 1;
						}

						// создание новых связей
						double N_ran1 = 0.01 * (rand() % 101);
						if (P_f > N_ran1)
						{
							n_bonds = n_bonds + 1;
						}

						if (n_bonds <= 0) // число связей не может быть отрицательным
						{
							n_bonds = 0;
						}
						resultPair.first = n_bonds; resultPair.second.x = ligand.x; resultPair.second.y = ligand.y; resultPair.second.z = ligand.z;
						return resultPair;
					}
					return resultPair;
				}
				return b;
			}
		}
	}
};


struct FindAdhesionForce
{
	double3 ligand;

	FindAdhesionForce(double3 _ligand) : ligand(_ligand){}

	__host__ __device__ double3 operator() (const double3* a, const std::pair<int, double3> b) const {
		double3 result;
		result.x = 0; result.y = 0; result.z = 0;
		if (a != NULL)
		{
			if (b.first != 0)
			{
				if (ligand.x == b.second.x && ligand.y == b.second.y && ligand.z == b.second.z)
				{
					double3 rec_lig;
					double3 norm_rec_lig;
					double distance;

					rec_lig.x = ligand.x - a->x;
					rec_lig.y = ligand.y - a->y;
					rec_lig.z = ligand.z - a->z;

					distance = sqrt(rec_lig.x*rec_lig.x + rec_lig.y*rec_lig.y + rec_lig.z*rec_lig.z);
					norm_rec_lig.x = rec_lig.x / distance;
					norm_rec_lig.y = rec_lig.y / distance;
					norm_rec_lig.z = rec_lig.z / distance;

					double adhesionForce = 10*0.001*b.first*(distance - 0.2);

					std::cout << "Сила адгезии с лигандом " << ligand.x << " " << ligand.y << " " << ligand.z << " = " << adhesionForce << ", дистанция = " << distance << ", толщина связи = " << b.first << std::endl;

					result.x = norm_rec_lig.x*adhesionForce;
					result.y = norm_rec_lig.y*adhesionForce;
					result.z = norm_rec_lig.z*adhesionForce;
				}
			}
		}
		return result;
	}
};
#endif
