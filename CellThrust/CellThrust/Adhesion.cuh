#ifndef __kerA__Adhesion__
#define __kerA__Adhesion__

#include <thrust\host_vector.h>
//функтор для поиска количества адгезионных связей
struct FindAdhesionContacts
{
	double3 ligand;

	FindAdhesionContacts(double3 _ligand) : ligand(_ligand){}

	__host__ __device__ int operator() (const double3* a, const int &b) const {

		int result = 0;
		if (a != NULL)
		{
			int n_bonds;
			double distance = sqrt((ligand.x - a->x)*(ligand.x - a->x) + (ligand.y - a->y)*(ligand.y - a->y) + (ligand.z - a->z)*(ligand.z - a->z));

			double k_f = 6 * exp(-0.0001*(distance - 0.1)*(distance - 0.1) / 2 / (1.38*pow(10, -8)) / 310); // k0_f*exp(-sigma_ts*(l - l0)*(l - l0) / 2 / K_B / T);
			double k_r = 1.6 * exp(-0.00005*(distance - 0.1)*(distance - 0.1) / 2 / (1.38*pow(10, -8)) / 310); // k0_r*exp((sigma_b-sigma_ts)*(l - l0)*(l - l0) / 2 / K_B / T)
			double P_f = 1 - exp(-k_f); //-k_f
			double P_r = 1 - exp(-k_r);

			if (b == 0) // связи отсутствуют
			{
				if (distance <= 1) // если лиганд достаточно близко
				{
					// создание новой связи
					double N_ran1 = 0.01 * (rand() % 101);
					if (P_f > N_ran1)
					{
						return 1;
					}
				}
				return 0; // нет лигандов поблизости
			}
			else // связи уже существуют
			{
				if (distance <= 1) // если лиганд достаточно близко
				{
					n_bonds = b;
					double P_f = 1 - exp(-k_f);
					double P_r = 1 - exp(-k_r);

					// разрыв существующих связей
					for (int i = 0; i < n_bonds; i++)
					{
						double N_ran2 = 0.01 * (rand() % 101);
						if (P_r > N_ran2)
						{
							n_bonds = b - 1;
						}
					}

					// создание новых связей
					double N_ran1 = 0.01 * (rand() % 101);
					if (P_f > N_ran1)
					{
						n_bonds = b + 1; 
					}

					if (n_bonds <= 0) // число связей не может быть отрицательным
					{
						n_bonds = 0;
					}
					return n_bonds;
				}
				return 0; // дистанция слишком велика
			}
		}
		return result;
	}
};


struct FindAdhesionForce
{
	double3 ligand;

	FindAdhesionForce(double3 _ligand) : ligand(_ligand){}

	__host__ __device__ double3 operator() (const double3* a, const int &b) const {
		double3 result;
		result.x = 0; result.y = 0; result.z = 0;
		if (a != NULL)
		{
			if (b != 0)
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

				double adhesionForce = 100*0.00095*b*(distance - 0.1);

				std::cout << "Сила адгезии с лигандом " << ligand.x << " " << ligand.y << " " << ligand.z << " = " << adhesionForce << ", толщина связи = " << b << std::endl;

				result.x = norm_rec_lig.x*adhesionForce;
				result.y = norm_rec_lig.y*adhesionForce;
				result.z = norm_rec_lig.z*adhesionForce;
			}
		}
		return result;
	}
};
#endif
