#ifndef __kerA__Actin__
#define __kerA__Actin__

#include <thrust\host_vector.h>

//функтор для поиска силы актина
struct FindActinForce
{
	double3 branch;
	double3 kernel;
	double actinRadius;
	double actinForce;
	double diameter;

	FindActinForce(double3 _branch, double3 _kernel, double _actinRadius, double _actinForce, double _diameter) : branch(_branch), kernel(_kernel), actinRadius(_actinRadius), actinForce(_actinForce), diameter(_diameter) {}

	__host__ __device__ double3 operator() (const double3 &a, const int &b) const {
		double3 result;
		result.x = 0; result.y = 0; result.z = 0;
		if (b == 1)
		{
			double3 ker_branch;
			double3 norm_ker_branch;
			double mod_ker_branch;

			ker_branch.x = branch.x - kernel.x;
			ker_branch.y = branch.y - kernel.y;
			ker_branch.z = branch.z - kernel.z;
			mod_ker_branch = sqrt(ker_branch.x*ker_branch.x + ker_branch.y*ker_branch.y + ker_branch.z*ker_branch.z);
			norm_ker_branch.x = ker_branch.x / mod_ker_branch;
			norm_ker_branch.y = ker_branch.y / mod_ker_branch;
			norm_ker_branch.z = ker_branch.z / mod_ker_branch;
			//вектор соединяющий ядро и мишень

			double3 ker_a;

			ker_a.x = a.x - kernel.x;
			ker_a.y = a.y - kernel.y;
			ker_a.z = a.z - kernel.z;
			//вектор соединяющий ядро и шар

			double pr_a_tau;
			double3 a_tau;

			pr_a_tau = (ker_a.x*ker_branch.x + ker_a.y*ker_branch.y + ker_a.z*ker_branch.z) / mod_ker_branch;
			a_tau.x = pr_a_tau*norm_ker_branch.x;
			a_tau.y = pr_a_tau*norm_ker_branch.y;
			a_tau.z = pr_a_tau*norm_ker_branch.z;

			double3 a_rad;
			double mod_a_rad;

			a_rad.x = ker_a.x - a_tau.x;
			a_rad.y = ker_a.y - a_tau.y;
			a_rad.z = ker_a.z - a_tau.z;
			mod_a_rad = sqrt(a_rad.x*a_rad.x + a_rad.y*a_rad.y + a_rad.z*a_rad.z);

			if ((pr_a_tau > diameter / 4)&(mod_a_rad < actinRadius)) //проекция на вектор направления должна быть положительна и больше d/4, радиальная проекция меньше радиуса вырезаемого цилиндра
			{
				result.x = norm_ker_branch.x*actinForce;
				result.y = norm_ker_branch.y*actinForce;
				result.z = norm_ker_branch.z*actinForce;
			}
		}
		return result;
	}
};

#endif

