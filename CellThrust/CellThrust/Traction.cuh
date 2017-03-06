#ifndef __kerA__Traction__
#define __kerA__Traction__

#include <thrust\host_vector.h>

//функтор дл€ поиска силы актина
struct FindTractionForce
{
	double3 branch;
	double3 kernel;
	double tractionForce;

	FindTractionForce(double3 _branch, double3 _kernel, double _tractionForce): branch(_branch), kernel(_kernel), tractionForce(_tractionForce) {}

	__host__ __device__ double3 operator() (const double3 &a, const int &b) const {
		double3 result;
		result.x = 0; result.y = 0; result.z = 0;
		if (b == 0)
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
			//вектор соедин€ющий €дро и мишень

			result.x = norm_ker_branch.x*tractionForce;
			result.y = norm_ker_branch.y*tractionForce;
			result.z = norm_ker_branch.z*tractionForce;
		}
		return result;
	}
};
#endif