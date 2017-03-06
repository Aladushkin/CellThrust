#include <thrust/host_vector.h>

//функтор силы ЛД для шаров цитоплазмы
struct FindForcesLJ
{
	double3 y;
	double maxForce;

	FindForcesLJ(double3 _y, double _maxForce) : y(_y), maxForce(_maxForce) {}

	__host__ __device__ double3 operator()(const double3& x) const {
		double3 normVec;
		normVec.x = x.x - y.x; normVec.y = x.y - y.y; normVec.z = x.z - y.z;

		double dist = sqrt(normVec.x * normVec.x + normVec.y * normVec.y + normVec.z * normVec.z);
		normVec.x = normVec.x / dist;
		normVec.y = normVec.y / dist;
		normVec.z = normVec.z / dist;

		double D = 0.16, a = 1.0;
		double temp = a / dist;

		temp = temp * temp * temp;
		temp = temp * temp;

		double force = 12 * D * (temp * temp - temp);

		//ограничиваем силу сверху
	
		if (abs(force) > maxForce)
		{
			if (force >= 0)
			{
				force = maxForce;
			}
			else
			{
				force = (-1)*maxForce;
			}
		}

		normVec.x = normVec.x * force;
		normVec.y = normVec.y * force;
		normVec.z = normVec.z * force;

		return normVec;
	}
};

struct FindForcesShellCytLJ
{
	double3 y;
	double maxForce;

	FindForcesShellCytLJ(double3 _y ,double _maxForce) : y(_y), maxForce(_maxForce) {}

	__host__ __device__ double3 operator()(const double3& x) const {
		double3 normVec;
		normVec.x = x.x - y.x; normVec.y = x.y - y.y; normVec.z = x.z - y.z;

		double dist = sqrt(normVec.x * normVec.x + normVec.y * normVec.y + normVec.z * normVec.z);
		normVec.x = normVec.x / dist;
		normVec.y = normVec.y / dist;
		normVec.z = normVec.z / dist;

		double D = 0.16, a = 1.0;
		double temp = a / dist;

		temp = temp * temp * temp;
		temp = temp * temp;

		double force = 12 * D * (temp * temp - temp);
		
		if (abs(force) > maxForce)
		{
			if (force >= 0)
			{
				force = maxForce;
			}
			else
			{
				force = (-1)*maxForce;
			}
		}
		

		normVec.x = normVec.x * force;
		normVec.y = normVec.y * force;
		normVec.z = normVec.z * force;

		return normVec;
	}
};

struct FindForcesShellLJ
{
	double a;
	double maxForce;

	FindForcesShellLJ(double _a, double _maxForce) : a(_a), maxForce(_maxForce) {}

	__host__ __device__ double3 operator() (const thrust::host_vector<double3*> &h_vector) const {
		double3 *root = h_vector.back();
		
		double3 result;
		result.x = 0; result.y = 0; result.z = 0;

		for (int i = 0; i < h_vector.size() - 1; i++) {
			double3 normVec, *node = h_vector[i];
			normVec.x = root->x - node->x; normVec.y = root->y - node->y; normVec.z = root->z - node->z;

			double dist = sqrt(normVec.x * normVec.x + normVec.y * normVec.y + normVec.z * normVec.z);
			normVec.x = normVec.x / dist;
			normVec.y = normVec.y / dist;
			normVec.z = normVec.z / dist;

			//пока а - статичен, нужно будет изменить
			/*double D = 0.42, a = 0.8;
			double temp = a / dist;

			temp = temp * temp * temp;
			temp = temp * temp;

			double force = 12 * D * (temp * temp - temp);
			*/

			double k = 1.5 ;
			double force = k * (a - dist);

			result.x += normVec.x * force;
			result.y += normVec.y * force;
			result.z += normVec.z * force;
		}

		return result;
	}
};

struct FindForceTest // вспомогательный функтор, для поиска максимальной силы взаимодействия меджу шарами 
{
	double3 y;

	FindForceTest(double3 _y) : y(_y) {}

	__host__ __device__ double operator()(const double3& x) const {
		double3 normVec;
		normVec.x = x.x - y.x; normVec.y = x.y - y.y; normVec.z = x.z - y.z;

		double dist = sqrt(normVec.x * normVec.x + normVec.y * normVec.y + normVec.z * normVec.z);
		normVec.x = normVec.x / dist;
		normVec.y = normVec.y / dist;
		normVec.z = normVec.z / dist;

		double D = 0.16, a = 1.0;
		double temp = a / dist;

		temp = temp * temp * temp;
		temp = temp * temp;

		double force = 12 * D * (temp * temp - temp);

		return force;
	}
};