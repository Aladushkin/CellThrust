#ifndef __MCA__Wall__
#define __MCA__Wall__

#include <thrust\host_vector.h>

//функтор для поиска силы коллизии со стенкой
struct FindCollision 
{
	double4 wall;

	FindCollision(double4 _wall) : wall(_wall) {}

	__host__ __device__ double3 operator() (const double3 &a) const {
		double3 result;
		result.x = 0; result.y = 0; result.z = 0;
		
		double norm = sqrt(wall.x * wall.x + wall.y * wall.y + wall.z * wall.z);

		//формула расстояния от точки до плоскости
		double dist = fabs(wall.x * a.x + wall.y * a.y + wall.z * a.z + wall.w) / norm;
		
		//рассояние, с которого начинает действовать сила отталкивания
		double interactionDist = 1.0;
		//множитель для силы отталкивания
		double eps = 0.1;

		if (dist <= interactionDist) {
			double force = eps / dist;
			result.x = (wall.x / norm) * force;
			result.y = (wall.y / norm) * force;
			result.z = (wall.z / norm) * force;
		}
		return result;
	}
};

#endif