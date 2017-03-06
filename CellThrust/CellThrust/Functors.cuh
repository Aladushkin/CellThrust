#include <thrust/host_vector.h>

// ������� ��� �������� ���� double3
struct SumDouble3
{
	__host__ __device__ double3 operator() (const double3 &a, const double3 &b) const {
		double3 c;
		c.x = a.x + b.x;
		c.y = a.y + b.y;
		c.z = a.z + b.z;
		return c;
	}
};


//������� ��� ���������� ����������
struct FindDistanse
{
	double3 a;

	FindDistanse(double3 _a) : a(_a) {}

	__host__ __device__ double operator() (const double3 &b) const {
		double3 c;
		c.x = a.x - b.x;
		c.y = a.y - b.y;
		c.z = a.z - b.z;

		double dist = sqrt(c.x * c.x + c.y * c.y + c.z * c.z);
		return dist;
	}
};

//������� �������� ���� ��������
struct PlusDouble3
{
	__host__ __device__ double3 operator() (const double3 &a, const double3 &b) const {
		double3 result;
		
		result.x = a.x + b.x;
		result.y = a.y + b.y;
		result.z = a.z + b.z;

		return result;
	}
};

//������� ���� �������
struct PlusG
{
	//����� ����� ���� �� �� ��������� ����� ����������, �� ��� �������� �����, ����� �� ������� �������� 
	//������ � ������ �������������, ������� ���� ���
	__host__ __device__ double3 operator() (const double3 &a) const {
		double3 result;
		result.x = a.x;
		result.y = a.y;
		result.z = a.z - 0.1;

		return result;
	}
};
