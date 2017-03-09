#include <thrust/device_vector.h>
#include "Wall.cuh"


class Figure
{
private:
	static int r; // ������ ����������
	static int iterationCount; // ���������� ��������

	static double actinRadius; // ������� ��������� ����
	static double actinForce;  // ���� ������

	static int receptorCount; // ���������� ����������
	static int ligandCount; // ���������� ��������
	static std::vector<int> shellType; //������, �� ����� �������, ����� ��� ����������� ����� ���������� � ��������
	static std::vector<int> ballType; //������, �� ����� �������, ����� ��� ����������� �������� � ������� �����

	static std::vector<std::pair<int, int>> v_ActinSchedule; //������ ��������� ���������� ��� ��������� ������
	static std::vector<std::pair<int, int>> v_TractionSchedule; //������ ��������� ���������� ��� ��������� ������������ ����

	/*
	static double l0; // 
	static double sigma_b; // spring constant (0.1 � 2.5 pN/nm)
	static double sigma_st; // transition spring constant (0.1 � 2.5 pN/nm)
	static double k0_f; // unstressed forward rate (2 � 200s^-1)
	static double k0_r; //unstressed reverse rate (2 � 200s^-1)
	*/

	static double tractionForce;  // ���� ����

	static double3 kernel; //���������� ����
	static double3 forceKernel;
	static double3 forceSumKernel;

	static thrust::host_vector<double3> h_balls; // ������ ����� � ������ CPU
	static thrust::host_vector<double3> h_forcesSumBall; // ������ ��� ���� �����
	static thrust::host_vector<double3> h_forcesBall; // ������ ��� ��� ���������� ����
	static thrust::host_vector<int> bondsCount; // ������ ��� �������� ������� ������ ����������� �����

	static thrust::host_vector<double3> h_ligands; 

	//static thrust::host_vector<double> h_forces_test; 
	//static thrust::host_vector<double> h_max_forces_test; 

	static thrust::host_vector<double3> h_shell;
	static thrust::host_vector<double3> h_forcesShell;
	static thrust::host_vector<double3> h_forcesSumShell;

	static double a; // ���������� ����� ������ ��������
	static double3 nullVec; //������� ������ ��� ������������ ���

	//������� � ������� �� ������� ��������
	static thrust::host_vector<thrust::host_vector<double3*>> neighbors;
	//������ ������ �� ���� ���������
	static thrust::host_vector<double3*> l_receptors;


	//������ ������������, ������������ ������ ������ �� ������, ������� double4
	//A - x, B - y, C - z, D - w
	static thrust::host_vector<double4> Figure::walls;

	//������ ������ ������
	static thrust::host_vector<double3> Figure::branches;

	//static thrust::host_vector<double3> h_balls; // ������ ����� � ������ CPU

	// ������������ � �������� ����������
	static void Config();
	
	static void CreateSphereCytoplasmAndKernel();
	static void CreateShell();

	static void FindActiveBalls();

	static void CreateReceptors();
	static void CreateLigands();

	static void CreateStepSchedule();
	
	static void InitializeVariable();
	static void FindDistanseShell();

	static void FindForcesBall(int iteration);
	static void FindForcesKernel(int iteration);
	static void FindForcesShell();

	static void FindRadiusVector();

	static void ToFile(int iteration);

public:
	static void Start();
	
};