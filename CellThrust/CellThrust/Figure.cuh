#include <thrust/device_vector.h>
#include "Wall.cuh"


class Figure
{
private:
	static int r; // радиус цитоплазмы
	static int iterationCount; // количество итераций

	static double actinRadius; // толщина актиновой нити
	static double actinForce;  // сила актина

	static int receptorCount; // количество рецепторов
	static int ligandCount; // количество лигандов
	static std::vector<int> shellType; //вектор, из нулей единицы, нужен для определения шаров рецепторов и оболочки
	static std::vector<int> ballType; //вектор, из нулей единицы, нужен для определения активных и обычных шаров

	static std::vector<std::pair<int, int>> v_ActinSchedule; //вектор временных интервалов для включения актина
	static std::vector<std::pair<int, int>> v_TractionSchedule; //вектор временных интервалов для включения транспортной силы

	/*
	static double l0; // 
	static double sigma_b; // spring constant (0.1 – 2.5 pN/nm)
	static double sigma_st; // transition spring constant (0.1 – 2.5 pN/nm)
	static double k0_f; // unstressed forward rate (2 – 200s^-1)
	static double k0_r; //unstressed reverse rate (2 – 200s^-1)
	*/

	static double tractionForce;  // сила тяги

	static double3 kernel; //координаты ядра
	static double3 forceKernel;
	static double3 forceSumKernel;

	static thrust::host_vector<double3> h_balls; // вектор шаров в памяти CPU
	static thrust::host_vector<double3> h_forcesSumBall; // вектор сил всех шаров
	static thrust::host_vector<double3> h_forcesBall; // вектор сил для отдельного шара
	static thrust::host_vector<int> bondsCount; // массив для хранения толищны каждой адгезионной связи

	static thrust::host_vector<double3> h_ligands; 

	//static thrust::host_vector<double> h_forces_test; 
	//static thrust::host_vector<double> h_max_forces_test; 

	static thrust::host_vector<double3> h_shell;
	static thrust::host_vector<double3> h_forcesShell;
	static thrust::host_vector<double3> h_forcesSumShell;

	static double a; // расстояние между шарами оболочки
	static double3 nullVec; //нулевой вектор для суммирования сил

	//матрица с линками на соседей оболочки
	static thrust::host_vector<thrust::host_vector<double3*>> neighbors;
	//массив линков на шары рецепторы
	static thrust::host_vector<double3*> l_receptors;


	//вектор поверхностей, произвольный объект класса не влезет, поэтому double4
	//A - x, B - y, C - z, D - w
	static thrust::host_vector<double4> Figure::walls;

	//вектор ветвей актина
	static thrust::host_vector<double3> Figure::branches;

	//static thrust::host_vector<double3> h_balls; // вектор шаров в памяти CPU

	// конфигурация и создание цитоплазмы
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