#include "Figure.cuh"
#include <iostream>
#include <fstream>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>

#include <locale>
#include <vector>

#include "InterfaceForShell.cuh"
#include "LennardJones.cuh"
#include "Functors.cuh"
#include "Wall.cuh"
#include "Actin.cuh"
#include "Adhesion.cuh"
#include "Traction.cuh"

double maxForceLJ = 0.5; // глобальная переменная для ограничения давления
int actinTime = 1000;
int tractionTime = 2500;
int stepTime = 1500;
int twoSteps = 0;

int Figure::r;
int Figure::iterationCount;

double Figure::actinRadius=0;
double Figure::actinForce=0;
std::vector<int> Figure::ballType;

int Figure::receptorCount=0;
int Figure::ligandCount=0;
std::vector<int> Figure::shellType;
thrust::host_vector <std::pair<int, double3>> Figure::bondsCount;

std::vector<std::pair<int, int>> Figure::v_ActinSchedule; 
std::vector<std::pair<int, int>> Figure::v_TractionSchedule;

double Figure::tractionForce = 0;

double3 Figure::kernel;
double3 Figure::forceKernel;
double3 Figure::forceSumKernel;

thrust::host_vector<double3> Figure::h_balls;
thrust::host_vector<double3> Figure::h_forcesBall;
thrust::host_vector<double3> Figure::h_forcesSumBall;

thrust::host_vector<double3> Figure::h_ligands;

//thrust::host_vector<double> Figure::h_forces_test;
//thrust::host_vector<double> Figure::h_max_forces_test;

thrust::host_vector<double3> Figure::h_shell;
thrust::host_vector<double3> Figure::h_forcesShell;
thrust::host_vector<double3> Figure::h_forcesSumShell;

double Figure::a;
double3 Figure::nullVec;

thrust::host_vector<thrust::host_vector<double3*>> Figure::neighbors;
thrust::host_vector<double3*> Figure::l_receptors;

thrust::host_vector<double4> Figure::walls;
thrust::host_vector<double3> Figure::branches;

void Figure::Start()
{
	Config();
	
	CreateSphereCytoplasmAndKernel();
	CreateShell();
	CreateLigands();
	
	InitializeVariable();
	FindDistanseShell();

	CreateReceptors();
	CreateStepSchedule();

	std::cout << "Количество шаров: " << h_balls.size() << std::endl;
	
	//сохраняем начальное положение
	ToFile(0);
	
	clock_t time;
	for (int i = 1; i < iterationCount; i++) {
		time = clock();
		//------------------------------------------------------------------------------------------
		if (i==actinTime)
		{
			FindActiveBalls();
		}

		FindForcesBall(i);
		FindForcesShell();
		FindForcesKernel(i);

		FindRadiusVector();
		ToFile(i);

		time = clock() - time;
		std::cout << "iteration = " << i << std::endl;
		std::cout << "runtime = " << time << std::endl;
	}
}

void Figure::Config()
{
	std::cout << "Введите диаметр клетки: "; // в микрометрах 
	std::cin >> r;

	int wallCount;
	std::cout << "Введите количество плоскостей\n";
	std::cin >> wallCount;

	walls = thrust::host_vector<double4>(wallCount);

	for (int i = 0; i < walls.size(); i++)
	{
		std::cout << "Введите коэффициенты уравнения плоскости " << i + 1 << "\n";
		
		std::cout << "A = ";
		std::cin >> walls[i].x;
		
		std::cout << "B = ";
		std::cin >> walls[i].y;
		
		std::cout << "C = ";
		std::cin >> walls[i].z;
		
		std::cout << "D = ";
		std::cin >> walls[i].w;
	}

	int actinCount;
	std::cout << "Введите количество актиновых нитей\n";
	std::cin >> actinCount;

	branches = thrust::host_vector<double3>(actinCount);

	for (int i = 0; i < branches.size(); i++)
	{
		std::cout << "Введите вектор роста актиновой нити " << i + 1 << "\n";
		
		std::cout << "x = ";
		std::cin >> branches[i].x;
		
		std::cout << "y = ";
		std::cin >> branches[i].y;
		
		std::cout << "z = ";
		std::cin >> branches[i].z;		
	}
	
	std::cout << "Введите число итераций\n";
	std::cin >> iterationCount;

	if (branches.size()!=0)
	{
	std::cout << "Введите радиус актиновой нити\n";
	std::cin >> actinRadius; // в микрометрах

	std::cout << "Введите силу актина \n";
	std::cin >> actinForce; // ???

	std::cout << "Введите силу тяги \n";
	std::cin >> tractionForce;
	}

	std::cout << "Введите количество лигандов \n";
	std::cin >> ligandCount;

	std::cout << "Введите количество рецепторов \n";
	std::cin >> receptorCount;

}

void Figure::CreateSphereCytoplasmAndKernel()
{
	//для работы с CUDA создаётся отдельный вектор на хосте, сейчас неактуален, но впоследствии пригодится
	thrust::host_vector<double3> tempBalls;
	double3 temp;

	double radius = r / 2;
	for (double i = -radius; i < radius; i += 1.8)
	{
		for (double j = -radius; j < radius; j += 1.8) 
		{
			for (double k = -radius; k < radius; k += 1.8) {
				if (i * i + j * j + k * k <= radius * radius &&
					i * i + j * j + k * k >= 0.81) {
					temp.x = i; temp.y = j; temp.z = k;
					tempBalls.push_back(temp);
				}
			}
		}
	}

	kernel.x = 0; kernel.y = 0; kernel.z = 0;

	h_balls = thrust::host_vector<double3>(tempBalls.begin(), tempBalls.end());
	ballType.resize(h_balls.size());
}

void Figure::CreateShell()
{
	std::vector<double3> shellVec;
	
	InterfaceForShell::Initialize(r / 2 + 1.8);
	InterfaceForShell::createShell(&shellVec);

	h_shell = thrust::host_vector<double3>(shellVec.begin(), shellVec.end());
}

void Figure::FindActiveBalls()
{
	int n = 0;
	for (int j = 0; j <= branches.size() - 1; j++)
	{
		for (int i = 0; i <= h_balls.size(); i++)
			{
				double3 ker_branch;
				double3 norm_ker_branch;
				double mod_ker_branch;

				ker_branch.x = branches[j].x - kernel.x;
				ker_branch.y = branches[j].y - kernel.y;
				ker_branch.z = branches[j].z - kernel.z;
				mod_ker_branch = sqrt(ker_branch.x*ker_branch.x + ker_branch.y*ker_branch.y + ker_branch.z*ker_branch.z);
				norm_ker_branch.x = ker_branch.x / mod_ker_branch;
				norm_ker_branch.y = ker_branch.y / mod_ker_branch;
				norm_ker_branch.z = ker_branch.z / mod_ker_branch;
				//вектор соединяющий ядро и мишень

				double3 ker_a;

				ker_a.x = h_balls[i].x - kernel.x;
				ker_a.y = h_balls[i].y - kernel.y;
				ker_a.z = h_balls[i].z - kernel.z;
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

				if ((pr_a_tau > r / 4)&(mod_a_rad < actinRadius)) //проекция на вектор направления должна быть положительна и больше d/4, радиальная проекция меньше радиуса вырезаемого цилиндра
				{
					ballType[i] = 1;
					n++;
				}
			}
		}
	std::cout << "Количество активных шаров: " << n << std::endl;
}

void Figure::CreateReceptors()
{
	shellType.resize(h_shell.size());

	int n = 0;
	while (n<receptorCount)
	{
		int x = rand() % (h_shell.size() - 1); // индекс для выбора случайного шара оболочки
		if (shellType[x] != 1)
			{
				shellType[x]=1; // если нет - добавляем (1 - шар рецептор)
				l_receptors[x]=&(h_shell[x]);
				n++;
			}
	} 
}

void Figure::CreateLigands()
{
	for (int j = 0; j <= walls.size()- 1; j++)
	{
		double3 n_vector; // нормальный вектор поверхности
		double h = 0.25; // h - высота расположения лиганда над поверхностью
		n_vector.x = walls[j].x; n_vector.y = walls[j].y; n_vector.z = walls[j].z - walls[j].w / walls[j].z;
		while (h_ligands.size() < ligandCount)
		{
			double x = (0.01 * (rand() % 101)) * 2 * r - r; // задаем круг, у которого радиус - двойной радиус цитоплазмы
			double y = (0.01 * (rand() % 101)) * 2 * r - r;
			if (pow(x, 2) + pow(y, 2) <= pow(r, 2)) // отсеиваем точки, что не входят в круг
			{
				double z = (walls[j].w - walls[j].x * x - walls[j].y * y) / -walls[j].z; // Пока работает только для плоскости где С!= 0 

				double3 n_vect; // нормальный вектор поверхности в точке расположения лиганда
				n_vect.x = n_vector.x + x; n_vect.y = n_vector.y + y; n_vect.z = n_vector.z + z; // переносим вектор к лиганду
				double dist = sqrt(n_vect.x*n_vect.x + n_vector.y*n_vector.y + n_vector.z*n_vector.z);
				n_vect.x = n_vect.x / dist; n_vect.y = n_vect.y / dist; n_vect.z = n_vect.z / dist; // нормируем

				x = x - n_vect.x*h; y = y - n_vect.y*h; z = z - n_vect.z*h; // приподнимаем лиганды над поверхностью (в связи с реализацией коллизий, при моделировании клетка "висит" над плоскостью)

				double3 ligand;
				ligand.x = x;
				ligand.y = y;
				ligand.z = z;

				h_ligands.push_back(ligand);
			}
		}
	}
}

void Figure::InitializeVariable()
{
	nullVec.x = 0; nullVec.y = 0; nullVec.z = 0;

	h_forcesBall = thrust::host_vector<double3>(h_balls.size());
	h_forcesSumBall = thrust::host_vector<double3>(h_balls.size());

	//h_forces_test = thrust::host_vector<double>(h_balls.size());

	h_forcesShell = thrust::host_vector<double3>(h_shell.size());
	h_forcesSumShell = thrust::host_vector<double3>(h_shell.size());

	neighbors = thrust::host_vector<thrust::host_vector<double3*>>(h_shell.size());

	double3 nullLigand; nullLigand.x = 0; nullLigand.y = 0; nullLigand.z = 0;
	std::pair<int, double3> nullPair; nullPair.first = 0; nullPair.second = nullLigand;
	bondsCount = thrust::host_vector<std::pair<int,double3>>(h_shell.size(), nullPair); // изначально адгезионных связей нет
	l_receptors = thrust::host_vector<double3*>(h_shell.size(), NULL); // массив рецепторов изначально пуст
}

void Figure::FindDistanseShell()
{
	thrust::host_vector<double> vecDist(h_shell.size());
	//здесь ссылку при переводе на CUDA надо будет поменять на device_reference

	double distShell = 0;
	const double max = r * 4;
	for (int i = 0; i < h_shell.size(); i++)
	{
		thrust::transform(h_shell.begin(), h_shell.end(), vecDist.begin(), FindDistanse(h_shell[i]));
		
		//переписать в отдельный метод, но пока так оставить
		vecDist[i] = max;

		for (int j = 0; j < 6; j++)
		{
			double dist = max;
			int it = -1;

			for (int k = 0; k < vecDist.size(); k++) {
				if (dist > vecDist[k]) {
					dist = vecDist[k]; 
					it = k;
				}
			}

			neighbors[i].push_back(&(h_shell[it]));
			vecDist[it] = max;
			distShell += dist;
		}

		//добавляем в конец шар относительно которого искали соседей 
		neighbors[i].push_back(&(h_shell[i]));
	}

	a = (distShell / (6 * h_shell.size()));
}

void Figure::CreateStepSchedule()
{
	int restOfTime = (iterationCount - actinTime) % stepTime;
	int stepCount = (iterationCount - actinTime) / stepTime;
	for (int i = 0; i < stepCount - 1; i = i + 2)
	{
		std::pair<int, int> tmpActin(i*stepTime + actinTime, (i + 1)*stepTime + actinTime);
		std::pair<int, int> tmpTraction((i + 1)*stepTime + actinTime, (i + 2)*stepTime + actinTime);
		v_ActinSchedule.push_back(tmpActin);
		v_TractionSchedule.push_back(tmpTraction);
	}

	if (stepTime % 2 == 1)
	{
		std::pair<int, int> lastStep(stepCount*stepTime + actinTime, (stepCount + 1)*stepTime + actinTime);
		v_TractionSchedule.push_back(lastStep);

		std::pair<int, int> lastStep2((stepCount + 1)*stepTime + actinTime, (stepCount + 2)*stepTime + actinTime);
		v_ActinSchedule.push_back(lastStep2);

	}
	else
	{
		std::pair<int, int> lastStep(stepCount*stepTime + actinTime, (stepCount + 1)*stepTime + actinTime);
		v_ActinSchedule.push_back(lastStep);

		std::pair<int, int> lastStep2((stepCount + 1)*stepTime + actinTime, (stepCount + 2)*stepTime + actinTime);
		v_TractionSchedule.push_back(lastStep2);
	}
	for (int i = 0; i < v_ActinSchedule.size(); i++)
	{
		std::cout << "Actin :" << v_ActinSchedule[i].first << " - " << v_ActinSchedule[i].second << std::endl;
	}
	for (int i = 0; i < v_TractionSchedule.size(); i++)
	{
		std::cout << "Traction :" << v_TractionSchedule[i].first << " - " << v_TractionSchedule[i].second << std::endl;
	}
}
//------------------------------------------------------------------------------------------
void Figure::FindForcesBall(int iteration)
{
	//Подсчёт сил ЛД для шаров цитоплазмы

	for (int i = 0; i < h_balls.size(); i++)
	{
		//подсчёт сил ЛД между шарами цитоплазмы
		thrust::transform(h_balls.begin(), h_balls.end(), h_forcesBall.begin(), FindForcesLJ(h_balls[i], maxForceLJ));
		/* это для поиска максимальной силы взаимодействия при нормальной работе (для ограни)
		//--------------------------------------------------------------------------------------
		thrust::transform(h_balls.begin(), h_balls.end(), h_forces_test.begin(), FindForceTest(h_balls[i]));
		double max = h_forces_test[0];
		for (int i = 1; i < h_forces_test.size(); i++)
		{
			if (abs(h_forces_test[i])>max)
			{
				max = h_forces_test[i];
			}
		}
		h_max_forces_test.push_back(max);
		*/
		//--------------------------------------------------------------------------------------


		h_forcesBall[i] = nullVec;
		//суммирование всех сил
		h_forcesSumBall[i] = thrust::reduce(h_forcesBall.begin(), h_forcesBall.end(), nullVec, SumDouble3());
	}

	//подсчёт сил ЛД между шарами цитоплазмы и оболочки
	for (int i = 0; i < h_balls.size(); i++)
	{
		thrust::transform(h_shell.begin(), h_shell.end(), h_forcesShell.begin(), FindForcesShellCytLJ(h_balls[i], maxForceLJ));
		//суммирование всех сил
		double3 force = thrust::reduce(h_forcesShell.begin(), h_forcesShell.end(), nullVec, SumDouble3());

		h_forcesSumBall[i].x += force.x;
		h_forcesSumBall[i].y += force.y;
		h_forcesSumBall[i].z += force.z;
	}

	//подсчёт сил ЛД между шарами цитоплазмы и ядром
	thrust::transform(h_balls.begin(), h_balls.end(), h_forcesBall.begin(), FindForcesLJ(kernel, maxForceLJ));
	thrust::transform(h_forcesSumBall.begin(), h_forcesSumBall.end(), h_forcesBall.begin(), h_forcesSumBall.begin(), PlusDouble3());

	//подсчёт коллизий от стен
	for (int i = 0; i < walls.size(); i++) {
		thrust::transform(h_balls.begin(), h_balls.end(), h_forcesBall.begin(), FindCollision(walls[i]));
		thrust::transform(h_forcesSumBall.begin(), h_forcesSumBall.end(), h_forcesBall.begin(), h_forcesSumBall.begin(), PlusDouble3());
	}

	//Подсчёт сил тяжести
	thrust::transform(h_forcesSumBall.begin(), h_forcesSumBall.end(), h_forcesSumBall.begin(), PlusG());


	if (iteration >= v_ActinSchedule[twoSteps].first && iteration < v_ActinSchedule[twoSteps].second)
	{
		for (int i = 0; i < branches.size(); i++)
		{
			//передаю лишние параметры в функтор - переделать----------------------------------------
			thrust::transform(h_balls.begin(), h_balls.end(), ballType.begin(), h_forcesBall.begin(), FindActinForce(branches[i], kernel, actinRadius, actinForce, r));
			thrust::transform(h_forcesSumBall.begin(), h_forcesSumBall.end(), h_forcesBall.begin(), h_forcesSumBall.begin(), PlusDouble3());
		}
	}
	//Подсчет силы подтягивания
	if (iteration >= v_TractionSchedule[twoSteps].first && iteration < v_TractionSchedule[twoSteps].second)
	{
		for (int i = 0; i < branches.size(); i++)
		{
			thrust::transform(h_balls.begin(), h_balls.end(), ballType.begin(), h_forcesBall.begin(), FindTractionForce(branches[i], kernel, tractionForce));
			thrust::transform(h_forcesSumBall.begin(), h_forcesSumBall.end(), h_forcesBall.begin(), h_forcesSumBall.begin(), PlusDouble3());
		}
	}
	if (iteration >= 0 && iteration < actinTime)
	{
	}
	else
	{
		if ((iteration >= v_ActinSchedule[twoSteps].first || iteration < v_TractionSchedule[twoSteps].second) == 0)
		{
			twoSteps++;
		}
	}
}

void Figure::FindForcesShell()
{
	//Подсчёт сил ЛД между шарами оболочки и цитоплазмы
	for (int i = 0; i < h_shell.size(); i++) {
		thrust::transform(h_balls.begin(), h_balls.end(), h_forcesBall.begin(), FindForcesShellCytLJ(h_shell[i], maxForceLJ));
		h_forcesSumShell[i] = thrust::reduce(h_forcesBall.begin(), h_forcesBall.end(), nullVec, SumDouble3());
	}

	//Подсчёт сил между шарами оболочки
	thrust::transform(neighbors.begin(), neighbors.end(), h_forcesShell.begin(), FindForcesShellLJ(a, maxForceLJ));
	thrust::transform(h_forcesSumShell.begin(), h_forcesSumShell.end(), h_forcesShell.begin(), h_forcesSumShell.begin(), PlusDouble3());

	//Подсчёт сил ЛД между шарами оболочки и ядром
	thrust::transform(h_shell.begin(), h_shell.end(), h_forcesShell.begin(), FindForcesShellCytLJ(kernel, maxForceLJ));
	thrust::transform(h_forcesSumShell.begin(), h_forcesSumShell.end(), h_forcesShell.begin(), h_forcesSumShell.begin(), PlusDouble3());

	//подсчёт коллизий от стен
	for (int i = 0; i < walls.size(); i++) {
		thrust::transform(h_shell.begin(), h_shell.end(), h_forcesShell.begin(), FindCollision(walls[i]));
		thrust::transform(h_forcesSumShell.begin(), h_forcesSumShell.end(), h_forcesShell.begin(), h_forcesSumShell.begin(), PlusDouble3());
	}

	//Подсчёт сил тяжести
	thrust::transform(h_forcesSumShell.begin(), h_forcesSumShell.end(), h_forcesSumShell.begin(), PlusG());
	//Подсчет силы адгезии
	for (int i = 0; i < h_ligands.size(); i++)
	{
		thrust::transform(l_receptors.begin(), l_receptors.end(), bondsCount.begin(), bondsCount.begin(), FindAdhesionContacts(h_ligands[i]));
		thrust::transform(l_receptors.begin(), l_receptors.end(), bondsCount.begin(), h_forcesShell.begin(), FindAdhesionForce(h_ligands[i]));
		thrust::transform(h_forcesSumShell.begin(), h_forcesSumShell.end(), h_forcesShell.begin(), h_forcesSumShell.begin(), PlusDouble3());
	}

}

void Figure::FindForcesKernel(int iteration)
{
	//Подсчёт сил ЛД между ядром и шарами цитоплазмы
	thrust::transform(h_balls.begin(), h_balls.end(), h_forcesBall.begin(), FindForcesLJ(kernel, maxForceLJ));
	forceSumKernel = thrust::reduce(h_forcesBall.begin(), h_forcesBall.end(), nullVec, SumDouble3());

	forceSumKernel.x = -forceSumKernel.x;
	forceSumKernel.y = -forceSumKernel.y;
	forceSumKernel.z = -forceSumKernel.z;

	//Подсчёт сил ЛД между ядром и шарами оболочки
	thrust::transform(h_shell.begin(), h_shell.end(), h_forcesShell.begin(), FindForcesShellCytLJ(kernel, maxForceLJ));
	forceKernel = thrust::reduce(h_forcesShell.begin(), h_forcesShell.end(), nullVec, SumDouble3());
	
	forceSumKernel.x += -forceKernel.x;
	forceSumKernel.y += -forceKernel.y;
	forceSumKernel.z += -forceKernel.z;

	//Подсчёт коллизий от стен
	for (int i = 0; i < walls.size(); i++) {
		forceKernel = FindCollision(walls[i])(kernel);

		forceSumKernel.x += forceKernel.x;
		forceSumKernel.y += forceKernel.y;
		forceSumKernel.z += forceKernel.z;
	}

	//Подсчёт сил тяжести
	forceSumKernel = PlusG()(forceSumKernel);

	//Подсчет силы подтягивания
	if (iteration > tractionTime)
	{
		for (int i = 0; i < branches.size(); i++)
		{
			forceKernel = FindTractionForce(branches[i],kernel, tractionForce)(kernel, 0);

			forceSumKernel.x += forceKernel.x;
			forceSumKernel.y += forceKernel.y;
			forceSumKernel.z += forceKernel.z;
		}
	}
}

void Figure::FindRadiusVector()
{
	//сделать через transform
	for (int i = 0; i < h_balls.size(); i++) {
		h_balls[i].x += h_forcesSumBall[i].x / 9;
		h_balls[i].y += h_forcesSumBall[i].y / 9;
		h_balls[i].z += h_forcesSumBall[i].z / 9;
	}

	for (int i = 0; i < h_shell.size(); i++) {
		h_shell[i].x += h_forcesSumShell[i].x / 9;
		h_shell[i].y += h_forcesSumShell[i].y / 9;
		h_shell[i].z += h_forcesSumShell[i].z / 9;
	}

	kernel.x += forceSumKernel.x / 9;
	kernel.y += forceSumKernel.y / 9;
	kernel.z += forceSumKernel.z / 9;
}

void Figure::ToFile(int iteration)
{
	//h_balls = d_balls;

	std::string path = "C:/Result/simple_balls/mca_out.csv." + std::to_string(iteration);
	std::string pathShell = "C:/Result/shell/mca_out.csv." + std::to_string(iteration);
	std::string pathKernel = "C:/Result/kernel/mca_out.csv." + std::to_string(iteration);
	std::string pathReceptor = "C:/Result/receptor/mca_out.csv." + std::to_string(iteration);
	std::string pathLigand = "C:/Result/ligand/mca_out.csv." + std::to_string(iteration);
	//std::string pathForce = "C:/Result/force/mca_out.csv.";

	std::ofstream out(path.c_str()); // каждую итерацию пишем в новый файл. В данном случае путь указан для обычных шаров
	std::ofstream outShell(pathShell.c_str());
	std::ofstream outKernel(pathKernel.c_str());
	std::ofstream outReceptor(pathReceptor.c_str());
	std::ofstream outLigand(pathLigand.c_str());
	//std::ofstream outForce(pathForce.c_str(), std::ios_base::app);

	out		    << "x,y,z\n";
	outShell    << "x,y,z\n";
	outKernel   << "x,y,z\n";
	outReceptor << "x,y,z\n";
	outLigand   << "x,y,z\n";

	for (int i = 0; i < h_balls.size(); i++) {
		out << h_balls[i].x << "," << h_balls[i].y << "," << h_balls[i].z << std::endl;
	}

	for (int i = 0; i < h_shell.size(); i++) {
		if (shellType.size() != 0) 
		{
			if (shellType[i] == 0)
			{
				outShell << h_shell[i].x << "," << h_shell[i].y << "," << h_shell[i].z << std::endl;
			}
			else {
				outReceptor << h_shell[i].x << "," << h_shell[i].y << "," << h_shell[i].z << std::endl;
			}
		}
		else
		{
			outShell << h_shell[i].x << "," << h_shell[i].y << "," << h_shell[i].z << std::endl;
		}
	}

	outKernel << kernel.x << "," << kernel.y << "," << kernel.z << std::endl;

	for (int i = 0; i < h_ligands.size(); i++) {
		outLigand << h_ligands[i].x << "," << h_ligands[i].y << "," << h_ligands[i].z << std::endl;
	}

	/*
	if (iteration > 0)
	{
		outForce << h_max_forces_test[iteration] << std::endl;
	}
	*/

	out.close();
	outShell.close();
	outKernel.close();
	outReceptor.close();
	outLigand.close();
	//outForce.close();
}