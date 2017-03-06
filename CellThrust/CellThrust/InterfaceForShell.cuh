#include <vector>
#include <thrust/host_vector.h>

class InterfaceForShell
{
private:
	static int id[20][3]; // массив индекса вершин икосаэдра
	static double V, W;		// параметры икосаэдра
	static double vertex[12][3];
	static double3 radiusVectorOfCenter; // center of cell in time of create
	static double radius;
	static void getNorm(double v1[3], double v2[3], double *out);	// норма вектора (единичный вектор нормали)
	static void Split(double *v1, double *v2, double *v3, long depth, std::vector<double3> *shell);
	static void setTria(double *v1, double *v2, double *v3, std::vector<double3> *shell);
	static void Scale(double *v); // нормирование вектора нормали ( любого вектора )

public:
	static void Initialize(double rad);
	static void createShell(std::vector<double3> *shell); // создание оболочки путем деления треугольников
};
