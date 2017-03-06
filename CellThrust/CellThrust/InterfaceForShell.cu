//
//  CApplication.h
//  MCA
//
//  Created by Dmitry Senushkin on 17.11.15.
//  Copyright (c) 2015 Dmitry Senushkin. All rights reserved.
//
#include "InterfaceForShell.cuh"

int InterfaceForShell::id[20][3] = { { 1,0,4 },{ 0,4,9 },{ 5, 4, 9 },{ 4, 8, 5 },{ 4, 1, 8 },
{ 8, 1, 10 },{ 8, 10, 3 },{ 5, 8, 3 },{ 5, 3, 2 },{ 2, 3, 7 },
{ 7, 3, 10 },{ 7, 10, 6 },{ 7, 6, 11 },{ 11, 6, 0 },{ 0, 6, 1 },
{ 6, 10, 1 },{ 9, 11, 0 },{ 9, 2, 11 },{ 9, 5, 2 },{ 7, 11, 2 } };

double3 InterfaceForShell::radiusVectorOfCenter;
double InterfaceForShell::radius;//5
double InterfaceForShell::V;
double InterfaceForShell::W;
double InterfaceForShell::vertex[12][3];

void InterfaceForShell::Initialize(double rad)
{
	radius = rad;

	radiusVectorOfCenter.x = 0;
	radiusVectorOfCenter.y = 0;
	radiusVectorOfCenter.z = 0;

	V = radius * cos(3. * atan(1.) / 2.5);
	W = radius * sin(3. * atan(1.) / 2.5);

	// красиво не получилось(( убогий костыль
	vertex[0][0] = -V, vertex[0][1] = 0., vertex[0][2] = W;
	vertex[1][0] = V, vertex[1][1] = 0., vertex[1][2] = W;
	vertex[2][0] = -V, vertex[2][1] = 0., vertex[2][2] = -W;
	vertex[3][0] = V, vertex[3][1] = 0., vertex[3][2] = -W;
	vertex[4][0] = 0., vertex[4][1] = W, vertex[4][2] = V;
	vertex[5][0] = 0., vertex[5][1] = W, vertex[5][2] = -V;
	vertex[6][0] = 0., vertex[6][1] = -W, vertex[6][2] = V;
	vertex[7][0] = 0., vertex[7][1] = -W, vertex[7][2] = -V;
	vertex[8][0] = W, vertex[8][1] = V, vertex[8][2] = 0.;
	vertex[9][0] = -W, vertex[9][1] = V, vertex[9][2] = 0.;
	vertex[10][0] = W, vertex[10][1] = -V, vertex[10][2] = 0.;
	vertex[11][0] = -W, vertex[11][1] = -V, vertex[11][2] = 0.;
}

void InterfaceForShell::Scale(double *v)
{
	double d = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

	if (d == 0.)
	{
		//MessageBox(0, L"Zero length vector", L"Error", MB_OK);
	}
	else
	{
		v[0] /= d;
		v[1] /= d;
		v[2] /= d;
		v[0] *= InterfaceForShell::radius;
		v[1] *= InterfaceForShell::radius;
		v[2] *= InterfaceForShell::radius;
	}
	return;
}

void InterfaceForShell::getNorm(double v1[3], double v2[3], double *out)
{
	out[0] = v1[1] * v2[2] - v1[2] * v2[1];
	out[1] = v1[2] * v2[0] - v1[0] * v2[2];
	out[2] = v1[0] * v2[1] - v1[1] * v2[0];
	Scale(out);

}

void InterfaceForShell::Split(double *v1, double *v2, double *v3, long depth, std::vector<double3> *shell)
{
	double v12[3], v23[3], v31[3];

	if (depth == 0)
	{
		setTria(v1, v2, v3, shell);

		return;
	}

	for (int i = 0; i< 3; i++)
	{
		v12[i] = v1[i] + v2[i];
		v23[i] = v2[i] + v3[i];
		v31[i] = v3[i] + v1[i];

	}

	Scale(v12);
	Scale(v23);
	Scale(v31);
	Split(v1, v12, v31, depth - 1, shell);
	Split(v2, v23, v12, depth - 1, shell);
	Split(v3, v31, v23, depth - 1, shell);
	Split(v12, v23, v31, depth - 1, shell);

}

void InterfaceForShell::setTria(double *v1, double *v2, double *v3, std::vector<double3> *shell)
{
	int counterVertex1 = 0;
	int counterVertex2 = 0;
	int counterVertex3 = 0;

	double3 vertex1, vertex2, vertex3;
	vertex1.x = v1[0] + radiusVectorOfCenter.x;
	vertex1.y = v1[1] + radiusVectorOfCenter.y;
	vertex1.z = v1[2] + radiusVectorOfCenter.z;

	vertex2.x = v2[0] + radiusVectorOfCenter.x;
	vertex2.y = v2[1] + radiusVectorOfCenter.y;
	vertex2.z = v2[2] + radiusVectorOfCenter.z;

	vertex3.x = v3[0] + radiusVectorOfCenter.x;
	vertex3.y = v3[1] + radiusVectorOfCenter.y;
	vertex3.z = v3[2] + radiusVectorOfCenter.z;

	if (shell->size() != 0)
	{
		for (auto i = 0; i < shell->size(); i++)
		{
			if (shell->at(i).x == vertex1.x && shell->at(i).y == vertex1.y && shell->at(i).z == vertex1.z)
				counterVertex1++;

			if (shell->at(i).x == vertex2.x && shell->at(i).y == vertex2.y && shell->at(i).z == vertex2.z)
				counterVertex2++;

			if (shell->at(i).x == vertex3.x && shell->at(i).y == vertex3.y && shell->at(i).z == vertex3.z)
				counterVertex3++;
		}
	}

	if (!counterVertex1)
		shell->push_back(vertex1);
	if (!counterVertex2)
		shell->push_back(vertex2);
	if (!counterVertex3)
		shell->push_back(vertex3);
}

void InterfaceForShell::createShell(std::vector<double3> *shell)
{
	for (int i = 0; i < 20; i++)
	{
		Split(InterfaceForShell::vertex[InterfaceForShell::id[i][0]], InterfaceForShell::vertex[InterfaceForShell::id[i][1]], InterfaceForShell::vertex[InterfaceForShell::id[i][2]], 3, shell);
	}
}