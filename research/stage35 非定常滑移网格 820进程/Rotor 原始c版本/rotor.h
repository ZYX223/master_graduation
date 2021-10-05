#include<iostream>
#include<stdlib.h>
using namespace std;
#ifndef ARRAY_H
#define ARRAY_H
template<class T>
T** malloc2D(int x_size, int y_size)
{
	T** ite;
	ite = (T**)malloc(y_size * sizeof(T*));
	ite[0] = (T*)malloc(y_size*x_size * sizeof(T));
	for (int i = 1; i < y_size; i++)
	{
		ite[i] = ite[i - 1] + x_size;
	}
	if (NULL == ite[0])
		return NULL;
	return ite;
}
template<class T>
T*** malloc3D(int x_size, int y_size, int z_size)
{
	T*** ite;
	int size1 = z_size*y_size;
	ite = (T***)malloc(z_size * sizeof(T**));
	ite[0] = (T**)malloc(size1 * sizeof(T*));
	ite[0][0] = (T*)malloc(size1*x_size * sizeof(T));
	for (int i = 1; i < z_size; i++)
		ite[i] = ite[i - 1] + y_size;
	for (int i = 1; i < size1; i++)
		ite[0][i] = ite[0][i - 1] + x_size;
	if (NULL == ite[0][0])
		return NULL;
	return ite;
}
template<class T>
T**** malloc4D(int x_size, int y_size, int z_size, int n_size)
{
	T**** ite;
	int size1 = n_size*z_size;
	int size2 = size1*y_size;
	ite = (T****)malloc(n_size * sizeof(T**));
	ite[0] = (T***)malloc(size1 * sizeof(T*));
	ite[0][0] = (T**)malloc(size2 * sizeof(T));
	ite[0][0][0] = (T*)malloc(size2 *x_size * sizeof(T));
	for (int i = 1; i < n_size; i++)
		ite[i] = ite[i - 1] + z_size;
	for (int i = 1; i < size1; i++)
		ite[0][i] = ite[0][i - 1] + y_size;
	for (int i = 1; i < size2; i++)
		ite[0][0][i] = ite[0][0][i - 1] + x_size;
	if (NULL == ite[0][0][0])
		return NULL;
	return ite;
}
template<class T>
bool free4D(T ****ite)
{
	free(ite[0][0][0]);
	free(ite[0][0]);
	free(ite[0]);
	free(ite);
	if (NULL == ite)
		return true;
	return false;
}
template<class T>
bool free3D(T ***ite)
{
	free(ite[0][0]);
	free(ite[0]);
	free(ite);
	if (NULL == ite)
		return true;
	return false;
}
template<class T>
bool free2D(T **ite)
{
	free(ite[0]);
	free(ite);
	if (NULL == ite)
		return true;
	return false;
}
#endif