#ifndef MAIN_H
#define MAIN_H

#include<iostream>
#include <cstdlib>
#include<fstream>
#include<math.h>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include "Global.h"
#include <algorithm>
#include <cstdlib>
#include<iomanip>
#include<ctime>
#include <cstdarg>
#include "Input_Class.h" 

using namespace std;


int Find_length(const std :: string& file_name, int rskip);
void Combine_Output(Input_Class* Input);
  
//void Read_PaQEvo(const std :: string& file_name, int n_t, double* t, double** P, double ** Q, double *V);

//void Calc_inter_dist(int n_t, double* t, double** Q, double** R);

//void Read_Traj_Tot(const std :: string& file_name, int nTraj, double** Traj);
//void Read_Traj_Tot(const std :: string& file_name, int nTraj, vector<vector<double>>& Traj);

//void Read_PaQ_Sol(const std :: string& file_name, int iProc, vector<vector<double>>& PaQSol);

//void Read_QEvo(const std :: string& file_name, int n_t, double* t, double ** Q);

double static get_diff(double a, double b, double c)
{
	double d_max = std::max(std::max(a,b),std::max(a,c));
	double d_min = std::min(std::min(a,b),std::min(a,c));
	
	return (d_max-d_min);
}

int static get_max(double a, double b, double c)
{
	double arr[3] = {a,b,c};
	double max = a;
	int j = 0;
	for(int i=0; i<3; i++)
	{
		if(arr[i] > max)
		{
			j = i;
			max = arr[i];
		}
	}
	
	return j;
	
}

int static get_min(double arr[], int n)
{
	//double arr[3] = {a,b,c};
	double min = arr[0];
	int j = 0;
	for(int i=0; i<n; i++)
	{
		if(arr[i] < min)
		{
			j = i;
			min = arr[i];
		}
	}
	
	return j;
	
}

double static dot(double a[3], double b[3])
{
	double res = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
	return res;
}

#endif /* MAIN_H */

