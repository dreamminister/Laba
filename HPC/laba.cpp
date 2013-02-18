// laba.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include <sstream>
//#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

double xmax = 8, ymax = 8, tmax = 10; 
double dx = 1000 / xmax, dy = 1000 / ymax, dt = 1 / tmax;
double k = 30, d = 30, f = 100;
//#pragma omp parallel


double detectCenter(int x, int y)
{
	if ( (x == int(xmax/2)) && (y == int(ymax/2)) ) {
		return 100;
	}
	else 
		return 0;
}
double*** alloc() {
double*** array = new double** [int(tmax)];
for (int i = 0; i < tmax; i++) {
array[i] = new double* [int(xmax)];
for (int j = 0; j < xmax; j++) {
array[i][j] = new double [int(ymax)];
}
}
return array;
}
double ***u = alloc();
	double *** utemp = alloc();
		double *** Res = alloc();

int _tmain(int argc, _TCHAR* argv[])
{

	
	for (int t = 0; t < tmax; t++){
		for (int i = 0; i < xmax; i++){
			for (int j = 0; j < ymax; j++)
				{
					u[t][i][j] = 0; utemp[t][i][j]=0; Res[t][i][j]=0;
			}}
	}
//////////INIT-END///////////////////////////////////////
	
	
				//u[0][i][j] = utemp[0][i][j] + dt * ( (k/(dx*dx))*(utemp[0][i+1][j]-2*utemp[0][i][j]+utemp[0][i-1][j])+(k/(dy*dy))*(utemp[0][i][j+1]-2*utemp[0][i][j]+utemp[0][i][j-1])+ (c/dx)*(utemp[0][i+1][j]-utemp[0][i][j])+ (c/dy)*(utemp[0][i][j+1]-utemp[0][i][j])+d*utemp[0][i][j]+f);
				Res[0][int(xmax/2)][int(ymax/2)] = f;//(u[t+1][i][j]-utemp[t][i][j])/dt;
				//utemp[t+1][i][j] = u[t+1][i][j];
			

	




	for (int t = 0; t < tmax-1; t++)
	{
		double Cx = 0.5*cos(t*dt);
		double Cy = 0.5*sin(t*dt);
		for ( int i = 1; i< xmax-1; i++)
		{
			for (int j=1; j<ymax-1; j++)
			{
				u[t+1][i][j] = utemp[t][i][j] + dt * ( (k/(dx*dx))*(utemp[t][i+1][j]-2*utemp[t][i][j]+utemp[t][i-1][j])+(k/(dy*dy))*(utemp[t][i][j+1]-2*utemp[t][i][j]+utemp[t][i][j-1])-(Cx/dx)*(utemp[t][i+1][j]-utemp[t][i][j])-(Cy/dy)*(utemp[t][i][j+1]-utemp[t][i][j])+d*utemp[t][i][j]+detectCenter(i, j));
				Res[t+1][i][j] = (u[t+1][i][j]-utemp[t][i][j])/dt;
				utemp[t+1][i][j] = u[t+1][i][j];
			}
		}

	}


	
		for ( int i = 0; i< xmax; i++)
		{
			for (int j=0; j<ymax; j++)
			{
				cout<<Res[2][i][j]<<"\t"; 
				
			}
			cout<<endl;
		}

	/////////////////////////////////////////////////
		string fname;
		fname = "data.txt";
	std::ofstream outF(fname.c_str() , ios::out);
      
      
     // outF << xmax << ' ' << ymax << endl;
		
    
      for (int i = 0; i <int( xmax); i++)
      {
         for(int j = 0; j < int(ymax); j++)
         {
            
            outF << Res[2][i][j] << " ";
         }
         outF << endl;
      }
      outF.close ();


	
	//cout << "Hello";
	cin.get();




	return 0;
}

