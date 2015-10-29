#ifndef GaussSeidel_HPP
#define GaussSeidel_HPP

//*************************
// Librerias necesarias
#include <iostream>
#include <cmath>
#include "CSR.hpp"
#include "Vector.hpp"
#include "Timer.hpp"
//*************************

class GaussSeidel
{
	private:
		// Variables necesarias
		int iterMax;
		int iter;
		double tol;
		double err;
		double etime;
		Timer timer;
	public:
		// Constructor de la clase
		GaussSeidel():iterMax(2000),tol(1e-6),iter(0),err(0),etime(0){};
		template<class T>
		void solve(T &A, Vector &x, Vector &b)
		{
			// Variables
			timer.tic(); // Comenzamos a medir el tiempo
			int j = 0;
			double diag = 0.0;
			double suma = 0.0;
			int n = x.size();
			// Comenzamos con el metodo de GaussSeidel
			for (int k = 0; k < iterMax; ++k)
			{
				for (int i = 0; i < n; ++i)
				{
					
					suma = 0.0;
					for (int l = A.irow[i]; l < A.irow[i+1]; ++l)
					{
						j = A.col[l];
						if (i != j)
						{
							suma = suma + A.data[l] * x[j];
						}
						else
						{
							diag = A.data[l];
						}
					}
					x[i] = (b[i] - suma) /diag;
				}
				// Iteraciones realizadas
				iter++;
				// Revisamos la tolerancia de paro
			}
			// Terminamos de medir el tiempo usado
			timer.toc();
      		etime = timer.etime();
		}
		// Devuelve el numero maximo de iteraciones realizadas
		int its(){return iter;};
		// Devuelve el error 
   		double error() {return err;};
   
   		void maxIts(int it){ iterMax = it;};
   		void tole(double t) {tol=t;};
   		void report()
   		{
      		std::cout << "GaussSeidel";
      		std::cout << "\n";
      		std::cout << "  Iteraciones  "<< iter << std::endl;
      		std::cout << "  Tiempo de calculo "<< etime << " msec" << std::endl;
   		};
};

#endif // GaussSeidel_HPP
