#ifndef SOR_HPP
#define SOR_HPP

//*************************
// Librerias necesarias
#include <iostream>
#include <cmath>
#include "CSR.hpp"
#include "Vector.hpp"
#include "Timer.hpp"
//*************************

class SOR
{
	private:
		// Variables necesarias
		int iterMax;
		int iter;
		double tol;
		double err;
		double etime;
		double omega;
		Timer timer;
	public:
		// Constructor de la clase
		SOR():iterMax(2000),tol(1e-6),iter(0),err(0),etime(0),omega(0){};
		template<class T>
		void solve(T &A, Vector &x, Vector &b)
		{
			// Variables
			timer.tic(); // Comenzamos a medir el tiempo
			int j = 0;
			double diag = 0.0;
			double suma = 0.0;
			double comp = 0.0;
			double omaga = 0.5;
			int n = x.size();
			// Comenzamos con el metodo de SOR
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
					comp = (b[i] - suma)/diag;
					x[i] = x[i]+omega*(comp-x[i]);
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
   		void setOmega(double o) {omega = o;}
   		void report()
   		{
      		std::cout << "SOR";
      		std::cout << "\n";
      		std::cout << "  Iteraciones  "<< iter << std::endl;
      		std::cout << "  Tiempo de calculo "<< etime << " msec" << std::endl;
   		};
};

#endif // SOR_HPP
