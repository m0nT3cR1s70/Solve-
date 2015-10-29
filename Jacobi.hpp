#ifndef JACOBI_HPP
#define JACOBI_HPP

//*************************
// Librerias necesarias
#include <iostream>
#include <cmath>
#include "CSR.hpp"
#include "Vector.hpp"
#include "Timer.hpp"
//*************************

class Jacobi
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
		Jacobi():iterMax(2000),tol(1e-6),iter(0),err(0),etime(0){};
		template<class T>
		void solve(T &A, Vector &x, Vector &b)
		{
			// Variables
			timer.tic(); // Comenzamos a medir el tiempo
			int j = 0;
			double diag = 0.0;
			int n = x.size();
			// Comenzamos con el metodo de Jacobi
			for (int k = 0; k < iterMax; ++k)
			{
				for (int i = 0; i < n; ++i)
				{
					x[i] = b[i];
					//std::cout  << i << ": " << x.asizeVect << std::endl;
					for (int l = A.irow[i]; l < A.irow[i+1]; ++l)
					{
						j = A.col[l];
						if( i != j) 
						{
							x[i] = x[i] - A.data[l]*x[j];
						}
						else	
						{
							diag = A.data[l];
						}
						/*std::cout<<"----------------"<<std::endl;
						std::cout<<"Iteracion: "<<i<<std::endl;
						x.print();
						std::cout<<"----------------"<<std::endl;*/
					}
					x[i] = x[i]/diag;
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
      		std::cout << "Jacobi";
      		std::cout << "\n";
      		std::cout << "  Iteraciones  "<< iter << std::endl;
      		std::cout << "  Tiempo de calculo "<< etime << " msec" << std::endl;
   		};
};

#endif // JACOBI_HPP
