#ifndef CSR_HPP
#define CSR_HPP

//**********************
// Librerias necesarias
#include <iostream>
#include "COO.hpp"
#include "Vector.hpp"
//**********************

class CSR
{
	public:
		// Arreglos necesarios
		double *data;
		int *col;
		int *irow;
		int *idiag;
		// Variables necesarias
		int nnzMax;
		int nnz;
		int n;
		// Metodos utiles
		void reserveMemory();
		void freeMemory();
	public:
		// Constructores
		CSR(){};
		CSR(int m):n(m),nnzMax(m*30),nnz(0){reserveMemory();zeros();};
		// Constructor copia
		CSR(CSR  const &Mtx):n(Mtx.n),nnzMax(Mtx.nnzMax),nnz(Mtx.nnz)
		{
			reserveMemory();
			// Copiamos los elementos de data y col
			for (int i = 0; i < nnz; ++i)
			{
				data[i] = Mtx.data[i];
				col[i] = Mtx.col[i];
			}
			// Copiamos los elementos de irow
			for (int i = 0; i < n; ++i)
			{
				irow[i] = Mtx.irow[i];
			}
		};
		// Libera la memoria
		//~CSR(){freeMemory();};
		// Funciones
		// Llenamos los arreglos con cero
		void zeros();
		// Convertimos una matriz COO a una CSR
		void convert(COO const &coo);
		// Busca un elemento usando busqueda binaria
		int binarySearch(int lmin,int lmax,int j);
		// Busca un elemento dentro de CSR
		double search(int i, int j);
		// Devuelve un valor de data
		inline double dat(int i){return data[i];};
		// Devuelve un valor de col
		inline int co(int i){return col[i];};
		// Devuelve un valor de irow
		inline int iro(int i){return irow[i];};
		// Sobrecargamos la multiplicacion Matriz-Vector
		Vector operator *(Vector const &x);
		// Iteracion de Jacobi
		void iterJacobi(int i, Vector &x);
		// Iteracion de Gauss-Seidel
		void iterGaussSeidel(int i, Vector &x, Vector const &b);
		// Sobrecargamos el operador = en CSR
		void operator =(CSR const &Mtx);
};

#endif // CSR_HPP
