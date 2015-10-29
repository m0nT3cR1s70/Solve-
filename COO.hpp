#ifndef COO_HPP
#define COO_HPP

//***********************
// Librerias necesarias
#include <iostream>
//***********************

class COO
{
	public:
		// Arreglos para almacenar
		double *data;
		int *row;
		int *col;
		// Variables utiles
		int nnz;
		int nnzMax;
		int n;
	private:
		// Metodos utiles
		void reserveMemory();
		void freeMemory();
	public:
		// Constructor vacio
		COO(){};
		// Constructor que almacena un elemento
		COO(int m):n(m),nnzMax(m*30),nnz(0){reserveMemory();zeros();};
		// El contenido de los arreglos es de cero
		void zeros();
		// Destructor de la clase
		~COO(){reserveMemory();};
		// Inserccion de elementos
		void insert(int i, int j, double val);
		// Impresion en formato matriz
		void impMatrix();
};

#endif // COO_HPP
