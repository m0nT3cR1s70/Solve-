#include "COO.hpp"

// Reservamos el espacio en memoria suficiente
void COO :: reserveMemory()
{
	// Solicitamos espacio para data
	data = new double[nnzMax];
	// Verificamos que el espacio fue asignado
	if (data == nullptr)
	{
		std::cout << "ERROR: ESPACIO NO ASIGNADO: COO: data" << std::endl;
		exit(0);
	}
	// Solicitamos espacio para row
	row = new int[nnzMax];
	// Verificamos que el espacio fue asignado
	if (row == nullptr)
	{
		std::cout << "ERROR: ESPACIO NO ASIGNADO: COO: row" << std::endl;
	}
	// Solicitamos espacio para col
	col = new int[nnzMax];
	// Verificamos que el espacio fue asignado
	if (col == nullptr)
	{
		std::cout << "ERROR: ESPACIO NO ASIGNADO: COO: row" << std::endl;
	}
}
// Liberamos el espacio asignado en memoria
void COO :: freeMemory()
{
	delete [] data;
	delete [] row;
	delete [] col;
}
// Llenamos de ceros los arreglos
void COO :: zeros()
{
	for (int i = 0; i < nnzMax; ++i)
	{
		data[i] = 0.0;
		row[i] = 0;
		col[i] = 0;
	}
}
// Insertamos los elementos en los arreglos
void COO :: insert(int i, int j, double val)
{
	// Verificamos que aun se tenga espacio para almacenar elementos
	if (nnz > nnzMax)
	{
		std::cout << "ERROR: MEMORIA ESTIMADA NO SUFICENTE" << std::endl;
		exit(0);
	}
	// Insertamos en caso de tener espacio
	data[nnz] = val;
	row[nnz] = i;
	col[nnz] = j;
	nnz++;
}
// Se imprime en forma de matriz para verficar 
void COO :: impMatrix()
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			bool x = false;
			for (int l = 0; l < nnz; ++l)
			{
				//std::cout << i << "," << j << "," << row[i] << "," << col[l] << std::endl;
				if ((row[l] == i) and (col[l] == j))
				{
					std::cout << data[l] << "\t";
					x = true;
					continue;
				}
			}
			if(x == false)std::cout << 0 << "\t";
		}
		std::cout << std::endl;
	}
}