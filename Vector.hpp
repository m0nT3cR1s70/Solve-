#ifndef VECTOR_HPP
#define VECTOR_HPP

//******************************
// Librerias necesarias
#include <iostream>
#include <cmath>
#include <fstream>
//******************************

class Vector
{
	public:
		// Variables Necesarias
		double *data; // Arreglo para almacenar los datos
		int n;  // Largo del vector
		// Metodos para asignar y liberar memoria
	private:
		void reserveMemory();
		void freeMemory();
	public:
		// Constructor Vacio
		//Vector():n(0),data(nullptr){};
		// Constructor para un vector lleno de ceros
		Vector(int o):n(o){reserveMemory();zeros();};
		// Constructor copia 
		Vector(Vector const &v):n(v.n){reserveMemory();for (int i = 0; i < n; ++i)data[i]=v[i];}
		// Devuelve el largo de un vector
		inline int size(){return n;};
		// Calcula la norma euclidiana de un vector
		double norm();
		// Inicializa un vector con ceros
		void zeros();
		// Almacena en un archivo el contenido de un vector
		void saveData(std::string name);
		// Imprime en pantalla el contenido de un vector
		void print();
		// Sobrecarga de operadores necesarias
		// Funcion get y funcion set del vector
		double operator  [](int i) const {return data[i];}; 
		double &operator [](int i) {return data[i];};
		// Operaciones vectoriales principales
		// Suma vectorial
		Vector operator +(Vector const &v);
		// Resta Vectorial
		Vector operator -(Vector const &v);
		// Escalar por vector
		friend Vector operator *(double escalar, Vector const &v);
		// Producto interno 
		double operator *(Vector const &v);
		// Sobrecarga del operador igual
		// Inicializa un vector con un valor dado
		Vector operator =(double escalar);
		// Iguala un vector con otro vector
		void operator =(Vector const &v);
		// Destructor
		~Vector(){freeMemory();};
	
};

#endif // VECTOR_HPP
