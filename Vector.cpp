#include "Vector.hpp"

// Asigna memoria al Vector
void Vector :: reserveMemory()
{
	// Comprobamos que sea un largo valido
	if (n < 1 )
	{
		std::cout << "ERROR: LARGO: " << n << " NO VALIDO: Vector reserveMemory" << std::endl;
		exit(0);
	}
	// Asignamos memoria a data
	data = new double[n];
	// Verificamos que la memoria fue asignada
	if (data == nullptr)
	{
		std::cout << "ERROR: MEMORIA INSUFICIENTE" << std::endl;
	}
}
// Liberamos la memoria asignada
void Vector :: freeMemory()
{
	delete [] data;
}
// Obtiene la norma de un vector
double Vector :: norm()
{
    double euclidiana = 0.0;
    for (int i = 0; i < n; ++i)
    {
        euclidiana += data[i]*data[i];
    }
    return sqrt(euclidiana);
}
// Funcion para llenar un vector de ceros
void Vector :: zeros()
{
	for (int i = 0; i < n; ++i)
	{
		data[i] = 0;
	}
}
// Almacena el contenido de un vector en un archivo
void Vector :: saveData(std::string name)
{
	// Almacenar en un archivo
	std::ofstream file;
	// Nombre del archivo
	file.open(name.c_str());
	// Recorremos el vector
	for (int i = 0; i < n; ++i)
	{
		file.width(9);
		file.fill(' ');
		file.setf( std::ios::fixed, std:: ios::floatfield );
		file.precision(12);
		file << data[i] << std::endl;
	}
}
// Imprime en pantalla el contenido de un vector
void Vector :: print()
{
	for (int i = 0; i < n; ++i)
	{
		std::cout<<data[i]<<std::endl;
	}
}
// Sobrecarga de operadores
// Suma vectorial
Vector Vector :: operator +(Vector const &v)
{
	// Largo del vector
	int nt = v.n;
	// Verificamos que ambos vectores sean del mismo orden
	if (nt != n)
	{
		std::cout << "ERROR: vectores de distinto largo" << std::endl;
		exit(0);
	}
	// Vector temporal para almacenar la solucion
	Vector tmp(nt);
	// Llevamos acabo la suma
	for (int i = 0; i < nt; ++i)
	{
		tmp[i] = data[i] + v[i];
	}
	// Devolvemos el vector con la solucion 
	return tmp; // INEFICIENTE(HACE UNA COPIA)
}
// Sustraccion vectorial
Vector Vector :: operator -(Vector const &v)
{
	// Largo del vector
	int nt = v.n;
	// Verificamos que ambos vectores sean del mismo largo
	if (nt != n)
	{
		std::cout << "ERROR: vectores de distinto largo" << std::endl;
		exit(0);
	}
	// Vector temporal para almacenar la solucion
	Vector tmp(nt);
	// Llevamos acabo la suma
	for (int i = 0; i < nt; ++i)
	{
		tmp[i] = data[i] - v[i];
	}
	// Devolvemos el vector con la solucion 
	return tmp; // INEFICIENTE(HACE UNA COPIA)
}
// Producto interno
double Vector :: operator *(Vector const &v)
{
	// Largo del vector
	int nt = v.n;
	double suma = 0.0;
	// Verificamos que ambos vectores sean del mismo largo
	if (nt != n)
	{
		std::cout << "ERROR: vectores de distinto largo" << std::endl;
		exit(0);
	}
	// Realizamos el producto interno
	for (int i = 0; i < nt; ++i)
	{
		suma += data[i]*v[i];
	}
	// Devolvemos el vector
	return suma;
}
// Producto de un escalar por un vector
Vector operator *(double escalar, Vector const &v)
{
	// Largo del vector
	int nt = v.n;
	// Vector temporal para almacenar los datos
	Vector a(nt);
	//Realizamos la multplicacion
	for (int i = 0; i < nt; ++i)
	{
		a[i] = escalar * v[i];
	}
	return a;
}
// Sobrecargamos el operador =
// Inicializa un vector 
Vector Vector :: operator =(double escalar)
{
	for (int i = 0; i < n; ++i)
	{
		data[i] = escalar;
	}
	// retornamos el vector
	return *this;
}
//Copia un vector en otro
void Vector :: operator =(Vector const &v)
{
	int nt = v.n;
	// Revisamos que los indices sean los mismos
	if(nt != n) 
	{
		std::cout<<"ERROR: vectores de distinto largo: =" << std::endl;
		exit(0);
	}
	for (int i = 0; i < n; ++i)
	{
		data[i] = v[i];
	}
} 