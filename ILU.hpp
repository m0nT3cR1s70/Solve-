#ifndef ILU_HPP
#define ILU_HPP

//-----------------------------
// Librerias Necesarias
#include <iostream>
#include "CSR.hpp"
//-----------------------------

template<class P>
class ILU
{
    private:
        // Objeto que contiene la matriz LU
        P LU;
        // Almacena del nombre del precondicionador
        std::string namep;
    public:
        // Constructor vacio
        ILU():namep("ILU"){};
        // Devuelve el nombre del precondicionador
        inline std::string name(){return namep;};
        // Calcula el precondicionador
        void calculate(P const &Mtx);
        // Resuelve un sistema de ecuaciones usando sustitucion hacia adelante y hacia atras
        void solve(Vector &x, Vector const &b);
        // Imprimir resultados
        void iSolve();
};
template<class P>
void ILU<P> :: calculate(P const &Mtx)
{
    // Copiamos el contenido de la matriz
    LU = Mtx;
    // Variables utiles
    int k = 0;
    double L = 0.0;
    int lli = 0;
    int llk = 0;
    int ji = 0;
    int jk = 0;
    // Comenzamos con la factorizacion
    for (int i = 0; i < LU.n; ++i)
    {
        for (int l = LU.irow[i]; l < LU.idiag[i]; ++l)
        {
            k = LU.col[l];
            L = LU.data[l]/LU.data[LU.idiag[k]];
            lli = LU.irow[i];
            llk = LU.irow[k];
            ji = LU.col[lli];
            jk = LU.col[llk];
            while(lli < LU.irow[i+1] && llk < LU.irow[k+1]) 
            {
                if (ji == jk)
                {
                    if (ji > k && jk > k)
                    {
                        LU.data[lli] = LU.data[lli] - L*LU.data[llk];
                    }
                    lli = lli + 1;
                    llk = llk + 1;
                    ji = LU.col[lli];
                    jk = LU.col[llk];
                }
                else if (ji < jk)
                {
                    lli = lli + 1;
                    ji = LU.col[lli];
                }
                else
                {
                    llk = llk + 1;
                    jk = LU.col[llk];
                }
            }
            LU.data[l]=L;
        }
    }
}
// Solve de ILU
template<class P>
void ILU<P> :: solve(Vector &x, Vector const &b)
{

    Vector aux(x.n);
    // Sustitucion hacia adelante
    //x = b;
    for (int i = 0; i < x.n; ++i)
    {
        aux[i] = b[i];
        for (int l = LU.irow[i]; l < LU.idiag[i]; ++l)
        {
            aux[i] = aux[i] - LU.data[l]*aux[LU.col[l]];
        }
    }
    // Sustitucion hacia atras
    //x[x.n - 1] = x[LU.n - 1] / LU.data[LU.idiag[LU.n]-1];
    for (int i = LU.n-1; i >= 0; --i)
    {
        x[i] = aux[i];
        for (int l = LU.idiag[i]+1; l <= LU.irow[i+1]-1; ++l)
        {
            x[i] = x[i] - LU.data[l]*x[LU.col[l]];
        }
        x[i] = x[i]/LU.data[LU.idiag[i]];
    }

}
template<class P>
void ILU<P> :: iSolve()
{
    for (int i = 0; i < LU.nnz; ++i)
    {
        std::cout<<LU.data[i]<< "\t\t\t\t" << LU.col[i] << "\t\t"<<LU.irow[i]<<std::endl;
    }
    std::cout<<std::endl;
}

#endif // ILU_HPP
