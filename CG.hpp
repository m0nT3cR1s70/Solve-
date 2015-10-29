
#ifndef CG_HPP
#define CG_HPP

#include <cmath>
#include <iostream>
#include <chrono>
#include <string>

//#include "Matrix.hpp"
#include "Vector.hpp"
#include "CSR.hpp"
#include "Timer.hpp"

class CG {
private:
   int mmaxIts;
   double mtol;
   int mits;
   double merror;
   double etime;
   bool precond;
   std::string namep;
   Timer timer;
public:
   CG()
   {
      mmaxIts = 2000;
      mtol = 1e-6;
      mits = 0;
      merror = 0;
      etime = 0;
   }
   template<class T>
   void solve(T &A, Vector &x, Vector &b)
   {

      timer.tic();
      double tol2 = mtol*mtol;
      int n = x.size();
      Vector r(n);
      Vector p(n);
      Vector Ap(n);
      double alpha;
      double beta;
      double rr0,rr;
      int k;
      
      r = b-A*x;
      rr0 = r*r;
      p=r;
      for(k =0;k<mmaxIts;++k)
      {
         Ap = A*p;
         alpha =  (rr0) / (Ap*p);
         x = x + alpha*p;
         r = r - alpha*Ap;
         rr = r*r;
         merror = rr;
         if(merror < tol2)
            break;
         beta = rr / rr0; 
         p = r + beta*p;
         rr0 = rr;
      }
      mits = k+1;
      merror = sqrt(merror);
      timer.toc();
      etime = timer.etime();
      precond = false;
   }
   
   template<class T, class U>
   void solve(T &A, Vector &x, Vector &b,U &M)  //T: tipo de matriz de entrada, 
                                                //U: tipo de precondicinador de entrada
   {

      timer.tic();                 //comienza a medir tiempo de ejecucion
      double tol2 = mtol*mtol;
      int n = x.size();
      Vector r(n),z(n),p(n),Ap(n);
      double alpha;
      double beta;
      double rr0,rr;
      int k;
      
      r = b-A*x;                    //primer residual operacion vector = vector - matriz*vector
      M.solve(z,r);                 //resuleve el sistema Mz=r
      p=z;                          //copia de vector
      rr0 = z*r;                    //producto interno de dos vectores
      for(k =0;k<mmaxIts;++k)
      { 
         Ap = A*p;                  //producto matriz vector
         alpha =  (rr0) / (Ap*p);   //producto interno en el divisor
         x = x + alpha*p;           //operacion vector = vector + escalar*vector (saxpy en BLAS)
         r = r - alpha*Ap;
         M.solve(z,r);              //resuleve el sistema Mz=r
         rr = z*r;                  //producto interno
         merror = r*r;              //producto interno
         if(merror < tol2)
            break;
         beta = rr / rr0; 
         p = z + beta*p;            //operacion vector = vector + escalar*vector
         rr0 = rr;
      }
      mits = k+1;
      merror = sqrt(merror);
      timer.toc();                 //termina de medir tiempo
      etime = timer.etime();       //guarda en etime el tiempo de ejecucion
      precond = true;
      namep = M.name();            //guarda el nombre del precondicionador
   }
   
   int its(){return mits;}
   double error() {return merror;}
   
   void maxIts(int it){ mmaxIts = it;}
   void tol(double t) {mtol=t;}
   void report()
   {
      std::cout << "Gradiente conjugado";
      if(precond)
          std::cout << " precondicionado con "<< namep ;
      std::cout << "\n";
      std::cout << "  Iteraciones  "<< mits<< std::endl;
      std::cout << "  Error "<< merror << std::endl;
      std::cout << "  Tiempo de calculo "<< etime << "msec" << std::endl;
   }
};

#endif
