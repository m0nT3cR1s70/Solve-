#ifndef BICGSTAB_HPP
#define BICGSTAB_HPP

#include <cmath>
#include <iostream>
#include <chrono>
#include <string>

//#include "Matrix.hpp"
#include "Vector.hpp"
#include "CSR.hpp"
#include "Timer.hpp"

class BICGSTAB {
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
   BICGSTAB()
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
      Vector ra(n);
      Vector p(n);
      Vector v(n);
      Vector Ap(n);
      Vector s(n);
      Vector t(n);
      Vector ro(n);
      double alpha;
      double beta;
      double omega;
      double ro0;
      double ss;
      int k;
      
      // Inicializamos lo pedido
      r = b - A*x;   //Calculamos el residual
      ra = r;      // r arbitraria
      p = 0;        //p0
      v = 0;      //v0
      ro[0] = 1;     //ro
      alpha = 1;  //alpha
      omega = 1;  //omega
      for(k =1; k<mmaxIts; ++k)
      {
         ro[k] = ra*r;
         beta = (ro[k]/ro[k-1]) * (alpha/omega);
         p = r + beta*(p-omega*v);
         v = A*p;
         alpha = ro[k]/(ra*v);
         s = r - alpha*v;
         t = A*s;
         omega = (t*s)/(t*t);
         x = x + alpha*p + omega*s;
         ss = s*s;
         merror = ss;
         if(merror < tol2)
            break;
         r = s - omega*t;
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
      timer.tic();
      double tol2 = mtol*mtol;
      int n = x.size();
      Vector r(n);
      Vector ra(n);
      Vector p(n);
      Vector v(n);
      Vector Ap(n);
      Vector s(n);
      Vector t(n);
      Vector ro(n);
      Vector o(n);
      Vector y(n);
      double alpha;
      double beta;
      double omega;
      double ro0;
      double ss;
      int k;
      
      // Inicializamos lo pedido
      r = b - A*x;   //Calculamos el residual
      ra = r;      // r arbitraria
      p = 0;        //p0
      v = 0;      //v0
      ro[0] = 1;     //ro
      alpha = 1;  //alpha
      omega = 1;  //omega

      for(k =1;k<mmaxIts;++k)
      { 
         // Se calcula ro de k
         ro[k] = ra*r;
         // Se calcula la beta                  
         beta = (ro[k]/ro[k-1]) * (alpha/omega);   
         p = r + beta*(p-omega*v);
         M.solve(y,p);              //resuleve el sistema My=p
         v = A*y;
         alpha = ro[k]/(ra*v);
         s = r - alpha*v;
         M.solve(o,s);               //resuleve el sistema Mo=p
         t = A*o;
         omega = (t*s)/(t*t);
         x = x + alpha*y + omega*o;
         ss = s*s;
         merror = ss;
         if(merror < tol2)
            break;
         r = s - omega*t;
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
      std::cout << "BICGSTAB";
      if(precond)
          std::cout << " precondicionado con "<< namep ;
      std::cout << "\n";
      std::cout << "  Iteracinoes  "<< mits<< std::endl;
      std::cout << "  Error "<< merror << std::endl;
      std::cout << "  Tempo de calculo "<< etime << " msec" << std::endl;
   }
};

#endif
