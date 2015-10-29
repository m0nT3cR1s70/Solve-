#include "COO.hpp"
#include "CSR.hpp"
#include "Vector.hpp"
#include "Timer.hpp"
#include "Jacobi.hpp"
#include "GaussSeidel.hpp"
#include "SOR.hpp"
#include "CG.hpp"
#include "BICGSTAB.hpp"
#include "JacobiP.hpp"
#include "ILU.hpp"
#include "MILU.hpp"
#include "ICHOL.hpp"
using namespace std;
int main(int argc, char const *argv[])
{
  int nx=1001;
  int ny=nx;
  double dx = 1./nx;
  double dy = 1./ny;
  int n = (nx - 1)*(ny - 1); // dimension del sistema a resolver

  //almacenamiento
  CSR A(n);                      //matriz en formato CSR
  Vector x(n),b(n),r(n);         //vectores: solucion, lado derecho, residual 
                                  //(no es de la biblioteca STL de C++)
    //x.print();

  //solvers: todos tienen por default 2000 iteraciones máximas y tolerancia de convergencia 1e-6
  Jacobi jac;                    //Metodo de Jacobi
  GaussSeidel gs;                //Metodo de Gauss Seidel
  SOR sor;                       //Metodo de sobre-relajacion
  CG cg;                         //Metodo de gradiente conjugado
  BICGSTAB bicg;                 //Metodo de gradiente biconjugado estabilizado

  //precondicionadores, reciben el formato de la matriz como un parametro template
  JacobiP<CSR> jacp;
  ILU<CSR> ilu;                   //Precondicionador LU incompleto
  MILU<CSR> milu;                //Precondicionador LU incompleto modificado
  ICHOL<CSR> ichol;              //Precondicionador Cholesky incompleto

  Timer timer;                   //mide tiempo de ejecucion
   
  //llena matriz COO
  cout<< endl << "Tamanio de problema " << n << "x" <<n<<endl<<endl;
  int l = 0;

  {
    COO Acoo(n);                 //matriz temporal en formato de coordenadas
    timer.tic();                 //comienza a medir tiempo
    for(int j=1;j<ny;++j)
    {
      for(int i=1;i<nx;++i)
      {
        // prototipo de funcion de insertar void insert(int i,int j,doble val);
        if(j>1)
          Acoo.insert(l, l-(nx-1),-1.);      
        if(i>1)
          Acoo.insert(l, l-1,-1.);
          Acoo.insert(l, l,4.);
        if(i<nx-1)
          Acoo.insert(l, l+1,-1.);
        if(j<ny-1)
          Acoo.insert(l, l+(ny-1),-1.);
          ++l;
      }
    }
    timer.toc();                 //termina de medir tiempo
    std::cout << "Tiempo de llenado de matriz   COO: " << timer.etime() << " ms" << std::endl;
    timer.tic();
    A.convert(Acoo);             //convierte matriz COO a formato CSR
    //Acoo.impMatrix();
    timer.toc();
    std::cout << "Tiempo de conversion de COO a CSR " <<timer.etime() << " ms" << std::endl;
  }//destruye Acoo

  b = 1.*dx*dx;                  //a todas las entradas del vector b les asigna el valor 1.*dx*dx
  /*    JACOBI     */
  cout<< endl << endl;
  x=0;                           //aproximacion inicial de la solucion
  //jac.maxIts(1);
  jac.solve(A,x,b);              //resuelve el sistema Ax=b, guarda el resultado en x
  jac.report();                  //reporta numero de iteraciones y tiempo de ejecucion
  x.saveData("jac.txt");
  r = b-A*x;                     //calcula el vector residual con la solucion x. 
                                  //Operacion vector = vector - Matroz*vector
  cout << "#Error ||b-A*x||: " << r.norm() << endl;  //imprime la norma del vector residual

  /*    Gauss Seidel     */
   
  cout<< endl << endl;
  x=0;
  gs.solve(A,x,b);
  gs.report();
  x.saveData("gauss.txt");
  r = b-A*x;
  cout << "#Error ||b-A*x||: " << r.norm() << endl;
  /*    SOR     */
   
  cout<< endl << endl;
  x=0;
  sor.setOmega(1.9);
  sor.solve(A,x,b);
  sor.report();
  r = b-A*x;
  cout << "#Error ||b-A*x||: " << r.norm() << endl;

   
  /**    Gradiente Conjugado   */
  cout<< endl << endl;
  x=0;
  cg.solve(A,x,b);
  cg.report();
  r = b-A*x;
  x.saveData("cg.txt");
  cout << "#Error ||b-A*x||: " << r.norm() << endl;

  /**    Gradiente Conjugado precondicionado con Jacobi incompleto   */
   
  cout<< endl << endl;
  jacp.calculate(A);            //calcula el precondicionador
  x=0;
  cg.solve(A,x,b,jacp);
  cg.report();
  r = b-A*x;
  x.saveData("cgpjacobi.txt");
  cout << "#Error ||b-A*x||: " << r.norm() << endl;

  /**    Gradiente Conjugado precondicionado con ILU incompleto   */

   cout<< endl << endl;
  ilu.calculate(A);            //calcula el precondicionador
  x=0;
  cg.solve(A,x,b,ilu);
  cg.report();
  r = b-A*x;
  x.saveData("cgpilu.txt");
  cout << "#Error ||b-A*x||: " << r.norm() << endl;

  /**    Gradiente Conjugado precondicionado con MILU incompleto   */

  cout<< endl << endl;
  milu.calculate(A);            //calcula el precondicionador
  x=0;
  cg.solve(A,x,b,milu);
  cg.report();
  r = b-A*x;
  x.saveData("cgpmilu.txt");
  cout << "#Error ||b-A*x||: " << r.norm() << endl;

  /**    Gradiente Conjugado precondicionado con ichol incompleto   */

  cout<< endl << endl;
  ichol.calculate(A);            //calcula el precondicionador
  x=0;
  cg.solve(A,x,b,ichol);
  cg.report();
  r = b-A*x;
  x.saveData("cgpichol.txt");
  cout << "#Error ||b-A*x||: " << r.norm() << endl;
   
  /**    Gradiente biconjugado estabilizado (BICGSTAB)      */
   
  cout<< endl << endl;
  x = 0;
  bicg.solve(A,x,b);
  bicg.report();
  r = b-A*x;
  x.saveData("bicg.txt");
  cout << "#Error ||b-A*x||: " << r.norm() << endl;
   
  /**    BICGSTAB  Precondicionado con LU incompleto   */
  
  cout<< endl << endl;
  ilu.calculate(A);              //calculo del precondicionador
  x=0;
  bicg.solve(A,x,b,ilu);
  bicg.report();
  r = b-A*x;
  x.saveData("bicgilu.txt");
  cout << "#Error ||b-A*x||: " << r.norm() << endl;
  
   
  /**    BICGSTAB  Precondicionado con LU incompleto modificado */
  
  cout<< endl << endl;
  milu.calculate(A);            //calculo del precondicionador
  x=0;
  bicg.solve(A,x,b,milu);
  bicg.report();
  r = b-A*x;
  x.saveData("bicgmilu.txt");
  cout << "#Error ||b-A*x||: " << r.norm() << endl;
  
   
  /**    BICGSTAB  Precondicionado con ichol */
  
  cout<< endl << endl;
  ichol.calculate(A);
  //ichol.imchol();
  cout<< endl ;
  x=0;
  bicg.solve(A,x,b,ichol);
  bicg.report();
  x.saveData("bicgIchol.txt");
  r = b-A*x;
  cout << "#Error ||b-A*x||: " << r.norm() << endl;
  
  return 0;
}
