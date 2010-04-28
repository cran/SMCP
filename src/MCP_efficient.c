#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#define Max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define Min( a, b ) ( ((a) < (b)) ? (a) : (b) )

double Sgn(double A)
{
  double value;
  if ( fabs(A) < 1E-12) 
     value = 0;
  if ( A > 0)
     value = 1;
  else if ( A < 0)
     value = -1;
  return value;
}

double Norm2(double *Beta_tilde, double *Old_beta, int P)
{
  int j;
  double value=0;
  
  for ( j = 0; j < P; j++ )
  {
    value = value + pow((*(Beta_tilde+j))-(*(Old_beta+j)),2);
  } 
  return value;
}

double S_func(double z , double gamma)  
{
  double value; 
  if (z > 0 && gamma < fabs(z)) 
      value = z - gamma;
  else if (z < 0 && gamma < fabs(z)) 
      value = z + gamma;
  else 
      value = 0;
  return value ;
}

double Solve_eq(double A, double B, double C)
{
  double solution=0;
  if ( C > 0 )
  {
    solution = S_func(-B/(2*A),C/(2*A));
  }
  else if ( C <= 0 )
  {
    solution = -B/(2*A)+Sgn(B)*C/(2*A); 
  }
  return solution;
}


double Coordinate_MCP(double *Y,double *X, int *N, double *Weight, double *W1, double *W2, 
                     double *Lambda, double *Beta_tilde, double *Gamma, int *Param, 
                     double *Epsilon, int *M)
// D is the value for upper and lower points;
// S is the active set;
{
  int i, j ,n, start, p= Param[0], iter=Param[1];
  double lambda1 = Lambda[0], lambda2 =Lambda[1], gamma = Gamma[0];  
//  int fusion_set_1st[Param[1]];  
  double /*beta_tilde[Param[1]],*/old_beta[Param[0]], diff; 
  double ak, bk, a, b, c, tmp, part1[Param[0]], part2[Param[0]];
  for ( j = 0; j < p; j++)
  {
    old_beta[j] = 0;
  }   
  for ( j = 0; j < p; j++)
  {  
    n = *(N + j); 
    start = 0;
    for ( i = 0; i < j; i++ )
    {
      start = start + *(N + i);
    }  
    part1[j] = 0;
    part2[j] = 0;
    for ( i = 0; i < n; i++ )
    {
      part1[j] = part1[j] + (*(Weight+start+i))*(*(X+start+i))*(*(X+start+i))/2/n;
      part2[j] = part2[j] - (*(Weight+start+i))*(*(Y+start+i))*(*(X+start+i))/n;
    }
  }  
  do 
  {
    diff = 0;
    a = 0;
    b = 0;
    for ( j = 0; j < p; j++)
    {  
      n = *(N + j); 
      start = 0;
      for ( i = 0; i < j; i++ )
      {
        start = start + *(N + i);
      }
      if ( j == 0 )
      {
          a = 0;
          b = 0;
          ak = fabs(Beta_tilde[j+1]);
          bk = 0;
          a = part1[j] + lambda2*W2[0]/2-W1[0]/gamma/2;
          b = part2[j];
          c = lambda1*W1[j] - lambda2*(W2[0]*ak);
          tmp = Solve_eq(a, b, c);
          if ( abs(tmp) < gamma*lambda1)
            Beta_tilde[j] = tmp;
          else 
          {
            a = a + W1[0]/gamma/2;
            c = - lambda2*(W2[0]*ak);
            tmp = Solve_eq(a, b, c);
            Beta_tilde[j] = tmp;
          }  
      }
      else if ( j == (p-1) )
      {
          a = 0;
          b = 0;
          ak = 0;
          bk = fabs(Beta_tilde[j-1]);
          a = part1[j] + lambda2*W2[p-2]/2-W1[p-1]/gamma/2;
          b = part2[j];
          c = lambda1*W1[p-1] - lambda2*(W2[p-2]*bk);
          tmp = Solve_eq(a, b, c); 
          if ( abs(tmp) < gamma*lambda1)
            Beta_tilde[j] = tmp;
          else 
          {
            a = a + W1[p-1]/gamma/2;
            c = - lambda2*(W2[p-2]*bk);
            tmp = Solve_eq(a, b, c);
            Beta_tilde[j] = tmp;
          }
      }
      else 
      {
          a = 0;
          b = 0;
          ak = fabs(Beta_tilde[j+1]);
          bk = fabs(Beta_tilde[j-1]);
          a = part1[j] + lambda2*(W2[j-1]+W2[j])/2-W1[j]/gamma/2;
          b = part2[j];
          c = lambda1*W1[j] - lambda2*(W2[j]*ak+W2[j-1]*bk);
          tmp = Solve_eq(a, b, c);           
          if ( abs(tmp) < gamma*lambda1)
            Beta_tilde[j] = tmp;
          else 
          {
            a = a + W1[j]/gamma/2;
            c = - lambda2*(W2[j]*ak+W2[j-1]*bk);
            tmp = Solve_eq(a, b, c);
            Beta_tilde[j] = tmp;    
          } 
      }  
     // if ( fabs(Beta_tilde[j]) < 1E-10 )
       //   Beta_tilde[j] = 0;
    }      
    M[0] = M[0] + 1;  
    diff = Norm2(Beta_tilde, old_beta, p);
    for ( j = 0; j < p; j++ )
    {
      old_beta[j] = Beta_tilde[j];
    }
  }
  while ( M[0] <= iter && diff > Epsilon[0]);
//  printf("m=%d\n",M[0]); 
  return 0;
}

/* This defines a data structure for registering the routine
   `pareto'. It records the number of arguments, which allows some
   error checking when the routine is called. */
static R_CMethodDef DotCEntries[] = {
  {"SMCP", (DL_FUNC) Coordinate_MCP, 12},
  {NULL}
};

/* This is called by the dynamic loader to register the routine. */
void R_init_pareto(DllInfo *info)
{
  R_registerRoutines(info, DotCEntries, NULL, NULL, NULL);
}
