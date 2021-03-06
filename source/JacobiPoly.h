#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <cmath>
#include <string.h>
#include "lapacke.h"

/*get generalized jacobi polynomials at points in array xi */
double* get_JACOBIPOLY(double* xi, int Q, int N, double alpha, double beta)
{
    double *LN = new double[Q];
    int n = N-1;
    double a1_n = 2.0*(n+1)*(1.0*n+alpha+beta+1.0)*(2.0*n+alpha + beta);
    double a2_n = (2.0*n + alpha + beta + 1.0)*(alpha*alpha - beta*beta);
    double a3_n = (2.0*n + alpha + beta)*(2.0*n + alpha + beta + 1.0)*(2.0*n + alpha + beta + 2.0);
    double a4_n = 2.0*(1.0*n + alpha)*(1.0*n + beta)*(2.0*n + alpha + beta + 2.0);
    if (N==0)
    {
        for (int i=0;i<Q;i++)
            LN[i] = 1.0;
    }
    else if (N==1)
        for (int i=0;i<Q;i++)
            LN[i] = 0.5*(alpha - beta + (alpha + beta + 2.0)*xi[i]);
    
    else
    {   
        double* LN_1 = get_JACOBIPOLY(xi, Q, N-1, alpha, beta);
        double* LN_2 = get_JACOBIPOLY(xi, Q, N-2, alpha, beta);
        for (int i=0; i<Q; i++)
            LN[i] = (a2_n +a3_n*xi[i])*LN_1[i]/a1_n - a4_n*LN_2[i]/a1_n;

    }
    return LN;

}

double* get_JACOBI_POLY_END(double alpha, double beta, int n)
{   
    double *P = 0;
    P = new double[2];
    int a1 = std::pow(2,n); 
    double num_1 = tgamma(alpha+1.0+1.0*n)/tgamma(alpha+1.0); 
    double den_1 = tgamma(alpha+beta+1.0*n+1.0)/tgamma(alpha+beta+2.0*n+1.0);
    P[1] = a1*num_1/den_1;

    double num_2 = std::pow(-1, n)*tgamma(beta+1.0+1.0*n)/tgamma(beta+1.0);
    P[0] = a1*num_2/den_1;
    return P;

}

double* get_JACOBI(double alpha, double beta, int Q)
{   
    double *J = 0;
    J = new double[Q*Q];
    for(int i=0 ; i<Q ; i++)
        for(int j=0 ; j<Q ; j++)
            J[i*Q+j] = 0.0;
    
    for(int i=0 ; i<Q-1 ; i++)
    {   
        /*Only for Legendre polynomial*/

        /*For jacobi polynomial of general order*/
        /*Lapack functions are column major...because FORTRAN*/
        int k = i+1; double P = 2.0*k + alpha + beta;
        J[i*Q + k] = sqrt(4.0*k*(1.0*k+alpha)*(1.0*k+beta)*(1.0*k+alpha+beta)/(P*P*(P+1.0)*(P-1.0)));
        J[i*Q + k-1] = (beta*beta - alpha*alpha)/(P*(P+2.0));

    }   
    return J;
}
double* get_JACOBI_RADAU(double alpha, double beta, int Q)
{
    double *JR = get_JACOBI(alpha, beta, Q);
    int k = Q-1; 
    /*double P = 2.0*k + alpha + beta;
    double beta_k = 4.0*k*(1.0*k+alpha)*(1.0*k+beta)*(1.0*k+alpha+beta)/(P*P*(P+1.0)*(P-1.0));*/
    
    double* P_q_2_ab = get_JACOBI_POLY_END(alpha, beta, Q-2);
    double* P_q_1_ab = get_JACOBI_POLY_END(alpha, beta, Q-1);

    //JR[Q*Q-1] = -1 - beta_k*P_q_2_ab[0]/P_q_1_ab[0];
    JR[Q*Q-1] = -1.0 +2.0*k*(1.0*k+alpha)/((2.0*k+alpha+beta)*(2.0*k+alpha+beta+1.0));
    return JR;
}

double* get_JACOBI_LOBATTO(double alpha, double beta, int Q)
{
    double *JL = get_JACOBI(alpha, beta, Q);
    int k = Q-1;
    double pk1_pk1 = 0.5*((alpha+beta+2.0*k)*(alpha+beta+2.0*k-1.0))/((1.0*k+alpha)*(1.0*k+alpha+beta));
    double pk1_pk0 = -0.5*((alpha+beta+2.0*k)*(alpha+beta+2.0*k-1.0))/((1.0*k+beta)*(1.0*k+alpha+beta));
    double beta_L = 2.0/(pk1_pk1 - pk1_pk0);
    double alpha_L = -1.0*(pk1_pk1 + pk1_pk0)/(pk1_pk1 - pk1_pk0);

    JL[Q*Q-Q-1] = sqrt(beta_L);
    JL[Q*Q-1] = alpha_L;

    return JL;
}



double* get_WEIGHTS(double* J, int Q, double alpha, double beta)
{
    double B_0 = std::pow(2, alpha+beta+1.0)*tgamma(alpha+1.0)*tgamma(beta+1.0)/tgamma(alpha + beta + 1.0);
    double* W = new double[Q]; 
    for (int i=0; i<Q; i++)
    {   
        double mag = 0.0;
        for (int j=0; j<Q; j++)
            mag += J[i*Q+j]*J[i*Q+j];
        
        W[i] = B_0*J[i*Q]*J[i*Q]/mag;
    }
    return W;
}
double* get_EIGS(double* J, int Q)
{
    char    TRANS = 'V'; char    DECOM = 'L';                                   //EGVEC flag, Upper-Lower flag
    int     INFO=3; int     LDA = Q; int     NRHS = 1;
    double* xi=0; xi = new double[Q];
    int LWORK = 3*Q;
    double     work[LWORK*sizeof(int)];
    //double xi[Q];
    dsyev_(&TRANS, &DECOM, &Q, J, &LDA, xi, work, &LWORK, &INFO);
    return xi;
}

/*get collocated differentiation matrix for Gauss-Lobatto*/
double** get_DIFFMATLEG(double* xi, int Q)
{   
    double* L_Q_1 = get_JACOBIPOLY(xi, Q, Q-1, 0.0 , 0.0);
    double** D = 0;
    D = new double*[Q];
    for (int i =0; i<Q; i++)
    {
        D[i] = new double[Q];
        for (int j = 0; j<Q; j++)
        {
            if(i==j)
            {   
                D[i][j] = 0.0;
                if (i==Q-1)
                    D[i][j] = 1.0*Q*(1.0*Q-1.0)/4.0;
                else if (i == 0)
                    D[i][j] = -1.0*Q*(1.0*Q-1.0)/4.0;
            }
            else
            {
                D[i][j] = L_Q_1[i]/(L_Q_1[j]*(xi[i]-xi[j]));
            }
        }
    }
    return D;
}
double* get_C0MODALBASIS_a(double* xi, int p, int P, int Q)
{
    double* phi_p = new double[Q];
    if (p==0)
    {
        for (int i = 0; i<Q; i++)
            phi_p[i] = 0.5*(1.0 - xi[i]);
    }
    else if (p==P)
    {
        for (int i = 0; i<Q; i++)
            phi_p[i] = 0.5*(1.0 + xi[i]);
    }
    else
    {
        phi_p = get_JACOBIPOLY(xi, Q, p-1, 1.0, 1.0);
        for (int i = 0; i<Q; i++)
            phi_p[i] = 0.25*(1.0 + xi[i])*(1.0 - xi[i])*phi_p[i];
    }
    return phi_p;
}   
double* get_ORTMODALBASIS_a(double* xi, int p, int P, int Q)
{
    double* phi_p = new double[Q];
    phi_p = get_JACOBIPOLY(xi, Q, p, 0.0, 0.0);
    
    return phi_p;
}
double* get_C0MODALBASIS_b(double* xi, int p, int q, int P, int Q)
{
    double* phi_pq = new double[Q];

    if ((p==0)||(p==P))
        phi_pq = get_C0MODALBASIS_a(xi, q, P, Q);

    else
    {
        if (q==0)
        {
            for (int i = 0; i<Q; i++)
                phi_pq[i] = std::pow((0.5*(1.0 - xi[i])), p+1);
        }
        else
        {
            phi_pq = get_JACOBIPOLY(xi, Q, q-1, 2.0*p+1.0,1.0);
            for (int i = 0; i<Q; i++)
                phi_pq[i] = 0.5*std::pow((0.5*(1.0 - xi[i])), p+1)*(1.0 + xi[i])*phi_pq[i];
        }
    } 

    return phi_pq;
}   
double* get_ORTMODALBASIS_b(double* xi, int p, int q, int P, int Q)
{
    double* phi_pq = new double[Q];
    phi_pq = get_JACOBIPOLY(xi, Q, q, 2.0*p+1, 0.0);
    for (int i = 0; i<Q; i++)
        phi_pq[i] = std::pow((0.5*(1.0 - xi[i])), p)*phi_pq[i];
    return phi_pq;
}
