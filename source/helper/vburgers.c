
//
// viscous Burgers' equation
//
//   u_t + u u_x = u_xx
//
// To compile:
//    gcc -Wall -std=c99 -o vburgers vburgers.c
//

#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

//
// Global Variables
//

int imax, imethod, n, nmax, nprint;
double dt, dx, t;
double *unew, *uold, *x;

const double xmin = -9.0, xmax = 9.0;
const double tmin = 0.1, tmax = 1.0;

bool first_time;
FILE *fp1, *fp2, *fp3;



double uanalytical(double x, double t) {
  return -2.*sinh(x) / (cosh(x) - exp(-t));
}



void tridiag(double *aa, double *dd, double *cc, double *bb, double *x, int imax) {
  //
  // Thomas's Tridiagonal Algorithm
  //
  // Description:
  //
  //   Solve a tridiagonal system of the form:
  //
  //   d[0] x[0] + c[0] x[1] = b[0]
  //   a[1] x[0] + d[1] x[1] + c[1] x[2] = b[1]
  //   a[2] x[1] + d[2] x[2] + c[2] x[3] = b[2]
  //   ...
  //   a[N-2] x[N-3] + d[N-2] x[N-2] + c[N-2] x[N-1] = b[N-2]
  //   a[N-1] x[N-2] + d[N-1] x[N-1] = b[N-1]
  //
  // Reference:
  //    Cheney and Kincaid (1994), Numerical Mathematics and Computing,
  //    3rd ed., Brooks/Cole Publishing Co., Pacific Grove CA, 1994,
  //    sec. 6.3, pp. 249-253
  //
  int i;
  double xmult;

  for (i = 1; i < imax; i++) {
    xmult = aa[i]/dd[i-1];
    dd[i] = dd[i] - xmult*cc[i-1];
    bb[i] = bb[i] - xmult*bb[i-1];
  }

  x[imax-1] = bb[imax-1]/dd[imax-1];

  for (i = imax-2; i >= 0; i--) {
    x[i] = (bb[i] - cc[i]*x[i+1])/dd[i];
  }
}



void setup(void) {
  first_time = true;

  imax = 91;
  nmax = 91;
  imethod = 1;

  printf("Enter imethod:\n");
  scanf("%d", &imethod);

  nprint = (nmax-1)/3;

  unew = malloc(imax*sizeof(double));
  uold = malloc(imax*sizeof(double));
  x = malloc(imax*sizeof(double));

  dx = (xmax - xmin)/(double)(imax-1);
  dt = (tmax - tmin)/(double)(nmax-1);

  for (int i = 0; i < imax; i++) {
    x[i] = xmin + (double)i*dx;
    unew[i] = uanalytical(x[i],tmin);
  }

  double umax = 0.;

  for (int i = 0; i < imax; i++) {
    double val = fabs(unew[i]);
    umax = (val > umax) ? val : umax;
  }

  double c = umax*dt/dx;
  double d = dt/(dx*dx);
  double Rec = umax*dx/1.0;

  switch (imethod) {
  case 1:
    fprintf(stderr,"FTCS\n");
    if (d > 0.5) fprintf(stderr, " *** d > 1/2\n");
    if (c > 1.0) fprintf(stderr, " *** c > 1\n");
    if (Rec > 2./c) fprintf(stderr, " *** Rec > 2/c\n");
    if (c*c > 2.*d) fprintf(stderr, " *** c^2 > 2d\n");
    break;
  case 2:
    fprintf(stderr,"MacCormack\n");
    if (c+2.*d > 1.) fprintf(stderr, " *** c + 2d > 1\n");
    break;
  case 3:
    fprintf(stderr,"BTCS\n");
    break;
  case 4:
   fprintf(stderr,"Roe\n");
    break;
  default:
    fprintf(stderr,"bad imethod\n");
    exit(EXIT_FAILURE);
    break;
  }

  fprintf(stderr,"dx, dt: %g %g c, d, Rec: %g %g %g\n", dx, dt, c, d, Rec);
}



void output(void) {
  if (first_time) {
    fp1 = fopen("solution.dat", "w");
    fp2 = fopen("error_rms.dat", "w");
    fp3 = fopen("error_profile.dat", "w");

    fprintf(fp1,"Title = \"viscous Burgers\' equation\"\n");
    fprintf(fp1,"Variables = \"x\",\"u(x,t)\",\"u<sub>a</sub>(x,t)\"\n");

    fprintf(fp2,"Title = \"viscous Burgers\' equation\"\n");
    fprintf(fp2,"Variables = \"t\",\"RMS Error\"\n");

    fprintf(fp3,"Title = \"viscous Burgers\' equation\"\n");
    fprintf(fp3,"Variables = \"x\",\"u<sup>n</sup><sub>i</sub> - u(x,t)\"\n");

    first_time = false;
  }

  if (n % nprint == 0) {
    fprintf(fp1, "Zone T = \"t = %8.3f\"\n", t);
    fprintf(fp3, "Zone T = \"t = %8.3f\"\n", t);

    for (int i = 0; i < imax; i++) {
      fprintf(fp1, "%15.8E %15.8E %15.8E\n", x[i], unew[i], uanalytical(x[i],t));
    }

    for (int i = 0; i < imax; i++) {
      fprintf(fp3, "%15.8E %15.8E\n", x[i], unew[i]-uanalytical(x[i],t));
    }

  }

  double err = 0.;

  for (int i = 0; i < imax; i++) {
    double tmp = unew[i] - uanalytical(x[i],t);
    err += tmp*tmp;
  }

  err = sqrt(err);

  fprintf(fp2, "%15.8E %15.8E\n", t, err);
}



void FTCS(void) {
  unew[0] = uold[0];
  unew[imax-1] = uold[imax-1];

  for (int i = 1; i < imax-1; i++) {
    double EP = 0.5*uold[i+1]*uold[i+1];
    double EM = 0.5*uold[i-1]*uold[i-1];
    unew[i] = uold[i] - 0.5*(dt/dx)*(EP - EM) + (dt/(dx*dx))*(uold[i+1] - 2.*uold[i] + uold[i-1]);
  }
}



void MacCormack(void) {
  double *ustar = malloc(imax*sizeof(double));

  ustar[0] = uold[0];
  ustar[imax-1] = uold[imax-1];

  for (int i = 1; i < imax-1; i++) {
    double EP = 0.5*uold[i+1]*uold[i+1];
    double EC = 0.5*uold[i]*uold[i];

    double du = - (dt/dx)*(EP - EC) + (dt/(dx*dx))*(uold[i+1] - 2.*uold[i] + uold[i-1]);
    ustar[i] = uold[i] + du;
  }

  unew[0] = uold[0];
  unew[imax-1] = uold[imax-1];

  for (int i = 1; i < imax-1; i++) {
    double EC = 0.5*ustar[i]*ustar[i];
    double EM = 0.5*ustar[i-1]*ustar[i-1];

    double du = - (dt/dx)*(EC - EM) + (dt/(dx*dx))*(ustar[i+1] - 2.*ustar[i] + ustar[i-1]);
    unew[i] = 0.5*(uold[i] + ustar[i] + du);
  }

  free(ustar);
}



void BTCS(void) {
  double *aa = malloc(imax*sizeof(double));
  double *bb = malloc(imax*sizeof(double));
  double *cc = malloc(imax*sizeof(double));
  double *dd = malloc(imax*sizeof(double));

  // unew[0] = uold[0]
  dd[0] = 1.;
  cc[0] = 0.;
  bb[0] = uold[0];

  // aa[i]*unew[i-1] + dd[i]*unew[i] + cc[i]*unew[i+1] = uold[i]


  for (int i = 1; i < imax-1; i++) {
    aa[i] = -dt/(dx*dx) - 0.5*uold[i]*dt/dx;
    dd[i] = 1. + 2.*dt/(dx*dx);
    cc[i] = -dt/(dx*dx) + 0.5*uold[i]*dt/dx;
    bb[i] = uold[i];
  }

  // unew[imax-1] = uold[imax-1]
  aa[imax-1] = 0.;
  dd[imax-1] = 1.;
  bb[imax-1] = uold[imax-1];

  tridiag(aa, dd, cc, bb, unew, imax);

  free(dd);
  free(cc);
  free(bb);
  free(aa);
}



void Roe(void) {
  fprintf(stderr, "*** PUT YOUR CODE FOR THE ROE SCHEME HERE ***\n");
  exit(EXIT_FAILURE);
}



int main(void) {
  n = 0;
  t = tmin;
  setup();
  output();

  for (n = 1; n < nmax; n++) {
    t = tmin + (double)n*dt;

    for (int i = 0; i < imax; i++) uold[i] = unew[i];

    switch (imethod) {
    case 1:
      FTCS();
      break;
    case (2):
      MacCormack();
      break;
    case 3:
      BTCS();
      break;
    case 4:
      Roe();
      break;
    default:
      exit(EXIT_FAILURE);
    }

    output();
  }

  return EXIT_SUCCESS;
}



