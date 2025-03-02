#include <R.h>
#include <R_ext/Print.h>

void cnlprt_C(char *msg, int *plen) {
  char buf[1000];
  int len = *plen;

  memmove(buf, msg, len);
  buf[len] = '\0';
  Rprintf("\n%s\n", buf);
}

/* 30   FORMAT(/10H   IT   NF,6X,1HF,7X,5HRELDF,3X,6HPRELDF,3X,5HRELDX,
   1       2X,13HMODEL  STPPAR) */
void F77_SUB(h30)(void) {
  Rprintf("\n    IT   NF      F       RELDF   PRELDF   RELDX  MODEL  STPPAR\n");
}
/* 40   FORMAT(/11H    IT   NF,7X,1HF,8X,5HRELDF,4X,6HPRELDF,4X,5HRELDX,
   1       3X,6HSTPPAR) */
void F77_SUB(h40)(void) {
  Rprintf("\n    IT   NF      F         RELDF    PRELDF    RELDX   STPPAR");
}

/* 70   FORMAT(/11H    IT   NF,6X,1HF,7X,5HRELDF,3X,6HPRELDF,3X,5HRELDX,
   1       2X,13HMODEL  STPPAR,2X,6HD*STEP,2X,7HNPRELDF) */
void F77_SUB(h70)(void) {
  Rprintf("\n    IT   NF      F       RELDF   PRELDF   RELDX  MODEL  STPPAR");
  Rprintf("   D*STEP   NPRELDF\n");
}

/* 80   FORMAT(/11H    IT   NF,7X,1HF,8X,5HRELDF,4X,6HPRELDF,4X,5HRELDX,
     1       3X,6HSTPPAR,3X,6HD*STEP,3X,7HNPRELDF) */
void F77_SUB(h80)(void) {
  Rprintf("\n    IT   NF      F         RELDF    PRELDF    RELDX   STPPAR");
  Rprintf("   D*STEP   NPRELDF\n");
}

/* 100  FORMAT(I6,I5,D10.3,2D9.2,D8.1,A3,A4,2D8.1,D9.2) */
void h100s_C(int *i1, int *i2, double *d1, double *d2, double *d3, double *d4,
             char *a1, char *a2, double *d5) {
  Rprintf("%6d%5d%10.3e%9.2e%9.2e%8.1e%3s%4s%8.1e\n", *i1, *i2, *d1, *d2, *d3,
          *d4, a1, a2, *d5);
}

void h100l_C(int *i1, int *i2, double *d1, double *d2, double *d3, double *d4,
             char *a1, char *a2, double *d5, double *d6, double *d7) {
  Rprintf("%6d%5d%10.3e%9.2e%9.2e%8.1e%3s%4s%8.1e%8.1e%e9.2\n", *i1, *i2, *d1,
          *d2, *d3, *d4, a1, a2, *d5, *d6, *d7);
}

/*  110  FORMAT(I6,I5,D11.3,2D10.2,3D9.1,D10.2) */
void F77_SUB(h110s)(int *i1, int *i2, double *d1, double *d2, double *d3,
                    double *d4, double *d5) {
  Rprintf("%6d%5d%11.3e%10.2e%10.2e%9.1e%9.1e\n", *i1, *i2, *d1, *d2, *d3, *d4,
          *d5);
}
void F77_SUB(h110l)(int *i1, int *i2, double *d1, double *d2, double *d3,
                    double *d4, double *d5, double *d6, double *d7) {
  Rprintf("%6d%5d%11.3e%10.2e%10.2e%9.1e%9.1e%9.1e%10.2e\n", *i1, *i2, *d1, *d2,
          *d3, *d4, *d5, *d6, *d7);
}

void F77_SUB(h380)(int *i) { Rprintf(" ***** IV(1) =%i5 *****\n", *i); }

void F77_SUB(h400)(int *p, double *x, double *d) {
  int i;

  Rprintf("\n     I     INITIAL X(I)        D(I)\n\n");
  for (i = 0; i < *p; i++) Rprintf(" %5i%17.6e%14.3e\n", i + 1, x[i], d[i]);
}

void F77_SUB(h410)(double *x) { Rprintf("     0    1%10.3e\n", *x); }

void F77_SUB(h420)(double *x) { Rprintf("     0    1%11.3e\n", *x); }

void F77_SUB(h450)(double *d1, double *d2, int *i1, int *i2, double *d3,
                   double *d4) {
  Rprintf("\n FUNCTION%17.6e   RELDX%17.3e\n", *d1, *d2);
  Rprintf(" FUNC. EVALS%8i         GRAD. EVALS%8u\n", *i1, *i2);
  Rprintf(" PRELDF%16.3e      NPRELDF%15.3e\n", *d3, *d4);
}

void F77_SUB(h460)(int *i) {
  Rprintf("\n %4d EXTRA FUNC. EVALS FOR COVARIANCE AND DIAGNOSTICS\n", *i);
}

void F77_SUB(h470)(int *i) {
  Rprintf("\n %4d EXTRA GRAD. EVALS FOR COVARIANCE AND DIAGNOSTICS\n", *i);
}

void F77_SUB(h500)(int *p, double *x, double *d, double *g) {
  int i;

  Rprintf("\n");
  for (i = 0; i < *p; i++)
    Rprintf(" %5i%16.6e%14.3e%14.3e\n", i + 1, x[i], d[i], g[i]);
}
