#define _GNU_SOURCE
#include <dlfcn.h>
#include <lapacke.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

/* Transparent audit interceptor: the actual LAPACKE call remains authoritative. */
static void write_matrix(const char *dir, int call, const char *name,
                         const double *a, lapack_int n, lapack_int ld, int layout) {
  char path[4096];
  snprintf(path, sizeof path, "%s/call_%02d_%s.mtx", dir, call, name);
  FILE *f = fopen(path, "w");
  if (!f) { fprintf(stderr,"capture open %s: %s\n",path,strerror(errno)); exit(97); }
  fprintf(f,"%%%%MatrixMarket matrix array real general\n%% Exact-double LAPACKE audit capture; column-major Matrix Market order.\n%d %d\n", (int)n, (int)n);
  for (lapack_int j=0;j<n;++j) for (lapack_int i=0;i<n;++i)
    fprintf(f,"%.17g\n",layout==LAPACK_COL_MAJOR?a[(size_t)j*ld+i]:a[(size_t)i*ld+j]);
  if (fclose(f)) exit(98);
}

lapack_int LAPACKE_dsygv(int layout, lapack_int itype, char jobz, char uplo,
                        lapack_int n, double *a, lapack_int lda,
                        double *b, lapack_int ldb, double *w) {
  typedef lapack_int (*real_fn)(int,lapack_int,char,char,lapack_int,double*,lapack_int,double*,lapack_int,double*);
  static int count=0;
  real_fn real = (real_fn)dlsym(RTLD_NEXT,"LAPACKE_dsygv");
  if (!real) { fprintf(stderr,"Cannot resolve real LAPACKE_dsygv\n"); exit(96); }
  const char *dir=getenv("M6A_LAPACK_CAPTURE_DIR");
  int call=++count;
  if (dir && *dir) {
    write_matrix(dir,call,"stiffness_matrix",a,n,lda,layout);
    write_matrix(dir,call,"mass_matrix",b,n,ldb,layout);
  }
  lapack_int info=real(layout,itype,jobz,uplo,n,a,lda,b,ldb,w);
  fprintf(stderr,"M6A_LAPACK_CALL call=%d layout=%d itype=%d jobz=%c uplo=%c n=%d info=%d\n",call,layout,(int)itype,jobz,uplo,(int)n,(int)info);
  if (dir && *dir && info==0) {
    write_matrix(dir,call,"eigenvectors",a,n,lda,layout);
    char path[4096]; snprintf(path,sizeof path,"%s/call_%02d_eigenvalues.csv",dir,call);
    FILE *f=fopen(path,"w"); if (!f) exit(97);
    fprintf(f,"index,lambda,kc\n");
    for(lapack_int i=0;i<n;++i) fprintf(f,"%d,%.17g,%.17g\n",(int)i,w[i],w[i]>0?sqrt(w[i]):0);
    if (fclose(f)) exit(98);
  }
  return info;
}
