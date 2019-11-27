#include <stdio.h>
#include <math.h>

#define N 729
#define reps 1000 
#include <omp.h> 

double a[N][N], b[N][N], c[N];
int jmax[N];

void init1(void);
void init2(void);
void runloop(int); 
void loop1chunk(int, int);
void loop2chunk(int, int);
void valid1(void);
void valid2(void);
void fred(int niterations, int nthreads);

int main(int argc, char *argv[]) {

  double start1,start2,end1,end2;
  int r;

  init1();

  start1 = omp_get_wtime();

  for (r=0; r<reps; r++){
    runloop(1);
  }

  end1  = omp_get_wtime();

  valid1();

  printf("Total time for %d reps of loop 1 = %f\n",reps, (float)(end1-start1));



    init2();

    start2 = omp_get_wtime();

    for (r=0; r<reps; r++){
      runloop(2);
    }

    end2  = omp_get_wtime();

    valid2();

    printf("Total time for %d reps of loop 2 = %f\n",reps, (float)(end2-start2));

} 

void init1(void){
  int i,j; 

  for (i=0; i<N; i++){ 
    for (j=0; j<N; j++){ 
      a[i][j] = 0.0; 
      b[i][j] = 3.142*(i+j); 
    }
  }
}

void init2(void){ 
  int i,j, expr;

  for (i=0; i<N; i++){
    expr =  i%( 3*(i/30) + 1);  // i=0, 0  i=1, 1.1 i=2, 1.2

    if (expr == 0) {  // 给jmax赋值1/729
      jmax[i] = N;
    }
    else {
      jmax[i] = 1;
    }
    c[i] = 0.0;     // 给c赋值0
  }

  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      b[i][j] = (double) (i*j+1) / (double) (N*N);  // 给b重新赋值
    }
  }
}

void runloop(int loopid)  {
    
  int *ipt;

#pragma omp parallel num_threads(8) default(none) shared(loopid, ipt)
  {
    int myid  = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
#pragma omp single
  {
    ipt = (int*)malloc(sizeof(int)*nthreads);
  }
    ipt[myid] = (int) ceil((double)N/(double)nthreads); // 729 / 4 = 183  假如10个进程，729 / 10 = 73

    int chunksize;
    int lo, hi;

    lo = myid*ipt[myid];
    hi = (myid+1)*ipt[myid];
    if (hi > N) hi = N;

    while(ipt[myid]>0) {
      chunksize = (int)ceil((double)ipt[myid]/(double)nthreads);

      hi = lo + chunksize;  // 183 366 549 732        // 73 146 219 292 365 438 511 584 657 730
      if (hi > N) hi = N;  // 183 366 549 729

      //      printf("lo on thread %d is %d\n", myid, lo);
      //      printf("hi on thread %d is %d\n", myid, hi);









      switch (loopid) {
        case 1: loop1chunk(lo,hi); break;
        case 2: loop2chunk(lo,hi); break;
      }
      ipt[myid] -= chunksize;
      lo += chunksize;
      //      printf("this is thread %d, ipt is %d\n", myid, ipt[myid]);
    }
  }
}

void loop1chunk(int lo, int hi) { 
  int i,j; 

  for (i=lo; i<hi; i++){    // i = 0 183 366 549; i < 183 366 549 729
    for (j=N-1; j>i; j--){  // j = 728; j 728 727, ..., 545(取值)     j = 728; j 728, 727, ..., 0
      a[i][j] += cos(b[i][j]);
    } 
  }
}

void loop2chunk(int lo, int hi) {
  int i,j,k;
  double rN2;

  rN2 = 1.0 / (double) (N*N);  // 0.00000188167

  for (i=lo; i<hi; i++){
    for (j=0; j<jmax[i]; j++){
      for (k=0; k<j; k++){
	c[i] += (k+1) * log (b[i][j]) * rN2;
      } 
    }
  }
}

void valid1(void) { 
  int i,j; 
  double suma; 
  
  suma= 0.0; 
  for (i=0; i<N; i++){ 
    for (j=0; j<N; j++){ 
      suma += a[i][j];
    }
  }
  printf("Loop 1 check: Sum of a is %lf\n", suma);
}

void valid2(void) { 
  int i; 
  double sumc; 
  
  sumc= 0.0; 
  for (i=0; i<N; i++){ 
    sumc += c[i];
  }
  printf("Loop 2 check: Sum of c is %f\n", sumc);
}
