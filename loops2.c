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
    expr =  i%( 3*(i/30) + 1);

    if (expr == 0) {
      jmax[i] = N;
    }
    else {
      jmax[i] = 1;
    }
    c[i] = 0.0;
  }

  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      b[i][j] = (double) (i*j+1) / (double) (N*N);
    }
  }
}

void runloop(int loopid)  {
    
  int *ipt;
  int *lo;
  int *hi;

#pragma omp parallel default(none) shared(loopid, ipt, lo, hi)
  {
    int myid  = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
#pragma omp single
  {
    ipt = (int*)malloc(sizeof(int)*nthreads);
    lo = (int*)malloc(sizeof(int)*nthreads);
    hi = (int*)malloc(sizeof(int)*nthreads);
  }
  ipt[myid] = (int) ceil((double)N/(double)nthreads); // 729 / 4 = 183
  lo[myid] = myid*ipt[myid];
  hi[myid] = (myid+1)*ipt[myid];
  if (hi[myid] > N) hi[myid] = N;

  int lo_exec = 0;     // temporary value
  int hi_exec = 0;     // temporary value
  int next_thread = 0; // loop condition, represents which thread will be executed next
  int chunksize = 0;   // chunksize
  int remaining = 0;   // remaining iterations in a thread
  /* 
   * These two variables are for after running out local set,
   * begining to search other sets.
   * max: the number of most iterations of some thread
   * maxthread: the id of thread with most iterations
   */
  int max;
  int maxthread;

  while(next_thread != -1) {

    switch (loopid) {
      case 1: loop1chunk(lo_exec,hi_exec); break;
      case 2: loop2chunk(lo_exec,hi_exec); break;
    }
#pragma omp critical
  {
    /* There is still remaining in local set, prioritize local tasks */
    if (hi[myid] - lo[myid] > 0) {
      next_thread = myid;

      remaining = hi[next_thread]-lo[next_thread]; 
      chunksize = (int)ceil((double)remaining/(double)nthreads);

      lo_exec = lo[next_thread];
      hi_exec = lo_exec + chunksize;
      lo[next_thread] = lo[next_thread] + chunksize;
    }
    /* There is no remaining in local set, finding in other threads */
    else {
      max = 0;
      maxthread = -1;

      int i;
      /* search thread with the most iterations */
      for (i=0; i<nthreads; i++) {
	if (myid == i)
	  continue;
	if (hi[i] - lo[i] > max) {
	  max = hi[i] - lo[i];  // max iterations
	  maxthread = i;  // thread index
	}
      }
      next_thread = maxthread;
      if (next_thread != -1) {
	int remaining = hi[next_thread]-lo[next_thread];  
	chunksize = (int)ceil((double)remaining/(double)nthreads);

	lo_exec = lo[next_thread];
	hi_exec = lo_exec + chunksize;
	lo[next_thread] = lo[next_thread] + chunksize;
      }
    }
  } // critical
  } // while
  } // parallel
  free(ipt);
  free(lo);
  free(hi);
}

void loop1chunk(int lo, int hi) { 
  int i,j; 

  for (i=lo; i<hi; i++){
    for (j=N-1; j>i; j--){
      a[i][j] += cos(b[i][j]);
    } 
  }
}

void loop2chunk(int lo, int hi) {
  int i,j,k;
  double rN2;

  rN2 = 1.0 / (double) (N*N);

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
