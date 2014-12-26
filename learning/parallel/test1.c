#include "omp.h"
int main()
{
  int X=0;
  int nthreads= omp_get_num_threads();	//here is 1
  printf("nthreads=%d\n", nthreads);
  printf("time start %f\n", omp_get_wtime());
  #pragma omp parallel
  {
	int nnthreads= omp_get_num_threads();	//here is 8 in omp threads
	printf(">>nnthreads=%d\n", nnthreads);
	int id = omp_get_thread_num();
	printf("hello(%d)", id);
	printf(" world(%d) ..\n",id);
  #pragma omp barrier
	printf("==>pass barrier:%d\n",id);

  #pragma omp critical
	printf(">>entering cr %d\n", id);

  #pragma omp atomic
	X++;
	printf("atomic X=%d by id=%d\n", X,id);
  }
  printf("time end %f\n", omp_get_wtime());
}