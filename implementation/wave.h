#include "matrix.h"



matrix_t *phi (const matrix_t* x,const matrix_t* y,const matrix_t* a,const matrix_t* b, const matrix_t* c, const matrix_t* d);
int* parite (int* parite_U, int* parite_V);
int* trapdoor (int lambda);
int* sign (int* sk, int m);
bool verify (int* pk, int m, int* signature);
matrix_t *syndrome (const matrix_t *e, const matrix_t *parite);
int inversion_of_f (int* parite_U, int* parite_V, int* inv);
int* invert_alg (int* sk, int* S);
int iteration_prange (int* parite, int syndrome);
int prange_alg (int* parite, int syndrome, int* info, int x);
