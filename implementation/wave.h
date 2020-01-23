#include "matrix.h"



matrix_t *phi (const matrix_t* x,const matrix_t* y);
matrix_t *syndrome (const matrix_t *e, const matrix_t *parite);
matrix_t *parite (const matrix_t *parite_U, const matrix_t *parite_V);
void key_gen (int lambda);

int* sign (int* sk, int m);
bool verify (int* pk, int m, int* signature);
int inversion_of_f (int* parite_U, int* parite_V, int* inv);
int* invert_alg (int* sk, int* S);
int iteration_prange (int* parite, int syndrome);
int prange_alg (int* parite, int syndrome, int* info, int x);
