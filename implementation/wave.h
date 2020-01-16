#include "matrix.h"

int* phi (int x,int y, int a, int b, int c, int d);
int* parite (int* parite_U, int* parite_V);
int* trapdoor (int lambda);
int* sign (int* sk, int m);
bool verify (int* pk, int m, int* signature);
int* syndrome (int e, int* parite_UV);
int inversion_of_f (int* parite_U, int* parite_V, int* inv);
int* invert_alg (int* sk, int* S);
int iteration_prange (int* parite, int syndrome);
int prange_alg (int* parite, int syndrome, int* info, int x);
