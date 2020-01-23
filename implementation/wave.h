#include "matrix.h"

/* Applique la fonction Phi sur x et y */
matrix_t *phi (const matrix_t* x,const matrix_t* y);

/* Renvoie s le syndrome de e pour le code associé à la matrice de parité */
matrix_t *syndrome (const matrix_t *e, const matrix_t *parite);

/* Renvoie parité du code UV en fonction des matrices de parité des codes U et V */
matrix_t *parite (const matrix_t *parite_U, const matrix_t *parite_V);

/* Génère les clés */
void key_gen (int lambda);

int* sign (int* sk, int m);
bool verify (int* pk, int m, int* signature);
int inversion_of_f (int* parite_U, int* parite_V, int* inv);
int* invert_alg (int* sk, int* S);
int iteration_prange (int* parite, int syndrome);
int prange_alg (int* parite, int syndrome, int* info, int x);
