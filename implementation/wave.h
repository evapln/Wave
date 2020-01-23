#include "matrix.h"

typedef struct sk_t sk_t;

/* Applique la fonction Phi sur x et y */
matrix_t *phi (const matrix_t* x,const matrix_t* y);

/* Renvoie s le syndrome de e pour le code associé à la matrice de parité */
matrix_t *syndrome (const matrix_t *e, const matrix_t *parite);

/* Renvoie parité du code UV en fonction des matrices de parité des codes U et V */
matrix_t *parite (const matrix_t *parite_U, const matrix_t *parite_V);

/* Si mode = 0, a=c=d=1 et b=0, si mode = 1 a,b,c,d random avec les bonnes propriétés */
void coeff_phi (int mode);

/* Génère les clés */
void key_gen (int lambda, matrix_t *pk, sk_t *sk, int mode);

int* sign (int* sk, int m);
bool verify (int* pk, int m, int* signature);
int inversion_of_f (int* parite_U, int* parite_V, int* inv);
int* invert_alg (int* sk, int* S);
int iteration_prange (int* parite, int syndrome);
int prange_alg (int* parite, int syndrome, int* info, int x);
