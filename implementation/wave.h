#include "matrix.h"

typedef struct sk_t sk_t;

typedef struct keys_t keys_t;

/* Alloue en mémoire l'espace d'une clé secrète */
sk_t *sk_alloc(int dim_U, int dim_V, int dim);

/* Libère l'espace alloué pour la clé secrète et ses composants */
void sk_free (sk_t *sk);

/* Alloue en mémoire l'espace des deux clés */
keys_t *key_alloc(int dim_U, int dim_V, int dim);

/* Libère l'espace alloué pour les clés */
void key_free (keys_t *keys);

/* Applique la fonction Phi sur x et y */
matrix_t *phi (const matrix_t* x,const matrix_t* y);

/* Renvoie s le syndrome de e pour le code associé à la matrice de parité */
matrix_t *syndrome (const matrix_t *e, const matrix_t *parite);

/* Renvoie la matrice de parité du code UV en fonction des matrices de parité des codes U et V */
matrix_t *parite (const matrix_t *parite_U, const matrix_t *parite_V);

/* Si mode = 0, a=c=d=1 et b=0, si mode = 1 a,b,c,d random avec les bonnes propriétés */
void coeff_phi (int mode);

/* Génère les clés */
keys_t *key_gen (int lambda, int mode);

// TO DO //
int* sign (int* sk, int m);
bool verify (int* pk, int m, int* signature);
int inversion_of_f (int* parite_U, int* parite_V, int* inv);
int* invert_alg (int* sk, int* S);
int iteration_prange (int* parite, int syndrome);
int prange_alg (int* parite, int syndrome, int* info, int x);
