#include "matrix.h"

typedef struct sk_t sk_t;

typedef struct keys_t keys_t;

typedef struct sign_t sign_t;

////////////////////////////////////////////////////////////////////////////////
//////////////////////////// gestion de la mémoire /////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Alloue en mémoire l'espace d'une clé secrète */
sk_t *sk_alloc(int dim_U, int dim_V, int dim);

/* Libère l'espace alloué pour la clé secrète et ses composants */
void sk_free (sk_t *sk);

/* Alloue en mémoire l'espace des deux clés */
keys_t *key_alloc(int dim_U, int dim_V, int dim);

/* Libère l'espace alloué pour les clés */
void key_free(keys_t *keys);

/* Alloue en mémoire l'espace de la signature */
sign_t *sign_alloc(void);

/* Libère l'espace alloué pour la signature */
void sign_free(sign_t *signature);

////////////////////////////////////////////////////////////////////////////////
/////////////////////////// fonctions secondaires //////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Applique la fonction Phi sur x et y */
matrix_t *phi(const matrix_t* x,const matrix_t* y);

/* Renvoie s le syndrome de e pour le code associé à la matrice de parité */
matrix_t *syndrome(const matrix_t *e, const matrix_t *parite);

/* Renvoie la matrice de parité du code UV en fonction des matrices de parité
   des codes U et V */
matrix_t *parite(const matrix_t *parite_U, const matrix_t *parite_V);

/* Si mode = 0, a=c=d=1 et b=0,
   Si mode = 1 a,b,c,d random avec les bonnes propriétés */
void coeff_phi(int mode);


////////////////////////////////////////////////////////////////////////////////
////////////////////////////// génération des clés /////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Génère les clés */
keys_t *key_gen(int mode);

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// décodage /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* decodage par ensemble d'information :
  G est la matrice génératrice du code V
  synd vecteur ligne est le syndrome cherché
  ev vecteur ligne est la sortie */
void decode_ev(matrix_t * ev, const matrix_t *G, const matrix_t *synd);

/* decodage par ensemble d'information :
  G est la matrice [-A|Id] avec A partie droite de H_U systematisée
  synd vecteur ligne est le syndrome cherché
  ev vecteur ligne estla sortie de decode_ev
  eu vecteur ligne est la sortie */
void decode_eu(matrix_t * eu, const sk_t *sk, const matrix_t *synd,
               const matrix_t *ev);

 /* decodage par ensemble d'information :
   synd vecteur ligne est le syndrome cherché
   alloue en mémoire et renvoie e vecteur ligne tq son syndrome vaille synd */
matrix_t *decode_uv(const sk_t *sk, const matrix_t *synd);

/* inverse la fonction syndrome */
matrix_t *invert_alg(const sk_t *sk, const matrix_t *synd);

/* alloue en mémoire en renvoie le vecteur x tel que Ax=y */
// matrix_t *resol_syst(const matrix_t *A, const matrix_t *y);

// int *freeset(const matrix_t *H, const matrix_t *ev, const int k);


////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// PRANGE ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Renvoie un tableau de longueur len aléatoire dont les coefficients sont
   inférieurs à n */
void infoset(int *info, const int n, const int len);

/* Calcul de e, appelé par iteration_prange */
matrix_t *prange_algebra(const matrix_t *parite, const matrix_t *syndrome,
                         const int *info, const int len_i, const matrix_t *x);

/* Renvoie e tel que eH_T = s */
matrix_t *iteration_prange(const matrix_t *parite, const matrix_t *syndrome);




////////////////////////////////////////////////////////////////////////////////
////////////////////////// signature et vérification ///////////////////////////
////////////////////////////////////////////////////////////////////////////////

matrix_t *hash(const matrix_t *m, const matrix_t *r, const matrix_t *pk);

sign_t *sign(const keys_t *keys, const matrix_t *m);

bool verify(const matrix_t *pk, const matrix_t *m, const sign_t *signature);
