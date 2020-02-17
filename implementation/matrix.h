#include "algebra.h"
#include <math.h>

typedef struct matrix_t matrix_t;


////////////////////////////////////////////////////////////////////////////////
////////////////////////// gestion de la mémoire ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Alloue en mémoire et renvoie une matrice de taille row x col vide */
matrix_t *matrix_alloc(const int row, const int col);

/* Libère la mémoire allouée pour matrix */
void matrix_free (matrix_t *matrix);




////////////////////////////////////////////////////////////////////////////////
////////////////////// gestion des paramètres de matrix_t //////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Renvoie le contenu de A[row_val][col_val] */
char matrix_get_cell(const matrix_t *A,const int row_val, const int col_val);

/* Renvoie le nombre de ligne de matrix */
int matrix_get_row(const matrix_t*matrix);

/* Renvoie le nombre de colonne de matrix */
int matrix_get_col(const matrix_t*matrix);

/* Mets  A[row][col] à val */
void matrix_set_cell(matrix_t *A, const int row, const int col, const char val);


////////////////////////////////////////////////////////////////////////////////
//////////////////////// création de matrices spécifiques //////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Alloue en mémoire une matrice de taille (row,col) et l'initialise à la valeur
   val partout */
matrix_t *matrix_init (const int row, const int col, const char val);

/* Alloue en mémoire et renvoie une copie de matrix */
matrix_t *matrix_copy(const matrix_t *matrix);

/* Copie matrix2 dans matrix1 sans toucher à aucune mémoire */
void matrix_copy2 (const matrix_t *matrix1, const matrix_t *matrix2);

/* Alloue en mémoire et renvoie une matrice identité de taille size */
matrix_t *matrix_identity (const int size);

/* Renvoie une matrice diagonale ayant val*vect en diagonale */
matrix_t *matrix_vect_to_diag (const matrix_t *vect, const char val);

/* Alloue en mémoire et renvoie une matrice aléatoire de taille row x col dont
   les coefficients sont dans Fq */
matrix_t *matrix_random(const int row, const int col);

/* Alloue et renvoie une matrice de permutation aléatoire de taille n */
matrix_t *matrix_perm_random(const int n);

/* Alloue en mémoire et renvoie une matrice de permutation aléatoire de taille n
  renvoyant info de longueur len_i sur les dernières coordonnées*/
matrix_t *matrix_perm_random_info(const int n, const int *info, const int len_i,
                                  const int dim);

/* Alloue en mémoire et renvoie la ransposée de matrix */
matrix_t *matrix_trans(const matrix_t *matrix);

/* Alloue en mémoire et renvoie la commatrice de mat */
matrix_t *matrix_com(const matrix_t *mat);

/* Alloue en mémoire et renvoie l'inverse de mat */
matrix_t *matrix_inv(const matrix_t *mat);

/* Alloue en mémoire et renvoie l'inverse de mat par la comatrice et le det */
matrix_t *matrix_inv_com(const matrix_t *A);

/* Alloue et renvoie la sous-matrice de A sans la ligne a et la colonne b */
matrix_t *matrix_sub (const matrix_t *A, const int a, const int b);

/* met dans A la partie gauche de matrix et dans B la partie droite */
void matrix_separate(const matrix_t *matrix, matrix_t *A, matrix_t *B);

/* Alloue et renvoie la concaténation de A et B, si mode = 0 ->(A|B),
   si mode != 0 ->  (A)
                    (B) */
matrix_t *matrix_concatenation(const matrix_t *A, const matrix_t *B, const int mode);

/* Alloue en mémoire et renvoie la matrice sum = matrix1 + matrix2 */
matrix_t *matrix_add(const matrix_t *matrix1, const matrix_t *matrix2);

/* met dest à dest + coef * src */
void matrix_add_modified(matrix_t *dest, const matrix_t *src, const char coef);

/* Alloue en mémoire et renvoie la matrice matrix dont tous les coefficients ont
   été multipliés par scal */
matrix_t *matrix_mul_by_scal(const matrix_t *matrix, const char scal);

/* Alloue en mémoire et renvoie la matrice prdouit de matrix1 et matrix2 */
matrix_t *matrix_prod(const matrix_t *matrix1, const matrix_t *matrix2);

/* alloue et renvoie la matrice de parité associée à la matrice génératrice */
matrix_t *matrix_parite(const matrix_t *gen);

/* alloue et renvoie une sous-matrice contenant les colonnes de A spécifées dans
   ind_col */
matrix_t *sub_col_matrix(const matrix_t *A, const int *ind_col, const int len);

/* systematise matrix */
void matrix_systematisation(matrix_t *matrix);



////////////////////////////////////////////////////////////////////////////////
///////////////////////////// calcul sur matrices //////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Renvoie le determinant de la matrice A */
char matrix_det(const matrix_t *A);

/* Renvoie le rang de la matrice A */
int matrix_rank(const matrix_t *A);


////////////////////////////////////////////////////////////////////////////////
//////////////////////////// test de la forme des matrices /////////////////////
////////////////////////////////////////////////////////////////////////////////

/* vérifie A==Id */
bool is_identity(const matrix_t *A);

/* Vérifie si la matrice est correctement systématisée */
bool matrix_is_syst(const matrix_t *matrix);

/* Vérifie si la matrice est triangulaire */
bool is_trigonalise (const matrix_t *A);



////////////////////////////////////////////////////////////////////////////////
///////////////////////////// utilisation de tableaux //////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* retourne true si val est dans array de longueur len, false sinon */
bool is_in_array(const int *array, const int len, const int val);

/* Mélange le vecteur array sur les n premières coordonnées */
void shuffle(int *array, const int n);

/* Mélange le vecteur array de taille len_a en mettant les indices contenus dans
   infos de taille len_i aux dim dernières coordonnées */
void shuffle_info(int *array, const int len_a, const int *info, const int len_i,
                  const int dim);



////////////////////////////////////////////////////////////////////////////////
////////////////////////////// gestion de vecteurs /////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* renvoie le poids du vecteur */
int weight(const matrix_t *vect);

/* renvoie le poids de vect sur les coordonnées de subset */
int sub_weight(const matrix_t *vect, const int *subset, const int len_s);

/* alloue en mémoire et renvoie un vecteur aléatoire de taille n */
matrix_t *vector_rand(const int n);

/* alloue en mémoire et renvoie un vecteur aléatoire de poids w et de taille n */
matrix_t *vector_rand_weight(const int n, const int w);

/* renvoie un vecteur aléatoire de taille n avec son poids sur les coordonnées
   de info égal à t */
matrix_t *vector_rand_sub_weight(const int n, const int *info, const int len_i,
                             const int t);

/* Alloue et renvoie le produit élément par élément de vect1 et vect2 */
matrix_t *vect_scal(const matrix_t *vect1, const matrix_t *vect2);



////////////////////////////////////////////////////////////////////////////////
////////////////// calculs et tests sur les lignes d'une matrice ///////////////
////////////////////////////////////////////////////////////////////////////////

/* met A[i][*] dans ligne */
void matrix_row(matrix_t *ligne, const matrix_t *A, const int row);

/* Supprime row1 de matrix. Alloue dans une nouvelle matrice et la renvoie */
matrix_t *matrix_del_row(const matrix_t *matrix, const int row1);

/* Regarde si row de matrix est entièrement nulle */
bool row_is_zero (const matrix_t *matrix, const int row);

/* Dans matrix, row1 prends row1 + coef*Row2 */
void matrix_add_row(matrix_t *matrix, const int row1, const int row2, const char coef);

/* Dans matrix, row1 prends coef*row1 */
void matrix_mul_row(matrix_t *matrix, const int row1, const char coef);

/* Echange row1 et row2 dans matrix */
void matrix_exchange_row(matrix_t *matrix, const int row1, const int row2);

/* Supprime les lignes nulles dans matrix. Renvoie une nouvelle matrice */
matrix_t *matrix_del_null_row (const matrix_t *matrix) ;



////////////////////////////////////////////////////////////////////////////////
/////////////////////////// gestion d'un code //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* met c à un mot aléatoire du code de matrice génératrice G */
void random_word(matrix_t *c, const matrix_t *G);



////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// affichage ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Ecrit matrix dans fd */
void matrix_print(const matrix_t *matrix, FILE *fd);
