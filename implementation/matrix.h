#include "algebra.h"
#include <math.h>

typedef struct matrix_t matrix_t;

/* Alloue en mémoire et renvoie une matrice de taille row x col vide */
matrix_t *matrix_alloc(int row, int col);

/* Libère la mémoire allouée pour matrix */
void matrix_free (matrix_t *matrix);

/* Renvoie le contenu de matrix[row_val][col_val] */
char matrix_get_cell(const matrix_t *matrix,const int row_val, const int col_val);

/* Renvoie le nombre de ligne de matrix */
int matrix_get_row(const matrix_t*matrix);

/* Renvoie le nombre de colonne de matrix */
int matrix_get_col(const matrix_t*matrix);

/* Mets  matrix[row_val][col_val] à val */
void matrix_set_cell(matrix_t *matrix, const int row_val, const int col_val, const char val);

/* Alloue en mémoire une matrice de taille (row,col) et l'initialise à la valeur val partout */
matrix_t *matrix_init (const int row, const int col, const char val);

/* Alloue en mémoire et renvoie une copie de matrix */
matrix_t *matrix_copy(const matrix_t *matrix);

/* Alloue en mémoire et renvoie une matrice identité de taille size */
matrix_t *matrix_identity (int size);

/* Renvoie une matrice diagonale ayant val*vect en diagonale */
matrix_t *matrix_vect_to_diag (const matrix_t *vect, const char val);

/* Alloue en mémoire et renvoie une matrice aléatoire de taille row x col dont les coefficients sont dans Fq */
matrix_t *matrix_random(const int row, const int col);

/* Mélange le vecteur array */
void shuffle(int *array, int n);

/* Alloue en mémoire et renvoie une matrice de permutation aléatoire */
matrix_t *matrix_perm_random (const int n);

/* Alloue en mémoire et renvoie la ransposée de matrix */
matrix_t *matrix_trans(const matrix_t *matrix);

/* Alloue en mémoire et renvoie la commatrice de mat */
matrix_t *matrix_com(matrix_t *mat);

/* Alloue en mémoire et renvoie l'inverse de mat */
matrix_t *matrix_inv(matrix_t *mat);

/* Alloue en mémoire et renvoie la sous-matrice de A sans la ligne a et la colonne b */
matrix_t *matrix_sub (const matrix_t *A, int a, int b);

/* Alloue en mémoire et renvoie la concaténation de A et B, si mode = 0 ->(A|B), si mode != 0 ->  (A)
                                                                                                  (B) */
matrix_t *matrix_concatenation(const matrix_t *A, const matrix_t *B, const int mode);

/* Alloue en mémoire et renvoie la matrice sum = matrix1 + matrix2 */
matrix_t *matrix_add(const matrix_t *matrix1, const matrix_t *matrix2);

/* Alloue en mémoire et renvoie le produit élément par élément de vect1 et vect2 */
matrix_t *vect_scal(const matrix_t *vect1, const matrix_t *vect2);

/* Alloue en mémoire et renvoie la matrice matrix dont tous les coefficients ont été multipliés par scal */
matrix_t *matrix_mul_by_scal(const matrix_t *matrix, const char scal);

/* Alloue en mémoire et renvoie la matrice prdouit de matrix1 et matrix2 */
matrix_t *matrix_prod(const matrix_t *matrix1, const matrix_t *matrix2);

/* Supprime row1 de matrix. Alloue dans une nouvelle matrice et la renvoie */
matrix_t *matrix_del_row(matrix_t *matrix, const int row1);

/* Regarde si row de matrix est entièrement nulle */
bool row_is_zero (matrix_t *matrix, const int row);

/* Dans matrix, row1 prends row1 + coef*Row2 */
void matrix_add_row(matrix_t *matrix, const int row1, const int row2, const char coef);

/* Dans matrix, row1 prends coef*row1 */
void matrix_mul_row(matrix_t *matrix, const int row1, const char coef);

/* Echange row1 et row2 dans matrix */
void matrix_exchange_row(matrix_t *matrix, const int row1, const int row2);

/* Supprime les lignes nulles dans matrix. Renvoie une nouvelle matrice */
matrix_t *matrix_del_null_row (matrix_t *matrix) ;

/* alloue en mémoire et renvoie la matrice de parité associée à la matrice génératrice */
matrix_t *matrix_parite(const matrix_t *gen);

/* Trigonalise matrix */
void matrix_systematisation(matrix_t *matrix);

/* Vérifie si la matrice est correctement systématisée */
bool matrix_is_syst(matrix_t *matrix);

/* Renvoie le determinant de la matrice A */
char matrix_det(matrix_t *A);

/* Ecrit matrix dans fd */
void matrix_print(matrix_t *matrix, FILE *fd);
