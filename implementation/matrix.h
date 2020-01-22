#include "algebra.h"
#include <math.h>

typedef struct matrix_t matrix_t;

/* alloue en mémoire et renvoie une matrice de taille row x col vide */
matrix_t *matrix_alloc(int row, int col);

/* libère la mémoire allouée pour matrix */
void matrix_free (matrix_t *matrix);

/* renvoie le contenu de matrix[row_val][col_val] */
char matrix_get_cell(const matrix_t *matrix,const int row_val, const int col_val);

/* renvoie le nombre de ligne de matrix */
int matrix_get_row(const matrix_t*matrix);

/* renvoie le nombre de colonne de matrix */
int matrix_get_col(const matrix_t*matrix);

/* mets  matrix[row_val][col_val] à val */
void matrix_set_cell(matrix_t *matrix, const int row_val, const int col_val, const char val);

/* alloue en mémoire et renvoie une copie de matrix */
matrix_t *matrix_copy(const matrix_t *matrix);

/* alloue en mémoire et renvoie une matrice identité de taille size */
matrix_t *matrix_identity (int size);

/* alloue en mémoire et renvoie une matrice aléatoire de taille row x col dont les coefficient sont dans Fq */
matrix_t *matrix_random(const int row, const int col);
void shuffle(int *array, int n);
matrix_t *matrix_perm_random (const int n);

/* alloue en mémoire et erenvoie la ransposée de matrix */
matrix_t *matrix_trans(const matrix_t *matrix);

/* alloue en mémoire et renvoie la commatrice de mat */
matrix_t *matrix_com(matrix_t *mat);

/* alloue en mémoire et renvoie l'inverse de mat */
matrix_t *matrix_inv(matrix_t *mat);

/* alloue en mémoire et renvoie la sous-matrice de A sans la lignaa et la colonne b */
matrix_t *matrix_sub (const matrix_t *A, int a, int b);

/* alloue en mémoire et renvoie la concaénation de A et B */
matrix_t *matrix_concatenation(const matrix_t *A, const matrix_t *B);

/* alloue en mémoire et renvoie la matrice sum = matrix1 + matrix2 */
matrix_t *matrix_add(const matrix_t *matrix1, const matrix_t *matrix2);

/* alloue en mémoire et renvoie le produit élément par élément de vect1 et vect2 */
matrix_t *vect_scal(const matrix_t *vect1, const matrix_t *vect2);

/* alloue en mémoire et renvoie la matrice matrix dont tous les coefficients ont été multipliés par scal */
matrix_t *matrix_mul_by_scal(const matrix_t *matrix, const char scal);

/* alloue en mémoire et renvoie la matrice prdouit de matrix1 et matrix2 */
matrix_t *matrix_prod(const matrix_t *matrix1, const matrix_t *matrix2);

matrix_t *matrix_del_row(matrix_t *matrix, const int row1);

bool row_is_zero (matrix_t *matrix, const int row);

void matrix_add_row(matrix_t *matrix, const int row1, const int row2, const char coef);
void matrix_mul_row(matrix_t *matrix, const int row1, const char coef);
void matrix_exchange_row(matrix_t *matrix, const int row1, const int row2);
matrix_t *matrix_del_null_row (matrix_t *matrix) ;
/*trigonalise matrix */
void matrix_trigonalisation(matrix_t *matrix);

/* renvoie le determinant de la matrice A */
char matrix_det(matrix_t *A);

/* écrit matrix dans fd */
void matrix_print(matrix_t *matrix, FILE *fd);
