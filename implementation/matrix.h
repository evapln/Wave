#include "algebra.h"
#include <math.h>

typedef struct matrix_t matrix_t;

matrix_t *matrix_alloc(int row, int col);
matrix_t *matrix_copy(const matrix_t *matrix);
matrix_t *matrix_identity (int size);
matrix_t *matrix_trans(const matrix_t *matrix);
void matrix_free (matrix_t *matrix);
char matrix_get_cell(const matrix_t *matrix,const int row_val, const int col_val);
int matrix_get_row(const matrix_t*matrix);
int matrix_get_col(const matrix_t*matrix);
void matrix_set_cell(matrix_t *matrix, const int row_val, const int col_val, const char val);
void matrix_print(matrix_t *matrix, FILE *fd);
void matrix_add(matrix_t *matrix, const matrix_t *matrix1, const matrix_t *matrix2);
void vect_scal(matrix_t *matrix, matrix_t *matrix1, matrix_t *matrix2);
void matrix_scal(matrix_t *matrix, matrix_t *matrix1, char scal);
void matrix_prod(matrix_t *matrix, matrix_t *matrix1, matrix_t *matrix2);
matrix_t *matrix_com(matrix_t *mat);
matrix_t *matrix_inv(matrix_t *mat);
char matrix_det(matrix_t *A);
matrix_t *matrix_sub (const matrix_t *A, int a, int b);
