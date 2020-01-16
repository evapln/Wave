#include "matrix.h"

struct matrix_t {
  int nb_col;
  int nb_row;
  char** mat;
};

matrix_t *matrix_alloc(int row, int col)
{
  matrix_t *matrix = malloc(sizeof(matrix_t));
  if (!matrix)
  {
    return NULL;
  }
  matrix->nb_col = col;
  matrix->nb_row = row;
  matrix->mat = calloc(row, sizeof(char*));
  if (!matrix->mat)
  {
    free(matrix);
    return NULL;
  }
  for (int i = 0; i < row; i++)
  {
    matrix->mat[i] = calloc(col, sizeof(char));
    if (!matrix->mat[i])
    {
      for (int j = 0; j < i; ++j)
      {
         free(matrix->mat[j]);
      }
      free(matrix->mat);
      free(matrix);
      return NULL;
    }
  }
  return matrix;
}

matrix_t *matrix_copy(const matrix_t *matrix)
{
  if (!matrix)
    return NULL;
  matrix_t *copy = matrix_alloc(matrix->nb_row, matrix->nb_col);
  if (!copy)
    return NULL;
  for (int i = 0; i < copy->nb_row; ++i)
    for (int j = 0; j < copy->nb_col; ++j)
      copy->mat[i][j] = matrix->mat[i][j];
  return copy;
}

void matrix_free (matrix_t *matrix)
{
  if (matrix != NULL)
  {
    for (int i = 0; i < matrix->nb_row; ++i)
    {
      free(matrix->mat[i]);
    }
    free(matrix->mat);
    free(matrix);
  }
}


char matrix_get_cell(const matrix_t *matrix,const int row_val, const int col_val)
{
  if (!matrix)
    return 2;
  int row = matrix->nb_row;
  int col = matrix->nb_col;
  if (row_val > row || col_val > col)
    return 2;
  return matrix->mat[row_val][col_val];
}

int matrix_get_row(const matrix_t *matrix)
{
  if (!matrix)
    return 0;
  return matrix->nb_row;
}

int matrix_get_col(const matrix_t *matrix)
{
  if (!matrix)
    return 0;
  return matrix->nb_col;
}

void matrix_set_cell(matrix_t *matrix, const int row_val, const int col_val,
                   const char val) {
  if (!matrix)
    return;
  int row = matrix->nb_row;
  int col = matrix->nb_col;
  if (row_val >= row || col_val >= col)
    return;
  matrix->mat[row_val][col_val] = val;
}

void identity (matrix_t *matrix){
  int size = matrix->nb_row;
  if (size != matrix->nb_col)
    return;
  for (int i = 0; i < size; i++){
    matrix->mat[i][i] = 1;
    for (int j = 0; j < size; j++){
      if (j != i)
        matrix->mat[i][j] = 0;
    }
  }
}

void matrix_print(matrix_t *matrix, FILE *fd)
{
  if (matrix != NULL && matrix->mat != NULL)
  {
    int col = matrix->nb_col;
    int row = matrix->nb_row;
    for (int i = 0; i < row; ++i)
    {
      char val;
      for (int j = 0; j < col; ++j)
      {
        val = matrix_get_cell(matrix,i,j);
        fprintf(fd,"%d ", val);
      }
      fputs("\n", fd);
    }
    fputs("\n\n", fd);
  }
}

void matrix_add(matrix_t *matrix, const matrix_t *matrix1, const matrix_t *matrix2) {
  int col = matrix->nb_col;
  if (matrix1->nb_col != col || matrix2->nb_col != col)
    return;
  int row = matrix->nb_row;
  if (matrix1->nb_row != row || matrix2->nb_col != col)
      return;
  for (int i = 0; i < row; ++i)
    for (int j = 0; j < col; ++j)
      matrix->mat[i][j] = matrix1->mat[i][j] + matrix2->mat[i][j];
}

void matrix_scal(matrix_t *matrix, matrix_t *matrix1, matrix_t *matrix2) {
  int col = matrix->nb_col;
  if (matrix1->nb_col != col || matrix2->nb_col != col)
    return;
  int row = matrix->nb_row;
  if (matrix1->nb_row != row || matrix2->nb_col != col)
      return;
  for (int i = 0; i < row; ++i)
    for (int j = 0; j < col; ++j)
      matrix->mat[i][j] = matrix1->mat[i][j] * matrix2->mat[i][j];
}

matrix_t *matrix_concatenation(matrix_t *A, matrix_t *B){
  int row = A->nb_row;
  if (row != B->nb_row)
    return NULL;
  int col_A = A->nb_col;
  int col_B = B->nb_col;
  int col = col_A + col_B;
  matrix_t *conc = NULL;
  conc = matrix_alloc(row,col);
  for (int i = 0; i < row; i++){
    for (int j = 0; j < col_A; j++)
      conc->mat[i][j] = A->mat[i][j];
    for (int j = col_A; j < col_B; j++)
      conc->mat[i][j] = B->mat[i][j-col_A];
  }
  return conc;
}

void matrix_prod(matrix_t *matrix, matrix_t *matrix1, matrix_t *matrix2) {
  int k = matrix1->nb_col;
  if (k != matrix2->nb_row)
    return;
  int col = matrix->nb_col;
  int row = matrix->nb_row;
  if (matrix1->nb_row != row || matrix2->nb_col != col)
    return;
  for (int i = 0; i < row; i++){
    for (int j = 0; j < col; j++){
      matrix->mat[i][j] = 0;
      for (int ind = 0; ind < k; ind++)
        matrix->mat[i][j] += matrix1->mat[i][ind]*matrix2->mat[ind][j];
    }
  }
}

void matrix_inv(matrix_t *mat_inv, matrix_t *A){
  matrix_t *I = NULL;
  I = matrix_copy(A);
  identity(I);
  char** mat = A->mat;
  matrix_t *B = NULL;
  B = matrix_concatenation(A,I);
  char** inv = mat_inv->mat;
  int size = B->nb_row;
  for(int i = 0; i < size; i++) {
		if(B->mat[i][i] == 0)
		  return;
		for(int j=0;j<size;j++) {
			if(i!=j) {
			  char ratio = mul_Fq(B->mat[j][i], inv_Fq(B->mat[i][i]));
			  for(int k=0;k<2*size;k++)
					inv[j][k] = add_Fq(mat[j][k], -mul_Fq(ratio,mat[i][k]));
			}
    }
	}
	/* Row Operation to Make Principal Diagonal to 1 */
	for(int i=0;i<size;i++)
		for(int j=0;j<2*size;j++)
		  inv[i][j] = mul_Fq(inv[i][j], inv_Fq(inv[i][i]));
  matrix_free(I);
  matrix_free(B);
}
