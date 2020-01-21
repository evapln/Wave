#include "matrix.h"

struct matrix_t {
  int nb_col;
  int nb_row;
  char** mat;
};

matrix_t *matrix_alloc(int row, int col) {
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

void matrix_free (matrix_t *matrix) {
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

char matrix_get_cell(const matrix_t *matrix,const int row_val, const int col_val) {
  if (!matrix)
    return 2;
  int row = matrix->nb_row;
  int col = matrix->nb_col;
  if (row_val > row || col_val > col)
    return 2;
  return matrix->mat[row_val][col_val];
}

int matrix_get_row(const matrix_t *matrix) {
  if (!matrix)
    return 0;
  return matrix->nb_row;
}

int matrix_get_col(const matrix_t *matrix) {
  if (!matrix)
    return 0;
  return matrix->nb_col;
}

void matrix_set_cell(matrix_t *matrix, const int row_val, const int col_val, const char val) {
  if (!matrix)
    return;
  int row = matrix->nb_row;
  int col = matrix->nb_col;
  if (row_val >= row || col_val >= col)
    return;
  char tmp = val % ORDER;
  if (tmp < 0)
    tmp += ORDER;
  matrix->mat[row_val][col_val] = tmp;
}

matrix_t *matrix_copy(const matrix_t *matrix) {
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

matrix_t *matrix_identity (int size) {
  matrix_t *id = NULL;
  id = matrix_alloc(size,size);
  for (int i = 0; i < size; i++){
    id->mat[i][i] = 1;
    for (int j = 0; j < size; j++){
      if (j != i)
        id->mat[i][j] = 0;
    }
  }
  return id;
}

matrix_t *matrix_random(const int row, const int col) {
  matrix_t *rand = NULL;
  rand = matrix_alloc(row,col);
  for (int i = 0; i < row; ++i)
    for (int j = 0; j < col; ++j)
      rand->mat[i][j] = rand_Fq();
  return rand;
}

matrix_t *matrix_trans(const matrix_t *matrix) {
  if (!matrix)
    return NULL;
  matrix_t *trans = matrix_alloc(matrix->nb_col, matrix->nb_row);
  if (!trans)
    return NULL;
  for (int i = 0; i < trans->nb_row; ++i)
    for (int j = 0; j < trans->nb_col; ++j)
      trans->mat[i][j] = matrix->mat[j][i];
  return trans;
}

matrix_t *matrix_com(matrix_t *A) {
  int col = A->nb_col;
  int row = A->nb_row;
  matrix_t *com = NULL;
  com = matrix_alloc(row,col);
  for (int i = 0; i < row; i++)
    for (int j = 0; j < col; j++){
      matrix_t *sub = NULL;
      sub = matrix_sub(A,i,j);
      com->mat[i][j] = mul_Fq(pow(-1,i+j),matrix_det(sub));
      matrix_free(sub);
    }
  return com;
}

matrix_t *matrix_inv(matrix_t *A) {
  // calcul du dterminant
  char det = matrix_det(A);
  if (det == 0) {
    puts("non inversible");
    return NULL;
  }
  // calcul de la comatrice
  matrix_t *com = NULL;
  com = matrix_com(A);
  // calcul de la transposé de la comatrice
  matrix_t *trans = NULL;
  trans = matrix_trans(com);
  // calcul final de l'inverse
  matrix_t *inv = NULL;
  inv = matrix_mul_by_scal(trans, inv_Fq(det));
  // clean up
  matrix_free(com);
  matrix_free(trans);
  return inv;
}

matrix_t *matrix_sub (const matrix_t *A, int a, int b) {
  matrix_t *sub = NULL;
  int row = A->nb_row;
  int col = A->nb_col;
  sub = matrix_alloc(row - 1, col - 1);
  int ind_row = 0;
  int ind_col = 0;
  for (int i = 0; i < row; i++){
    ind_col = 0;
    if (i != a){
      for (int j = 0; j < col; j++){
        if (j != b){
          sub->mat[ind_row][ind_col] = A->mat[i][j];
          ind_col++;
        }
      }
      ind_row++;
    }
  }
  return sub;
}

matrix_t *matrix_concatenation(const matrix_t *A, const matrix_t *B) {
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

matrix_t *matrix_add(const matrix_t *matrix1, const matrix_t *matrix2) {
  int row = matrix1->nb_row;
  if (matrix2->nb_row != row)
    return NULL;
  int col = matrix1->nb_col;
  if (matrix2->nb_col != col)
    return NULL;
  matrix_t *sum = NULL;
  sum = matrix_alloc(row, col);
  for (int i = 0; i < row; ++i)
    for (int j = 0; j < col; ++j)
      sum->mat[i][j] = add_Fq(matrix1->mat[i][j], matrix2->mat[i][j]);
  return sum;
}

matrix_t *vect_scal(const matrix_t *vect1, const matrix_t *vect2) {
  int row = vect1->nb_row;
  if (vect2->nb_row != row)
    return NULL;
  int col = vect1->nb_col;
  if (vect2->nb_col != col)
    return NULL;
  matrix_t *mul = NULL;
  mul = matrix_alloc(row,col);
  for (int i = 0; i < row; ++i)
    for (int j = 0; j < col; ++j)
      mul->mat[i][j] = mul_Fq(vect1->mat[i][j], vect2->mat[i][j]);
  return mul;
}

matrix_t *matrix_mul_by_scal(const matrix_t *matrix, const char scal) {
  int row = matrix->nb_row;
  int col = matrix->nb_col;
  matrix_t *mul = NULL;
  mul = matrix_alloc(row,col);
  for (int i = 0; i < row; ++i)
    for (int j = 0; j < col; ++j)
      mul->mat[i][j] = mul_Fq(matrix->mat[i][j], scal);
  return mul;
}

matrix_t *matrix_prod(const matrix_t *matrix1, const matrix_t *matrix2) {
  int k = matrix1->nb_col;
  if (k != matrix2->nb_row)
    return NULL;
  int row = matrix1->nb_row;
  int col = matrix2->nb_col;
  matrix_t *prod = NULL;
  prod = matrix_alloc(row,col);
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      prod->mat[i][j] = 0;
      for (int ind = 0; ind < k; ind++)
        prod->mat[i][j] = add_Fq(mul_Fq(matrix1->mat[i][ind], matrix2->mat[ind][j]), prod->mat[i][j]);
    }
  }
  return prod;
}

void matrix_add_row(matrix_t *matrix, const int row1, const int row2, const char coef) {
  int row = matrix->nb_row;
  if (row1 > row || row2 > row)
    return;
  for (int i = 0; i < matrix->nb_col; ++i)
    matrix->mat[row1][i] = add_Fq(matrix->mat[row1][i], mul_Fq(coef,matrix->mat[row2][i]));
}

void matrix_mul_row(matrix_t *matrix, const int row1, const char coef) {
  if (row1 > matrix->nb_row)
    return;
  for (int i = 0; i < matrix->nb_col; ++i)
    matrix->mat[row1][i] = mul_Fq(coef,matrix->mat[row1][i]);
}

void matrix_exchange_row(matrix_t *matrix, const int row1, const int row2) {
  int row = matrix->nb_row;
  if (row1 > row || row2 > row)
    return;
  char tmp;
  for (int i = 0; i < matrix->nb_col; ++i) {
    tmp = matrix->mat[row1][i];
    matrix->mat[row1][i] = matrix->mat[row2][i];
    matrix->mat[row2][i] = tmp;
  }
}

void matrix_trigonalisation(matrix_t *matrix) {
  int row = matrix->nb_row;
  int col = matrix->nb_col;
  // étape 1
  // pour toutes les colonnes
  for (int j = 0; j < col; ++j) {
    // parcours de chaque colonne
    for (int i = j; i < row; ++i) {
      // on met un coefficient non nul en haut de chaque colonne
      if (i != j && matrix->mat[i][j] % ORDER != 0) {
        matrix_exchange_row(matrix, j, i);
        break;
      }
    }
    // annulation du reste de la colonne
    for (int i = j + 1; i < row; ++i) {
      char coef = mul_Fq(matrix->mat[i][j], inv_Fq(matrix->mat[j][j]));
      if (coef != 0) {
        matrix_add_row(matrix, i, j, -coef);
      }
    }
  }
  // on met tous les pivits à 1
  // pour toutes les lignes
  for (int i = 0; i < row; ++i) {
    // on parcours les lignes
    for (int j = 0; j < col; ++j) {
      // on trouve le premier coef non nul, le pivot
      if (matrix->mat[i][j] != 0) {
        matrix_mul_row(matrix, i, inv_Fq(matrix->mat[i][j]));
        break;
      }
    }
  }
  // reduction à revoir
  // pour toutes le lignes
  for (int i = row - 1; i >= 0; --i) {
    // parcours de la ligne
    for (int j = 0; j < col; ++j) {
      // on trouve le premier coef non nul
      if (matrix->mat[i][j] != 0) {
        break;
      }
      for (int k = 0; k < i; ++k) {
        if (matrix->mat[k][j] != 0) {
          matrix_add_row(matrix, k, i, -1);
        }
      }
    }
  }


  // ancienne version
  // annulation des début de ligne
  // for (int i = 0; i < row; i++) {
  //   for (int j = 0; j < i && j < col; j++) {
  //     if (matrix->mat[i][j] % ORDER == 1)
  //       matrix_add_row(matrix, i, j, -1);
  //     if (matrix->mat[i][j] % ORDER == 2)
  //       matrix_add_row(matrix, i, j, 1);
  //   }
  //   // si A[i][i] == 0
  //   if (i < col) {
  //     if (matrix->mat[i][i] % ORDER == 0) {
  //       for (int k = i + 1; k < row; ++k) {
  //         if (matrix->mat[k][i] % ORDER != 0) {
  //           matrix_exchange_row(matrix, i, k);
  //           if (i > 0)
  //             i -= 1;
  //         }
  //       }
  //     }
  //   }
  //   //
  //   for (int j = 0; j < col; ++j) {
  //     if (matrix->mat[i][j] != 0) {
  //       matrix_mul_row(matrix, i, inv_Fq(matrix->mat[i][i]));
  //       break;
  //     }
  //   }
  // }
  // for (int i = 0; i < row; i++) {
  //   if (i < col && matrix->mat[i][i] % ORDER == 0) {
  //     for (int k = i + 1; k < row; ++k) {
  //       if (matrix->mat[k][i] % ORDER != 0) {
  //         matrix_exchange_row(matrix, i, k);
  //       }
  //     }
  //   }
  // }
  // matrix_print(matrix,stdout);
  // rajouter les addition et soustraction de ligne pour diagonaliser en remontant
  // for (int i = row - 1; i >= 0; --i) {
  //
  // }
}

char matrix_det(matrix_t *A) {
  int size = A->nb_col;
  if (size == 2)
    return (add_Fq(mul_Fq(A->mat[0][0],A->mat[1][1]),-mul_Fq(A->mat[0][1],A->mat[1][0])));
  char det = 0;
  matrix_t *sub = NULL;
  for (int i = 0; i < size; ++i) {
    sub = matrix_sub(A,0,i);
    char sub_det = matrix_det(sub);
    char val1 = mul_Fq(pow(-1,i),A->mat[0][i]);
    det = add_Fq(det,mul_Fq(val1,sub_det));
    matrix_free(sub);
  }
  return det;
}

void matrix_print(matrix_t *matrix, FILE *fd) {
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
        if (val == 2 && ORDER == 3)
          fprintf(fd,"-1 ");
        else
          fprintf(fd," %d ", val);
      }
      fputs("\n", fd);
    }
    fputs("\n\n", fd);
  }
}
