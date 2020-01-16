#include <lapacke.h>
#include "wave.h"



int main(int argc, char **argv) {
  int col = 2;
  int row = 3;
  int k = 2;

  matrix_t *A = NULL;
  A = matrix_alloc(row,col);
  if (!A)
    return EXIT_FAILURE;
  for (int i = 0; i < row; ++i)
    for (int j = 0; j < col; ++j) {
      short v = i+j;
      matrix_set_cell(A, i, j, v);
    }
  puts("A :");
  matrix_print(A,stdout);

  matrix_t *B = NULL;
  B = matrix_alloc(col,k);
  matrix_set_cell(B, 0, 0, 3);
  matrix_set_cell(B, 0, 1, 2);
  matrix_set_cell(B, 1, 0, 4);
  matrix_set_cell(B, 1, 1, 1);
  matrix_print(B,stdout);

  matrix_t *C = NULL;
  C = matrix_alloc (row,k);
  matrix_prod(C,A,B);
  matrix_print(C,stdout);
  matrix_free(C);
  matrix_free(A);
  matrix_free(B);
  return EXIT_SUCCESS;
}
