#include <lapacke.h>
#include "wave.h"



int main(int argc, char **argv) {
  int col = 2;
  int row = 3;
  matrix_t *matrix = NULL;
  matrix = matrix_alloc(row,col);
  if (!matrix)
    return EXIT_FAILURE;
  for (int i = 0; i < row; ++i)
    for (int j = 0; j < col; ++j) {
      short v = i+j;
      matrix_set_cell(matrix, i, j, v);
    }
  puts("A :");
  matrix_print(matrix,stdout);
  matrix_t *sum = NULL;
  sum = matrix_copy(matrix);
  matrix_scal(sum, matrix, matrix);
  puts("sum :");
  matrix_print(sum,stdout);
  matrix_free(matrix);
  matrix_free(sum);
  return EXIT_SUCCESS;
}
