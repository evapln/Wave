#include <lapacke.h>
#include "wave.h"

// matrix_t* phi (matrix_t* x,matrix_t* y,matrix_t* a,matrix_t* b, matrix_t* c, matrix_t* d) {
//   matrix_t *res = NULL;
//   res = matrix_alloc ();
//   return res;
// }


int main(int argc, char **argv) {
  // int col = 2;
  // int row = 3;
  // int k = 2;

  matrix_t *A = NULL;
  A = matrix_random(5,3);
  matrix_print(A, stdout);
  matrix_trigonalisation(A);
  matrix_print(A, stdout);
  // for (int i = 0; i < 5; i++){
  //   printf("%d \n", row_is_zero(A,i) );
  // }
  matrix_t *B = NULL;
  B = matrix_del_null_row(A);

  matrix_print(B, stdout);
  matrix_free(B);

  matrix_t *P = NULL;
  P = matrix_perm_random(10);
  matrix_print(P, stdout);
  matrix_free(P);
  // for (int i = 3; i < 7; ++i) {
  //   printf("row : %d\t col : %d\n",i,4);
  //   A = matrix_random(i,4);
  //   matrix_print(A, stdout);
  //   matrix_trigonalisation(A);
  //   matrix_print(A, stdout);
  //   matrix_free(A);
  // }
  // matrix_t *A = NULL;
  // A = matrix_alloc(row,col);
  // if (!A)
  //   return EXIT_FAILURE;
  // for (int i = 0; i < row; ++i)
  //   for (int j = 0; j < col; ++j) {
  //     char v = i+j;
  //     matrix_set_cell(A, i, j, v);
  //   }
  // puts("A :");
  // matrix_print(A,stdout);
  //
  // matrix_t *B = NULL;
  // B = matrix_random(3,3);
  // matrix_set_cell(B, 0, 0, (char)0);
  // matrix_set_cell(B, 0, 1, (char)1);
  // matrix_set_cell(B, 0, 2, (char)1);
  // matrix_set_cell(B, 1, 0, (char)1);
  // matrix_set_cell(B, 1, 1, (char)2);
  // matrix_set_cell(B, 1, 2, (char)1);
  // matrix_set_cell(B, 2, 0, (char)0);
  // matrix_set_cell(B, 2, 1, (char)0);
  // matrix_set_cell(B, 2, 2, (char)2);
  // puts("B :");
  // matrix_print(B,stdout);
  //
  // matrix_t *C = NULL;
  // C = matrix_alloc (row,k);
  // matrix_prod(C,A,B);
  // matrix_print(C,stdout);

  // matrix_t *sub = NULL;
  // sub = matrix_sub(B,1,1);
  // puts("sub:");
  // matrix_print(sub,stdout);

  // char det = matrix_det(B);
  // char det = (char)(-2)%3;
  // printf("det = %d \n", det);

  // matrix_t *C = NULL;
  // C = matrix_inv(B);
  // puts("B^-1 :");
  // matrix_print(C,stdout);
  //
  // matrix_t *id = NULL;
  // id = matrix_prod(B,C);
  // matrix_print(id,stdout);

  // matrix_t *inv = NULL;
  // inv = matrix_alloc(col,k);
  // matrix_inv(inv, B);
  // matrix_print(C,stdout);
  // matrix_free(inv);
  // matrix_free(C);
  // matrix_free(A);
  // matrix_free(B);
  // matrix_free(id);
  // matrix_free(C);
  return EXIT_SUCCESS;
}
