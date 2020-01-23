#include <lapacke.h>
#include "wave.h"

const int SIZE;
const int OMEGA;
matrix_t* A = NULL;
matrix_t* B = NULL;
matrix_t* C = NULL;
matrix_t* D = NULL;

matrix_t* phi (const matrix_t* x,const matrix_t* y) {
  matrix_t *res = NULL;
  matrix_t *res_g = NULL;
  matrix_t *res_d = NULL;
  matrix_t *ax = NULL;
  matrix_t *by = NULL;
  matrix_t *cx = NULL;
  matrix_t *dy = NULL;
  ax = vect_scal(A,x);
  by = vect_scal(B,y);
  cx = vect_scal(C,x);
  dy = vect_scal(D,y);
  res_g = matrix_add(ax,by);
  res_d = matrix_add(cx,dy);
  res = matrix_concatenation (res_g, res_d, 0);
  matrix_free(ax);
  matrix_free(by);
  matrix_free(cx);
  matrix_free(dy);
  matrix_free(res_d);
  matrix_free(res_g);
  return res;
}

matrix_t* syndrome (const matrix_t *e, const matrix_t *parite){
  matrix_t *trans = NULL;
  trans = matrix_trans(parite);
  matrix_t *res = NULL;
  res = matrix_prod (e, trans);
  matrix_free(trans);
  return res;
}

matrix_t* parite (const matrix_t* parite_U,const matrix_t* parite_V) {
  matrix_t *parite_UV = NULL;
  matrix_t *A_diag = NULL;
  matrix_t *B_diag = NULL;
  matrix_t *C_diag = NULL;
  matrix_t *D_diag = NULL;
  A_diag = matrix_vect_to_diag (A,(char)1);
  B_diag = matrix_vect_to_diag (B,opp_Fq());
  C_diag = matrix_vect_to_diag (C,opp_Fq());
  D_diag = matrix_vect_to_diag (D,(char)1);
  matrix_t *HU_D = NULL;
  matrix_t *HU_B = NULL;
  matrix_t *HV_C = NULL;
  matrix_t *HV_A = NULL;
  HU_D = matrix_prod (parite_U, D_diag);
  HU_B = matrix_prod (parite_U, B_diag);
  HV_C = matrix_prod (parite_V, C_diag);
  HV_A = matrix_prod (parite_V, A_diag);
  matrix_t *haut = NULL;
  matrix_t *bas = NULL;
  haut = matrix_concatenation (HU_D, HU_B, 0);
  bas = matrix_concatenation (HV_C, HV_A, 0);
  parite_UV = matrix_concatenation (haut,bas, 1);
  matrix_free (haut);
  matrix_free (bas);
  matrix_free (HU_D);
  matrix_free (HU_B);
  matrix_free (HV_C);
  matrix_free (HV_A);
  matrix_free (A_diag);
  matrix_free (B_diag);
  matrix_free (C_diag);
  matrix_free (D_diag);
  return parite_UV;
}

int main(int argc, char **argv) {

  // NE PAS TOUCHER !! INIT DES vecteurs a,b,c,d de phi //
  A = matrix_init (1,4,(char)1);
  puts("A : \n");
  matrix_print(A, stdout);
  B = matrix_init (1,4,(char)0);
  puts("B : \n");
  matrix_print(B, stdout);
  C = matrix_init (1,4,(char)1);
  puts("C : \n");
  matrix_print(C, stdout);
  D = matrix_init (1,4,(char)1);
  puts("D : \n");
  matrix_print(D, stdout);
  // NE PAS TOUCHER !! INIT DES vecteurs a,b,c,d de phi //


  matrix_t *X = NULL;
  X = matrix_random(3,4);
  puts("X : \n");
  matrix_print(X, stdout);
  matrix_t *Y = NULL;
  Y = matrix_random(2,4);
  puts("Y : \n");
  matrix_print(Y, stdout);
  matrix_t *pariteUV = NULL;
  pariteUV = parite (X,Y);
  puts("parite : \n");
  matrix_print(pariteUV, stdout);



  matrix_free(A);
  matrix_free(B);
  matrix_free(C);
  matrix_free(D);
  matrix_free(X);
  matrix_free(Y);
  matrix_free(pariteUV);
  return EXIT_SUCCESS;
}
