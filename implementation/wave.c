#include <lapacke.h>
#include "wave.h"

const int SIZE = 16;
const int DIM = 10;
const int OMEGA;
matrix_t* A = NULL;
matrix_t* B = NULL;
matrix_t* C = NULL;
matrix_t* D = NULL;

struct sk_t {
  matrix_t *parite_U;
  matrix_t *parite_V;
  matrix_t *S;
  matrix_t *permut;
};

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

void coeff_phi (int mode) {
  char a; char b; char c; char d;
  if (mode == 0) {
    A = matrix_init (1,SIZE/2,(char)1);
    B = matrix_init (1,SIZE/2,(char)0);
    C = matrix_init (1,SIZE/2,(char)1);
    D = matrix_init (1,SIZE/2,(char)1);
  }
  else {
      A = matrix_alloc (1,SIZE/2);
      B = matrix_alloc (1,SIZE/2);
      C = matrix_alloc (1,SIZE/2);
      D = matrix_alloc (1,SIZE/2);
      for (int i = 0; i < SIZE/2; i++) {
        a = 0;
        b = 0;
        c = 0;
        d = 0;
        while (add_Fq(mul_Fq(a,d),mul_Fq(-b,c)) == 0) {
          a = 0;
          while (a == 0)
            a = rand_Fq();
          b = rand_Fq();
          c = 0;
          while (c == 0)
            c = rand_Fq();
          d = rand_Fq();
        }
        matrix_set_cell(A,0,i,a);
        matrix_set_cell(B,0,i,b);
        matrix_set_cell(C,0,i,c);
        matrix_set_cell(D,0,i,d);
      }
  }
  puts("A : \n");
  matrix_print(A, stdout);
  puts("B : \n");
  matrix_print(B, stdout);
  puts("C : \n");
  matrix_print(C, stdout);
  puts("D : \n");
  matrix_print(D, stdout);
  return;
}

void key_gen (int lambda, matrix_t *pk, sk_t *sk, int mode) {
  // initialise a,b,c,d :
  coeff_phi (mode);

  // crée les matrices génératrices de U et V et vérifie qu'elles sont ok:
  bool answer = 0;
  matrix_t *gen_U_temp = NULL;
  matrix_t *gen_U = NULL;
  while (answer == 0){
    gen_U_temp = matrix_random(DIM/2,SIZE/2);
    matrix_systematisation(gen_U_temp);
    gen_U = matrix_del_null_row(gen_U_temp);
    answer = matrix_is_syst(gen_U);
    if (answer == 0)
      matrix_free (gen_U);
  }

  answer = 0;
  matrix_t *gen_V_temp = NULL;
  matrix_t *gen_V = NULL;
  while (answer == 0){
    gen_V_temp = matrix_random(DIM/2,SIZE/2);
    matrix_systematisation(gen_V_temp);
    gen_V = matrix_del_null_row(gen_V_temp);
    answer = matrix_is_syst(gen_V);
    if (answer == 0)
      matrix_free (gen_U);
    }

  // crée H_U H_V
  // crée H
  // S, P aléatoire
  // pk = SHP
  // sk = H_U,H_V,S,P
  matrix_free(gen_U);
  matrix_free(gen_V);
  return;
}

int main(int argc, char **argv) {


  // matrix_t *X = NULL;
  // X = matrix_random(3,4);
  // puts("X : \n");
  // matrix_print(X, stdout);
  // matrix_t *Y = NULL;
  // Y = matrix_random(2,4);
  // puts("Y : \n");
  // matrix_print(Y, stdout);
  // matrix_t *pariteUV = NULL;
  // pariteUV = parite (X,Y);
  // puts("parite : \n");
  // matrix_print(pariteUV, stdout);
  matrix_t *pk = NULL;
  sk_t *sk = NULL;
  key_gen(3,pk,sk,0);

  matrix_free(A);
  matrix_free(B);
  matrix_free(C);
  matrix_free(D);
  // matrix_free(X);
  // matrix_free(Y);
  // matrix_free(pariteUV);
  return EXIT_SUCCESS;
}
