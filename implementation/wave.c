#include <lapacke.h>
#include "wave.h"
// #include "matrix.c"

const int SIZE = 16;
const int DIM = 10;
const int OMEGA;
matrix_t* A = NULL;
matrix_t* B = NULL;
matrix_t* C = NULL;
matrix_t* D = NULL;
matrix_t *gen_U = NULL;
matrix_t *gen_V = NULL;

static bool coef_init = false;

struct sk_t {
  matrix_t *parite_U;
  matrix_t *parite_V;
  matrix_t *S;
  matrix_t *permut;
};

sk_t *sk_alloc(void) {
  // allocation sk
  sk_t *sk = malloc(sizeof(sk_t));
  if (!sk)
    return NULL;
  // allocation matrive de parité U
  sk->parite_U = malloc(sizeof(matrix_t));
  if (!sk->parite_U) {
    free(sk);
    return NULL;
  }
  // sk->parite_U = par_U;
  // allocation matrice de parité V
  sk->parite_V = malloc(sizeof(matrix_t));
  if (!sk->parite_V) {
    free(sk->parite_U);
    free(sk);
    return NULL;
  }
  // sk->parite_V = par_V;
  // allocation matrice S
  sk->S = malloc(sizeof(matrix_t));
  if (!sk->S) {
    free(sk->parite_V);
    free(sk->parite_U);
    free(sk);
    return NULL;
  }
  // sk->S = S;
  // allocation matrice de permutation P
  sk->permut = matrix_alloc(SIZE,SIZE);
  if (!sk->permut) {
    free(sk->S);
    free(sk->parite_V);
    free(sk->parite_U);
    free(sk);
    return NULL;
  }
  // sk->permut = permut;
  // retour
  return sk;
}

// sk_t *sk_alloc(int par_U_row, int par_U_col, int par_V_row, int par_V_col, int S_size) {
//   // allocation sk
//   sk_t *sk = malloc(sizeof(sk_t));
//   if (!sk)
//     return NULL;
//   // allocation matrive de parité U
//   matrix_t *par_U = matrix_alloc(par_U_row, par_U_col);
//   if (!par_U) {
//     free(sk);
//     return NULL;
//   }
//   sk->parite_U = par_U;
//   // allocation matrice de parité V
//   matrix_t *par_V = matrix_alloc(par_V_row, par_V_col);
//   if (!par_V) {
//     matrix_free(sk->parite_U);
//     free(sk);
//     return NULL;
//   }
//   sk->parite_V = par_V;
//   // allocation matrice S
//   matrix_t *S = matrix_alloc(S_size, S_size);
//   if (!S) {
//     matrix_free(sk->parite_V);
//     matrix_free(sk->parite_U);
//     free(sk);
//     return NULL;
//   }
//   sk-> S = S;
//   // allocation matrice de permutation P
//   matrix_t *permut = matrix_alloc(SIZE,SIZE);
//   if (!permut) {
//     matrix_free(sk->S);
//     matrix_free(sk->parite_V);
//     matrix_free(sk->parite_U);
//     free(sk);
//     return NULL;
//   }
//   sk->permut = permut;
//   // retour
//   return sk;
// }

// sk_t *sk_alloc(void) {
//   sk_t *sk = malloc(sizeof(sk));
//   return sk;
// }

void sk_free (sk_t *sk) {
  if (sk != NULL)
  {
    matrix_free(sk->parite_U);
    matrix_free(sk->parite_V);
    matrix_free(sk->S);
    matrix_free(sk->permut);
  }
  free(sk);
}

matrix_t* phi (const matrix_t* x,const matrix_t* y) {
  if (!coef_init)
    return NULL;
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
  char a, b, c, d;
  if (mode == 0) {
    A = matrix_init(1,SIZE/2,(char)1);
    B = matrix_init(1,SIZE/2,(char)0);
    C = matrix_init(1,SIZE/2,(char)1);
    D = matrix_init(1,SIZE/2,(char)1);
  }
  else {
      A = matrix_alloc(1,SIZE/2);
      B = matrix_alloc(1,SIZE/2);
      C = matrix_alloc(1,SIZE/2);
      D = matrix_alloc(1,SIZE/2);
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
  coef_init = true;
}

void key_gen (int lambda, matrix_t *pk, sk_t *sk, int mode) {
  // initialise a,b,c,d :
  coeff_phi (mode);

  // crée les matrices génératrices de U et V et vérifie qu'elles sont ok:
  bool answer = false;
  matrix_t *gen_U_temp = NULL;
  // matrix_t *gen_U = NULL;
  while (!answer){
    gen_U_temp = matrix_random(DIM/2,SIZE/2);
    matrix_systematisation(gen_U_temp);
    gen_U = matrix_del_null_row(gen_U_temp);
    answer = matrix_is_syst(gen_U);
  }
  matrix_free (gen_U_temp);

  answer = false;
  matrix_t *gen_V_temp = NULL;
  // matrix_t *gen_V = NULL;
  while (!answer){
    gen_V_temp = matrix_random(DIM/2,SIZE/2);
    matrix_systematisation(gen_V_temp);
    gen_V = matrix_del_null_row(gen_V_temp);
    answer = matrix_is_syst(gen_V);
  }
  matrix_free (gen_V_temp);

  // crée H_U H_V
  matrix_t *parite_U = matrix_parite(gen_U);
  matrix_t *parite_V = matrix_parite(gen_V);

  // dimension des codes
  int dim_U = matrix_get_col(parite_U) - matrix_get_row(parite_U);
  int dim_V = matrix_get_col(parite_V) - matrix_get_row(parite_V);
  int dim_parite = dim_U + dim_V;

  // crée H
  matrix_t *H = parite(parite_U, parite_V);

  // S, P aléatoire
  matrix_t *S = matrix_random(SIZE-dim_parite,SIZE-dim_parite);
  matrix_t *P = matrix_perm_random(SIZE);

  // pk = SHP
  matrix_t *SH = matrix_prod(S,H);
  pk = matrix_prod(SH,P);

  // sk = H_U,H_V,S,P
  sk = sk_alloc();
  puts("1");
  sk->parite_U = parite_U;
  puts("2");
  sk->parite_V = parite_V;
  puts("3");
  sk->S = S;
  puts("4");
  sk->permut = P;
  puts("5");


  // clean up
  matrix_free(SH);
  matrix_free(H);

  // matrix_free(gen_U);
  // matrix_free(gen_V);
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
  //
  // matrix_t *A = NULL;
  // A = matrix_random(4,6);
  // matrix_print(A, stdout);
  // matrix_systematisation(A);
  // matrix_print(A, stdout);
  // if (matrix_is_syst(A))
  //   puts("A syst");
  // matrix_t *B = NULL;
  // B = matrix_del_null_row(A);
  // matrix_print(B, stdout);
  // if (matrix_is_syst(B))
  //   puts("B syst");
  // matrix_free(A);
  // matrix_free(B);

  // matrix_t *A = NULL;
  // A = matrix_alloc(4,6);
  // matrix_set_cell(A,0,0,1);
  // matrix_set_cell(A,0,1,1);
  // matrix_set_cell(A,0,2,1);
  // matrix_set_cell(A,0,3,0);
  // matrix_set_cell(A,0,4,1);
  // matrix_set_cell(A,0,5,1);
  // matrix_set_cell(A,1,0,0);
  // matrix_set_cell(A,1,1,1);
  // matrix_set_cell(A,1,2,0);
  // matrix_set_cell(A,1,3,0);
  // matrix_set_cell(A,1,4,1);
  // matrix_set_cell(A,1,5,1);
  // matrix_set_cell(A,2,0,1);
  // matrix_set_cell(A,2,1,0);
  // matrix_set_cell(A,2,2,1);
  // matrix_set_cell(A,2,3,1);
  // matrix_set_cell(A,2,4,0);
  // matrix_set_cell(A,2,5,1);
  // matrix_set_cell(A,3,0,0);
  // matrix_set_cell(A,3,1,1);
  // matrix_set_cell(A,3,2,1);
  // matrix_set_cell(A,3,3,1);
  // matrix_set_cell(A,3,4,0);
  // matrix_set_cell(A,3,5,1);
  // matrix_print(A,stdout);
  // matrix_systematisation(A);
  // matrix_print(A,stdout);
  // // matrix_free(A);
  // matrix_t *B = NULL;
  // B = matrix_parite(A);
  // if (!B)
  //   puts("error");
  // else
  //   matrix_print(B, stdout);
  // matrix_free(B);

  // matrix_t *pk = NULL;
  sk_t *sk = NULL;
  sk = sk_alloc();
  sk->parite_U = matrix_random(3,4);
  sk->parite_V = matrix_random(3,4);
  sk->S = matrix_random(3,4);
  sk->permut = matrix_random(3,4);
  // key_gen(3,pk,sk,1);
  // matrix_free(pk);
  sk_free(sk);

  // matrix_free(A);
  // matrix_free(B);
  // matrix_free(C);
  // matrix_free(D);

  // matrix_free(X);
  // matrix_free(Y);
  // matrix_free(pariteUV);
  return EXIT_SUCCESS;
}
