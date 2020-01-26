#include <lapacke.h>
#include "wave.h"

const int SIZE = 256;
const int OMEGA;
matrix_t* A = NULL;
matrix_t* B = NULL;
matrix_t* C = NULL;
matrix_t* D = NULL;
matrix_t* H = NULL;
int DIM;
static bool coef_init = false;

struct keys_t {
  sk_t *sk;
  matrix_t *pk;
};

struct sk_t {
  matrix_t *parite_U;
  matrix_t *parite_V;
  matrix_t *S;
  matrix_t *permut;
};

sk_t *sk_alloc(int dim_U, int dim_V, int dim) {
  sk_t *sk = malloc(sizeof(sk_t));
  if (!sk)
    return NULL;
  sk->parite_U = matrix_alloc ((SIZE/2)-dim_U, SIZE/2);
  if (!sk->parite_U){
    free(sk);
    return NULL;
  }
  sk->parite_V = matrix_alloc ((SIZE/2)-dim_V, SIZE/2);
  if (!sk->parite_V){
    matrix_free(sk->parite_U);
    free(sk);
    return NULL;
  }
  sk->S = matrix_alloc (SIZE-dim, SIZE-dim);
  if (!sk->S) {
    matrix_free(sk->parite_V);
    matrix_free(sk->parite_U);
    free(sk);
    return NULL;
  }
  sk->permut = matrix_alloc (SIZE,SIZE);
  if (!sk->permut) {
    matrix_free(sk->S);
    matrix_free(sk->parite_V);
    matrix_free(sk->parite_U);
    free(sk);
    return NULL;
  }
  return sk;
}

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

keys_t *key_alloc(int dim_U, int dim_V, int dim){
  keys_t *keys = malloc(sizeof(keys_t));
  if (!keys)
    return NULL;
  keys->sk = sk_alloc (dim_U, dim_V, dim);
  if (!keys->sk){
    free (keys);
    return NULL;
  }
  keys->pk = matrix_alloc (SIZE-dim, SIZE);
  if (!keys->pk){
    sk_free(keys->sk);
    free (keys);
    return NULL;
  }
  return keys;
}

void key_free (keys_t *keys) {
  if (keys != NULL)
  {
    matrix_free(keys->pk);
    sk_free(keys->sk);
  }
  free(keys);
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
  coef_init = true;
}

keys_t *key_gen (int lambda, int mode) {
  // initialise a,b,c,d :
  coeff_phi(mode);

  matrix_t *gen_U = NULL;
  matrix_t *gen_V = NULL;

  // crée les matrices génératrices de U et V et vérifie qu'elles sont ok:
  bool answer = false;
  matrix_t *gen_U_temp = NULL;
  while (!answer) {
    gen_U_temp = matrix_random(SIZE/2-2,SIZE/2);
    matrix_systematisation(gen_U_temp);
    gen_U = matrix_del_null_row(gen_U_temp);
    answer = matrix_is_syst(gen_U);
    if (!answer) {
      matrix_free(gen_U_temp);
      matrix_free(gen_U);
    }
  }
  matrix_free(gen_U_temp);

  answer = false;
  matrix_t *gen_V_temp = NULL;
  // matrix_t *gen_V = NULL;
  while (!answer) {
    gen_V_temp = matrix_random(SIZE/2-2,SIZE/2);
    matrix_systematisation(gen_V_temp);
    gen_V = matrix_del_null_row(gen_V_temp);
    answer = matrix_is_syst(gen_V);
    if (!answer) {
      matrix_free(gen_V_temp);
      matrix_free(gen_V);
    }
  }
  matrix_free(gen_V_temp);

  // dimension des codes
  int dim_U = matrix_get_row(gen_U);
  int dim_V = matrix_get_row(gen_V);
  DIM = dim_U + dim_V;

  // crée H_U H_V
  matrix_t *H_U = matrix_parite(gen_U);
  matrix_t *H_V = matrix_parite(gen_V);

  // alloue la structure clés et remplie avec les matrices de parités
  keys_t *keys = NULL;
  keys = key_alloc(dim_U,dim_V, DIM);
  matrix_copy2(keys->sk->parite_U,H_U);
  matrix_free(H_U);
  matrix_copy2(keys->sk->parite_V,H_V);
  matrix_free(H_V);

  // crée H
  H = parite(keys->sk->parite_U, keys->sk->parite_V);

  // S, P aléatoire
  matrix_t *S = matrix_random(SIZE-DIM,SIZE-DIM);
  matrix_copy2(keys->sk->S, S);
  matrix_free(S);
  matrix_t *permut = matrix_perm_random(SIZE);
  matrix_copy2(keys->sk->permut, permut);
  matrix_free(permut);

  // crée la clé publique SHP et la met dans la structure keys
  matrix_t *SH = matrix_prod(keys->sk->S,H);
  matrix_t *SHP = matrix_prod(SH,keys->sk->permut);
  matrix_copy2(keys->pk, SHP);

  // clean up
  matrix_free(SH);
  matrix_free(SHP);
  // matrix_free(H);
  matrix_free(gen_U);
  matrix_free(gen_V);
  return keys;
}

void infoset(int *info,const int n, const int len) {
  prng_init(time(NULL) + getpid());
  int i = 0;
  int r;
  while (i < len) {
    r = rand() % n;
    if (!is_in_array(info, i, r)) {
      info[i] = r;
      ++i;
    }
  }
}

matrix_t *prange_algebra(const matrix_t *parite, const matrix_t *syndrome, const int *info, const int len_i, const matrix_t *x) {
  // matrice de permutation envoyant info sur les DIM dernière coordonnées
  puts("A");
  matrix_t *P = matrix_perm_random_info(SIZE, info, len_i, DIM);
  puts("B");
  // (left | right) <- HP
  matrix_t *HP = matrix_prod(parite,P);
  puts("C");
  int row_HP = matrix_get_row(HP);
  matrix_t *left = matrix_alloc(row_HP,SIZE-DIM);
  matrix_t *right = matrix_alloc(row_HP,DIM);
  matrix_separate(HP, left, right);
  // (zero | ep) <- x
  matrix_t *zero = matrix_alloc(1, SIZE-DIM);
  matrix_t *ep = matrix_alloc(1,DIM);
  matrix_separate(x, zero, ep);
  // calcul de e
  matrix_t *right_T = matrix_trans(right);
  matrix_t *epright_T = matrix_prod(ep,right_T);
  matrix_t *epright_Tless = matrix_mul_by_scal(epright_T, -1);
  matrix_t *s = matrix_add(syndrome, epright_Tless);
  matrix_t *left_inv = matrix_inv(left);
  matrix_t *left_inv_T = matrix_trans(left_inv);
  matrix_t *couple_left = matrix_prod(s, left_inv_T);
  matrix_t *couple = matrix_concatenation(couple_left, ep, 0);
  matrix_t *P_T = matrix_trans(P);
  matrix_t *e = matrix_prod(couple, P_T);

  // clean up
  matrix_free(P);
  matrix_free(HP);
  matrix_free(left);
  matrix_free(right);
  matrix_free(zero);
  matrix_free(ep);
  matrix_free(right_T);
  matrix_free(epright_T);
  matrix_free(epright_Tless);
  matrix_free(s);
  matrix_free(left_inv);
  matrix_free(left_inv_T);
  matrix_free(couple_left);
  matrix_free(couple);
  matrix_free(P_T);
  // matric_free(P);
  return e;
}

matrix_t *iteration_prange(const matrix_t *parite, const matrix_t *syndrome) {
  puts("1");
  // initialisation
  prng_init(time(NULL) - getpid());
  puts("a");
  int t = rand() % DIM;
  puts("b");
  int col_H = matrix_get_col(parite);
  puts("c");
  int len_i = 0;
  puts("2");
  while(len_i == 0)
    len_i = rand() % col_H;
  // choix de l'ensemble d'information info de H
  int info[len_i];
  puts("3");
  infoset(info, col_H, len_i);
  puts("4");
  // tirage de x de poids t sur les coordonnées info
  matrix_t *x = vector_rand_weight(SIZE, info, len_i, t);
  puts("5");
  // recuperation de la solution
  matrix_t *e = prange_algebra(parite, syndrome, info, len_i, x);
  puts("6");
  //clean up
  matrix_free(x);
  puts("7");
  return e;
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

  // // test de sub_weight
  // matrix_t *vect = matrix_random(1,13);
  // for (int i = 0; i <13; ++i) printf(" %d ", i);
  // puts("");
  // matrix_print(vect, stdout);
  // int len_s = 6;
  // int subset[6] = {2,5,6,8,10,12};
  // printf("%d\n", sub_weight(vect, subset, len_s));
  // matrix_free(vect);


  // // test de infoset
  // int len_i = 10;
  // int info[len_i];
  // infoset(info, 20, len_i);
  // for (int i = 0; i < len_i; ++i)
  //   printf("%d ", info[i]);
  // puts("");
  //
  //
  // // test de vector_rand_weight
  // int t = 7;
  // matrix_t *vect = vector_rand_weight(20, info, len_i, t);
  // for (int i = 0; i <20; ++i) printf(" %d ", i);
  // puts("");
  // matrix_print(vect, stdout);
  // matrix_free(vect);

  // // test matrix_perm_random_info et matrix_perm_random
  // matrix_t *matrix = matrix_perm_random(20);
  // puts("permutation normale :");
  // matrix_print(matrix, stdout);
  //
  // int tab[5] = {0,1,2,3,4};
  // DIM = 7;
  // matrix_t *matr = matrix_perm_random_info(20, tab, 5, DIM);
  // puts("permutation avec les 5 premières coordonnées à la fin :");
  // matrix_print(matr, stdout);
  //
  // matrix_free(matrix);
  // matrix_free(matr);

  //// test matrix_separate
  // matrix_t *matrix = matrix_random(10,17);
  // matrix_t *left = matrix_alloc(10,6);
  // matrix_t *right = matrix_alloc(10,11);
  // matrix_separate(matrix, left, right);
  // matrix_print(matrix, stdout);
  // matrix_print(left, stdout);
  // matrix_print(right, stdout);
  // matrix_free(matrix);
  // matrix_free(left);
  // matrix_free(right);

  //// test key_gen
  // keys_t *keys = key_gen(0,0);
  // matrix_print(keys->pk,stdout);
  // key_free(keys);

  // test iteration_prange
  keys_t *keys = key_gen(0,0);
  printf("%d\n", DIM);
  matrix_t *e = matrix_random(1,SIZE);
  matrix_t *synd = syndrome(e, H);
  matrix_t *ep = iteration_prange(H, synd);
  puts("e");
  matrix_print(e, stdout);
  puts("ep");
  matrix_print(ep, stdout);
  matrix_free(e);
  matrix_free(synd);
  matrix_free(ep);
  key_free(keys);

  // clean up
  matrix_free(A);
  matrix_free(B);
  matrix_free(C);
  matrix_free(D);
  matrix_free(H);

  return EXIT_SUCCESS;
}
