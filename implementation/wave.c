#include <lapacke.h>
#include "wave.h"

const float LAMBDA = 3.08;
const int SIZE = 26;
const int OMEGA = 0.9261*SIZE;
const int d = 0;
const int K_U = 0.7978*SIZE/2;
const int K_V = 0.4201*SIZE/2;
// const float LAMBDA = 3.08;
// const int SIZE = 16;
// const int OMEGA = 13;
// const int d = 0;
// const int K_U = SIZE/4;
// const int K_V = SIZE/4;
matrix_t* A = NULL;
matrix_t* B = NULL;
matrix_t* C = NULL;
matrix_t* D = NULL;
matrix_t* H = NULL;
matrix_t *gen_U = NULL;
matrix_t *gen_V = NULL;
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
  int dim_U;
  int dim_V;
};


////////////////////////////////////////////////////////////////////////////////
//////////////////////////// gestion de la mémoire /////////////////////////////
////////////////////////////////////////////////////////////////////////////////

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
  if (sk) {
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
  keys->sk = sk_alloc(dim_U, dim_V, dim);
  if (!keys->sk){
    free (keys);
    return NULL;
  }
  keys->pk = matrix_alloc(SIZE-dim, SIZE);
  if (!keys->pk){
    sk_free(keys->sk);
    free (keys);
    return NULL;
  }
  return keys;
}

void key_free (keys_t *keys) {
  if (keys) {
    matrix_free(keys->pk);
    sk_free(keys->sk);
  }
  free(keys);
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////// fonctions secondaires //////////////////////////////
////////////////////////////////////////////////////////////////////////////////

matrix_t* phi (const matrix_t* x,const matrix_t* y) {
  if (!coef_init || !x || !y)
    return NULL;
  matrix_t *ax = vect_scal(A,x);
  matrix_t *by = vect_scal(B,y);
  matrix_t *cx = vect_scal(C,x);
  matrix_t *dy = vect_scal(D,y);
  if (!ax || !by || !cx || !dy) {
    matrix_free(ax);
    matrix_free(by);
    matrix_free(cx);
    matrix_free(dy);
    return NULL;
  }
  matrix_t *res_g = matrix_add(ax,by);
  matrix_t *res_d = matrix_add(cx,dy);
  if (!res_g || !res_d) {
    matrix_free(ax);
    matrix_free(by);
    matrix_free(cx);
    matrix_free(dy);
    matrix_free(res_d);
    matrix_free(res_g);
    return NULL;
  }
  matrix_t *res = matrix_concatenation (res_g, res_d, 0);
  matrix_free(ax);
  matrix_free(by);
  matrix_free(cx);
  matrix_free(dy);
  matrix_free(res_d);
  matrix_free(res_g);
  return res;
}

matrix_t* syndrome (const matrix_t *e, const matrix_t *parite){
  if (!e || !parite)
    return NULL;
  matrix_t *trans = matrix_trans(parite);
  if (!trans)
    return NULL;
  matrix_t *res = matrix_prod (e, trans);
  matrix_free(trans);
  return res;
}

matrix_t* parite (const matrix_t* parite_U,const matrix_t* parite_V) {
  if (!parite_U || !parite_V)
    return NULL;
  matrix_t *A_diag = matrix_vect_to_diag (A,(char)1);
  matrix_t *B_diag = matrix_vect_to_diag (B,opp_Fq());
  matrix_t *C_diag = matrix_vect_to_diag (C,opp_Fq());
  matrix_t *D_diag = matrix_vect_to_diag (D,(char)1);
  if (!A_diag || !B_diag || !C_diag || !D_diag) {
    matrix_free(A_diag);
    matrix_free(B_diag);
    matrix_free(C_diag);
    matrix_free(D_diag);
    return NULL;
  }
  matrix_t *HU_D = matrix_prod (parite_U, D_diag);
  matrix_t *HU_B = matrix_prod (parite_U, B_diag);
  matrix_t *HV_C = matrix_prod (parite_V, C_diag);
  matrix_t *HV_A = matrix_prod (parite_V, A_diag);
  matrix_free(A_diag);
  matrix_free(B_diag);
  matrix_free(C_diag);
  matrix_free(D_diag);
  if (!HU_D || !HU_B || !HV_C || !HV_A) {
    matrix_free(HU_D);
    matrix_free(HU_B);
    matrix_free(HV_C);
    matrix_free(HV_A);
    return NULL;
  }
  matrix_t *haut = matrix_concatenation (HU_D, HU_B, 0);
  matrix_t *bas = matrix_concatenation (HV_C, HV_A, 0);
  matrix_free(HU_D);
  matrix_free(HU_B);
  matrix_free(HV_C);
  matrix_free(HV_A);
  if (!haut || !bas) {
    matrix_free(haut);
    matrix_free(bas);
    return NULL;
  }
  matrix_t *parite_UV = matrix_concatenation (haut,bas, 1);
  matrix_free(haut);
  matrix_free(bas);
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
        while (add_Fq(mul_Fq(a,d),mul_Fq(-b,c)) != 1) {
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
  if (A && B && C && D)
    coef_init = true;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////// génération des clés /////////////////////////////
////////////////////////////////////////////////////////////////////////////////

keys_t *key_gen (int mode) {
  // initialise a,b,c,d :
  puts("\tchoix des coefficients de la fonction phi...");
  coeff_phi(mode);
  if (!coef_init)
    return NULL;

  // crée les matrices génératrices de U et V et vérifie qu'elles sont ok:
  puts("\tchoix des matrices génératrices des codes U et V...");
  bool answer = false;
  matrix_t *gen_U_temp = NULL;
  while (!answer) {
    matrix_free(gen_U_temp);
    matrix_free(gen_U);
    // on impose la dimensionde H_V à K_V
    gen_U_temp = matrix_random(K_U,SIZE/2);
    if (!gen_U_temp)
      continue;
    matrix_systematisation(gen_U_temp);
    gen_U = matrix_del_null_row(gen_U_temp);
    if (!gen_U) {
      matrix_free(gen_U_temp);
      continue;
    }
    answer = matrix_is_syst(gen_U);
  }
  matrix_free(gen_U_temp);

  answer = false;
  matrix_t *gen_V_temp = NULL;
  // matrix_t *gen_V = NULL;
  while (!answer) {
    matrix_free(gen_V_temp);
    matrix_free(gen_V);
    // on impose la dimensionde H_V à K_V
    gen_V_temp = matrix_random(K_V,SIZE/2);
    if (!gen_V_temp)
      continue;
    // puts("gen_V_temp");
    // matrix_print(gen_V_temp, stdout);
    matrix_systematisation(gen_V_temp);
    // puts("gen_V_temp"); matrix_print(gen_V_temp, stdout);
    // puts("gen_V_temp syst");
    // matrix_print(gen_V_temp, stdout);
    gen_V = matrix_del_null_row(gen_V_temp);
    // puts("gen_V"); matrix_print(gen_V, stdout);
    if (!gen_V) {
      matrix_free(gen_V_temp);
      continue;
    }
    // puts("gen_V");
    // matrix_print(gen_V, stdout)
    answer = matrix_is_syst(gen_V);
  }
  matrix_free(gen_V_temp);

  // dimension des codes
  int dim_U = matrix_get_row(gen_U);
  int dim_V = matrix_get_row(gen_V);
  // printf("dim_U : %d\n dim_V : %d\n", dim_U, dim_V);
  // printf("\tc'est bon on les a !\tla dimension de U est %d et celle de V est %d\n", dim_U, dim_V);
  DIM = dim_U + dim_V;

  // crée H_U H_V
  puts("\ton crée les matrices de parités associées...");
  matrix_t *H_U = matrix_parite(gen_U);
  matrix_t *H_V = matrix_parite(gen_V);
  // puts("parité V"); matrix_print(H_V,stdout);
  printf("\tHu : %d x %d\n\tHv : %d x %d\n", matrix_get_row(H_U),matrix_get_col(H_U),matrix_get_row(H_V),matrix_get_col(H_V));
  if (!H_U || !H_V) {
    matrix_free(H_U);
    matrix_free(H_V);
    return NULL;
  }

  // alloue la structure clés et remplie avec les matrices de parités
  keys_t *keys = key_alloc(dim_U,dim_V, DIM);
  if (!keys) {
    matrix_free(H_U);
    matrix_free(H_V);
    return NULL;
  }
  matrix_copy2(keys->sk->parite_U,H_U);
  matrix_copy2(keys->sk->parite_V,H_V);
  matrix_free(H_U);
  matrix_free(H_V);
  if (!keys->sk->parite_U ||! keys->sk->parite_V) {
    key_free(keys);
    return NULL;
  }

  // crée H
  puts("\tpuis la matrice de parité du code (U,U+V)...");
  H = parite(keys->sk->parite_U, keys->sk->parite_V);
  if (!H) {
    key_free(keys);
    return NULL;
  }
  // S, P aléatoire
  matrix_t *S = matrix_random(SIZE-DIM,SIZE-DIM);
  if (!S) {
    key_free(keys);
    return NULL;
  }
  matrix_copy2(keys->sk->S, S);
  matrix_free(S);
  if (!keys->sk->S) {
    key_free(keys);
    return NULL;
  }
  matrix_t *permut = matrix_perm_random(SIZE);
  if (!permut) {
    key_free(keys);
    return NULL;
  }
  matrix_copy2(keys->sk->permut, permut);
  matrix_free(permut);
  if (!keys->sk->permut) {
    key_free(keys);
    return NULL;
  }
  keys->sk->dim_U = dim_U;
  keys->sk->dim_V = dim_V;

  // crée la clé publique SHP et la met dans la structure keys
  puts("\tet enfin notre clé publique...");
  // printf("H : %d x %d\nS : %d x %d\n", matrix_get_row(H),matrix_get_col(H),matrix_get_row(keys->sk->S),matrix_get_col(keys->sk->S));
  matrix_t *SH = matrix_prod(keys->sk->S,H);
  if (!SH) {
    key_free(keys);
    return NULL;
  }
  // puts("a");
  matrix_t *SHP = matrix_prod(SH,keys->sk->permut);
  matrix_free(SH);
  if (!SHP) {
    key_free(keys);
    return NULL;
  }
  // puts("b");
  matrix_copy2(keys->pk, SHP);
  matrix_free(SHP);
  if (!keys->pk) {
    key_free(keys);
    return NULL;
  }
  // puts("c");
  return keys;
}

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// décodage /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void decode_ev(matrix_t * ev,const matrix_t *G, const matrix_t *synd, keys_t *keys) {
  if (!ev || !G || !synd)
    return;
  // printf("\tsyndrome cherché : "); matrix_print(synd,stdout);
  int col = matrix_get_col(G);
  if (matrix_get_col(ev) != col || matrix_get_row(ev) != 1)
    return;
  random_word(ev, G);
  // matrix_t *verif = syndrome(ev, keys->sk->parite_V);
  // printf("\tdoit être à 0 : "); matrix_print(verif, stdout);
  // matrix_free(verif);
  matrix_t *s = matrix_alloc(1,col);
  int col_synd = matrix_get_col(synd);
  // printf("\tlongueur G : %d\n\tlongueur syndrome cherché : %d\n",col, col_synd);
  for (int i = 0; i < col - col_synd; ++i)
    matrix_set_cell(s, 0, i, 0);
  // printf("\tsynd cherché : "); matrix_print(synd, stdout);
  for (int i = col - col_synd; i < col; ++i) {
    // printf("%d: %d\n", i - col + col_synd, matrix_get_cell(synd, 0, i - col + col_synd));
    matrix_set_cell(s, 0, i, matrix_get_cell(synd, 0, i - col + col_synd));
  }
  // printf("\ts : "); matrix_print(s,stdout);
  matrix_add_modified(ev, s, 1);
  matrix_free(s);
}

// matrix_t *resol_syst(const matrix_t *ev, int *ind) {
//   if (!ev)
//     return NULL;
//   int n = A->nb_row;
//   int m = A ->nb_col;
//   if (y->nb_col != 1 || y->nb_row != m)
//     return NULL;
// }

void decode_eu(matrix_t * eu, const keys_t *keys, const matrix_t *synd,
               const matrix_t *ev) {
  // Vérification des entrées
  if (!eu || !synd || !ev)
    return;
  // int col = matrix_get_col(G);
  if (matrix_get_row(ev) != 1 || matrix_get_row(eu) != 1)
    return;
  // initialistaions
  int dim_U = keys->sk->dim_U;
  int a, b, c, d, w, r;
  char no1, no2, evi, eui, eui_not;
  matrix_t *eu_int = NULL;
  matrix_t *x = NULL;
  matrix_t *e = NULL;
  matrix_t * eu_not = NULL;
  // matrix_print(eu, stdout);
  // création de I ensemble d'information de H : [1,k] k la dimension de H
  // int ens_I[DIM];
  // for (int i = 0; i < DIM; ++i)
  //   ens_I[i] = i;

  int ens[SIZE/2];
  for (int i = 0; i < SIZE/2; ++i)
    ens[i] = i;
  do {
    matrix_free(eu_int);
    // création de I aléatoirement
    shuffle(ens,SIZE/2);
    int len_I = rand() % (SIZE/2 - dim_U) + dim_U; // aléa entre dim_U et SIZE/2
    // printf("\tlongueur de I : %d\n",len_I);
    int ens_I[len_I];
    for (int i = 0; i < len_I; ++i)
      ens_I[i] = ens[i];
    // printf("I = {");
    // for (int i = 0; i < len_I; ++i)
    //   printf(" %d", ens_I[i]);
    // printf(" }\n");
    // if(eu_int)
    //   puts("existe");
    // while (!eu_int) {
    do {
      // puts("ah");
      eu_not = matrix_alloc(1,SIZE/2);
      if (!eu_not)
        return;
      for (int i = 0; i < SIZE/2; ++i) {
        matrix_set_cell(eu, 0, i, '*');
        matrix_set_cell(eu_not, 0, i, '*');
      }
      // création du tableau J contenant les k_u positions choisies parmi I
      shuffle(ens_I, len_I);
      int ens_J[dim_U];
      for (int i = 0; i < dim_U; ++i)
        ens_J[i] = ens_I[i];
      //   ens_J[i] = i;
      // tableau random
      // int r;
      // prng_init(time(NULL) + getpid());
      // int i = 0;
      // while (i < dim_U) {
      //   r = rand() % col;
      //   if (!is_in_array(ind, dim_U, r)) {
      //     ind[i] = r;
      //     ++i;
      //   }
      // }
      // printf("J = {");
      // for (int i = 0; i < dim_U; ++i)
      //   printf(" %d", ens_J[i]);
      // printf(" }\n");

      //// reste à résoudre le système linéaire

      // sur les ku positions choisies (dans J), on fixe eu de telle sorte que
      // A(i)*eu(i) + B(i)*ev(i) != 0 et  C(i)*eu(i) + D(i)*ev(i) != 0
      // quand on a une solution exacte, on met dans eu,
      // quand on a juste une impossibilité, on met dans eu_not
      for (int i = 0; i < dim_U; ++i) {
        // printf("abcd %d\n", ens_J[i]);
        a = matrix_get_cell(A,0,ens_J[i]);
        b = matrix_get_cell(B,0,ens_J[i]);
        c = matrix_get_cell(C,0,ens_J[i]);
        d = matrix_get_cell(D,0,ens_J[i]);
        // puts("evi");
        evi = matrix_get_cell(ev, 0, ens_J[i]);
        // puts("no");
        no1 = mul_Fq(inv_Fq(a), mul_Fq(evi, -b));
        no2 = mul_Fq(inv_Fq(c), mul_Fq(evi, -d));
        // printf("ev(%d) doit être différent de %d et de %d\n", ens_J[i], no1, no2);
        if (no1 != no2) {
          if (no1 != 0 && no2 !=0)
            eui = 0;
          if (no1 != 1 && no2 != 1)
            eui = 1;
          if (no1 != 2 && no2 != 2)
            eui = 2;
          matrix_set_cell(eu, 0, ens_J[i], eui);
        }
        else
          matrix_set_cell(eu_not, 0, ens_J[i], no1);
      }
      // puts("eu"); matrix_print(eu, stdout);
      // puts("eu_not"); matrix_print(eu_not, stdout);


      // on choisit les ensembldes libres I et J
      // int *ens_I = freeset(keys->sk->parite_U,ev,8);

      // on choisi un x aléatoire avec les contraintes précédentes
      x = matrix_copy(eu);
      for (int i = 0; i < SIZE/2; ++i) {
        // puts("*");
        if(matrix_get_cell(x,0,i) == '*') {
          eui_not = matrix_get_cell(eu_not,0,i);
          do {
            r = rand_Fq();
          } while (r == eui_not);
          matrix_set_cell(x, 0, i, r);
        }
      }
      matrix_free(eu_not);
      // puts("x"); matrix_print(x, stdout);

      // on résoud le système en appelant prange_algebra
      // puts("Hu:"); matrix_print(keys->sk->parite_U,stdout);
      // puts("synd:"); matrix_print(synd,stdout);
      // puts("Hu:"); matrix_print(keys->sk->parite_U,stdout);
      eu_int = prange_algebra(keys->sk->parite_U, synd, ens_I, len_I, x);
      matrix_free(x);
      // printf("\teu_int : "); matrix_print(eu_int, stdout);
      // puts("non");
    } while(!eu_int);
    // puts("oui");
    // puts("eu"); matrix_print(eu, stdout);
    // puts("eu_not"); matrix_print(eu_not, stdout);
    // puts("x"); matrix_print(x, stdout);
    // if(!eu_int) {
    //   matrix_free(x);
    //   matrix_free(eu_not);
    //   puts("pas de eu int");
    //   return;
    // }
    e = phi(eu_int, ev);
    w = weight(e);
    matrix_free(e);
    // printf("\twe = %d\n", w);
    // puts("e"); matrix_print(e, stdout); printf("poids : %d\n", w);
  } while (w != OMEGA);
  matrix_copy2(eu, eu_int);
  // on clean
  matrix_free(eu_int);
  // matrix_free(e);
  // matrix_free(x);
  // matrix_free(eu_not);
}

int *freeset(const matrix_t *H, const matrix_t *ev, const int k) {
  if (!H || !ev)
    return NULL;
  // création su supp de ev
  int len_supp = 0;
  int col = matrix_get_col(ev);
  for (int j = 0; j < col; ++j) {
    if (matrix_get_cell(ev,0,j) != 0 % ORDER) {
      ++len_supp;
    }
  }
  int supp[len_supp];
  int i = 0;
  for (int j = 0; j < col; ++j) {
    if (matrix_get_cell(ev,0,j) != 0 % ORDER) {
      supp[i] = j;
      ++i;
    }
  }
  // création de [1,n]\supp(ev)
  int len_inv_supp = SIZE - len_supp;
  int inv_supp[len_inv_supp];
  // int i = 0;
  for (int j = 0; j < SIZE; ++j) {
    if (!is_in_array(supp, len_supp, j)) {
      inv_supp[i] = j;
      ++i;
    }
  }
  // int J1[k];
  // int J2[DIM - SIZE - k];
  int J[DIM-d];
  int rank = 0;
  while (rank != SIZE - DIM) {
    // création de J = J1 U J2
    // ajout de J1 dans la première partie de J
    shuffle(supp, len_supp);
    for (int j = 0; j < k; ++j)
      J[j] = supp[j];
    // ajout de J2 dans la deuxième partie de J
    shuffle(inv_supp, len_inv_supp);
    for (int j = 0; j < DIM-d-k; ++j)
      J[j+k] = supp[j];
    // créatino de [1,SIZE]\J : inv_J
    int len_inv_J = SIZE - DIM - d;
    int inv_J[len_inv_J];
    int i = 0;
    for (int j = 0; j < SIZE; ++j) {
      if (!is_in_array(J, DIM-d, j)) {
        inv_J[i] = j;
        ++i;
      }
    }
    matrix_t *sub_H = sub_col_matrix(H, inv_J, len_inv_J);
    puts("sub_H : ");
    matrix_print(sub_H,stdout);
    rank = matrix_rank(sub_H);
    printf("rang de sub_H : %d\nn-k = %d\nd = %d\n", rank, SIZE-DIM, d);
    matrix_free(sub_H);
  }
  return NULL;
}


////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// PRANGE ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void infoset(int *info, const int n, const int len) {
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

matrix_t *prange_algebra(const matrix_t *parite, const matrix_t *syndrome,
                         const int *info, const int len_i, const matrix_t *x) {
  if (!parite || !syndrome || !x)
    return NULL;
  // matrice de permutation envoyant info sur les DIM dernière coordonnées
  matrix_t *P = NULL;
  matrix_t *HP = NULL;
  matrix_t *left = NULL;
  matrix_t *right = NULL;
  matrix_t *left_inv = NULL;
  int row_HP, col_HP;
  int cpt = 0;
  int n = matrix_get_col(parite);
  // printf("H_U : \n"); matrix_print(parite, stdout);
  // printf("H: %d x %d\n", matrix_get_row(parite), matrix_get_col(parite));
  // printf("n : %d\n",n);
  // choix de P pour avoir A inversible
  while (!left_inv) {
    // puts("left non inversible");
    matrix_free(left);
    matrix_free(right);
    matrix_free(P);
    P = matrix_perm_random_info(n, info, len_i, n);
    // printf("P: %d x %d\n", matrix_get_row(P), matrix_get_col(P));
    if (!P)
      return NULL;
    // puts("ok1");
    // (left | right) <- HP
    // puts("H : "); matrix_print(parite,stdout);
    // puts("P : "); matrix_print(P,stdout);
    HP = matrix_prod(parite,P);
    if (!HP) {
      // puts("pb HP");
      matrix_free(P);
      return NULL;
    }
    // puts("ok2");
    row_HP = matrix_get_row(HP);
    col_HP = matrix_get_col(HP);
    // printf("HP: %d x %d\n", matrix_get_row(HP), matrix_get_col(HP));
    // matrix_print(HP,stdout);
    left = matrix_alloc(row_HP,row_HP);
    right = matrix_alloc(row_HP,col_HP-row_HP);
    if (!left || !right) {
      matrix_free(left);
      matrix_free(right);
      matrix_free(P);
      matrix_free(HP);
      return NULL;
    }
    // puts("ok3");
    matrix_separate(HP, left, right);
    // matrix_free(HP);
    if (!left || !right) {
      matrix_free(left);
      matrix_free(right);
      matrix_free(P);
      matrix_free(HP);
      return NULL;
    }
    // puts("ok4");
    // printf("A: %d x %d\n", matrix_get_row(left), matrix_get_col(left));
    left_inv = matrix_inv(left);
    // if (!left_inv)
      // puts("non inve");
      // matrix_print(left_inv, stdout);
    ++cpt;
    if (cpt > 50) {
      // puts("trop d'essais");
      return NULL;
    }
    // puts("ok5");
  }
  // printf("P : \n"); matrix_print(P, stdout);
  // printf("H*P : \n"); matrix_print(HP, stdout);
  // printf("A : \n"); matrix_print(left, stdout);
  // printf("B : \n"); matrix_print(right, stdout);
  // printf("A-1 : \n"); matrix_print(left_inv, stdout);
  matrix_free(left);
  matrix_free(HP);
  // puts("ok");
  // (zero | ep) <- x
  matrix_t *zero = matrix_alloc(1, SIZE/2-K_U);
  matrix_t *ep = matrix_alloc(1,K_U);
  if (!zero || !ep) {
    matrix_free(zero);
    matrix_free(ep);
    matrix_free(P);
    // matrix_free(left);
    matrix_free(right);
    return NULL;
  }
  // puts("ok6");
  matrix_separate(x, zero, ep);
  // printf("doit être 0 : "); matrix_print(zero, stdout);
  // printf("e' : "); matrix_print(ep, stdout);
  matrix_free(zero);
  if (!ep) {
    matrix_free(ep);
    matrix_free(P);
    // matrix_free(left);
    matrix_free(right);
    return NULL;
  }
  // puts("ok7");
  // calcul de e
  matrix_t *right_T = matrix_trans(right);
  // printf("BT transposée: \n"); matrix_print(right_T, stdout);
  matrix_free(right);
  if (!right_T) {
    matrix_free(ep);
    matrix_free(P);
    // matrix_free(left);
    return NULL;
  }
  // puts("ok8");
  matrix_t *epright_T = matrix_prod(ep,right_T);
  // printf("ep: %d x %d\n", matrix_get_row(ep), matrix_get_col(ep));
  // printf("right_T: %d x %d\n", matrix_get_row(right_T), matrix_get_col(right_T));
  matrix_free(right_T);
  if (!epright_T) {
    matrix_free(ep);
    matrix_free(P);
    // matrix_free(left);
    return NULL;
  }
  // printf("e'*Bt : \n"); matrix_print(epright_T, stdout);
  // puts("ok9");
  matrix_t *epright_Tless = matrix_mul_by_scal(epright_T, -1);
  matrix_free(epright_T);
  if (!epright_Tless) {
    matrix_free(ep);
    matrix_free(P);
    // matrix_free(left);
    return NULL;
  }
  // printf("-e'*Bt : "); matrix_print(epright_Tless, stdout);
  matrix_t *s = matrix_add(syndrome, epright_Tless);
  matrix_free(epright_Tless);
  if (!s) {
    matrix_free(ep);
    matrix_free(P);
    // matrix_free(left);
    return NULL;
  }
  // printf("s-e'*Bt : "); matrix_print(s, stdout);
  // puts("ok10");
  // matrix_t *left_inv = matrix_inv(left);
  // matrix_free(left);
  // if (!left_inv) {
  //   matrix_free(s);
  //   matrix_free(ep);
  //   matrix_free(P);
  //   return NULL;
  // }
  matrix_t *left_inv_T = matrix_trans(left_inv);
  matrix_free(left_inv);
  if (!left_inv_T) {
    matrix_free(s);
    matrix_free(ep);
    matrix_free(P);
    return NULL;
  }
  // printf("(A-1)t : \n"); matrix_print(left_inv_T, stdout);
  // puts("ok11");
  matrix_t *couple_left = matrix_prod(s, left_inv_T);
  matrix_free(s);
  matrix_free(left_inv_T);
  if (!couple_left) {
    matrix_free(ep);
    matrix_free(P);
    return NULL;
  }
  // printf("(s-e'Bt)*(A-1)t : \n"); matrix_print(couple_left, stdout);
  matrix_t *couple = matrix_concatenation(couple_left, ep, 0);
  matrix_free(couple_left);
  matrix_free(ep);
  if (!couple) {
    matrix_free(P);
    return NULL;
  }
  // printf("((s-e'Bt)*(A-1)t,e') : \n"); matrix_print(couple, stdout);
  // puts("ok12");
  matrix_t *P_T = matrix_trans(P);
  matrix_free(P);
  if (!P_T)
    return NULL;
  // printf("P transposée : \n"); matrix_print(P_T, stdout);
  matrix_t *e = matrix_prod(couple, P_T);
  // printf("e : \n"); matrix_print(e, stdout);
  matrix_free(P_T);
  matrix_free(couple);
  return e;
}

matrix_t *iteration_prange(const matrix_t *parite, const matrix_t *syndrome) {
  if (!parite || !syndrome)
    return NULL;
  // initialisation
  prng_init(time(NULL) - getpid());
  int t = rand() % DIM;
  int col_H = matrix_get_col(parite);
  int len_i = 0;
  while(len_i == 0){
    // puts("...");
    len_i = rand() % col_H;
  }
  // choix de l'ensemble d'information info de H
  int info[len_i];
  infoset(info, col_H, len_i);
  // tirage de x de poids t sur les coordonnées info
  matrix_t *x = vector_rand_sub_weight(SIZE, info, len_i, t);
  if (!x)
    return NULL;
  // recuperation de la solution
  matrix_t *e = prange_algebra(parite, syndrome, info, len_i, x);
  matrix_free(x);
  return e;
}


////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// fonction main ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// main (à mettre dans un autre fichier peut-être)
int main(void) {


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
  // int len_i = 6;
  // int info[len_i];
  // infoset(info, 10, len_i);
  // printf("info : ");
  // for (int i = 0; i < len_i; ++i)
  //   printf("%d ", info[i]);
  // puts("");

  //
  // // test de vector_rand_sub_weight
  // int t = 4;
  // matrix_t *vect = vector_rand_sub_weight(10, info, len_i, t);
  // puts("vect");
  // for (int i = 0; i <10; ++i) printf(" %d ", i);
  // puts("");
  // matrix_print(vect, stdout);
  // matrix_free(vect);

  // // test de vector_rand_sub_weight
  // int w = 7;
  // matrix_t *vect = vector_rand_weight(10, w);
  // puts("vect");
  // for (int i = 0; i <10; ++i) printf(" %d ", i);
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

  // test key_gen
  // keys_t *keys = key_gen(0);
  // matrix_print(keys->pk,stdout);
  // key_free(keys);

  ///////////////////////////////// GENERATION DE CLES
  puts("génération des clés de chiffrement avec les paramètres :");
  printf("lambda = %f, n = %d, w = %d, Ku = %d, Kv = %d\n",LAMBDA, SIZE, OMEGA, K_U, K_V);
  keys_t *keys = key_gen(1);
  if (keys == NULL) {
    puts ("Key_gen revoit NULL\n");
    return EXIT_FAILURE;
  }
  puts("\tgénération des clés terminée !!");


  // ////////////////////////////// test decode ev seul
  // matrix_t *ev = matrix_alloc(1,SIZE/2);
  // matrix_t *ev_base = vector_rand(SIZE/2);
  // matrix_t *synd_v = syndrome(ev_base, keys->sk->parite_V);
  // decode_ev(ev, gen_V, synd_v, keys);
  // matrix_t *verif = syndrome(ev, keys->sk->parite_V);
  // printf("\tsyn v : "); matrix_print(synd_v, stdout);
  // printf("\tverif : "); matrix_print(verif, stdout);
  // matrix_free(ev);
  // matrix_free(ev_base);
  // matrix_free(synd_v);
  // matrix_free(verif);

  ///////////////////////////////////////test decode_ev

  matrix_t *ev = matrix_alloc(1,SIZE/2);
  matrix_t *eu = matrix_alloc(1,SIZE/2);
  matrix_t *e = vector_rand_weight(SIZE, OMEGA);
  printf("\te : "); matrix_print(e, stdout);
  matrix_t *synd = syndrome(e, keys->pk);
  printf("\tsyndrome : "); matrix_print(synd, stdout);
  matrix_t *synd_U = matrix_alloc(1,SIZE/2 - K_U);
  matrix_t *synd_V = matrix_alloc(1,SIZE/2 - K_V);
  matrix_separate(synd, synd_U, synd_V);
  if (!ev || !e || !synd || !synd_U || !synd_V) {
    puts("error 1");
    return EXIT_FAILURE;
  }
  puts("test de decode_ev...");
  // printf("dim V %d\n", keys->sk->dim_V);
  decode_ev(ev, gen_V, synd_V, keys);
  if (!ev) {
    puts("error 2");
    return EXIT_FAILURE;
  }
  printf("\tcol : %d, row : %d\n",matrix_get_col(synd_V), matrix_get_row(synd_V));
  // matrix_print(synd_V, stdout);
  matrix_t *verif = syndrome(ev, keys->sk->parite_V);
  if (!verif) {
    puts("\terror 3");
    return EXIT_FAILURE;
  }
  printf("\tsv :    ");
  matrix_print(synd_V, stdout);
  printf("\tverif : ");
  matrix_print(verif, stdout);
  // printf("ev obtenu :   "); matrix_print(ev, stdout);
  // matrix_t *supp_ev = vect_supp(ev);
  // printf("suppp de ev : "); matrix_print(supp_ev, stdout);
  // printf("taille du supp : %d\nKu: %d\n", matrix_get_col(supp_ev), keys->sk->dim_U);
  // matrix_free(supp_ev);
  /////////////////////////////////////////test decode_ev
  puts("test de decove_eu...");
  puts("\tmatrice de parité de U"); matrix_print(keys->sk->parite_U, stdout);
  // matrix_t *synd_U = syndrome(e, keys->sk->parite_U);
  decode_eu(eu, keys, synd_U, ev);
  if (!eu) {
    puts("\terror 2");
    return EXIT_FAILURE;
  }
  // printf("\tcol : %d, row : %d\n",matrix_get_col(synd_U), matrix_get_row(synd_U));
  // matrix_print(synd_U, stdout);
  matrix_t *verif2 = syndrome(eu, keys->sk->parite_U);
  if (!verif2) {
    puts("\terror 3");
    return EXIT_FAILURE;
  }
  printf("eu sortant : "); matrix_print(eu, stdout);
  printf("\tsu :    "); matrix_print(synd_U, stdout);
  printf("\tverif : "); matrix_print(verif2, stdout);

  printf("ev obtenu : "); matrix_print(ev, stdout);
  matrix_t *e_ob = phi(eu,ev);
  printf("e obtenu : "); matrix_print(e_ob, stdout);
  matrix_t *sy = syndrome(e_ob, keys->pk);
  printf("syndrome voulu :  "); matrix_print(synd, stdout);
  printf("syndrome obtenu : "); matrix_print(sy, stdout);
  // printf("a : "); matrix_print(A, stdout);
  // printf("b : "); matrix_print(B, stdout);
  // printf("c : "); matrix_print(C, stdout);
  // printf("d : "); matrix_print(D, stdout);
  matrix_free(e_ob);
  matrix_free(sy);

  // matrix_t *gen_U_T = matrix_trans(gen_U);
  // puts("matrice génératrice de U transposée"); matrix_print(gen_U_T, stdout);
  // matrix_t *gen_U_T_inv = matrix_inv(gen_U_T);
  // puts("inverse de la matrice génératrice de U transposée");
  // matrix_print(gen_U_T_inv, stdout);
  // matrix_t *ver = matrix_prod(gen_U_T, gen_U_T_inv);
  // puts("doit être identité"); matrix_print(ver, stdout);
  //
  // matrix_free(gen_U_T);
  // matrix_free(gen_U_T_inv);
  // matrix_free(ver);
  matrix_free(ev);
  matrix_free(eu);
  matrix_free(verif);
  matrix_free(verif2);
  matrix_free(e);
  matrix_free(synd);
  matrix_free(synd_V);
  matrix_free(synd_U);
  key_free(keys);


  ///////////////////////////////////////////// test matrix_rank
  // matrix_t *m = NULL;
  // m = matrix_random(160,207);
  // int rang = matrix_rank(m);
  // //puts("A"); matrix_print(m,stdout);
  // printf("rang : %d\n",rang);
  // matrix_free(m);

  ///////////////////////////////////////////// test sub_col_matrix
  // matrix_t *m = matrix_random(10,12);
  // int ind[6] = {0,3,4,6,7,11};
  // for (int i = 0; i < 6; ++i)
  //   printf("%d\t",ind[i]);
  // puts("");
  // matrix_t *sub = sub_col_matrix(m,ind,6);
  // puts("m"); matrix_print(m, stdout);
  // printf("rang de m : %d\n",matrix_rank(m));
  // puts("sub"); matrix_print(sub, stdout);
  // printf("rang de sub : %d\n", matrix_rank(sub));
  // matrix_free(m);
  // matrix_free(sub);


  ///////////////////////////////////////////// test random_word
  // matrix_t *G = matrix_random(4,6);
  // if (!G)
  //   return EXIT_FAILURE;
  // matrix_systematisation(G);
  // if (!G)
  //   return EXIT_FAILURE;
  // while (!matrix_is_syst(G)) {
  //   matrix_free(G);
  //   G = matrix_random(4,6);
  //   if (!G)
  //     return EXIT_FAILURE;
  //   matrix_systematisation(G);
  //   if (!G)
  //     return EXIT_FAILURE;
  // }
  // matrix_t *c = matrix_alloc(1,6);
  // if (!G || !c) {
  //   matrix_free(G);
  //   matrix_free(c);
  //   return EXIT_FAILURE;
  // }
  // random_word(c,G);
  // if (!c){
  //   puts("erreur random_word");
  //   matrix_free(G);
  //   return EXIT_FAILURE;
  // }
  // matrix_t *H = matrix_parite(G);
  // if (!H) {
  //   puts("erreur matrix_parite");
  //   matrix_free(G);
  //   matrix_free(c);
  //   return EXIT_FAILURE;
  // }
  // matrix_t *syndrome_c = syndrome(c,H);
  // if (!syndrome_c) {
  //   puts("erreur syndrome");
  //   matrix_free(G);
  //   matrix_free(H);
  //   matrix_free(c);
  //   return EXIT_FAILURE;
  // }
  // puts("G :");
  // matrix_print(G, stdout);
  // puts("H :");
  // matrix_print(H, stdout);
  // puts("c :");
  // matrix_print(c, stdout);
  // puts("syndrome :");
  // matrix_print(syndrome_c, stdout);
  // matrix_free(G);
  // matrix_free(H);
  // matrix_free(c);
  // matrix_free(syndrome_c);



  // printf("%d\n", DIM);

  // // TEST PRANGE
  //
  // // Défini e
  // int ind[SIZE];
  // for (int i = 0; i < SIZE; i++)
  //   ind[i] = i;
  // matrix_t *e = vector_rand_weight(SIZE, OMEGA);
  // // matrix_t *e = matrix_random(1,SIZE);
  //
  // // calcule s le syndrome de e
  // matrix_t *synd = syndrome(e, H);

  // // test inv
  // matrix_t *matrix = matrix_random(20,20);
  // matrix_t *inv = matrix_inv(matrix);
  // matrix_t *prod = matrix_prod(matrix,inv);
  // printf("res = %d\n", is_identity(prod));
  // matrix_free (matrix);
  // matrix_free (inv);
  // matrix_free (prod);

  // // // calcule d'un vecteur erreur associé au syndrome s avec prange iteration
  // matrix_t *ep = NULL;
  // int sb = 0;
  // while (!ep || sb != OMEGA) {
  //   matrix_free(ep);
  //   // puts("ok1");
  //   ep = iteration_prange(H, synd);
  //   if(!ep)
  //     puts("pas de ep");
  //   // puts("ok2");
  //   sb = sub_weight(ep, ind, SIZE);
  //   printf(" %d ",sb);
  // }
  // puts("e");
  // matrix_print(e, stdout);
  // puts("ep");
  // matrix_print(ep, stdout);
  //
  // // Vérifie s = syndrome(ep)
  // matrix_t *verif = syndrome(ep, H);
  // puts("syndrome :");
  // matrix_print(synd, stdout);
  // puts("verif :");
  // matrix_print(verif, stdout);
  // matrix_free(verif);
  // matrix_free(e);
  // matrix_free(synd);
  // matrix_free(ep);
  // key_free(keys);

  // clean up
  puts("Maintenant on nettoie notre bazard...");
  matrix_free(A);
  matrix_free(B);
  matrix_free(C);
  matrix_free(D);
  matrix_free(H);
  matrix_free(gen_U);
  matrix_free(gen_V);
  return EXIT_SUCCESS;
}
