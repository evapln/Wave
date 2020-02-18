#include "scheme.h"

const float LAMBDA = 3.08;
const int SIZE = 50;
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
  matrix_t *P;
  int dim_U;
  int dim_V;
};

struct sign_t {
  matrix_t *e;
  matrix_t *r;
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
  sk->P = matrix_alloc (SIZE,SIZE);
  if (!sk->P) {
    matrix_free(sk->S);
    matrix_free(sk->parite_V);
    matrix_free(sk->parite_U);
    free(sk);
    return NULL;
  }
  return sk;
}

void sk_free(sk_t *sk) {
  if (sk) {
    matrix_free(sk->parite_U);
    matrix_free(sk->parite_V);
    matrix_free(sk->S);
    matrix_free(sk->P);
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

void key_free(keys_t *keys) {
  if (keys) {
    matrix_free(keys->pk);
    sk_free(keys->sk);
  }
  free(keys);
}

sign_t *sign_alloc(void) {
  sign_t *signature = malloc(sizeof(sign_t));
  if (!signature)
    return NULL;
  signature->e = matrix_alloc(1, SIZE);
  if (!signature->e){
    free(signature);
    return NULL;
  }
  signature->r = matrix_alloc(1, SIZE);
  if (!signature->r){
    matrix_free(signature->e);
    free(signature);
    return NULL;
  }
  return signature;
}

void sign_free(sign_t *signature) {
  if (signature) {
    matrix_free(signature->e);
    matrix_free(signature->r);
  }
  free(signature);
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
  // initialisation de a,b,c et d
  coeff_phi(mode);
  if (!coef_init)
    return NULL;

  // création de la matrice génératrice de U
  bool answer = false;
  matrix_t *gen_U_temp = NULL;
  while (!answer) {
    matrix_free(gen_U_temp);
    matrix_free(gen_U);
    // on impose la dimension à K_U au code U
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

  // création de la matrice génératrice de V
  answer = false;
  matrix_t *gen_V_temp = NULL;
  while (!answer) {
    matrix_free(gen_V_temp);
    matrix_free(gen_V);
    // on impose la dimension K_V au code V
    gen_V_temp = matrix_random(K_V,SIZE/2);
    if (!gen_V_temp)
      continue;
    matrix_systematisation(gen_V_temp);
    gen_V = matrix_del_null_row(gen_V_temp);
    if (!gen_V) {
      matrix_free(gen_V_temp);
      continue;
    }
    answer = matrix_is_syst(gen_V);
  }
  matrix_free(gen_V_temp);

  // dimension des codes
  int dim_U = matrix_get_row(gen_U);
  int dim_V = matrix_get_row(gen_V);
  DIM = dim_U + dim_V;

  // création de H_U et H_V
  matrix_t *H_U = matrix_parite(gen_U);
  matrix_t *H_V = matrix_parite(gen_V);
  if (!H_U || !H_V) {
    matrix_free(H_U);
    matrix_free(H_V);
    return NULL;
  }

  // alloue la structure keys_t et remplie avec les matrices de parités
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

  // création de H
  H = parite(keys->sk->parite_U, keys->sk->parite_V);
  if (!H) {
    key_free(keys);
    return NULL;
  }
  // création de S aléaoire inversible
  matrix_t *S = NULL;
  matrix_t *S_inv = NULL;
  while (!S || !S_inv) {
    matrix_free(S);
    matrix_free(S_inv);
    S = matrix_random(SIZE-DIM,SIZE-DIM);
    S_inv = matrix_inv(S);
  }
  matrix_free(S_inv);
  matrix_copy2(keys->sk->S, S);
  matrix_free(S);
  if (!keys->sk->S) {
    key_free(keys);
    return NULL;
  }
  // créationde P matrice de permutation aléatoire
  matrix_t *permut = matrix_perm_random(SIZE);
  if (!permut) {
    key_free(keys);
    return NULL;
  }
  matrix_copy2(keys->sk->P, permut);
  matrix_free(permut);
  if (!keys->sk->P) {
    key_free(keys);
    return NULL;
  }
  keys->sk->dim_U = dim_U;
  keys->sk->dim_V = dim_V;

  // création de la clé publique SHP
  puts("\tet enfin notre clé publique...");
  matrix_t *SH = matrix_prod(keys->sk->S,H);
  if (!SH) {
    key_free(keys);
    return NULL;
  }
  matrix_t *SHP = matrix_prod(SH,keys->sk->P);
  matrix_free(SH);
  if (!SHP) {
    key_free(keys);
    return NULL;
  }
  matrix_copy2(keys->pk, SHP);
  matrix_free(SHP);
  if (!keys->pk) {
    key_free(keys);
    return NULL;
  }
  return keys;
}

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// décodage /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void decode_ev(matrix_t * ev,const matrix_t *G, const matrix_t *synd) {
  if (!ev || !G || !synd)
    return;
  int col = matrix_get_col(G);
  if (matrix_get_col(ev) != col || matrix_get_row(ev) != 1)
    return;
  random_word(ev, G);
  matrix_t *s = matrix_alloc(1,col);
  int col_synd = matrix_get_col(synd);
  for (int i = 0; i < col - col_synd; ++i)
    matrix_set_cell(s, 0, i, 0);
  for (int i = col - col_synd; i < col; ++i)
    matrix_set_cell(s, 0, i, matrix_get_cell(synd, 0, i - col + col_synd));
  matrix_add_modified(ev, s, 1);
  matrix_free(s);
}

void decode_eu(matrix_t * eu, const sk_t *sk, const matrix_t *synd,
               const matrix_t *ev) {
  // Vérification des entrées
  if (!eu || !synd || !ev || !sk)
    return;
  if (matrix_get_row(ev) != 1 || matrix_get_row(eu) != 1)
    return;
  // déclarations
  int dim_U = sk->dim_U;
  int a, b, c, d, w, r;
  char no1, no2, evi, eui, eui_not;
  matrix_t *eu_int = NULL;
  matrix_t *x = NULL;
  matrix_t *e = NULL;
  matrix_t * eu_not = NULL;

  // création des ensembles I et J
  int ens[SIZE/2];
  for (int i = 0; i < SIZE/2; ++i)
    ens[i] = i;
  do {
    matrix_free(eu_int); // en cas de repassage dans la boucle
    // création de I aléatoirement
    shuffle(ens,SIZE/2);
    int len_I = rand() % (SIZE/2 - dim_U) + dim_U; // aléa entre dim_U et SIZE/2
    int ens_I[len_I];
    for (int i = 0; i < len_I; ++i)
      ens_I[i] = ens[i];
    do {
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
      /* sur les ku positions choisies (dans J), on fixe eu de telle sorte que
      A(i)*eu(i) + B(i)*ev(i) != 0 et  C(i)*eu(i) + D(i)*ev(i) != 0
      quand on a une solution exacte, on met dans eu,
      quand on a juste une impossibilité, on met dans eu_not */
      for (int i = 0; i < dim_U; ++i) {
        a = matrix_get_cell(A,0,ens_J[i]);
        b = matrix_get_cell(B,0,ens_J[i]);
        c = matrix_get_cell(C,0,ens_J[i]);
        d = matrix_get_cell(D,0,ens_J[i]);
        evi = matrix_get_cell(ev, 0, ens_J[i]);
        no1 = mul_Fq(inv_Fq(a), mul_Fq(evi, -b));
        no2 = mul_Fq(inv_Fq(c), mul_Fq(evi, -d));
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

      // on choisi un x aléatoire avec les contraintes précédentes
      x = matrix_copy(eu);
      for (int i = 0; i < SIZE/2; ++i) {
        if(matrix_get_cell(x,0,i) == '*') {
          eui_not = matrix_get_cell(eu_not,0,i);
          do {
            r = rand_Fq();
          } while (r == eui_not);
          matrix_set_cell(x, 0, i, r);
        }
      }
      matrix_free(eu_not);

      // on résoud le système en appelant prange_algebra
      eu_int = prange_algebra(sk->parite_U, synd, ens_I, len_I, x);
      matrix_free(x);
    } while(!eu_int);
    e = phi(eu_int, ev);
    w = weight(e);
    matrix_free(e);
  } while (w != OMEGA);
  matrix_copy2(eu, eu_int);
  matrix_free(eu_int);
}

matrix_t *decode_uv(const sk_t *sk, const matrix_t *synd) {
  // vérifiaction des entrées
  if(!sk || !synd)
    return NULL;

  //initialisation
  matrix_t *su = matrix_alloc(1,SIZE/2 - K_U);
  matrix_t *sv = matrix_alloc(1,SIZE/2 - K_V);
  matrix_separate(synd, su, sv);
  if (!su || !sv)
    return NULL;

  // decodage de ev
  matrix_t *ev = matrix_alloc(1,SIZE/2);
  decode_ev(ev, gen_V, sv);
  if (!ev)
    return NULL;
  matrix_t *verifv = syndrome(ev, sk->parite_V);
  if (!verifv)
    return  NULL;

  // decodage de eu
  matrix_t *eu = matrix_alloc(1,SIZE/2);
  decode_eu(eu, sk, su, ev);
  if (!eu)
    return NULL;
  matrix_t *verifu = syndrome(eu, sk->parite_U);
  if (!verifu)
    return NULL;
  matrix_t *e = phi(eu,ev);
  matrix_free(ev);
  matrix_free(eu);
  matrix_free(verifv);
  matrix_free(verifu);
  matrix_free(su);
  matrix_free(sv);
  return e;
}

matrix_t *invert_alg(const sk_t *sk, const matrix_t *synd) {
  if (!sk || !synd)
    return NULL;
  matrix_t *S_inv = matrix_inv(sk->S);
  if (!S_inv)
    return NULL;
  matrix_t *S_inv_T = matrix_trans(S_inv);
  matrix_free(S_inv);
  if (!S_inv_T)
    return NULL;
  matrix_t *entree = matrix_prod(synd, S_inv_T);
  matrix_free(S_inv_T);
  if (!entree)
    return NULL;
  matrix_t *e = decode_uv(sk, entree);
  matrix_free(entree);
  if (!e)
    return NULL;
  matrix_t *eP = matrix_prod(e,sk->P);
  matrix_free(e);
  return eP;
}

// int *freeset(const matrix_t *H, const matrix_t *ev, const int k) {
//   if (!H || !ev)
//     return NULL;
//   // création su supp de ev
//   int len_supp = 0;
//   int col = matrix_get_col(ev);
//   for (int j = 0; j < col; ++j) {
//     if (matrix_get_cell(ev,0,j) != 0 % ORDER) {
//       ++len_supp;
//     }
//   }
//   int supp[len_supp];
//   int i = 0;
//   for (int j = 0; j < col; ++j) {
//     if (matrix_get_cell(ev,0,j) != 0 % ORDER) {
//       supp[i] = j;
//       ++i;
//     }
//   }
//   // création de [1,n]\supp(ev)
//   int len_inv_supp = SIZE - len_supp;
//   int inv_supp[len_inv_supp];
//   // int i = 0;
//   for (int j = 0; j < SIZE; ++j) {
//     if (!is_in_array(supp, len_supp, j)) {
//       inv_supp[i] = j;
//       ++i;
//     }
//   }
//   // int J1[k];
//   // int J2[DIM - SIZE - k];
//   int J[DIM-d];
//   int rank = 0;
//   while (rank != SIZE - DIM) {
//     // création de J = J1 U J2
//     // ajout de J1 dans la première partie de J
//     shuffle(supp, len_supp);
//     for (int j = 0; j < k; ++j)
//       J[j] = supp[j];
//     // ajout de J2 dans la deuxième partie de J
//     shuffle(inv_supp, len_inv_supp);
//     for (int j = 0; j < DIM-d-k; ++j)
//       J[j+k] = supp[j];
//     // créatino de [1,SIZE]\J : inv_J
//     int len_inv_J = SIZE - DIM - d;
//     int inv_J[len_inv_J];
//     int i = 0;
//     for (int j = 0; j < SIZE; ++j) {
//       if (!is_in_array(J, DIM-d, j)) {
//         inv_J[i] = j;
//         ++i;
//       }
//     }
//     matrix_t *sub_H = sub_col_matrix(H, inv_J, len_inv_J);
//     puts("sub_H : ");
//     matrix_print(sub_H,stdout);
//     rank = matrix_rank(sub_H);
//     printf("rang de sub_H : %d\nn-k = %d\nd = %d\n", rank, SIZE-DIM, d);
//     matrix_free(sub_H);
//   }
//   return NULL;
// }


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
  // déclarations
  matrix_t *P = NULL; // matrice de permutation envoyant info sur les DIM dernières coordonnées
  matrix_t *HP = NULL; // produit parite * P
  matrix_t *left = NULL; // partie gauche de HP (nommée A dans le rapport)
  matrix_t *right = NULL; // parti droite de HP (nommée B dans le rapport)
  matrix_t *left_inv = NULL; //inverse de la partie gauche de HP
  int row_HP, col_HP; // nombre de lignes et de colonnes de HP
  int cpt = 0; // compteur de tour de boucle
  int n = matrix_get_col(parite); // nombre de colonne de parite
  // choix de P pour avoir A inversible
  while (!left_inv) {
    // liberation de la mémoire en cas de repassage dans la boucle
    matrix_free(right);
    matrix_free(P);
    // P aléatoire
    P = matrix_perm_random_info(n, info, len_i, n);
    if (!P)
      return NULL;
    // produit parite * P
    HP = matrix_prod(parite,P);
    if (!HP) {
      matrix_free(P);
      return NULL;
    }
    // on rempli left et right telles que (left | right) <- HP
    row_HP = matrix_get_row(HP);
    col_HP = matrix_get_col(HP);
    left = matrix_alloc(row_HP,row_HP);
    right = matrix_alloc(row_HP,col_HP-row_HP);
    if (!left || !right) {
      matrix_free(left);
      matrix_free(right);
      matrix_free(P);
      matrix_free(HP);
      return NULL;
    }
    matrix_separate(HP, left, right);
    matrix_free(HP);
    if (!left || !right) {
      matrix_free(left);
      matrix_free(right);
      matrix_free(P);
      return NULL;
    }
    // inversion de left
    left_inv = matrix_inv(left);
    matrix_free(left);
    ++cpt;
    if (cpt > 10) {
      return NULL;
    }
  }
  // declaration et remplissage de zero et ep tels que (zero | ep) <- x
  matrix_t *zero = matrix_alloc(1, SIZE/2-K_U);
  matrix_t *ep = matrix_alloc(1,K_U); // nommé e' dans le rapport
  if (!zero || !ep) {
    matrix_free(zero);
    matrix_free(ep);
    matrix_free(P);
    matrix_free(right);
    return NULL;
  }
  matrix_separate(x, zero, ep);
  matrix_free(zero);
  if (!ep) {
    matrix_free(ep);
    matrix_free(P);
    matrix_free(right);
    return NULL;
  }
  // calcul de e
  matrix_t *right_T = matrix_trans(right); // transposée de B
  matrix_free(right);
  if (!right_T) {
    matrix_free(ep);
    matrix_free(P);
    return NULL;
  }
  matrix_t *epright_T = matrix_prod(ep,right_T); // e'*Bt
  matrix_free(right_T);
  if (!epright_T) {
    matrix_free(ep);
    matrix_free(P);
    return NULL;
  }
  matrix_t *s = matrix_copy(syndrome);
  if (!s) {
    matrix_free(ep);
    matrix_free(P);
    return NULL;
  }
  matrix_add_modified(s, epright_T, -1); // s - e'*Bt
  matrix_free(epright_T);
  if (!s) {
    matrix_free(ep);
    matrix_free(P);
    return NULL;
  }
  matrix_t *left_inv_T = matrix_trans(left_inv); // (A-1)t
  matrix_free(left_inv);
  if (!left_inv_T) {
    matrix_free(s);
    matrix_free(ep);
    matrix_free(P);
    return NULL;
  }
  matrix_t *couple_left = matrix_prod(s, left_inv_T); // (s - e'*Bt)*(A-1)t
  matrix_free(s);
  matrix_free(left_inv_T);
  if (!couple_left) {
    matrix_free(ep);
    matrix_free(P);
    return NULL;
  }
  matrix_t *couple = matrix_concatenation(couple_left, ep, 0); // [(s - e'*Bt)*(A-1)t, e']
  matrix_free(couple_left);
  matrix_free(ep);
  if (!couple) {
    matrix_free(P);
    return NULL;
  }
  matrix_t *P_T = matrix_trans(P); // Pt
  matrix_free(P);
  if (!P_T)
    return NULL;
  matrix_t *e = matrix_prod(couple, P_T);
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
  while(len_i == 0)
    len_i = rand() % col_H;
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
////////////////////////// signature et vérification ///////////////////////////
////////////////////////////////////////////////////////////////////////////////

matrix_t *hash(const matrix_t *m, const matrix_t *r, const matrix_t *pk) {
  if (!m || !r)
    return NULL;
  int col = matrix_get_col(r);
  matrix_t *m_high = matrix_init(1, col, 0);
  if (!m_high)
    return NULL;
  for (int i = 0; i < matrix_get_col(m); ++i)
    matrix_set_cell(m_high, 0, i, matrix_get_cell(m,0,i));
  matrix_add_modified(m_high, r, 1);
  if (!m_high)
    return NULL;
  matrix_t *s = syndrome(m_high, pk);
  return s;
}

sign_t *sign(const keys_t *keys, const matrix_t *m) {
  // vérifiations des entrées
  if (!keys || !m)
    return NULL;
  // tirage de r aléatoirement de taille SIZE
  matrix_t *r = matrix_random(1,SIZE);
  if (!r)
    return NULL;
  // hachage de (m,r)
  matrix_t *s = hash(m,r, keys->pk);
  if (!s) {
    matrix_free(r);
    return NULL;
  }
  // inversion de f
  matrix_t *e = invert_alg(keys->sk, s);
  matrix_free(s);
  if (!e) {
    matrix_free(r);
    return NULL;
  }
  // allocation de la structure sign_t
  sign_t *signature = sign_alloc();
  if (!signature) {
    matrix_free(r);
    matrix_free(e);
    return NULL;
  }
  matrix_copy2(signature->e, e);
  matrix_copy2(signature->r, r);
  matrix_free(e);
  matrix_free(r);
  return signature;
}

bool verify(const matrix_t *pk, const matrix_t *m, const sign_t *signature) {
  // vérifiations des entrées
  if (!pk || !m || !signature || !signature->r || !signature->e)
    return false;
  // calcul du poids de e
  if (weight(signature->e) != OMEGA)
    return false;
  // hachage de (m,r)
  matrix_t *s = hash(m,signature->r, pk);
  if (!s)
    return false;
  // calcul du syndrome de e
  matrix_t *synd = syndrome(signature->e, pk);
  if (!synd) {
    matrix_free(s);
    return false;
  }
  // vérification de l'égalité s = e*Ht
  if (!matrix_is_equal(s, synd)) {
    matrix_free(s);
    matrix_free(synd);
    return false;
  }
  matrix_free(s);
  matrix_free(synd);
  return true;
}

// main (à mettre dans un autre fichier peut-être)
// int main(void) {
//
//
//   // matrix_t *X = NULL;
//   // X = matrix_random(3,4);
//   // puts("X : \n");
//   // matrix_print(X, stdout);
//   // matrix_t *Y = NULL;
//   // Y = matrix_random(2,4);
//   // puts("Y : \n");
//   // matrix_print(Y, stdout);
//   // matrix_t *pariteUV = NULL;
//   // pariteUV = parite (X,Y);
//   // puts("parite : \n");
//   // matrix_print(pariteUV, stdout);
//   //
//   // matrix_t *A = NULL;
//   // A = matrix_random(4,6);
//   // matrix_print(A, stdout);
//   // matrix_systematisation(A);
//   // matrix_print(A, stdout);
//   // if (matrix_is_syst(A))
//   //   puts("A syst");
//   // matrix_t *B = NULL;
//   // B = matrix_del_null_row(A);
//   // matrix_print(B, stdout);
//   // if (matrix_is_syst(B))
//   //   puts("B syst");
//   // matrix_free(A);
//   // matrix_free(B);
//
//   // matrix_t *A = NULL;
//   // A = matrix_alloc(4,6);
//   // matrix_set_cell(A,0,0,1);
//   // matrix_set_cell(A,0,1,1);
//   // matrix_set_cell(A,0,2,1);
//   // matrix_set_cell(A,0,3,0);
//   // matrix_set_cell(A,0,4,1);
//   // matrix_set_cell(A,0,5,1);
//   // matrix_set_cell(A,1,0,0);
//   // matrix_set_cell(A,1,1,1);
//   // matrix_set_cell(A,1,2,0);
//   // matrix_set_cell(A,1,3,0);
//   // matrix_set_cell(A,1,4,1);
//   // matrix_set_cell(A,1,5,1);
//   // matrix_set_cell(A,2,0,1);
//   // matrix_set_cell(A,2,1,0);
//   // matrix_set_cell(A,2,2,1);
//   // matrix_set_cell(A,2,3,1);
//   // matrix_set_cell(A,2,4,0);
//   // matrix_set_cell(A,2,5,1);
//   // matrix_set_cell(A,3,0,0);
//   // matrix_set_cell(A,3,1,1);
//   // matrix_set_cell(A,3,2,1);
//   // matrix_set_cell(A,3,3,1);
//   // matrix_set_cell(A,3,4,0);
//   // matrix_set_cell(A,3,5,1);
//   // matrix_print(A,stdout);
//   // matrix_systematisation(A);
//   // matrix_print(A,stdout);
//   // // matrix_free(A);
//   // matrix_t *B = NULL;
//   // B = matrix_parite(A);
//   // if (!B)
//   //   puts("error");
//   // else
//   //   matrix_print(B, stdout);
//   // matrix_free(B);
//
//   // // test de sub_weight
//   // matrix_t *vect = matrix_random(1,13);
//   // for (int i = 0; i <13; ++i) printf(" %d ", i);
//   // puts("");
//   // matrix_print(vect, stdout);
//   // int len_s = 6;
//   // int subset[6] = {2,5,6,8,10,12};
//   // printf("%d\n", sub_weight(vect, subset, len_s));
//   // matrix_free(vect);
//
//
//   // // test de infoset
//   // int len_i = 6;
//   // int info[len_i];
//   // infoset(info, 10, len_i);
//   // printf("info : ");
//   // for (int i = 0; i < len_i; ++i)
//   //   printf("%d ", info[i]);
//   // puts("");
//
//   //
//   // // test de vector_rand_sub_weight
//   // int t = 4;
//   // matrix_t *vect = vector_rand_sub_weight(10, info, len_i, t);
//   // puts("vect");
//   // for (int i = 0; i <10; ++i) printf(" %d ", i);
//   // puts("");
//   // matrix_print(vect, stdout);
//   // matrix_free(vect);
//
//   // // test de vector_rand_sub_weight
//   // int w = 7;
//   // matrix_t *vect = vector_rand_weight(10, w);
//   // puts("vect");
//   // for (int i = 0; i <10; ++i) printf(" %d ", i);
//   // puts("");
//   // matrix_print(vect, stdout);
//   // matrix_free(vect);
//
//   // // test matrix_perm_random_info et matrix_perm_random
//   // matrix_t *matrix = matrix_perm_random(20);
//   // puts("permutation normale :");
//   // matrix_print(matrix, stdout);
//   //
//   // int tab[5] = {0,1,2,3,4};
//   // DIM = 7;
//   // matrix_t *matr = matrix_perm_random_info(20, tab, 5, DIM);
//   // puts("permutation avec les 5 premières coordonnées à la fin :");
//   // matrix_print(matr, stdout);
//   //
//   // matrix_free(matrix);
//   // matrix_free(matr);
//
//   //// test matrix_separate
//   // matrix_t *matrix = matrix_random(10,17);
//   // matrix_t *left = matrix_alloc(10,6);
//   // matrix_t *right = matrix_alloc(10,11);
//   // matrix_separate(matrix, left, right);
//   // matrix_print(matrix, stdout);
//   // matrix_print(left, stdout);
//   // matrix_print(right, stdout);
//   // matrix_free(matrix);
//   // matrix_free(left);
//   // matrix_free(right);
//
//   // test key_gen
//   // keys_t *keys = key_gen(0);
//   // matrix_print(keys->pk,stdout);
//   // key_free(keys);
//
//   ///////////////////////////////// GENERATION DE CLES
//   puts("génération des clés de chiffrement avec les paramètres :");
//   printf("lambda = %f, n = %d, w = %d, Ku = %d, Kv = %d\n",LAMBDA, SIZE, OMEGA, K_U, K_V);
//   keys_t *keys = key_gen(0);
//   if (keys == NULL) {
//     puts ("Key_gen revoit NULL\n");
//     return EXIT_FAILURE;
//   }
//   puts("\tgénération des clés terminée !!");
//
//   matrix_t *m = vector_rand(10);
//   printf("m : "); matrix_print(m, stdout);
//   sign_t *signature = sign(keys, m);
//   if(!signature) {
//     puts("ouch");
//     return EXIT_FAILURE;
//   }
//   printf("signature de m : \ne :"); matrix_print(signature->e, stdout);
//   printf("r :"); matrix_print(signature->r, stdout);
//   printf("la signature est correct : %d\n", verify(keys->pk, m, signature));
//   matrix_free(m);
//   sign_free(signature);
//   /////////////////////////////////// test decode_uv
//   // // puts("initialisation de e_test :");
//   // matrix_t *e_test = vector_rand_weight(SIZE, OMEGA);
//   // if (!e_test)
//   //   return EXIT_FAILURE;
//   // // puts("initialisation de synd_test :");
//   // matrix_t *synd_test = syndrome(e_test, keys->pk);
//   // if (!synd_test)
//   //   return EXIT_FAILURE;
//   // // printf("\tS : \n"); matrix_print(keys->sk->S, stdout);
//   // // printf("\tS-1 : \n"); matrix_print(S_inv, stdout);
//   // // printf("\tS-1t : \n"); matrix_print(S_inv_T, stdout);
//   // // printf("\tentre : \n"); matrix_print(entree, stdout);
//   // matrix_t *e = invert_alg(keys, synd_test);
//   // matrix_t *synd = syndrome(e, keys->pk);
//   // if (!synd)
//   //   return EXIT_FAILURE;
//   // printf("\te de base : "); matrix_print(e_test, stdout);
//   // printf("\te obtenu  : "); matrix_print(eP, stdout);
//   // printf("\tsyndrome de base : "); matrix_print(synd_test, stdout);
//   // printf("\tsyndrome obtenu  : "); matrix_print(synd, stdout);
//   // printf("vraie %d\n", true);
//   // printf("syndromes égaux : %d \n", matrix_is_equal(synd, synd_test));
//   // matrix_free(e_test);
//   // matrix_free(e);
//   // matrix_free(synd_test);
//   // matrix_free(synd);
//
//   // long long mm = 1911191110;
//   // matrix_t *m = int_to_vect(mm);
//   // long long mmm = vect_to_int(m);
//   // printf("m = %lld = ", mm); matrix_print(m, stdout);
//   // printf("= %lld\n", mmm);
//   // matrix_free(m);
//   // printf("a : "); matrix_print(A, stdout);
//   // printf("b : "); matrix_print(B, stdout);
//   // printf("c : "); matrix_print(C, stdout);
//   // printf("d : "); matrix_print(D, stdout);
//
//   // ////////////////////////////// test decode ev seul
//   // matrix_t *ev = matrix_alloc(1,SIZE/2);
//   // matrix_t *ev_base = vector_rand(SIZE/2);
//   // matrix_t *synd_v = syndrome(ev_base, keys->sk->parite_V);
//   // decode_ev(ev, gen_V, synd_v, keys);
//   // matrix_t *verif = syndrome(ev, keys->sk->parite_V);
//   // printf("\tsyn v : "); matrix_print(synd_v, stdout);
//   // printf("\tverif : "); matrix_print(verif, stdout);
//   // matrix_free(ev);
//   // matrix_free(ev_base);
//   // matrix_free(synd_v);
//   // matrix_free(verif);
//
//   ///////////////////////////////////////test decode_ev
//   //
//   // matrix_t *ev = matrix_alloc(1,SIZE/2);
//   // matrix_t *eu = matrix_alloc(1,SIZE/2);
//   // matrix_t *e = vector_rand_weight(SIZE, OMEGA);
//   // printf("\te : "); matrix_print(e, stdout);
//   // matrix_t *synd = syndrome(e, keys->pk);
//   // printf("\tsyndrome : "); matrix_print(synd, stdout);
//   // matrix_t *synd_U = matrix_alloc(1,SIZE/2 - K_U);
//   // matrix_t *synd_V = matrix_alloc(1,SIZE/2 - K_V);
//   // matrix_separate(synd, synd_U, synd_V);
//   // if (!ev || !e || !synd || !synd_U || !synd_V) {
//   //   puts("error 1");
//   //   return EXIT_FAILURE;
//   // }
//   // puts("test de decode_ev...");
//   // // printf("dim V %d\n", keys->sk->dim_V);
//   // decode_ev(ev, gen_V, synd_V, keys);
//   // if (!ev) {
//   //   puts("error 2");
//   //   return EXIT_FAILURE;
//   // }
//   // printf("\tcol : %d, row : %d\n",matrix_get_col(synd_V), matrix_get_row(synd_V));
//   // // matrix_print(synd_V, stdout);
//   // matrix_t *verif = syndrome(ev, keys->sk->parite_V);
//   // if (!verif) {
//   //   puts("\terror 3");
//   //   return EXIT_FAILURE;
//   // }
//   // printf("\tsv :    ");
//   // matrix_print(synd_V, stdout);
//   // printf("\tverif : ");
//   // matrix_print(verif, stdout);
//   // // printf("ev obtenu :   "); matrix_print(ev, stdout);
//   // // matrix_t *supp_ev = vect_supp(ev);
//   // // printf("suppp de ev : "); matrix_print(supp_ev, stdout);
//   // // printf("taille du supp : %d\nKu: %d\n", matrix_get_col(supp_ev), keys->sk->dim_U);
//   // // matrix_free(supp_ev);
//   // /////////////////////////////////////////test decode_ev
//   // puts("test de decove_eu...");
//   // puts("\tmatrice de parité de U"); matrix_print(keys->sk->parite_U, stdout);
//   // // matrix_t *synd_U = syndrome(e, keys->sk->parite_U);
//   // decode_eu(eu, keys, synd_U, ev);
//   // if (!eu) {
//   //   puts("\terror 2");
//   //   return EXIT_FAILURE;
//   // }
//   // // printf("\tcol : %d, row : %d\n",matrix_get_col(synd_U), matrix_get_row(synd_U));
//   // // matrix_print(synd_U, stdout);
//   // matrix_t *verif2 = syndrome(eu, keys->sk->parite_U);
//   // if (!verif2) {
//   //   puts("\terror 3");
//   //   return EXIT_FAILURE;
//   // }
//   // printf("eu sortant : "); matrix_print(eu, stdout);
//   // printf("\tsu :    "); matrix_print(synd_U, stdout);
//   // printf("\tverif : "); matrix_print(verif2, stdout);
//   //
//   // printf("ev obtenu : "); matrix_print(ev, stdout);
//   // matrix_t *e_ob = phi(eu,ev);
//   // printf("e obtenu : "); matrix_print(e_ob, stdout);
//   // matrix_t *sy = syndrome(e_ob, keys->pk);
//   // printf("syndrome voulu :  "); matrix_print(synd, stdout);
//   // printf("syndrome obtenu : "); matrix_print(sy, stdout);
//   // // printf("a : "); matrix_print(A, stdout);
//   // // printf("b : "); matrix_print(B, stdout);
//   // // printf("c : "); matrix_print(C, stdout);
//   // // printf("d : "); matrix_print(D, stdout);
//   // matrix_free(e_ob);
//   // matrix_free(sy);
//   //
//   // // matrix_t *gen_U_T = matrix_trans(gen_U);
//   // // puts("matrice génératrice de U transposée"); matrix_print(gen_U_T, stdout);
//   // // matrix_t *gen_U_T_inv = matrix_inv(gen_U_T);
//   // // puts("inverse de la matrice génératrice de U transposée");
//   // // matrix_print(gen_U_T_inv, stdout);
//   // // matrix_t *ver = matrix_prod(gen_U_T, gen_U_T_inv);
//   // // puts("doit être identité"); matrix_print(ver, stdout);
//   // //
//   // // matrix_free(gen_U_T);
//   // // matrix_free(gen_U_T_inv);
//   // // matrix_free(ver);
//   // matrix_free(ev);
//   // matrix_free(eu);
//   // matrix_free(verif);
//   // matrix_free(verif2);
//   // matrix_free(e);
//   // matrix_free(synd);
//   // matrix_free(synd_V);
//   // matrix_free(synd_U);
//   // key_free(keys);
//
//
//   ///////////////////////////////////////////// test matrix_rank
//   // matrix_t *m = NULL;
//   // m = matrix_random(160,207);
//   // int rang = matrix_rank(m);
//   // //puts("A"); matrix_print(m,stdout);
//   // printf("rang : %d\n",rang);
//   // matrix_free(m);
//
//   ///////////////////////////////////////////// test sub_col_matrix
//   // matrix_t *m = matrix_random(10,12);
//   // int ind[6] = {0,3,4,6,7,11};
//   // for (int i = 0; i < 6; ++i)
//   //   printf("%d\t",ind[i]);
//   // puts("");
//   // matrix_t *sub = sub_col_matrix(m,ind,6);
//   // puts("m"); matrix_print(m, stdout);
//   // printf("rang de m : %d\n",matrix_rank(m));
//   // puts("sub"); matrix_print(sub, stdout);
//   // printf("rang de sub : %d\n", matrix_rank(sub));
//   // matrix_free(m);
//   // matrix_free(sub);
//
//
//   ///////////////////////////////////////////// test random_word
//   // matrix_t *G = matrix_random(4,6);
//   // if (!G)
//   //   return EXIT_FAILURE;
//   // matrix_systematisation(G);
//   // if (!G)
//   //   return EXIT_FAILURE;
//   // while (!matrix_is_syst(G)) {
//   //   matrix_free(G);
//   //   G = matrix_random(4,6);
//   //   if (!G)
//   //     return EXIT_FAILURE;
//   //   matrix_systematisation(G);
//   //   if (!G)
//   //     return EXIT_FAILURE;
//   // }
//   // matrix_t *c = matrix_alloc(1,6);
//   // if (!G || !c) {
//   //   matrix_free(G);
//   //   matrix_free(c);
//   //   return EXIT_FAILURE;
//   // }
//   // random_word(c,G);
//   // if (!c){
//   //   puts("erreur random_word");
//   //   matrix_free(G);
//   //   return EXIT_FAILURE;
//   // }
//   // matrix_t *H = matrix_parite(G);
//   // if (!H) {
//   //   puts("erreur matrix_parite");
//   //   matrix_free(G);
//   //   matrix_free(c);
//   //   return EXIT_FAILURE;
//   // }
//   // matrix_t *syndrome_c = syndrome(c,H);
//   // if (!syndrome_c) {
//   //   puts("erreur syndrome");
//   //   matrix_free(G);
//   //   matrix_free(H);
//   //   matrix_free(c);
//   //   return EXIT_FAILURE;
//   // }
//   // puts("G :");
//   // matrix_print(G, stdout);
//   // puts("H :");
//   // matrix_print(H, stdout);
//   // puts("c :");
//   // matrix_print(c, stdout);
//   // puts("syndrome :");
//   // matrix_print(syndrome_c, stdout);
//   // matrix_free(G);
//   // matrix_free(H);
//   // matrix_free(c);
//   // matrix_free(syndrome_c);
//
//
//
//   // printf("%d\n", DIM);
//
//   // // TEST PRANGE
//   //
//   // // Défini e
//   // int ind[SIZE];
//   // for (int i = 0; i < SIZE; i++)
//   //   ind[i] = i;
//   // matrix_t *e = vector_rand_weight(SIZE, OMEGA);
//   // // matrix_t *e = matrix_random(1,SIZE);
//   //
//   // // calcule s le syndrome de e
//   // matrix_t *synd = syndrome(e, H);
//
//   // // test inv
//   // matrix_t *matrix = matrix_random(20,20);
//   // matrix_t *inv = matrix_inv(matrix);
//   // matrix_t *prod = matrix_prod(matrix,inv);
//   // printf("res = %d\n", matrix_is_identity(prod));
//   // matrix_free (matrix);
//   // matrix_free (inv);
//   // matrix_free (prod);
//
//   // // // calcule d'un vecteur erreur associé au syndrome s avec prange iteration
//   // matrix_t *ep = NULL;
//   // int sb = 0;
//   // while (!ep || sb != OMEGA) {
//   //   matrix_free(ep);
//   //   // puts("ok1");
//   //   ep = iteration_prange(H, synd);
//   //   if(!ep)
//   //     puts("pas de ep");
//   //   // puts("ok2");
//   //   sb = sub_weight(ep, ind, SIZE);
//   //   printf(" %d ",sb);
//   // }
//   // puts("e");
//   // matrix_print(e, stdout);
//   // puts("ep");
//   // matrix_print(ep, stdout);
//   //
//   // // Vérifie s = syndrome(ep)
//   // matrix_t *verif = syndrome(ep, H);
//   // puts("syndrome :");
//   // matrix_print(synd, stdout);
//   // puts("verif :");
//   // matrix_print(verif, stdout);
//   // matrix_free(verif);
//   // matrix_free(e);
//   // matrix_free(synd);
//   // matrix_free(ep);
//   // key_free(keys);
//
//   // clean up
//   puts("Maintenant on nettoie notre bazard...");
//   matrix_free(A);
//   matrix_free(B);
//   matrix_free(C);
//   matrix_free(D);
//   matrix_free(H);
//   matrix_free(gen_U);
//   matrix_free(gen_V);
//   key_free(keys);
//   return EXIT_SUCCESS;
// }
