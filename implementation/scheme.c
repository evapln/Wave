#include "scheme.h"

const int SIZE = floor(LAMBDA/0.0154)+((int)(floor(LAMBDA/0.0154))%2);
const int OMEGA = 0.9261*SIZE;
const int d = 3;
const int K_U = 0.7978*SIZE/2;
const int K_V = 0.4201*SIZE/2;
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

keys_t *keys_alloc(int dim_U, int dim_V, int dim){
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

void keys_free(keys_t *keys) {
  if (keys) {
    matrix_free(keys->pk);
    sk_free(keys->sk);
  }
  free(keys);
  matrix_free(A);
  matrix_free(B);
  matrix_free(C);
  matrix_free(D);
  matrix_free(H);
  matrix_free(gen_U);
  matrix_free(gen_V);
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

void sk_copy(sk_t *dest,const sk_t *src) {
  if(!src || !dest)
    return;
  matrix_copy2(dest->parite_U, src->parite_U);
  matrix_copy2(dest->parite_V, src->parite_V);
  matrix_copy2(dest->S, src->S);
  matrix_copy2(dest->P, src->P);
  dest->dim_U = src->dim_U;
  dest->dim_V = src->dim_V;
}

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// affichage ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* Ecrit keys->sk dans secret et keys->pk dans public */
void keys_print(const keys_t *keys, FILE *secret, FILE *public) {
  if (keys && keys->pk && keys->sk && keys->sk->parite_U && keys->sk->parite_V && keys->sk->S && keys->sk->P) {
    puts("génération des clés de chiffrement avec les paramètres :");
    printf("lambda = %f, n = %d, w = %d, Ku = %d, Kv = %d, d = %d\n", LAMBDA, SIZE, OMEGA, K_U, K_V, d);
    fputs("parité U :\n", secret); matrix_print(keys->sk->parite_U, secret);
    fputs("parité V :\n", secret); matrix_print(keys->sk->parite_V, secret);
    fputs("S :\n", secret); matrix_print(keys->sk->S, secret);
    fputs("P :\n", secret); matrix_print(keys->sk->P, secret);
    fprintf(secret, "dim U : %d\n", keys->sk->dim_U);
    fprintf(secret, "dim V : %d\n\n", keys->sk->dim_V);
    fputs("Pk :\n", public); matrix_print(keys->pk, public);
  }
}

void sign_print(const sign_t *signature, FILE *file) {
  if (signature && signature->e && signature->r) {
    fputs("signature :",file);
    fprintf(file, "e :\n"); matrix_print(signature->e, file);
    fprintf(file, "r :\n"); matrix_print(signature->r, file);
  }
}


////////////////////////////////////////////////////////////////////////////////
///////////////////// gestion des paramètres de keys_t et sk_t /////////////////
////////////////////////////////////////////////////////////////////////////////

void keys_get_sk(sk_t *sk, const keys_t *keys) {
  if (!keys)
    return;
  sk_copy(sk,keys->sk);
}

void keys_get_pk(matrix_t *pk, const keys_t *keys) {
  if (!keys)
    return;
  matrix_copy2(pk, keys->pk);
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

int m1(matrix_t *x) {
  int m1 = 0;
  if (!x || matrix_get_row(x) != 1 || matrix_get_col(x) != SIZE)
    return -1;
  for (int i = 0; i < SIZE/2; ++i) {
    char xi = matrix_get_cell(x,0,i), xin = matrix_get_cell(x,0,i+SIZE/2);
    if ((xi == 0 && xin != 0) || (xi != 0 && xin == 0))
      ++m1;
  }
  return m1;
}

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// calcul du rejet ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


int binom(int k, int n) {
  if (n < k)
    return 0;
  if (n == k)
    return 1;
  long long num = 1;
  long long den = 1;
  if (k < n-k) {
    for (int i = 0; i < k; ++i) {
      num *= (n-i);
      den *= (k-i);
    }
  } else {
    for (int i = 1; i <= n-k; ++i) {
      num *= k+i;
      den *= i;
    }
  }
  return num/den;
}

double reject(int s, int t) {
  double pu = proba_unif(s,t);
  double p = proba(s,t);
  double m = max_proba(t);
  if (p == 0 || m == 0)
    return 0;
  return ((pu/p)/m);
}

double proba_unif(int s, int t) {
  if ((OMEGA + s) % 2 == 1)
    return 0;
  // printf("s: %d, t : %d\ts parmi t :%d %d\n",s,t,binom(s,t), binom((OMEGA+s)/2-t,SIZE/2-t));
  double num = binom(s,t)*binom((OMEGA+s)/2-t,SIZE/2-t)*2^(3*s/2);
  double den = 0;
  for (int p = 0; p <= t; ++p)
    den += binom(p,t)*binom((OMEGA+p)/2-t,SIZE/2-t)*2^(3*p/2);
  // printf("num = %lf, den = %lf\n", num , den);
  double div = num/den;
  // printf("div = %lf\n", div);
  return div;
}

double proba(int s, int t) {
  if ((OMEGA - s) % 2 == 1)
    return 0;
  double prob = 0, num, den;
  int k;
  for (int i = t + K_U - d -SIZE/2; i <= t; ++i) {
    k = K_U - d - i;
    num = binom(s,t-i)*binom((OMEGA+s)/2-t-k,SIZE/2-t-k)*2^(3*s/2);
    den = 0;
    for (int p = 0; p <= t-i; ++p)
      den += binom(p,t-i)*binom((OMEGA+p)/2-t-k,SIZE/2-t-k)*2^(3*p/2);
    prob += num/(den*K_U);
  }
  return prob;
}

double max_proba(int t) {
  double max = 0, quo;
  for (int s = 0; s <= t; ++s) {
    quo = proba_unif(s,t)/proba(s,t);
    if (quo > max)
      max = quo;
  }
  return max;
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
  keys_t *keys = keys_alloc(dim_U,dim_V, DIM);
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
    keys_free(keys);
    return NULL;
  }

  // création de H
  H = parite(keys->sk->parite_U, keys->sk->parite_V);
  if (!H) {
    keys_free(keys);
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
    keys_free(keys);
    return NULL;
  }
  // créationde P matrice de permutation aléatoire
  matrix_t *permut = matrix_perm_random(SIZE);
  if (!permut) {
    keys_free(keys);
    return NULL;
  }
  matrix_copy2(keys->sk->P, permut);
  matrix_free(permut);
  if (!keys->sk->P) {
    keys_free(keys);
    return NULL;
  }
  keys->sk->dim_U = dim_U;
  keys->sk->dim_V = dim_V;

  // création de la clé publique SHP
  matrix_t *SH = matrix_prod(keys->sk->S,H);
  if (!SH) {
    keys_free(keys);
    return NULL;
  }
  matrix_t *SHP = matrix_prod(SH,keys->sk->P);
  matrix_free(SH);
  if (!SHP) {
    keys_free(keys);
    return NULL;
  }
  matrix_copy2(keys->pk, SHP);
  matrix_free(SHP);
  if (!keys->pk) {
    keys_free(keys);
    return NULL;
  }
  return keys;
}

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// décodage /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void decode_ev(matrix_t * ev,const matrix_t *G, const matrix_t *synd) {
  // vérifications des entreés
  if (!ev || !G || !synd)
    return;
  int col = matrix_get_col(G);
  if (matrix_get_col(ev) != col || matrix_get_row(ev) != 1)
    return;
  // choix d'un mot prix au hasard dans le code
  random_word(ev, G);
  // padage du syndrome avec des 0
  matrix_t *s = matrix_alloc(1,col);
  int col_synd = matrix_get_col(synd);
  for (int i = 0; i < col - col_synd; ++i)
    matrix_set_cell(s, 0, i, 0);
  for (int i = col - col_synd; i < col; ++i)
    matrix_set_cell(s, 0, i, matrix_get_cell(synd, 0, i - col + col_synd));
  // addition du mot aléatoire du code avec le syndrome paddé
  matrix_add_modified(ev, s, 1);
  matrix_free(s);
}

int decode_eu(matrix_t * eu, const sk_t *sk, const matrix_t *synd,
               const matrix_t *ev) {
  // Vérification des entrées
  if (!eu || !synd || !ev || !sk)
    return 1;
  if (matrix_get_row(ev) != 1 || matrix_get_row(eu) != 1)
    return 1;
  // déclarations
  int dim_U = sk->dim_U;
  int a, b, c, d, w, r;
  int cpt_boucle;
  char no1, no2, evi, eui = 0, eui_not;
  matrix_t *eu_int = NULL;
  matrix_t *x = NULL;
  matrix_t *e = NULL;
  matrix_t *eu_not = NULL;
  // création des ensembles I et J
  int ens[SIZE/2];
  for (int i = 0; i < SIZE/2; ++i)
    ens[i] = i;
  int cpt = 0;
  do {
    ++cpt;
    matrix_free(eu_int); // en cas de repassage dans la boucle
    // création de I aléatoirement
    shuffle(ens,SIZE/2);
    int len_I = rand() % (SIZE/2 - dim_U) + dim_U; // aléa entre dim_U et SIZE/2
    int ens_I[len_I];
    for (int i = 0; i < len_I; ++i)
      ens_I[i] = ens[i];
    cpt_boucle = 0;
    do {
      eu_not = matrix_alloc(1,SIZE/2);
      if (!eu_not)
        return 1;
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
      ++cpt_boucle;
      if (cpt_boucle  > 20) {
        break;
      }
    } while(!eu_int);
    if (!eu_int)
      continue;
    e = phi(eu_int, ev);
    w = weight(e);
    matrix_free(e);
    if (cpt > 100000) {
      matrix_free(eu_int);
      return 1;
    }
  } while (w != OMEGA);
  matrix_copy2(eu, eu_int);
  matrix_free(eu_int);
  return 0;
}

int decode_uv(matrix_t *er, const sk_t *sk, const matrix_t *synd) {
  // vérifiaction des entrées
  if(!sk || !synd || !er)
    return -1;

  //initialisation
  matrix_t *su = matrix_alloc(1,SIZE/2 - K_U);
  matrix_t *sv = matrix_alloc(1,SIZE/2 - K_V);
  matrix_separate(synd, su, sv);
  if (!su || !sv)
    return -1;

  // decodage de ev
  matrix_t *ev = NULL;
  decodeV:
  ev = matrix_alloc(1,SIZE/2);
  decode_ev(ev, gen_V, sv);
  if (!ev)
    return -1;
  matrix_t *verifv = syndrome(ev, sk->parite_V);
  if (!verifv)
    return -1;
  matrix_free(verifv);
  int w = weight(ev);

  // decodage de eu
  matrix_t *eu = NULL;
  matrix_t *verifu = NULL;
  matrix_t *e = NULL;
  int rej = -1;
  double rejet;
  double r;
  prng_init(time(NULL) + getpid());
  int tem;
  do {
    ++ rej;
    matrix_free(eu);
    matrix_free(verifu);
    matrix_free(e);
    eu = matrix_alloc(1,SIZE/2);
    tem = decode_eu(eu, sk, su, ev);
    if (tem == 1) {
      matrix_free(eu);
      matrix_free(ev);
      goto decodeV;
    }
    if (!eu)
      return -1;
    verifu = syndrome(eu, sk->parite_U);
    if (!verifu)
      return -1;
    e = phi(eu,ev);
    printf("eu : "); matrix_print(eu,stdout);
    printf("e  : "); matrix_print(e,stdout);
    rejet = reject(m1(e),w);
    printf("valeur rejet : %f\n",rejet);
    r = (double)rand() / (double)RAND_MAX;
    printf("rand = %f\n", r);
  } while (r > rejet);
  matrix_free(ev);
  matrix_free(eu);
  matrix_free(verifu);
  matrix_free(su);
  matrix_free(sv);
  matrix_copy2(er, e);
  matrix_free(e);
  return rej;
}

int invert_alg(matrix_t *er, const sk_t *sk, const matrix_t *synd) {
  if (!sk || !synd || !er)
    return -1;
  matrix_t *S_inv = matrix_inv(sk->S);
  if (!S_inv)
    return -1;
  matrix_t *S_inv_T = matrix_trans(S_inv);
  matrix_free(S_inv);
  if (!S_inv_T)
    return -1;
  matrix_t *entree = matrix_prod(synd, S_inv_T);
  matrix_free(S_inv_T);
  if (!entree)
    return -1;
  matrix_t *e = matrix_alloc(1,SIZE);
  int rej = decode_uv(e, sk, entree);
  matrix_free(entree);
  if (!e)
    return -1;
  matrix_t *eP = matrix_prod(e,sk->P);
  matrix_free(e);
  matrix_copy2(er,eP);
  matrix_free(eP);
  return rej;
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
      matrix_free(left_inv);
      matrix_free(P);
      matrix_free(right);
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
    matrix_free(left_inv);
    return NULL;
  }
  matrix_separate(x, zero, ep);
  matrix_free(zero);
  if (!ep) {
    matrix_free(ep);
    matrix_free(P);
    matrix_free(right);
    matrix_free(left_inv);
    return NULL;
  }
  // calcul de e
  matrix_t *right_T = matrix_trans(right); // transposée de B
  matrix_free(right);
  if (!right_T) {
    matrix_free(ep);
    matrix_free(P);
    matrix_free(left_inv);
    return NULL;
  }
  matrix_t *epright_T = matrix_prod(ep,right_T); // e'*Bt
  matrix_free(right_T);
  if (!epright_T) {
    matrix_free(ep);
    matrix_free(P);
    matrix_free(left_inv);
    return NULL;
  }
  matrix_t *s = matrix_copy(syndrome);
  if (!s) {
    matrix_free(ep);
    matrix_free(P);
    matrix_free(left_inv);
    return NULL;
  }
  matrix_add_modified(s, epright_T, -1); // s - e'*Bt
  matrix_free(epright_T);
  if (!s) {
    matrix_free(ep);
    matrix_free(P);
    matrix_free(left_inv);
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
  matrix_free(m_high);
  return s;
}

int sign(sign_t *signature, const keys_t *keys, const matrix_t *m) {
  // vérifiations des entrées
  if (!keys || !m || !signature || !signature->e || !signature->r)
    return -1;
  // tirage de r aléatoirement de taille SIZE
  matrix_t *r = matrix_random(1,SIZE);
  if (!r)
    return -1;
  // hachage de (m,r)
  matrix_t *s = hash(m,r, keys->pk);
  if (!s) {
    matrix_free(r);
    return -1;
  }
  // inversion de f
  matrix_t *e = matrix_alloc(1,SIZE);
  int rej = invert_alg(e, keys->sk, s);
  matrix_free(s);
  if (!e) {
    matrix_free(r);
    return -1;
  }
  // allocation de la structure sign_t
  // if (!signature) {
  //   matrix_free(r);
  //   matrix_free(e);
  //   return NULL;
  // }
  matrix_copy2(signature->e, e);
  matrix_copy2(signature->r, r);
  matrix_free(e);
  matrix_free(r);
  return rej;
}

bool verify(const keys_t *keys, const matrix_t *m, const sign_t *signature) {
  // vérifiations des entrées
  if (!keys || !keys->pk || !m || !signature || !signature->r || !signature->e)
    return false;
  // calcul du poids de e
  if (weight(signature->e) != OMEGA)
    return false;
  // hachage de (m,r)
  matrix_t *s = hash(m,signature->r, keys->pk);
  if (!s)
    return false;
  // calcul du syndrome de e
  matrix_t *synd = syndrome(signature->e, keys->pk);
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
