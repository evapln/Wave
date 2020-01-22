#include <lapacke.h>
#include "wave.h"

int SIZE;
int OMEGA;


matrix_t* phi (const matrix_t* x,const matrix_t* y,const matrix_t* a,const matrix_t* b, const matrix_t* c, const matrix_t* d) {
  matrix_t *res = NULL;
  matrix_t *res_g = NULL;
  matrix_t *res_d = NULL;
  matrix_t *ax = NULL;
  matrix_t *by = NULL;
  matrix_t *cx = NULL;
  matrix_t *dy = NULL;
  ax = vect_scal(a,x);
  by = vect_scal(b,y);
  cx = vect_scal(c,x);
  dy = vect_scal(d,y);
  res_g = matrix_add(ax,by);
  res_d = matrix_add(cx,dy);
  res = matrix_concatenation (res_g, res_d);
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

int main(int argc, char **argv) {

  // NE PAS TOUCHER !! INIT DES vecteurs a,b,c,d de phi //
  matrix_t *A = NULL;
  A = matrix_init (1,4,(char)1);
  puts("A : \n");
  matrix_print(A, stdout);
  matrix_t *B = NULL;
  B = matrix_init (1,4,(char)0);
  puts("B : \n");
  matrix_print(B, stdout);
  matrix_t *C = NULL;
  C = matrix_init (1,4,(char)1);
  puts("C : \n");
  matrix_print(C, stdout);
  matrix_t *D = NULL;
  D = matrix_init (1,4,(char)1);
  puts("D : \n");
  matrix_print(D, stdout);
  // NE PAS TOUCHER !! INIT DES vecteurs a,b,c,d de phi //


  matrix_t *X = NULL;
  X = matrix_random(1,4);
  puts("X : \n");
  matrix_print(X, stdout);
  matrix_t *Y = NULL;
  Y = matrix_random(1,4);
  puts("Y : \n");
  matrix_print(Y, stdout);
  matrix_t *Phi = NULL;
  Phi = phi(X,Y,A,B,C,D);
  puts("Phi : \n");
  matrix_print(Phi, stdout);



  matrix_free(A);
  matrix_free(B);
  matrix_free(C);
  matrix_free(D);
  matrix_free(X);
  matrix_free(Y);
  matrix_free(Phi);
  return EXIT_SUCCESS;
}
