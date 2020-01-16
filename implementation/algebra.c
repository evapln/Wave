#include "algebra.h"

static char tab_inv[3] = {0,1,2};

char add_Fq(char a, char b) {
  char sum = (a + b) % ORDER;
  if (sum < 0)
    sum += ORDER;
  return sum;
}

char mul_Fq(char a, char b) {
  char prod = (a * b) % ORDER;
  if (prod < 0)
    prod += ORDER;
  return prod;
}

char inv_Fq(char a) {
  return tab_inv[(int)a];
}
