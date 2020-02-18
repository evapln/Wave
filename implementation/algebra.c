#include "algebra.h"


/* tableau des inverses pour q = 3 */
static char tab_inv[3] = {0,1,2};

void prng_init(unsigned int seed) {
  static bool seed_init = false;
  if (!seed_init)
  {
    srand(seed);
    seed_init = true;
  }
}

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

char opp_Fq(void) {
  return (ORDER - 1);
}

char rand_Fq(void) {
  prng_init(time(NULL) + getpid());
  char rd = rand() % ORDER;
  return rd;
}

int mod_Fq(char a) {
  if (a  == -1)
    return 2;
  return a % ORDER;
}
