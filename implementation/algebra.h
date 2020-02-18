#define ORDER 3

#include <time.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

void prng_init(unsigned int seed);

/* Addition dans le corps Fq */
char add_Fq(char a, char b);

/* Multiplication dans le corps Fq */
char mul_Fq(char a, char b);

/* Inversion dans le corps Fq */
char inv_Fq(char a);

/* Opposé dans le corps Fq */
char opp_Fq(void);

/* Renvoi un nombre aléatoire modulo q */
char rand_Fq(void);

int mod_Fq(char a);
