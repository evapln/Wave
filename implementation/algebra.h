#define ORDER 3

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

void prng_init(unsigned int seed);

/* addition dans le corps Fq */
char add_Fq(char a, char b);

/* multiplication dans le corps Fq */
char mul_Fq(char a, char b);

/* inversion dans le corps Fq */
char inv_Fq(char a);

/* renvoi un nombre aléatoire modulo q */
char rand_Fq(void);
