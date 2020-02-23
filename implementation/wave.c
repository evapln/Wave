#include "scheme.h"
#include "wave.h"

#include <err.h>
#include <getopt.h>

int NB_IT = 100;

int main(int argc, char *argv[]) {
  int optc;
  int iteration = 0;

  /* long options struct */
  const struct option long_opts[] =
  {
    {"generate",      no_argument, NULL, 'g'},
    {"reject",        no_argument, NULL, 'r'},
    {"help",          no_argument, NULL, 'h'},
    {NULL,            0,           NULL,  0 },
  };

  char *opt = "grh";
  while ((optc = getopt_long(argc, argv, opt, long_opts, NULL)) != -1) {
    switch (optc) {
      case 'g':
        iteration = 1;
        break;
      case 'r':
        iteration = NB_IT;
        break;
      case 'h':
        printf("Usage :\twave [-g][-h]\n"
                "\twave -r[-h]\n\n"
                "A signature scheme based on Codes with max secutiry \n\n"
                "-g, --generate\t\tgenerate wave's keys then a signature "
                "on a random message\n"
                "-r, --reject\t\tgenerate wave's keys, %d signatures on"
                " random messages and print the number of rejects\n"
                "-h, --help\t\tdisplay this help\n\n",NB_IT);
        return EXIT_SUCCESS;
    }
  }


  if (iteration == 0)
    errx(EXIT_FAILURE, "error: choose an option\nFor help tap ./wave -h\n");

  ///////////////////////////////// GENERATION DE CLES
  keys_t *keys = key_gen(1);
  if (keys == NULL)
    errx(EXIT_FAILURE, "error: failed to generate keys\n");
  // keys_print(keys, stdout, stdout);

  if (iteration == 1) {
    ////////////////////////////// choix du message et signature
    matrix_t *m = vector_rand(50);
    printf("m : "); matrix_print(m, stdout);
    sign_t *signature = sign_alloc();
    if(!signature)
      errx(EXIT_FAILURE, "error: failed to allocate signature\n");
    int rej = sign(signature, keys, m);
    if(!signature)
      errx(EXIT_FAILURE, "error: failed to generate signature\n");
    sign_print(signature, stdout);
    printf("Nous savons eu besoin de %d rejet(s)\nLa signature est correcte :"
           " %d\n\n", rej, verify(keys, m, signature));

    //////////////////////////////// CLEAN UP
    matrix_free(m);
    sign_free(signature);
    keys_free(keys);
    return EXIT_SUCCESS;
  }

  int reject[iteration];
  for (int i = 0; i < iteration; ++i) {
    printf("\n%d\n",i);
    //////////////////////////////// choix du message et signature
    matrix_t *m = vector_rand(50);
    printf("m : "); matrix_print(m, stdout);
    sign_t *signature = sign_alloc();
    if(!signature)
      errx(EXIT_FAILURE, "error: failed to allocate signature\n");
    int rej = sign(signature, keys, m);
    reject[i] = rej;
    if(!signature)
      errx(EXIT_FAILURE, "error: failed to generate signature\n");
    // sign_print(signature, stdout);
    printf("la signature est correct : %d\n\n", verify(keys, m, signature));

    //////////////////////////////// CLEAN UP
    matrix_free(m);
    sign_free(signature);
  }
  printf("nombre de rejets : ");
  for (int i = 0; i < iteration; ++i)
    printf("%d ", reject[i]);
  puts("");
  keys_free(keys);
  return EXIT_SUCCESS;
}
