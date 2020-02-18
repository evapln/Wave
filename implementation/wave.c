#include "scheme.h"

#include <err.h>
#include <getopt.h>

int main(int argc, char *argv[]) {
  // à faire autrement, lecture dans le lancement de l'executable
  int optc;

  /* long options struct */
  const struct option long_opts[] =
  {
    {"help",     no_argument,       NULL, 'h'},
    {NULL,       0,                 NULL,  0 },
  };

  char *opt = "h";
  while ((optc = getopt_long(argc, argv, opt, long_opts, NULL)) != -1) {
    switch (optc) {
      case 'h':
        puts("Usage :\twave [-h]\n"
            "A signature scheme based on Codes\n"
            "\t-h,--help\t\tdisplay this help\n\n");
        return EXIT_SUCCESS;

      default:
        errx(EXIT_FAILURE, "error: invalid option '%s'!\n",
              argv[1]);
    }
  }
  return EXIT_SUCCESS;
}
