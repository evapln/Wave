#include "scheme.h"
#include "wave.h"

#include <err.h>
#include <getopt.h>

int main(void) {
  printf("Wave : A signature scheme based on Codes with secutiry %f\n\n", LAMBDA);
  ///////////////////////////////// GENERATION DE CLES
  keys_t *keys = key_gen(1);
  if (keys == NULL)
    errx(EXIT_FAILURE, "error: failed to generate keys\n");
  keys_print(keys, stdout, stdout);
  sk_t * sk = keys_get_sk(keys);
  matrix_t *pk = keys_get_pk(keys);
  if (!pk || !sk)
    errx(EXIT_FAILURE, "error: failed to generate keys\n");

  //////////////////////////////// choix du message et signature
  matrix_t *m = vector_rand(10);
  printf("m : "); matrix_print(m, stdout);
  sign_t *signature = sign(keys, m);
  if(!signature)
    errx(EXIT_FAILURE, "error: failed to generate signature\n");
  sign_print(signature, stdout);
  printf("la signature est correct : %d\n", verify(pk, m, signature));

  //////////////////////////////// CLEAN UP
  matrix_free(m);
  sign_free(signature);
  // matrix_free(pk);
  // sk_free(sk);
  keys_free(keys);
  return EXIT_SUCCESS;
}
// int main(int argc, char *argv[]) {
//   int optc;
//   FILE *public = stdout;
//   FILE *secret = stdout;
//   bool generate = false;
//   bool sign = false;
//   bool verify = false;
//   float lambda = 3;
//   int mode = 0;
//   /* long options struct */
//   const struct option long_opts[] =
//   {
//     {"generate",      optional_argument, NULL, 'g'},
//     {"mode",          required_argument, NULL, 'm'},
//     {"public_output", required_argument, NULL, 'P'},
//     {"secret_output", required_argument, NULL, 'S'},
//     {"sign",          required_argument, NULL, 's'},
//     {"verify",        required_argument, NULL, 'v'},
//     {"help",          no_argument,       NULL, 'h'},
//     {NULL,            0,                 NULL,  0 },
//   };
//
//   char *opt = "g::m:po:so:s:v:h";
//   while ((optc = getopt_long(argc, argv, opt, long_opts, NULL)) != -1) {
//     switch (optc) {
//       case 'P':
//         if (optarg != NULL) {
//           fclose(public);
//           public = fopen(optarg, "w");
//           /* opening test in write mode */
//           if (public == NULL)
//             errx(EXIT_FAILURE, "error: invalid file : check the permission to access\n");
//         }
//         break;
//       case 'S':
//         if (optarg != NULL) {
//           fclose(secret);
//           secret = fopen(optarg, "w");
//           /* opening test in write mode */
//           if (secret == NULL)
//             errx(EXIT_FAILURE, "error: invalid file : check the permission to access\n");
//         }
//         break;
//       case 'g':
//         generate = true;
//         if (optarg != 0) {
//           lambda = atof(optarg);
//           if (lambda > SEC_MAX)
//             errx(EXIT_FAILURE, "error : I can only generate keys with a security of %d\n", SEC_MAX);
//         }
//         break;
//       case 'm':
//         if (optarg != 0) {
//           mode = atoi(optarg);
//           if (mode != 0 && mode != 1)
//             errx(EXIT_FAILURE, "error : the mode have to be 0 or 1");
//         }
//         break;
//       case 's':
//         sign = true;
//         break;
//       case 'v':
//         verify = true;
//         break;
//       case 'h':
//         printf("Usage :\twave -g[SECURITY] [-P FILE |-S FILE|-h]\n"
//             "\twave -s FILE [-h]\n"
//             "\twave -v FILE [-h]\n\n"
//             "A signature scheme based on Codes with max secutiry %d\n\n"
//             "-g[K], --generate[=K]\t\t\tgenerate wave's public and "
//             "private keys with K bits of security\n"
//             "-S FILE, --secret_output FILE\t\twrite private key in FILE\n"
//             "-P FILE, --public_output FILE\t\twrite public key in FILE\n\n"
//             "-m[M], --mode[=M]\t\t\tswitch mode in key generation : \n"
//             "\t\t\t\t\tM = 0 -> a=c=d=1, b = 0\n"
//             "\t\t\t\t\tM = 1 -> a,b,c,d random\n"
//             "-s PRIVATE_KEY, --sign PRIVVATE_KEY\tsign messages with the "
//             "private key entered\n"
//             "-v PUBLIC_KEY, --verify PUBLIC_KEY\tverfy messages with the "
//             "public key entered\n\n"
//             "-h,--help\t\t\t\tdisplay this help\n\n", SEC_MAX);
//         return EXIT_SUCCESS;
//
//       default:
//         errx(EXIT_FAILURE, "error: invalid option '%s'!\n",
//               argv[1]);
//     }
//   }
//   if (!generate && !sign && !verify)
//     errx(EXIT_FAILURE, "error: you have to choose a mode : generate, sign or verify!\n");
//   if (generate) {
//     keys_t *keys = key_gen(mode);
//     if (!keys)
//       errx(EXIT_FAILURE, "error: failed to generate keys\n");
//     keys_print(keys, secret, public);
//     keys_free(keys);
//     return EXIT_SUCCESS;
//   }
//   if (sign) {
//     keys_t *keys = key_gen(mode);
//     if (!keys)
//       errx(EXIT_FAILURE, "error: failed to generate keys\n");
//     keys_print(keys, secret, public);
//     matrix_t *m = vector_rand(10);
//     printf("m : "); matrix_print(m, stdout);
//     sign_t *signature = sign(keys, m);
//     if(!signature)
//       errx(EXIT_FAILURE, "error: failed to generate signature\n");
//     sign_print(signature, stdout);
//     printf("la signature est correct : %d\n", verify(keys->pk, m, signature));
//     matrix_free(m);
//     sign_free(signature);
//     return EXIT_SUCCESS;
//   }
//   if (verify) {
//     return EXIT_SUCCESS;
//   }
//   return EXIT_SUCCESS;
// }
