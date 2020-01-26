#include "matrix.h"

struct matrix_t {
  int nb_col;
  int nb_row;
  char** mat;
};

matrix_t *matrix_alloc(int row, int col) {
  matrix_t *matrix = malloc(sizeof(matrix_t));
  if (!matrix)
    return NULL;
  matrix->nb_col = col;
  matrix->nb_row = row;
  matrix->mat = calloc(row, sizeof(char*));
  if (!matrix->mat)
  {
    free(matrix);
    return NULL;
  }
  for (int i = 0; i < row; i++)
  {
    matrix->mat[i] = calloc(col, sizeof(char));
    if (!matrix->mat[i])
    {
      for (int j = 0; j < i; ++j)
      {
         free(matrix->mat[j]);
      }
      free(matrix->mat);
      free(matrix);
      return NULL;
    }
  }
  return matrix;
}

void matrix_free (matrix_t *matrix) {
  if (matrix != NULL)
  {
    for (int i = 0; i < matrix->nb_row; ++i)
    {
      free(matrix->mat[i]);
    }
    free(matrix->mat);
    free(matrix);
  }
}

char matrix_get_cell(const matrix_t *matrix,const int row_val, const int col_val) {
  if (!matrix)
    return 2;
  int row = matrix->nb_row;
  int col = matrix->nb_col;
  if (row_val > row || col_val > col)
    return 2;
  return matrix->mat[row_val][col_val];
}

int matrix_get_row(const matrix_t *matrix) {
  if (!matrix)
    return 0;
  return matrix->nb_row;
}

int matrix_get_col(const matrix_t *matrix) {
  if (!matrix)
    return 0;
  return matrix->nb_col;
}

void matrix_set_cell(matrix_t *matrix, const int row_val, const int col_val, const char val) {
  if (!matrix)
    return;
  int row = matrix->nb_row;
  int col = matrix->nb_col;
  if (row_val >= row || col_val >= col)
    return;
  char tmp = val % ORDER;
  if (tmp < 0)
    tmp += ORDER;
  matrix->mat[row_val][col_val] = tmp;
}

matrix_t *matrix_init (const int row, const int col, const char val) {
  matrix_t *matrix = NULL;
  matrix = matrix_alloc (row, col);
  for (int i = 0; i < row; ++i)
    for (int j = 0; j < col; ++j)
      matrix->mat[i][j] = val;
  return matrix;
}

matrix_t *matrix_copy(const matrix_t *matrix) {
  if (!matrix)
    return NULL;
  matrix_t *copy = matrix_alloc(matrix->nb_row, matrix->nb_col);
  if (!copy)
    return NULL;
  for (int i = 0; i < copy->nb_row; ++i)
    for (int j = 0; j < copy->nb_col; ++j)
      copy->mat[i][j] = matrix->mat[i][j];
  return copy;
}

void matrix_copy2 (const matrix_t *matrix1,const matrix_t *matrix2) {
  if (!matrix1 || !matrix2)
    return;
  if (matrix1->nb_row != matrix2->nb_row || matrix1->nb_col != matrix2->nb_col)
    return;
  for (int i = 0; i < matrix1->nb_row; ++i)
    for (int j = 0; j < matrix1->nb_col; ++j)
      matrix1->mat[i][j] = matrix2->mat[i][j];
  return;
}

matrix_t *matrix_identity (int size) {
  matrix_t *id = NULL;
  id = matrix_alloc(size,size);
  for (int i = 0; i < size; i++){
    id->mat[i][i] = 1;
    for (int j = 0; j < size; j++){
      if (j != i)
        id->mat[i][j] = 0;
    }
  }
  return id;
}

matrix_t *matrix_vect_to_diag (const matrix_t *vect, const char val) {
  int size = vect->nb_col;
  matrix_t *diag = NULL;
  diag = matrix_alloc (size, size);
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++) {
      if (j == i)
        diag->mat[i][j] = mul_Fq(vect->mat[0][i], val);
      else
        diag->mat[i][j] = 0;
    }
  return diag;
}

matrix_t *matrix_random(const int row, const int col) {
  matrix_t *rand = NULL;
  rand = matrix_alloc(row,col);
  for (int i = 0; i < row; ++i)
    for (int j = 0; j < col; ++j)
      rand->mat[i][j] = rand_Fq();
  return rand;
}

void shuffle(int *array, const int n) {
  if (n > 1)
  {
    for (int i = 0; i < n - 1; i++)
    {
      prng_init(time(NULL) + getpid());
      int k = rand() % (n - i);
      int j = i + k;
      int t = array[j];
      array[j] = array[i];
      array[i] = t;
    }
  }
}

bool is_in_array(const int *array, const int len, const int val) {
    for (int i = 0; i < len; ++i)
      if (array[i] == val)
        return true;
    return false;
}

void shuffle_info(int *array, const int len_a, const int *info, const int len_i, const int dim) {
  // test longueurs correctes
  if (len_i > len_a || dim > len_a || dim < len_i)
    return;
  prng_init(time(NULL) + getpid());
  // initialisation du tableau contenant les dim dernières coordonnées pas encore mélangées
  int fin[dim];
  int i;
  for (i  = 0; i < len_i; ++i)
    fin[i] = info[i];
  int random;
  while (i < dim) {
    random = rand() % len_a;
    if (!is_in_array(fin, i, array[random])) {
      fin[i] = array[random];
      ++i;
    }
  }
  // initialisation du tableau contenant les indices de début de array privés de ceux de fin
  int len_deb = len_a - dim;
  int debut[len_deb];
  int j = 0;
  for (int i = 0; i < len_a; ++i) {
    if (!is_in_array(fin, dim, array[i])) {
      debut[j] = array[i];
      ++j;
    }
  }
  // mélange
  shuffle(fin, dim);
  shuffle(debut, len_deb);
  // concaténation
  for (int i = 0; i < len_deb; ++i)
      array[i] = debut[i];
  for (int i  = 0; i < dim; ++i)
      array[i + len_deb] = fin[i];
}

matrix_t *matrix_perm_random(const int n) {
  matrix_t *matrix = NULL;
  matrix = matrix_alloc(n,n);
  int indices[n];
  for (int i = 0; i < n; i++)
    indices[i] = i;
  for (int i = 0; i < 5; i++)
    shuffle(indices, n);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
      matrix->mat[j][i] = 0;
    matrix->mat[indices[i]][i] = 1;
  }
  return matrix;
}

matrix_t *matrix_perm_random_info(const int n, const int *info, const int len_i, const int dim) {
  matrix_t *perm = matrix_alloc(n,n);
  int indices[n];
  for (int i = 0; i < n; i++)
    indices[i] = i;
  for (int i = 0; i < 5; i++)
    shuffle_info(indices, n, info, len_i, dim);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
      perm->mat[j][i] = 0;
    perm->mat[indices[i]][i] = 1;
  }
  return perm;
}

matrix_t *matrix_trans(const matrix_t *matrix) {
  if (!matrix)
    return NULL;
  matrix_t *trans = matrix_alloc(matrix->nb_col, matrix->nb_row);
  if (!trans)
    return NULL;
  for (int i = 0; i < trans->nb_row; ++i)
    for (int j = 0; j < trans->nb_col; ++j)
      trans->mat[i][j] = matrix->mat[j][i];
  return trans;
}

matrix_t *matrix_com(matrix_t *A) {
  int col = A->nb_col;
  int row = A->nb_row;
  matrix_t *com = NULL;
  com = matrix_alloc(row,col);
  for (int i = 0; i < row; i++)
    for (int j = 0; j < col; j++){
      matrix_t *sub = NULL;
      sub = matrix_sub(A,i,j);
      com->mat[i][j] = mul_Fq(pow(-1,i+j),matrix_det(sub));
      matrix_free(sub);
    }
  return com;
}

matrix_t *matrix_inv(matrix_t *A) {
  // calcul du determinant
  char det = matrix_det(A);
  if (det == 0) {
    puts("non inversible");
    return NULL;
  }
  // calcul de la comatrice
  matrix_t *com = NULL;
  com = matrix_com(A);
  // calcul de la transposé de la comatrice
  matrix_t *trans = NULL;
  trans = matrix_trans(com);
  // calcul final de l'inverse
  matrix_t *inv = NULL;
  inv = matrix_mul_by_scal(trans, inv_Fq(det));
  // clean up
  matrix_free(com);
  matrix_free(trans);
  return inv;
}

matrix_t *matrix_sub (const matrix_t *A, int a, int b) {
  matrix_t *sub = NULL;
  int row = A->nb_row;
  int col = A->nb_col;
  sub = matrix_alloc(row - 1, col - 1);
  int ind_row = 0;
  int ind_col = 0;
  for (int i = 0; i < row; i++){
    ind_col = 0;
    if (i != a){
      for (int j = 0; j < col; j++){
        if (j != b){
          sub->mat[ind_row][ind_col] = A->mat[i][j];
          ind_col++;
        }
      }
      ind_row++;
    }
  }
  return sub;
}

void matrix_separate(const matrix_t *matrix, matrix_t *A, matrix_t *B) {
  int row = matrix->nb_row;
  int col_A = A->nb_col;
  int col_B = B->nb_col;
  // erreur dimensions
  if (A->nb_row != row || B->nb_row != row || col_A + col_B != matrix->nb_col) {
    A = NULL;
    B = NULL;
    return;
  }
  // initialisation des matrices
  for (int i = 0; i < row; ++i) {
    // matrice A
    for (int j = 0; j < col_A; ++j)
      A->mat[i][j] = matrix->mat[i][j];
    // matrice B
    for (int j = 0; j < col_B; ++j)
      B->mat[i][j] = matrix->mat[i][j+col_A];
  }
}

matrix_t *matrix_concatenation(const matrix_t *A, const matrix_t *B, const int mode) {
  matrix_t *conc = NULL;
  if (mode == 0){
    int row = A->nb_row;
    if (row != B->nb_row)
      return NULL;
    int col_A = A->nb_col;
    int col_B = B->nb_col;
    int col = col_A + col_B;
    conc = matrix_alloc(row,col);
    for (int i = 0; i < row; i++){
      for (int j = 0; j < col_A; j++)
        conc->mat[i][j] = A->mat[i][j];
      for (int j = col_A; j < col; j++)
        conc->mat[i][j] = B->mat[i][j-col_A];
    }
  }
  else{
    int col = A->nb_col;
    if (col != B->nb_col)
      return NULL;
    int row_A = A->nb_row;
    int row_B = B->nb_row;
    int row = row_A + row_B;
    conc = matrix_alloc(row,col);
    for (int i = 0; i < row_A; i++)
      for (int j = 0; j < col; j++)
        conc->mat[i][j] = A->mat[i][j];
    for (int i = row_A; i < row; i++)
      for (int j = 0; j < col; j++)
        conc->mat[i][j] = B->mat[i - row_A][j];
  }
  return conc;
}

matrix_t *matrix_add(const matrix_t *matrix1, const matrix_t *matrix2) {
  int row = matrix1->nb_row;
  if (matrix2->nb_row != row)
    return NULL;
  int col = matrix1->nb_col;
  if (matrix2->nb_col != col)
    return NULL;
  matrix_t *sum = NULL;
  sum = matrix_alloc(row, col);
  for (int i = 0; i < row; ++i)
    for (int j = 0; j < col; ++j)
      sum->mat[i][j] = add_Fq(matrix1->mat[i][j], matrix2->mat[i][j]);
  return sum;
}

int sub_weight(const matrix_t *vect, const int *subset, const int len_s) {
  int w = 0;
  for (int i = 0; i < len_s; ++i)
    if (vect->mat[0][subset[i]] % ORDER != 0)
      ++w;
  return w;
}

matrix_t *vector_rand_weight(const int n, const int *info, const int len_i, const int t) {
  // test valeur rentées
  if (t > len_i)
    return NULL;
  // allocation vecteur
  matrix_t *vect = matrix_alloc(1,n);
  // remplissage des coordonnées sans contraintes à rand_Fq et celles avec contraintes à 0
  for (int i = 0; i < n; ++i) {
    if (!is_in_array(info, len_i, i))
      vect->mat[0][i] = rand_Fq();
    else
      vect->mat[0][i] = 0;
  }
  // remplissage des coordonnées avec contraintes avec des rand_Fq
  int w;
  char r;
  for (int i = 0; i < len_i; ++i) {
    w = sub_weight(vect, info, len_i);
    // si le poids est déjà atteint, remplissage avec des 0
    if (w == t) {
      vect->mat[0][info[i]] = 0;
    } else {
      // si il reste juste assez de cases vides pour obtenir le poids souhaité, remplissage avec des valeurs non nulles
      if (len_i - i <= t - w) {
        r = 0;
        while (r == 0)
          r = rand_Fq();
        vect->mat[0][info[i]] = r;
      } else {
        // aucun des deux cas précédents, on remplie aléatoirement
        vect->mat[0][info[i]] = rand_Fq();
      }
    }
  }
  return vect;
}

matrix_t *vect_scal(const matrix_t *vect1, const matrix_t *vect2) {
  int row = vect1->nb_row;
  if (vect2->nb_row != row)
    return NULL;
  int col = vect1->nb_col;
  if (vect2->nb_col != col)
    return NULL;
  matrix_t *mul = NULL;
  mul = matrix_alloc(row,col);
  for (int i = 0; i < row; ++i)
    for (int j = 0; j < col; ++j)
      mul->mat[i][j] = mul_Fq(vect1->mat[i][j], vect2->mat[i][j]);
  return mul;
}

matrix_t *matrix_mul_by_scal(const matrix_t *matrix, const char scal) {
  int row = matrix->nb_row;
  int col = matrix->nb_col;
  matrix_t *mul = NULL;
  mul = matrix_alloc(row,col);
  for (int i = 0; i < row; ++i)
    for (int j = 0; j < col; ++j)
      mul->mat[i][j] = mul_Fq(matrix->mat[i][j], scal);
  return mul;
}

matrix_t *matrix_prod(const matrix_t *matrix1, const matrix_t *matrix2) {
  int k = matrix1->nb_col;
  if (k != matrix2->nb_row)
    return NULL;
  int row = matrix1->nb_row;
  int col = matrix2->nb_col;
  matrix_t *prod = NULL;
  prod = matrix_alloc(row,col);
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      prod->mat[i][j] = 0;
      for (int ind = 0; ind < k; ind++)
        prod->mat[i][j] = add_Fq(mul_Fq(matrix1->mat[i][ind], matrix2->mat[ind][j]), prod->mat[i][j]);
    }
  }
  return prod;
}

matrix_t *matrix_del_row(matrix_t *matrix, const int row1) {
  int row = matrix->nb_row;
  int col = matrix->nb_col;
  if (row1 > row)
    return NULL;
  matrix_t *matrix2 = NULL;
  matrix2 = matrix_alloc(row - 1, col);
  int id_row = 0;
  for (int i = 0; i < row; ++i){
    if (i != row1){
      for (int j = 0; j < col; ++j){
        matrix2->mat[id_row][j] = matrix->mat[i][j];
      }
      id_row++;
    }
  }
  // matrix_free(matrix);
  return matrix2;
}

bool row_is_zero (matrix_t *matrix, const int row) {
  for (int col = 0; col < matrix->nb_col; col++){
    if (matrix->mat[row][col] != 0) {
      return false;
    }
  }
  return true;
}

void matrix_add_row(matrix_t *matrix, const int row1, const int row2, const char coef) {
  int row = matrix->nb_row;
  if (row1 > row || row2 > row)
    return;
  for (int i = 0; i < matrix->nb_col; ++i)
    matrix->mat[row1][i] = add_Fq(matrix->mat[row1][i], mul_Fq(coef,matrix->mat[row2][i]));
}

void matrix_mul_row(matrix_t *matrix, const int row1, const char coef) {
  if (row1 > matrix->nb_row)
    return;
  for (int i = 0; i < matrix->nb_col; ++i)
    matrix->mat[row1][i] = mul_Fq(coef,matrix->mat[row1][i]);
}

void matrix_exchange_row(matrix_t *matrix, const int row1, const int row2) {
  int row = matrix->nb_row;
  if (row1 > row || row2 > row)
    return;
  char tmp;
  for (int i = 0; i < matrix->nb_col; ++i) {
    tmp = matrix->mat[row1][i];
    matrix->mat[row1][i] = matrix->mat[row2][i];
    matrix->mat[row2][i] = tmp;
  }
}

matrix_t *matrix_del_null_row (matrix_t *matrix) {
  int row = matrix->nb_row;
  int tab_row[row];
  for (int i = 0; i < row; ++i)
    tab_row[i] = 0;
  int nb_row = 0;
  for (int i = 0; i < matrix->nb_row; ++i) {
    if (!row_is_zero(matrix, i)) {
      tab_row[nb_row] = i;
      ++nb_row;
    }
  }
  int col = matrix->nb_col;
  matrix_t *matrix2 = NULL;
  matrix2 = matrix_alloc(nb_row, col);
  for (int i = 0; i < nb_row; ++i) {
    for (int j = 0; j < col; ++j) {
      matrix2->mat[i][j] = matrix->mat[tab_row[i]][j];
    }
  }
  return matrix2;
}

matrix_t *matrix_parite(const matrix_t *gen) {
  matrix_t *copy_gen = NULL;
  copy_gen = matrix_copy(gen);
  if (!copy_gen)
    return NULL;
  matrix_systematisation(copy_gen);
  matrix_t *syst = NULL;
  syst = matrix_del_null_row(copy_gen);
  matrix_free(copy_gen);
  if (!matrix_is_syst(syst)) {
    matrix_free(syst);
    return NULL;
  }
  int row_syst = syst->nb_row;
  int col_syst = syst->nb_col;
  int row_p = col_syst - row_syst;
  int col_p = row_syst + row_p;
  matrix_t *parite = NULL;
  parite = matrix_alloc(row_p,col_p);
  for (int i = 0; i < row_p; ++i) {
    for (int j = 0; j < row_syst; ++j)
      parite->mat[i][j] = syst->mat[j][row_syst + i];
    for (int j = row_syst; j < col_p; ++j) {
      if (j - row_syst == i)
        parite->mat[i][j] = 1;
      else
        parite->mat[i][j] = 0;
    }
  }
  matrix_free(syst);
  return parite;
}

void matrix_systematisation(matrix_t *matrix) {
  int row = matrix->nb_row;
  int col = matrix->nb_col;
  // étape 1
  // pour toutes les colonnes
  for (int j = 0; j < col; ++j) {
    // parcours de chaque colonne
    for (int i = j; i < row; ++i) {
      // on met un coefficient non nul en haut de chaque colonne
      if (i != j && matrix->mat[i][j] % ORDER != 0) {
        matrix_exchange_row(matrix, j, i);
        break;
      }
    }
    // annulation du reste de la colonne
    for (int i = j + 1; i < row; ++i) {
      char coef = mul_Fq(matrix->mat[i][j], inv_Fq(matrix->mat[j][j]));
      if (coef != 0) {
        matrix_add_row(matrix, i, j, -coef);
      }
    }
  }

  // on met tous les pivots à 1
  // pour toutes les lignes
  for (int i = 0; i < row; ++i) {
    // on parcours les lignes
    for (int j = 0; j < col; ++j) {
      // on trouve le premier coef non nul, le pivot
      if (matrix->mat[i][j] != 0) {
        matrix_mul_row(matrix, i, inv_Fq(matrix->mat[i][j]));
        break;
      }
    }
  }
  for (int i = row - 1; i >= 0; --i) {
    // parcours de la ligne
    int j = 0;
    while ((matrix->mat[i][j] == 0) && (j != col - 1))
      j += 1;
    for (int k = 0; k < i; ++k) {
      if (matrix->mat[k][j] != 0)
        matrix_add_row(matrix, k, i, -inv_Fq(matrix->mat[k][j]));
    }
  }
}

bool matrix_is_syst (matrix_t *matrix) {
  int row = matrix->nb_row;
  int col = matrix->nb_col;
  for (int i = 0; i < row; i++)
    for (int j = 0; j < col && j < row; j++){
      // char mat = matrix2->mat[i][j];
      if (i == j && (matrix->mat[i][j] != 1)){
        // matrix_free(matrix2);
        return false;
      }
      if (i != j && (matrix->mat[i][j] != 0)){
        // matrix_free(matrix2);
        return false;
      }
    }
  // matrix_free(matrix2);
  return true;
}

char matrix_det(matrix_t *A) {
  int size = A->nb_col;
  if (size == 2)
    return (add_Fq(mul_Fq(A->mat[0][0],A->mat[1][1]),-mul_Fq(A->mat[0][1],A->mat[1][0])));
  char det = 0;
  matrix_t *sub = NULL;
  for (int i = 0; i < size; ++i) {
    sub = matrix_sub(A,0,i);
    char sub_det = matrix_det(sub);
    char val1 = mul_Fq(pow(-1,i),A->mat[0][i]);
    det = add_Fq(det,mul_Fq(val1,sub_det));
    matrix_free(sub);
  }
  return det;
}

void matrix_print(matrix_t *matrix, FILE *fd) {
  if (matrix != NULL && matrix->mat != NULL)
  {
    int col = matrix->nb_col;
    int row = matrix->nb_row;
    for (int i = 0; i < row; ++i)
    {
      char val;
      for (int j = 0; j < col; ++j)
      {
        val = matrix_get_cell(matrix,i,j);
        if (val == 2 && ORDER == 3)
          fprintf(fd,"-1 ");
        else
          fprintf(fd," %d ", val);
      }
      fputs("\n", fd);
    }
    fputs("\n\n", fd);
  }
}
