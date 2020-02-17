#include "matrix.h"

struct matrix_t {
  int nb_col;
  int nb_row;
  char** mat;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////// gestion de la mémoire ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////

matrix_t *matrix_alloc(const int row, const int col) {
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
  if (matrix)
  {
    for (int i = 0; i < matrix->nb_row; ++i)
    {
      free(matrix->mat[i]);
    }
    free(matrix->mat);
    free(matrix);
  }
}



////////////////////////////////////////////////////////////////////////////////
////////////////////// gestion des paramètres de matrix_t //////////////////////
////////////////////////////////////////////////////////////////////////////////

char matrix_get_cell(const matrix_t *matrix,const int row_val, const int col_val) {
  if (!matrix)
    return 'e';
  int row = matrix->nb_row;
  int col = matrix->nb_col;
  if (row_val > row || col_val > col)
    return 'e';
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
  char tmp = val;
  if (tmp == '*') {
    char tmp = val % ORDER;
    if (tmp < 0)
      tmp += ORDER;
  }
  matrix->mat[row_val][col_val] = tmp;
}


////////////////////////////////////////////////////////////////////////////////
//////////////////////// création de matrices spécifiques //////////////////////
////////////////////////////////////////////////////////////////////////////////

matrix_t *matrix_init (const int row, const int col, const char val) {
  matrix_t *matrix = matrix_alloc (row, col);
  if (!matrix)
    return NULL;
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

void matrix_copy2 (const matrix_t *matrix1, const matrix_t *matrix2) {
  if (!matrix1 || !matrix2)
    return;
  if (matrix1->nb_row != matrix2->nb_row || matrix1->nb_col != matrix2->nb_col)
    return;
  for (int i = 0; i < matrix1->nb_row; ++i)
    for (int j = 0; j < matrix1->nb_col; ++j)
      matrix1->mat[i][j] = matrix2->mat[i][j];
  return;
}

matrix_t *matrix_identity (const int size) {
  matrix_t *id = matrix_alloc(size,size);
  if (!id)
    return NULL;
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
  if (!vect)
    return NULL;
  int size = vect->nb_col;
  matrix_t *diag = matrix_alloc (size, size);
  if (!diag)
    return NULL;
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
  matrix_t *rand = matrix_alloc(row,col);
  if (!rand)
    return NULL;
  for (int i = 0; i < row; ++i)
    for (int j = 0; j < col; ++j)
      rand->mat[i][j] = rand_Fq();
  return rand;
}

matrix_t *matrix_perm_random(const int n) {
  matrix_t *matrix = matrix_alloc(n,n);
  if (!matrix)
    return NULL;
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
  if (!perm)
    return NULL;
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

matrix_t *matrix_com(const matrix_t *A) {
  if (!A)
    return NULL;
  int col = A->nb_col;
  int row = A->nb_row;
  matrix_t *com = NULL;
  com = matrix_alloc(row,col);
  if (!com)
    return NULL;
  for (int i = 0; i < row; i++)
    for (int j = 0; j < col; j++){
      matrix_t *sub = matrix_sub(A,i,j);
      if (!sub)
        return NULL;
      com->mat[i][j] = mul_Fq(pow(-1,i+j),matrix_det(sub));
      matrix_free(sub);
    }
  return com;
}

matrix_t *matrix_inv(const matrix_t *A){
  if (!A)
    return NULL;
  int row = A->nb_row;
  int col = A->nb_col;
  if (col != row)
    return NULL;
  matrix_t *id = matrix_identity(col);
  if (!id)
    return NULL;
  matrix_t *conc = matrix_concatenation(A,id,0);
  if (!conc) {
    matrix_free(id);
    return NULL;
  }
  matrix_systematisation(conc);
  matrix_t *inv = matrix_alloc(row,col);
  if (!inv) {
    matrix_free(conc);
    matrix_free(id);
    return NULL;
  }
  matrix_separate(conc,id,inv);
  if (!inv || !conc || !id) {
    matrix_free(inv);
    matrix_free(conc);
    matrix_free(id);
    return NULL;
  }
  if (!is_identity(id)) {
    matrix_free(inv);
    matrix_free(conc);
    matrix_free(id);
    return NULL;
  }
  matrix_free(conc);
  matrix_free(id);
  return inv;
}

matrix_t *matrix_inv_com(const matrix_t *A) {
  if (!A)
    return NULL;
  // calcul du determinant
  char det = matrix_det(A);
  if (det == 0)
    return NULL;
  // calcul de la comatrice
  matrix_t *com = matrix_com(A);
  if (!com)
    return NULL;
  // calcul de la transposé de la comatrice
  matrix_t *trans = matrix_trans(com);
  if (!trans)
    return NULL;
  // calcul final de l'inverse
  matrix_t *inv = matrix_mul_by_scal(trans, inv_Fq(det));
  if (!inv)
    return NULL;
  // clean up
  matrix_free(com);
  matrix_free(trans);
  return inv;
}

matrix_t *matrix_sub (const matrix_t *A, const int a, const int b) {
  if (!A)
    return NULL;
  int row = A->nb_row;
  int col = A->nb_col;
  matrix_t *sub = matrix_alloc(row - 1, col - 1);
  if (!sub)
    return NULL;
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
  if (!matrix)
    return;
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
  if (!A || !B)
    return NULL;
  matrix_t *conc = NULL;
  if (mode == 0){
    int row = A->nb_row;
    if (row != B->nb_row)
      return NULL;
    int col_A = A->nb_col;
    int col_B = B->nb_col;
    int col = col_A + col_B;
    conc = matrix_alloc(row,col);
    if (!conc)
      return NULL;
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
    if (!conc)
      return NULL;
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
  if (!matrix1 || !matrix2)
    return NULL;
  int row = matrix1->nb_row;
  if (matrix2->nb_row != row)
    return NULL;
  int col = matrix1->nb_col;
  if (matrix2->nb_col != col)
    return NULL;
  matrix_t *sum = matrix_alloc(row, col);
  if (!sum)
    return NULL;
  for (int i = 0; i < row; ++i)
    for (int j = 0; j < col; ++j)
      sum->mat[i][j] = add_Fq(matrix1->mat[i][j], matrix2->mat[i][j]);
  return sum;
}

void matrix_add_modified(matrix_t *dest, const matrix_t *src, const char coef) {
  if (!dest || !src)
    return;
  int row = dest->nb_row;
  if (src->nb_row != row)
    return;
  int col = dest->nb_col;
  if (src->nb_col != col)
    return;
  for (int i = 0; i < row; ++i)
    for (int j = 0; j < col; ++j)
      dest->mat[i][j] = add_Fq(dest->mat[i][j], mul_Fq(src->mat[i][j], coef));
}

matrix_t *matrix_mul_by_scal(const matrix_t *matrix, const char scal) {
  if (!matrix)
    return NULL;
  int row = matrix->nb_row;
  int col = matrix->nb_col;
  matrix_t *mul = matrix_alloc(row,col);
  if (!mul)
    return NULL;
  for (int i = 0; i < row; ++i)
    for (int j = 0; j < col; ++j)
      mul->mat[i][j] = mul_Fq(matrix->mat[i][j], scal);
  return mul;
}

matrix_t *matrix_prod(const matrix_t *matrix1, const matrix_t *matrix2) {
  if (!matrix1 || !matrix2)
    return NULL;
  int k = matrix1->nb_col;
  if (k != matrix2->nb_row)
    return NULL;
  int row = matrix1->nb_row;
  int col = matrix2->nb_col;
  matrix_t *prod = matrix_alloc(row,col);
  if (!prod)
    return NULL;
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      prod->mat[i][j] = 0;
      for (int ind = 0; ind < k; ind++)
        prod->mat[i][j] = add_Fq(mul_Fq(matrix1->mat[i][ind], matrix2->mat[ind][j]), prod->mat[i][j]);
    }
  }
  return prod;
}

matrix_t *matrix_parite(const matrix_t *gen) {
  if (!gen)
    return NULL;
  matrix_t *copy_gen = matrix_copy(gen);
  if (!copy_gen)
    return NULL;
  matrix_systematisation(copy_gen);
  matrix_t *syst = matrix_del_null_row(copy_gen);
  if (!syst)
    return NULL;
  matrix_free(copy_gen);
  if (!matrix_is_syst(syst)) {
    matrix_free(syst);
    return NULL;
  }
  int row_syst = syst->nb_row;
  int col_syst = syst->nb_col;
  int row_p = col_syst - row_syst;
  int col_p = row_syst + row_p;
  matrix_t *parite = matrix_alloc(row_p,col_p);
  if (!parite)
    return NULL;
  for (int i = 0; i < row_p; ++i) {
    for (int j = 0; j < row_syst; ++j)
      parite->mat[i][j] = mul_Fq(-1,syst->mat[j][row_syst + i]);
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

matrix_t *sub_col_matrix(const matrix_t *A, const int *ind_col, const int len) {
  if (!A)
    return NULL;
  int row = A->nb_row;
  matrix_t *B = matrix_alloc(row,len);
  if (!B)
    return NULL;
  for (int i = 0; i < row; ++i) {
    for (int j = 0; j < len; ++j) {
      B->mat[i][j] = A->mat[i][ind_col[j]];
    }
  }
  return B;
}

void matrix_systematisation(matrix_t *matrix) {
  if (!matrix)
    return;
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


////////////////////////////////////////////////////////////////////////////////
///////////////////////////// calcul de determinant ////////////////////////////
////////////////////////////////////////////////////////////////////////////////

char matrix_det(const matrix_t *A) {
  if (!A)
    return 0;
  int size = A->nb_col;
  if (is_trigonalise(A)){
    // puts("ok");
    char det = 1;
    for (int i = 0; i < size; i++){
      det *= A->mat[i][i];
      if (det == 0)
        return 0;
    }
    return det;
  }
  // printf("size = %d",size);
  if (size == 2)
    return (add_Fq(mul_Fq(A->mat[0][0],A->mat[1][1]),-mul_Fq(A->mat[0][1],A->mat[1][0])));
  char det = 0;
  matrix_t *sub = NULL;
  for (int i = 0; i < size; ++i) {
    sub = matrix_sub(A,0,i);
    if (!sub)
      return 0;
    if (A->mat[0][i] == 0)
      continue;
    char sub_det = matrix_det(sub);
    char val1 = mul_Fq(pow(-1,i),A->mat[0][i]);
    det = add_Fq(det,mul_Fq(val1,sub_det));
    matrix_free(sub);
  }
  return det;
}

int matrix_rank(const matrix_t *A) {
  if (!A)
    return 0;
  matrix_t *B = matrix_copy(A);
  if (!B)
    return 0;
  matrix_systematisation(B);
  if (!B)
    return 0;
  matrix_t *C = matrix_del_null_row(B);
  if (!C)
    return 0;
  matrix_free(B);
  int rank = C->nb_row;
  matrix_free(C);
  return rank;
}

////////////////////////////////////////////////////////////////////////////////
//////////////////////////// test de la forme des matrices /////////////////////
////////////////////////////////////////////////////////////////////////////////

bool is_identity(const matrix_t *A){
  if (!A)
    return false;
  int row = A->nb_row;
  int col = A->nb_col;
  if (col != row)
    return false;
  for (int i = 0; i < row; i++)
    for (int j = 0; j < col; j++){
      if (i == j && A->mat[i][j] != 1)
        return false;
      if (i != j && A->mat[i][j] != 0)
        return false;
    }
  return true;
}

bool matrix_is_syst (const matrix_t *matrix) {
  if (!matrix)
    return false;
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

bool is_trigonalise (const matrix_t *A) {
  int size = A->nb_row;
  for (int i = 0; i < size; i++){
    for (int j = 0; j < i; j++){
      if (A->mat[i][j]%ORDER != 0)
        return false;
    }
  }
  return true;
}


////////////////////////////////////////////////////////////////////////////////
///////////////////////////// utilisation de tableaux //////////////////////////
////////////////////////////////////////////////////////////////////////////////

bool is_in_array(const int *array, const int len, const int val) {
  for (int i = 0; i < len; ++i)
    if (array[i] == val)
      return true;
  return false;
}

void shuffle(int *array, const int n) {
  if (n > 1) {
    for (int i = 0; i < n - 1; i++) {
      prng_init(time(NULL) + getpid());
      int k = rand() % (n - i);
      int j = i + k;
      int t = array[j];
      array[j] = array[i];
      array[i] = t;
    }
  }
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



////////////////////////////////////////////////////////////////////////////////
////////////////////////////// gestion de vecteurs /////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int weight(const matrix_t *vect) {
  if (!vect)
    return -1;
  int row = vect->nb_row;
  int col = vect->nb_col;
  if (row != 1 && col != 1)
    return -1; // ce n'est pas un vecteur
  int w = 0;
  if (col != 1) {
    for (int i = 0; i < col; ++i)
      if (vect->mat[0][i] % ORDER != 0)
        ++w;
  } else {
    for (int i = 0; i < row; ++i)
      if (vect->mat[i][0] % ORDER != 0)
        ++w;
  }
  return w;
}

int sub_weight(const matrix_t *vect, const int *subset, const int len_s) {
  if (!vect)
    return -1;
  int w = 0;
  for (int i = 0; i < len_s; ++i)
    if (vect->mat[0][subset[i]] % ORDER != 0)
      ++w;
  return w;
}

matrix_t *vector_rand(const int n) {
  matrix_t *vect = matrix_alloc(1,n);
  if (!vect)
    return NULL;
  for (int i = 0; i < n; ++i)
    vect->mat[0][i] = rand_Fq();
  return vect;
}

matrix_t *vector_rand_weight(const int n, const int w) {
  matrix_t *vect = NULL;
  int w_vect = 0;
  while (!vect || w_vect != w) {
    matrix_free(vect);
    vect = vector_rand(n);
    w_vect = weight(vect);
  }
  return vect;
}

matrix_t *vector_rand_sub_weight(const int n, const int *info, const int len_i, const int t) {
  matrix_t *vect = NULL;
  int w_vect = 0;
  while (!vect || w_vect != t) {
    matrix_free(vect);
    vect = vector_rand(n);
    w_vect = sub_weight(vect, info, len_i);
  }
  return vect;
}

matrix_t *vect_scal(const matrix_t *vect1, const matrix_t *vect2) {
  if (!vect1 || !vect2)
    return NULL;
  int row = vect1->nb_row;
  if (vect2->nb_row != row)
    return NULL;
  int col = vect1->nb_col;
  if (vect2->nb_col != col)
    return NULL;
  matrix_t *mul = NULL;
  mul = matrix_alloc(row,col);
  if (!mul)
    return NULL;
  for (int i = 0; i < row; ++i)
    for (int j = 0; j < col; ++j)
      mul->mat[i][j] = mul_Fq(vect1->mat[i][j], vect2->mat[i][j]);
  return mul;
}


////////////////////////////////////////////////////////////////////////////////
////////////////// calculs et tests sur les lignes d'une matrice ///////////////
////////////////////////////////////////////////////////////////////////////////

void matrix_row(matrix_t *ligne, const matrix_t *A, const int row) {
  if (!A)
    return;
  for (int i = 0; i < A->nb_col; ++i)
    ligne->mat[0][i] = A->mat[row][i];
}

matrix_t *matrix_del_row(const matrix_t *matrix, const int row1) {
  if (!matrix)
    return NULL;
  int row = matrix->nb_row;
  int col = matrix->nb_col;
  if (row1 > row)
    return NULL;
  matrix_t *matrix2 = matrix_alloc(row - 1, col);
  if (!matrix2)
    return NULL;
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

bool row_is_zero (const matrix_t *matrix, const int row) {
  if (!matrix)
    return false;
  for (int col = 0; col < matrix->nb_col; col++){
    if (matrix->mat[row][col] != 0) {
      return false;
    }
  }
  return true;
}

void matrix_add_row(matrix_t *matrix, const int row1, const int row2, const char coef) {
  if (!matrix)
    return;
  int row = matrix->nb_row;
  if (row1 > row || row2 > row)
    return;
  for (int i = 0; i < matrix->nb_col; ++i)
    matrix->mat[row1][i] = add_Fq(matrix->mat[row1][i], mul_Fq(coef,matrix->mat[row2][i]));
}

void matrix_mul_row(matrix_t *matrix, const int row1, const char coef) {
  if (!matrix)
    return;
  if (row1 > matrix->nb_row)
    return;
  for (int i = 0; i < matrix->nb_col; ++i)
    matrix->mat[row1][i] = mul_Fq(coef,matrix->mat[row1][i]);
}

void matrix_exchange_row(matrix_t *matrix, const int row1, const int row2) {
  if (!matrix)
    return;
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

matrix_t *matrix_del_null_row (const matrix_t *matrix) {
  if (!matrix)
    return NULL;
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
  matrix_t *matrix2 = matrix_alloc(nb_row, col);
  if (!matrix2)
    return NULL;
  for (int i = 0; i < nb_row; ++i) {
    for (int j = 0; j < col; ++j) {
      matrix2->mat[i][j] = matrix->mat[tab_row[i]][j];
    }
  }
  return matrix2;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////// gestion d'un code //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void random_word(matrix_t *c, const matrix_t *G) {
  if (!G)
    return;
  int row = matrix_get_row(G);
  int col = matrix_get_col(G);
  if (matrix_get_col(c) != col || matrix_get_row(c) != 1)
    return;
  matrix_t *ligne = matrix_alloc(1,col);
  matrix_t *lambda = matrix_random(1,col);
  if (!ligne || !lambda) {
    matrix_free(ligne);
    matrix_free(lambda);
    return;
  }
  for (int i = 0; i < col; ++i)
    matrix_set_cell(c, 0, i, 0);
  for (int i = 0; i < row; ++i) {
    matrix_row(ligne, G, i);
    matrix_add_modified(c, ligne, lambda->mat[0][i]);
  }
  matrix_free(ligne);
  matrix_free(lambda);
}


////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// affichage ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void matrix_print(const matrix_t *matrix, FILE *fd) {
  if (matrix && matrix->mat)
  {
    int col = matrix->nb_col;
    int row = matrix->nb_row;
    for (int i = 0; i < row; ++i)
    {
      char val;
      for (int j = 0; j < col; ++j)
      {
        val = matrix_get_cell(matrix,i,j);
        if (val == '*')
          fprintf(fd, " * ");
        else if (val == 2 && ORDER == 3)
          fprintf(fd,"-1 ");
        else
          fprintf(fd," %d ", val);
      }
      fputs("\n", fd);
    }
    fputs("\n\n", fd);
  }
}
