/* ------- file: -------------------------- matrix.c ----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Dec 29 14:15:27 1999 --

       --------------------------                      ----------RH-- */

/* --- Routines to create 2-dimensional arrays that are contiguous in
       memory. Therefore the whole array pointed to by pointer **matrix_??
       can be addressed as either 1-dimensional array with pointer
       matrix_??[0], or as 2-dimensional array as matrix_??[i][j].

 Note: space is initialized to zeros by using calloc rather than malloc.

       --                                              -------------- */
 
#include <stdlib.h>

#include "rh.h"
#include "error.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern char   messageStr[];


/* ------- begin -------------------------- matrix_char.c ----------- */

char **matrix_char(int Nrow, int Ncol)
{
  register int i;

  char *theMatrix, **Matrix;
  int   typeSize = sizeof(char), pointerSize = sizeof(char *);

  theMatrix = (char *)  calloc(Nrow * Ncol, typeSize);
  Matrix    = (char **) malloc(Nrow * pointerSize);
  for (i = 0;  i < Nrow;  i++, theMatrix += Ncol)
    Matrix[i] = theMatrix;

  return Matrix;
}
/* ------- end ---------------------------- matrix_char.c ----------- */

/* ------- begin -------------------------- matrix_int.c ------------ */

int **matrix_int(int Nrow, int Ncol)
{
  register int i;

  int *theMatrix, **Matrix, typeSize = sizeof(int),
       pointerSize = sizeof(int *);

  theMatrix = (int *)  calloc(Nrow * Ncol, typeSize);
  Matrix    = (int **) malloc(Nrow * pointerSize);
  for (i = 0;  i < Nrow;  i++, theMatrix += Ncol)
    Matrix[i] = theMatrix;

  return Matrix;
}
/* ------- end ---------------------------- matrix_int.c ------------ */

/* ------- begin -------------------------- matrix_double.c --------- */

double **matrix_double(int Nrow, int Ncol)
{
  register int i;

  int     typeSize = sizeof(double), pointerSize = sizeof(double *);
  double *theMatrix, **Matrix;

  theMatrix = (double *)  calloc(Nrow * Ncol, typeSize);
  Matrix    = (double **) malloc(Nrow * pointerSize);
  for (i = 0;  i < Nrow;  i++, theMatrix += Ncol)
    Matrix[i] = theMatrix;

  return Matrix;
}
/* ------- end ---------------------------- matrix_double.c --------- */

/* ------- begin -------------------------- matrix3d_double.c --------- */
double **matrix3d_double(int Nrow, int Ncol, int Ndep)
{
  register int i, j;

  int     typeSize = sizeof(double), pointerSize = sizeof(double *), doublepointerSize = sizeof(double **);
  double  ***Matrix3d;

  Matrix3d = calloc(Nrow, sizeof(*Matrix3d));
  Matrix3d[0] = malloc(Nrow*Ncol * sizeof(*Matrix3d[0]));
  Matrix3d[0][0] = malloc(Nrow*Ncol*Ndep * sizeof(Matrix3d[0][0])); // Contiguous
  
  for (i = 1; i < Nrow; i++) {
    Matrix3d[i] = &Matrix3d[0][i*Ncol];
  }
  
  for (j = 1; j < Nrow*Ncol; j++) {
    Matrix3d[0][j] = &Matrix3d[0][0][j*Ndep];
  }

  return Matrix3d;
}
/* ------- end ---------------------------- matrix3d_double.c --------- */


/* ------- begin -------------------------- freeMatrix.c ------------ */

void freeMatrix(void **matrix)
{
  const char routineName[] = "freeMatrix";

  if (matrix == NULL  ||  matrix[0] == NULL) {
    sprintf(messageStr, "Trying to free NULL pointer");
    Error(ERROR_LEVEL_2, routineName, messageStr);
  } else {
    /* --- Free the memory allocated for matrix matrix -- ----------- */

    free(matrix[0]);
    free(matrix);
  }
}
/* ------- end ---------------------------- freeMatrix.c ------------ */

void freeMatrix3d(void ***matrix3d)
{
  const char routineName[] = "freeMatrix";

  if (matrix3d == NULL  ||  matrix3d[0] == NULL || matrix3d[0,0] == NULL) {
    sprintf(messageStr, "Trying to free NULL pointer");
    Error(ERROR_LEVEL_2, routineName, messageStr);
  } else {
    /* --- Free the memory allocated for matrix matrix3d -- ----------- */

    free(matrix3d[0][0]);
    free(matrix3d[0]);
    free(matrix3d);
  }
}
