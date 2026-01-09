/* ------- file: -------------------------- bezier_aux.c ------------

   Auxiliary routines for Cubic DELO-Bezier (polarized) and
   cubic short-char Bezier solvers.
   
   References: de la Cruz Rodriguez & Piskunov (2013), Auer (2003)
               (Derivatives) Fritsch & Butland (1984),
	       
   Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

   Modifications:
           2017-03-12, JdlCR: Created!

       Last modified: Thu May 24 16:47:42 2018 --

       --------------------------                      ----------RH-- */


#include <math.h>
#include <string.h>    // memcpy, memset
#include <stdio.h>

#ifdef SIMDON
#include <x86intrin.h> // Intrinsic SSE instructions
#endif

#include "bezier.h"
#include "headers/xdr.h"
#include "rh.h"
#include "inputs.h"

/* --- Global variables --                             -------------- */

extern InputData input;

/* ------- begin -------------------------- cent_deriv.c ------------ */

inline double cent_deriv(double dsup,double dsdn, 
			 double chiup,double chic, double chidn)
{
  /* --- Derivative Fritsch & Butland (1984) --        -------------- */

  double fim1, fi, alpha, wprime;
  
  fim1=(chic - chiup) / dsup;
  fi=(chidn - chic) / dsdn;

  if (fim1*fi > 0) {
    alpha = 0.333333333333333333333333 * ( 1.0 + dsdn / (dsdn+dsup) );
    wprime = (fim1 * fi) / ( (1.0 - alpha) * fim1 + alpha*fi );
  } else {
    wprime=0.0;
  }
  return wprime;
}
/* ------- end ---------------------------- cent_deriv.c ------------ */


/* ------- begin -------------------------- cent_deriv_mat.c -------- */

inline void cent_deriv_mat(double wprime[4][4], double dsup, double dsdn,
			   double chiup[4][4], double chic[4][4],
			   double chidn[4][4])
{
  register int i,j;
  
  for(j = 0;  j<4;  j++)
    for(i = 0;  i < 4;  i++)
      wprime[j][i] = cent_deriv(dsup, dsdn, chiup[j][i], chic[j][i], chidn[j][i]);
}

/* ------- end ---------------------------- cent_deriv_mat.c -------- */


/* ------- begin -------------------------- cent_deriv_vec.c -------- */

inline void cent_deriv_vec(double wprime[4], double dsup, double dsdn,
		    double chiup[4], double chic[4], double chidn[4])
{
  register int i;
  
  for(i = 0;  i < 4;  i++)
    wprime[i] = cent_deriv(dsup, dsdn, chiup[i], chic[i], chidn[i]);
  
}
/* ------- end ---------------------------- cent_deriv_vec.c -------- */


/* ------- begin -------------------------- m4m.c ------------------- */


inline void m4m(double a[4][4], double b[4][4], double c[4][4])
{

  /* --- Matrix multiplication --                      -------------- */

  register int i, j, k;
  memset(&c[0][0],0,sizeof(double)*16);
  for(j = 0; j<4; j++)
    for(i = 0; i<4; i++)
      for(k = 0; k<4; k++)
	c[j][i] += a[k][i]*b[j][k]; 
}
/* ------- end ---------------------------- m4m.c ------------------- */


/* ------- begin -------------------------- m4v.c ------------------- */

/* --- Matrix/vector multiplication.
       We use matrix as float to be able to use Intel's
       matrix inversion --                         ------------------ */

inline void m4v(float a[4][4], double b[4], double c[4])
{
  register int k, i;
  memset(&c[0],0,sizeof(double)*4);
  for(i = 0; i<4; i++)
    for(k = 0; k<4; k++)
      c[i] += ((double)a[i][k]) * b[k];
}

inline void m4v_double(double a[4][4], double *b, double c[4])
{
  register int k, i;
  memset(&c[0],0,sizeof(double)*4);
  for(i = 0; i<4; i++)
    for(k = 0; k<4; k++)
      c[i] += a[i][k] * b[k];
}

/* ------- end ---------------------------- m4v.c ------------------- */


/* ------- begin -------------------------- Svec.c ------------------ */

inline void Svec(int k, double **S, double *Sf)
{
  /* --- Extracts the Source vector at depth-points k -- ------------ */

  Sf[0] = S[0][k], Sf[1] = S[1][k], Sf[2] = S[2][k], Sf[3] = S[3][k];
}

/* ------- end ---------------------------- Svec.c ------------------ */


/* ------- begin -------------------------- SIMD_MatInv.c ----------- */
#ifdef SIMDON
void MatInv(float* mat)
{
  /* --- 

     Very fast in-place 4x4 Matrix inversion using SIMD instrutions
     Only works with 32-bits floats. It uses Cramer's rule.
     
     Provided by Intel

     Requires SSE instructions but all x86 machines since 
     Pentium III have them.
     
     --                                            ------------------ */
  
  __m128 minor0, minor1, minor2, minor3;
  __m128 row0, row1, row2, row3;
  __m128 det, tmp1;
  
  // -----------------------------------------------
  tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(mat)), (__m64*)(mat+ 4));
  row1 = _mm_loadh_pi(_mm_loadl_pi(row1, (__m64*)(mat+8)), (__m64*)(mat+12));
  row0 = _mm_shuffle_ps(tmp1, row1, 0x88);
  row1 = _mm_shuffle_ps(row1, tmp1, 0xDD);
  tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(mat+ 2)), (__m64*)(mat+ 6));
  row3 = _mm_loadh_pi(_mm_loadl_pi(row3, (__m64*)(mat+10)), (__m64*)(mat+14));
  row2 = _mm_shuffle_ps(tmp1, row3, 0x88);
  row3 = _mm_shuffle_ps(row3, tmp1, 0xDD);
  // -----------------------------------------------
  tmp1 = _mm_mul_ps(row2, row3);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
  minor0 = _mm_mul_ps(row1, tmp1);
  minor1 = _mm_mul_ps(row0, tmp1);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
  minor0 = _mm_sub_ps(_mm_mul_ps(row1, tmp1), minor0);
  minor1 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor1);
  minor1 = _mm_shuffle_ps(minor1, minor1, 0x4E);
  // -----------------------------------------------
  tmp1 = _mm_mul_ps(row1, row2);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
  minor0 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor0);
  minor3 = _mm_mul_ps(row0, tmp1);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
  minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row3, tmp1));
  minor3 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor3);
  minor3 = _mm_shuffle_ps(minor3, minor3, 0x4E);
  // -----------------------------------------------
  tmp1 = _mm_mul_ps(_mm_shuffle_ps(row1, row1, 0x4E), row3);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
  row2 = _mm_shuffle_ps(row2, row2, 0x4E);
  minor0 = _mm_add_ps(_mm_mul_ps(row2, tmp1), minor0);
  minor2 = _mm_mul_ps(row0, tmp1);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
  minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row2, tmp1));
  minor2 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor2);
  minor2 = _mm_shuffle_ps(minor2, minor2, 0x4E);
  // -----------------------------------------------
  tmp1 = _mm_mul_ps(row0, row1);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
  minor2 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor2);
  minor3 = _mm_sub_ps(_mm_mul_ps(row2, tmp1), minor3);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
  minor2 = _mm_sub_ps(_mm_mul_ps(row3, tmp1), minor2);
  minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row2, tmp1));
  // -----------------------------------------------
  tmp1 = _mm_mul_ps(row0, row3);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
  minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row2, tmp1));
  minor2 = _mm_add_ps(_mm_mul_ps(row1, tmp1), minor2);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
  minor1 = _mm_add_ps(_mm_mul_ps(row2, tmp1), minor1);
  minor2 = _mm_sub_ps(minor2, _mm_mul_ps(row1, tmp1));
  // -----------------------------------------------
  tmp1 = _mm_mul_ps(row0, row2);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
  minor1 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor1);
  minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row1, tmp1));
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
  minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row3, tmp1));
  minor3 = _mm_add_ps(_mm_mul_ps(row1, tmp1), minor3);
  // -----------------------------------------------
  det = _mm_mul_ps(row0, minor0);
  det = _mm_add_ps(_mm_shuffle_ps(det, det, 0x4E), det);
  det = _mm_add_ss(_mm_shuffle_ps(det, det, 0xB1), det);
  tmp1 = _mm_rcp_ss(det);
  det = _mm_sub_ss(_mm_add_ss(tmp1, tmp1),
		   _mm_mul_ss(det, _mm_mul_ss(tmp1, tmp1)));
  det = _mm_shuffle_ps(det, det, 0x00);
  minor0 = _mm_mul_ps(det, minor0);
  _mm_storel_pi((__m64*)(mat), minor0);
  _mm_storeh_pi((__m64*)(mat+2), minor0);
  minor1 = _mm_mul_ps(det, minor1);
  _mm_storel_pi((__m64*)(mat+4), minor1);
  _mm_storeh_pi((__m64*)(mat+6), minor1);
  minor2 = _mm_mul_ps(det, minor2);
  _mm_storel_pi((__m64*)(mat+ 8), minor2);
  _mm_storeh_pi((__m64*)(mat+10), minor2);
  minor3 = _mm_mul_ps(det, minor3);
  _mm_storel_pi((__m64*)(mat+12), minor3);
  _mm_storeh_pi((__m64*)(mat+14), minor3);
}
/* ------- end ---------------------------- SIMD_MatInv.c ----------- */
#else
void MatInv(float *mat)
{
  float tmp[12]; /* temp array for pairs*/
  float src[16]; /* array of transpose source matrix */
  float det;     /* determinant*/
  
  // output matrix -- will be copied into 'mat' at the end
  float *dst;
  dst = (float *) malloc(16 * sizeof(float));

  /* transpose matrix */
  for (int i = 0; i < 4; i++) {
    src[i] = mat[i*4];
    src[i + 4] = mat[i*4 + 1];
    src[i + 8] = mat[i*4 + 2];
    src[i + 12] = mat[i*4 + 3];
    }
  
  /* calculate pairs for first 8 elements (cofactors) */
  tmp[0] = src[10] * src[15];
  tmp[1] = src[11] * src[14];
  tmp[2] = src[9] * src[15];
  tmp[3] = src[11] * src[13];
  tmp[4] = src[9] * src[14];
  tmp[5] = src[10] * src[13];
  tmp[6] = src[8] * src[15];
  tmp[7] = src[11] * src[12];
  tmp[8] = src[8] * src[14];
  tmp[9] = src[10] * src[12];
  tmp[10] = src[8] * src[13];
  tmp[11] = src[9] * src[12];

  /* calculate first 8 elements (cofactors) */
  dst[0] = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7];
  dst[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];
  dst[1] = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7];
  dst[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7];
  dst[2] = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];
  dst[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7];
  dst[3] = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6];
  dst[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6];
  dst[4] = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3];
  dst[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3];
  dst[5] = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3];
  dst[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];
  dst[6] = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3];
  dst[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];
  dst[7] = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];
  dst[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];

  /* calculate pairs for second 8 elements (cofactors) */
  tmp[0] = src[2]*src[7];
  tmp[1] = src[3]*src[6];
  tmp[2] = src[1]*src[7];
  tmp[3] = src[3]*src[5];
  tmp[4] = src[1]*src[6];
  tmp[5] = src[2]*src[5];
  tmp[6] = src[0]*src[7];
  tmp[7] = src[3]*src[4];
  tmp[8] = src[0]*src[6];
  tmp[9] = src[2]*src[4];
  tmp[10] = src[0]*src[5];
  tmp[11] = src[1]*src[4];

  /* calculate second 8 elements (cofactors) */
  dst[8] = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15];
  dst[8] -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];
  dst[9] = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15];
  dst[9] -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15];
  dst[10] = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];
  dst[10]-= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15];
  dst[11] = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];
  dst[11]-= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14];
  dst[12] = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];
  dst[12]-= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10];
  dst[13] = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10];
  dst[13]-= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8];
  dst[14] = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8];
  dst[14]-= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9];
  dst[15] = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9];
  dst[15]-= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8];

  /* calculate determinant */
  det = src[0]*dst[0]+src[1]*dst[1]+src[2]*dst[2]+src[3]*dst[3];

  /* calculate matrix inverse */
  det = 1.0/det;
  for (int j = 0; j < 16; j++){
    dst[j] *= det;
  }

  // get back the result to input matrix
  memcpy(mat, dst, sizeof(float)*16);
}
#endif

/* ------- begin -------------------------- Bezier3_coeffs ---------- */

inline void Bezier3_coeffs(double dt, double *alpha, double *beta,
		    double *gamma, double *theta, double *eps)
{
  /* --- Integration coeffs. for cubic Bezier interpolants
         Use Taylor expansion if dtau is small --  ------------------ */
  
  double dt2 = dt*dt, dt3 = dt2 * dt,  dt4;
    
  if(dt >= 5.e-2){

    *eps = exp(-dt);

    *alpha = (-6.0 + 6.0 * dt - 3.0 * dt2 + dt3 + 6.0 * eps[0])        / dt3;
    dt3 = 1.0/dt3;
    *beta  = (6.0 + (-6.0 - dt * (6.0 + dt * (3.0 + dt))) * eps[0])    * dt3;
    *gamma = 3.0 * (6.0 + (-4.0 + dt)*dt - 2.0 * (3.0 + dt) * eps[0])  * dt3;
    *theta = 3.0 * ( eps[0] * (6.0 + dt2 + 4.0 * dt) + 2.0 * dt - 6.0) * dt3;
  } else{
    dt4 = dt2*dt2;
    *eps = 1.0 - dt + 0.5 * dt2 - dt3 / 6.0 + dt4 / 24.0;

    *alpha = 0.25 * dt - 0.05 * dt2 + dt3 / 120.0 - dt4 / 840.0;
    *beta  = 0.25 * dt - 0.20 * dt2 + dt3 / 12.0  - dt4 / 42.0; 
    *gamma = 0.25 * dt - 0.10 * dt2 + dt3 * 0.025 - dt4 / 210.0; 
    *theta = 0.25 * dt - 0.15 * dt2 + dt3 * 0.05  - dt4 / 84.0; 
  }
}
/* ------- end ---------------------------- Bezier3_coeffs ---------- */


/* ------- begin -------------------------- m4inv.c ----------------- */

void m4inv(double MI[4][4])
{

  /* --- In-place Shipley-Coleman matrix inversion
         Fast, but ... how accurate??
         Pivoting is always done in the diagonal.
         Copied here just in case the SIMD matrix inversion 
         gives troubles. --                            -------------- */
  
  register int k, i, j;
  
  for (k = 0;  k < 4;  k++){

    /* --- The pivot element --                        -------------- */
    
    MI[k][k] = -1.0 / MI[k][k];

    /* --- The pivot column --                         -------------- */
    
    for(i = 0;  i < 4;  ++i) if(i != k) MI[i][k]*=MI[k][k];

    /* --- Elements not in a pivot row or column --    -------------- */
    
    for(i = 0;  i < 4;  ++i) {
      if(i != k)
	for(j = 0;  j < 4;  ++j)
	  if(j != k)
	    MI[i][j] += MI[i][k] * MI[k][j];
    }
    /* --- Elements in a pivot row --                  -------------- */
    
    for(i = 0;  i < 4;  ++i) {
      if(i != k)
	MI[k][i] *= MI[k][k];
    }
  }
  
  for(i = 0;  i < 4;  ++i) {
    for(j = 0;  j < 4;  ++j) MI[i][j] = -MI[i][j];
  }
  return;
}
/* ------- end ---------------------------- m4inv.c ----------------- */
