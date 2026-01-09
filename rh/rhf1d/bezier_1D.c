/* ------- file: -------------------------- bezier_1D.c -------------

   Cubic DELO-Bezier (polarized) and cubic short-char Bezier solvers.
   
   References: de la Cruz Rodriguez & Piskunov (2013), Auer (2003)
               (Derivatives) Fritsch & Butland (1984),
	       
   Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

   Modifications:
           2017-03-12, JdlCR: Created!

       Last modified: Thu May 31 13:34:48 2018 --

       --------------------------                      ----------RH-- */


#include <math.h>
#include <string.h>

#include "../rh.h"
#include "../error.h"
#include "../atom.h"
#include "../atmos.h"
#include "geometry.h"
#include "../spectrum.h"
#include "../bezier.h"

#include "../inputs.h"


/* --- Identity matrix --                              -------------- */

static const double ident[4][4] =
  {{1.0, 0.0, 0.0, 0.0},
   {0.0, 1.0, 0.0, 0.0},
   {0.0, 0.0, 1.0, 0.0},
   {0.0, 0.0, 0.0, 1.0}};


/* --- Global variables --                             -------------- */

extern Geometry geometry;
extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- Piece_Stokes_Bezier3_1D.c */

void Piece_Stokes_Bezier3_1D(int nspect, int mu, bool_t to_obs,
			     double *chi, double **S, double **I,
			     double *Psi)
{
  /* --- Cubic DELO-Bezier solver for polarized light
         Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

         Reference(s):
         J. de la Cruz Rodriguez & N. Piskunov (2013)
         --                                        ------------------ */
  
  const char routineName[] = "Piece_Stokes_Bezier3_1D";
  register int k, n, m, i, j;
  
  int    Ndep = geometry.Ndep, k_start, k_end, dk;
  double dtau_uw, dtau_dw = 0.0, c1, c2, w[3], dsdn2, dchi_dn,
         I_upw[4], Bnu[2];
  double dchi_up,dchi_c,dt03;
  double dsup,dsdn,dt,eps=0,alpha=0,beta=0,gamma=0,theta=0;
  double Ku[4][4], K0[4][4], Kd[4][4], dKu[4][4], dK0[4][4];
  double Su[4], S0[4], Sd[4], dSu[4], dS0[4];
  double A[4][4], Ma[4][4], Mb[4][4], Mc[4][4], V0[4], V1[4];
  double imu = 1.0 / geometry.muz[mu];
  float Md[4][4];
  double *z = geometry.height;

  FILE *fptr;
  fptr = fopen("I.txt", "a");
  
  if (to_obs) {
    dk      = -1;
    k_start = Ndep-1;
    k_end   = 0;
  } else {
    dk      = 1;
    k_start = 0;
    k_end   = Ndep-1;
  }
  dtau_uw = 0.5 * imu * (chi[k_start] + chi[k_start+dk]) *
    fabs(z[k_start] - z[k_start+dk]);
  
  /* --- Boundary conditions --                        -------------- */

  if (to_obs) {
    switch (geometry.vboundary[BOTTOM]) {
    case ZERO:
      for (n = 0;  n < 4;  n++) I_upw[n] = 0.0;
      break;
    case THERMALIZED:
      Planck(2, &atmos.T[Ndep-2], spectrum.lambda[nspect], Bnu, -1);
      I_upw[0] = Bnu[1] - (Bnu[0] - Bnu[1]) / dtau_uw;
      for (n = 1;  n < 4;  n++) I_upw[n] = 0.0;
      break;
    case IRRADIATED:
      I_upw[0] = geometry.Ibottom[nspect][mu];
      for (n = 1;  n < 4;  n++) I_upw[n] = 0.0;
      break;
    case REFLECTIVE:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[BOTTOM]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  } else {
    switch (geometry.vboundary[TOP]) {
    case ZERO:
      for (n = 0;  n < 4;  n++) I_upw[n] = 0.0;
      break;
    case IRRADIATED:
      I_upw[0] = geometry.Itop[nspect][mu];
      for (n = 1;  n < 4;  n++) I_upw[n] = 0.0;
      break;
    default:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[TOP]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  }

  for (n = 0;  n < 4;  n++) I[n][k_start] = I_upw[n];
  if (Psi) Psi[k_start] = 0.0;

  k=k_start+dk;
  dsup = fabs(z[k] - z[k-dk]) * imu;
  dsdn = fabs(z[k+dk] - z[k]) * imu;
  dchi_up= (chi[k] - chi[k-dk])/dsup;

  
  /* ---  dchi/ds at central point--               ------------------ */
  
  dchi_c = cent_deriv(dsup,dsdn,chi[k-dk],chi[k],chi[k+dk]);
  
  /* --- Upwind path_length (BEzier3 integration) -- ---------------- */

  c2 = MAX(chi[k]    - (dsup/3.0) * dchi_c,  0.0);
  c1 = MAX(chi[k-dk] + (dsup/3.0) * dchi_up, 0.0);
  
  dtau_uw = 0.25 * dsup * (chi[k] + chi[k-dk] + c1 + c2);
  
  /* --- Ku, K0 and dKu, dSu -                     ------------------ */
  
  StokesK(nspect, k_start,    chi[k_start],    Ku);
  StokesK(nspect, k_start+dk, chi[k_start+dk], K0);

  Svec(k_start,    S, Su);
  Svec(k_start+dk, S, S0);

  /* --- Assume side derivative in the first interval -- ------------ */
  
  for(n = 0;  n < 4;  n++){
    dSu[n] = (S0[n] - Su[n]) / dtau_uw;
    
    for(m = 0;  m < 4;  m++)
      dKu[n][m] = (K0[n][m] - Ku[n][m]) / dtau_uw;
  }
  
  /* --- Solve transfer along ray --                   -------------- */

  if (to_obs) fprintf(fptr, "%2.8e %2.8e %2.8e %2.8e\n", I[0][k_start], I[1][k_start], I[2][k_start], I[3][k_start]);

  for (k = k_start+dk;  k != k_end;  k += dk) { 

    // if (nspect==49 && k>50 && to_obs){
    //   printf("k = %d | chi = (%2.8e, %2.8e, %2.8e, %2.8e)\n", k, spectrum.chi_c_lam[nspect][k], spectrum.chi_c_lam[nspect][atmos.Nspace+k], spectrum.chi_c_lam[nspect][2*atmos.Nspace+k], spectrum.chi_c_lam[nspect][3*atmos.Nspace+k]);
    // }
      
    /* --- dchi/ds at downwind point --                -------------- */
      
    dsdn = fabs(z[k+dk] - z[k]) * imu;
      
    if(fabs(k - k_end) > 1){
      dsdn2 = fabs(z[k+2*dk] - z[k+dk]) * imu;
      dchi_dn = cent_deriv(dsdn, dsdn2, chi[k], chi[k+dk], chi[k+2*dk]);       
    } else
      dchi_dn = (chi[k+dk] - chi[k])/dsdn;      
      
    /* --- Make sure that c1 and c2 don't do below zero -- ---------- */
      
    c2 = MAX(chi[k]    + (dsdn/3.0) * dchi_c , 0.0);
    c1 = MAX(chi[k+dk] - (dsdn/3.0) * dchi_dn, 0.0);
          
    /* --- Bezier3 integrated dtau --              ------------------ */
      
    dtau_dw = 0.25 * dsdn * (chi[k] + chi[k+dk] + c1 + c2);
    dt = dtau_uw, dt03 = dt / 3.0;
  
    /* --- Bezier3 coeffs. --                      ------------------ */
      
    Bezier3_coeffs(dt, &alpha, &beta, &gamma, &theta, &eps);
   
    /* --- Diagonal operator --                    ------------------ */
      
    if(Psi) Psi[k] = alpha + gamma;
   
    /* ---- get algebra in place --                ------------------ */
      
    StokesK(nspect, k+dk, chi[k+dk], Kd);
    Svec(k+dk, S, Sd);

    cent_deriv_mat(dK0, dtau_uw, dtau_dw, Ku, K0, Kd);
    cent_deriv_vec(dS0, dtau_uw, dtau_dw, Su, S0, Sd);

    m4m(Ku, Ku, Ma); // Ku # Ku
    m4m(K0, K0, A ); // K0 # K0

    for(j = 0;  j < 4;  j++){
      for(i = 0;  i < 4;  i++){
        Md[j][i] = ident[j][i] + alpha * K0[j][i] + gamma*(-dt03 * (A[j][i] + dK0[j][i] + K0[j][i]) + K0[j][i]);
          
        Ma[j][i] = eps * ident[j][i] - beta * Ku[j][i] - theta*(dt03 * (Ma[j][i] + dKu[j][i] + Ku[j][i]) + Ku[j][i]);
          
        Mc[j][i] = alpha* ident[j][i] + gamma * (ident[j][i] - dt03 * K0[j][i]);
        Mb[j][i] = beta * ident[j][i] + theta * (ident[j][i] + dt03 * Ku[j][i]);
      }
    }
      
    /* --- Here I am doing Ma*stk + Mb * Su + Mc * S0 + 
           (gam * dS0 - theta * dSu) * dtau / 3.0 to compute the 
           right-hand term
           --                                      ------------------ */

    memset(V0, 0, 4*sizeof(double));

    for(i = 0;  i < 4;  i++){
      for(j = 0;  j < 4;  j++){
      	V0[i] += Ma[i][j] * I[j][k-dk] + Mb[i][j] * Su[j] + Mc[i][j] * S0[j];
      }
      V0[i] += dt03 * (-gamma * dS0[i] + theta * dSu[i]);
    }
    /* --- Solve linear system to get the intensity -- -------------- */
      
    // SIMD_MatInv(Md[0]);   // Invert Md
    MatInv(Md[0]);
    m4v(Md, V0, V1);      // Multiply Md^-1 * V0

    for(i=0;i<4;i++){
      I[i][k] = V1[i];
      if (to_obs) fprintf(fptr, "%2.8e ", I[i][k]);
    }
    if (to_obs) fprintf(fptr, "\n");
      
    /* --- Shift values for next depth --          ------------------ */
      
    memcpy(Su,   S0, 4*sizeof(double));
    memcpy(S0,   Sd, 4*sizeof(double));
    memcpy(dSu, dS0, 4*sizeof(double));
      
    memcpy(Ku[0],   K0[0], 16*sizeof(double));
    memcpy(K0[0],   Kd[0], 16*sizeof(double));
    memcpy(dKu[0], dK0[0], 16*sizeof(double));
      
    dtau_uw = dtau_dw;
    dsup    = dsdn;
    dchi_up = dchi_c;
    dchi_c  = dchi_dn;    
  }
      
  /* --- Linear integration in the last interval -- ----------------- */
  
  k = k_end;
  dtau_uw = 0.5*imu * (chi[k] + chi[k-dk]) *
    fabs(geometry.height[k] - geometry.height[k-dk]);
  w3(dtau_uw, w);

  /* --- dSu is defined negative in Han's implementation ------------ */
  
  for (n = 0;  n < 4;  n++)
    V0[n] = w[0]*S[n][k] + w[1] * -dSu[n];
  
  if (Psi) Psi[k] = w[0] - w[1] / dtau_uw;
      
  for (n = 0;  n < 4;  n++) {
    for (m = 0;  m < 4;  m++) {
      A[n][m]  = -w[1]/dtau_uw * Ku[n][m];
      Md[n][m] = (w[0] - w[1]/dtau_uw) * K0[n][m];
    }
    A[n][n]  = 1.0 - w[0];
    Md[n][n] = 1.0;
  }
      
  for (n = 0;  n < 4;  n++) 
    for (m = 0;  m < 4;  m++) 
      V0[n] += A[n][m] * I[m][k-dk];

  /* --- Solve linear system --                    ------------------ */
  
  // SIMD_MatInv(Md[0]); // Invert Md
  MatInv(Md[0]);
  m4v(Md,V0,V1);      // Multiply Md^-1 * V0
  
  for (n = 0;  n < 4;  n++){
    I[n][k] = V1[n];
    if (to_obs) fprintf(fptr, "%2.8e ", I[n][k]);
  }
  if (to_obs) fprintf(fptr, "\n");

  fclose(fptr);
}
/* ------- end ------------------------- Piece_Stokes_Bezier3_1D.c -- */


void Piece_Stokes_Bezier3_1D_RFs(int nspect, int mu, bool_t to_obs,
			     double *chi, double **S, double **I, double ***dI)
{
  /* --- Cubic DELO-Bezier solver for polarized light perturbation
         Copied from: Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

         Reference(s):
         J. de la Cruz Rodriguez & N. Piskunov (2013)
         --                                        ------------------ */
  
  const char routineName[] = "Piece_Stokes_Bezier3_1D";
  register int k, n, m, i, j;
  
  int    Ndep = geometry.Ndep, k_start, k_end, dk;
  double dtau_uw, dtau_dw = 0.0, c1, c2, w[3], dsdn2, dchi_dn,
         I_upw[4], Bnu[2];
  double dchi_up,dchi_c,dt03;
  double dsup,dsdn,dt,eps=0,alpha=0,beta=0,gamma=0,theta=0;
  double Ku[4][4], K0[4][4], Kd[4][4], dKu[4][4], dK0[4][4];
  double Su[4], S0[4], Sd[4], dSu[4], dS0[4];
  double A[4][4], Ma[4][4], Mb[4][4], Mc[4][4], V0[4], V1[4];
  double imu = 1.0 / geometry.muz[mu];
  float Md[4][4];
  double *z = geometry.height;

  // RFs variables
  int idp;
  double **Gu, **G0, **Gd;
  double **dGu, **dG0;
  double Ml[4][4];

  FILE *fptr;
  // Open a file in writing mode
  fptr = fopen("dI.txt", "a");
 
  Gu      = matrix_double(input.n_atomic_pars,4);
  G0      = matrix_double(input.n_atomic_pars,4);
  Gd      = matrix_double(input.n_atomic_pars,4);
  dGu     = matrix_double(input.n_atomic_pars,4);
  dG0     = matrix_double(input.n_atomic_pars,4);

  if (to_obs) {
    dk      = -1;
    k_start = Ndep-1;
    k_end   = 0;
  } else {
    dk      = 1;
    k_start = 0;
    k_end   = Ndep-1;
  }
  dtau_uw = 0.5 * imu * (chi[k_start] + chi[k_start+dk]) * fabs(z[k_start] - z[k_start+dk]);
  
  /* --- Boundary conditions --                        -------------- */
  
  for (n=0; n<4; n++){
    for (idp=0; idp<input.n_atomic_pars; idp++){
      dI[n][k_start][idp] = 0.0;
    }
  }

  k=k_start+dk;
  dsup = fabs(z[k] - z[k-dk]) * imu;
  dsdn = fabs(z[k+dk] - z[k]) * imu;
  dchi_up= (chi[k] - chi[k-dk])/dsup;

  /* ---  dchi/ds at central point--               ------------------ */
  
  dchi_c = cent_deriv(dsup,dsdn,chi[k-dk],chi[k],chi[k+dk]);
  
  /* --- Upwind path_length (Bezier3 integration) -- ---------------- */

  c2 = MAX(chi[k]    - (dsup/3.0) * dchi_c,  0.0);
  c1 = MAX(chi[k-dk] + (dsup/3.0) * dchi_up, 0.0);
  
  dtau_uw = 0.25 * dsup * (chi[k] + chi[k-dk] + c1 + c2);

  /* --- Ku, K0 and dKu, dSu -                     ------------------ */
  
  StokesK(nspect, k_start,    chi[k_start],    Ku);
  StokesK(nspect, k_start+dk, chi[k_start+dk], K0);
  
  Gvec(nspect, k_start, chi[k_start], I, Gu);
  Gvec(nspect, k_start+dk, chi[k_start+dk], I, G0);

  /* --- Assume side derivative in the first interval -- ------------ */

  for(n = 0;  n < 4;  n++){
    for (idp=0; idp<input.n_atomic_pars; idp++){
      dGu[idp][n] = (G0[idp][n] - Gu[idp][n]) / dtau_uw;
    }
    for (m=0; m<4; m++){
      dKu[n][m] = (K0[n][m] - Ku[n][m]) / dtau_uw;
    }
  }

  fprintf(fptr, "%2.8e %2.8e %2.8e %2.8e\n", dI[0][k_start][idp], dI[1][k_start][idp], dI[2][k_start][idp], dI[3][k_start][idp]);

  /* --- Solve transfer along ray --                   -------------- */
  for (k = k_start+dk;  k != k_end;  k += dk) {     

    // if (nspect==49 && k>50 && to_obs){
    //   printf("k = %d | dchi = (%2.8e, %2.8e, %2.8e, %2.8e)\n", k, spectrum.dchi_c_lam[nspect][k][0], spectrum.dchi_Q[nspect][k][0], spectrum.dchi_U[nspect][k][0], spectrum.dchi_V[nspect][k][0]);
    // }
    
    /* --- dchi/ds at downwind point --                -------------- */
      
    dsdn = fabs(z[k+dk] - z[k]) * imu;
      
    if(fabs(k - k_end) > 1){
      dsdn2 = fabs(z[k+2*dk] - z[k+dk]) * imu;
      dchi_dn = cent_deriv(dsdn, dsdn2, chi[k], chi[k+dk], chi[k+2*dk]);       
    } else
      dchi_dn = (chi[k+dk] - chi[k])/dsdn;      
      
    /* --- Make sure that c1 and c2 don't do below zero -- ---------- */
      
    c1 = MAX(chi[k+dk] - (dsdn/3.0) * dchi_dn, 0.0);
    c2 = MAX(chi[k]    + (dsdn/3.0) * dchi_c , 0.0);
          
    /* --- Bezier3 integrated dtau --              ------------------ */
      
    dtau_dw = 0.25 * dsdn * (chi[k] + chi[k+dk] + c1 + c2);
    dt = dtau_uw, dt03 = dt / 3.0;
  
    /* --- Bezier3 coeffs. --                      ------------------ */
      
    Bezier3_coeffs(dt, &alpha, &beta, &gamma, &theta, &eps);
    w3(dtau_uw, w);
   
    /* ---- get algebra in place --                ------------------ */
    
    StokesK(nspect, k+dk, chi[k+dk], Kd);
    Gvec(nspect, k+dk, chi[k+dk], I, Gd);

    cent_deriv_mat(dK0, dtau_uw, dtau_dw, Ku, K0, Kd);
    for (idp=0; idp<input.n_atomic_pars; idp++){
      cent_deriv_vec(dG0[idp], dtau_uw, dtau_dw, Gu[idp], G0[idp], Gd[idp]);
    }
    
    m4m(K0, K0, Ma); // K x K
    m4m(Ku, Ku, Ml);
    for (j = 0;  j < 4;  j++){
      for (i = 0;  i < 4;  i++){
        Ma[j][i] += dK0[j][i] + K0[j][i];
        Md[j][i] = ident[j][i] + alpha*K0[j][i] + gamma*(K0[j][i] - dt03*Ma[j][i]); // delta_I_k coeff
        
        Ml[j][i] += dKu[j][i] + Ku[j][i];
        Ml[j][i] = eps*ident[j][i] - beta*Ku[j][i] - theta*(dt03*Ml[j][i] + Ku[j][i]); // delta_I_k+1 coeff --> the one from the previous iteration

        Mc[j][i] = alpha* ident[j][i] + gamma * (ident[j][i] - dt03 * K0[j][i]);
        Mb[j][i] = beta * ident[j][i] + theta * (ident[j][i] + dt03 * Ku[j][i]);
      }
    }
    MatInv(Md[0]);
      
    for (idp=0; idp<input.n_atomic_pars; idp++){
        
      memset(V0, 0, 4*sizeof(double));
      
      for(i = 0;  i < 4;  i++){
        for(j = 0;  j < 4;  j++){
          V0[i] += Ml[i][j]*dI[j][k-dk][idp] + Mb[i][j]*Gu[idp][j] + Mc[i][j]*G0[idp][j];
        }
        V0[i] += dt03*(-gamma*dG0[idp][i] + theta*dGu[idp][i]);
      }
      
      /* --- Solve linear system to get the intensity perturbation -- -------------- */
        
      m4v(Md, V0, V1);      // Multiply Md^-1 * V0

      for(i=0;i<4;i++)
      {
        dI[i][k][idp] = V1[i];
        fprintf(fptr, "%2.8e ", dI[i][k][idp]);
      }
      fprintf(fptr, "\n");
    }
    
    /* --- Shift values for next depth --          ------------------ */
      
    memcpy(Ku[0],   K0[0], 16*sizeof(double));
    memcpy(K0[0],   Kd[0], 16*sizeof(double));
    memcpy(dKu[0], dK0[0], 16*sizeof(double));

    memcpy(Gu[0],   G0[0], input.n_atomic_pars*4*sizeof(double));
    memcpy(G0[0],   Gd[0], input.n_atomic_pars*4*sizeof(double));
    memcpy(dGu[0], dG0[0], input.n_atomic_pars*4*sizeof(double));

    dtau_uw = dtau_dw;
    dsup    = dsdn;
    dchi_up = dchi_c;
    dchi_c  = dchi_dn;      
  }
      
  /* --- Linear integration in the last interval -- ----------------- */
  
  k = k_end;
  dtau_uw = 0.5*imu * (chi[k] + chi[k-dk]) *
    fabs(geometry.height[k] - geometry.height[k-dk]);
  w3(dtau_uw, w);

  /* --- dSu is defined negative in Han's implementation ------------ */

  for (idp=0; idp<input.n_atomic_pars; idp++){

    for (j=0; j<4; j++){
      for (n=0; n<4; n++){
        Ma[j][n] = ident[j][n]*(1.0 - w[0]) - w[1]/dtau_uw*K0[j][n];
        Md[j][n] = ident[j][n] + w[0]*Kd[j][n] - w[1]/dtau_uw*Kd[j][n];
      }
    }

    memset(V0, 0, 4*sizeof(double));
    for (n = 0;  n < 4;  n++){
      for (j=0; j<4; j++){
        V0[n] += Ma[n][j] * dI[n][k-dk][idp];
      }
      V0[n] += w[0]*Gd[idp][n] - w[1]*dG0[idp][n];  
    }
  
    /* --- Solve linear system --                    ------------------ */
    
    MatInv(Md[0]);
    m4v(Md,V0,V1);      // Multiply Md^-1 * V0
    
    for(i=0;i<4;i++) {
      dI[i][k][idp] = V1[i]; 
      fprintf(fptr, "%2.8e ", dI[i][k][idp]);
    }
    fprintf(fptr, "\n");
  }

  fclose(fptr);

  freeMatrix(Gu);
  freeMatrix(G0);
  freeMatrix(Gd);
  freeMatrix(dGu);
  freeMatrix(dG0);
}
/* ------- end ------------------------- Piece_Stokes_Bezier3_1D.c -- */

/* ------- begin ----------------------- Piecewise_Bezier3_1D.c ----- */

void Piecewise_Bezier3_1D(int nspect, int mu, bool_t to_obs,
			  double *chi, double *S, double *I, double *Psi, double **dI)
{
  
  /* --- Cubic Bezier solver for unpolarized light
         Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

         Reference:
         J. de la Cruz Rodriguez & N. Piskunov (2013), Auer (2003)

         Comments: 
           JdlCR: We only check that the control points of the opacity
	          and source function are always above zero to avoid having
	          a negative interpolant.
         --                                            -------------- */
  
  register int k;
  const char routineName[] = "Piecewise_Bezier3_1D";

  int    k_start, k_end, dk, Ndep = geometry.Ndep;
  double dtau_uw, dtau_dw, dS_uw, I_upw, c1, c2, w[3], zmu, Bnu[2];
  double dsup, dsdn, dt03, eps=0, alpha=0, beta=0, gamma=0, theta=0;
  double dS_up, dS_c, dchi_up, dchi_c, dchi_dn, dsdn2;
  int idp;
  double Zk, Zkm1, Zkp1, *dZk, dZkm1, *dZup, *dI_upw;

  zmu = 1.0 / geometry.muz[mu];

  /* --- Distinguish between rays going from BOTTOM to TOP
         (to_obs == TRUE), and vice versa --           -------------- */

  if (to_obs) {
    dk      = -1;
    k_start = Ndep-1;
    k_end   = 0;
  } else {
    dk      = 1;
    k_start = 0;
    k_end   = Ndep-1;
  }
  
  dtau_uw = 0.5 * zmu * (chi[k_start] + chi[k_start+dk]) *
    fabs(geometry.height[k_start] - geometry.height[k_start+dk]);

  /* --- Boundary conditions --                        -------------- */

  if (to_obs) {
    switch (geometry.vboundary[BOTTOM]) {
    case ZERO:
      I_upw = 0.0;
      break;
    case THERMALIZED:
      Planck(2, &atmos.T[Ndep-2], spectrum.lambda[nspect], Bnu, -1);
      I_upw = Bnu[1] - (Bnu[0] - Bnu[1]) / dtau_uw;
      break;
    case IRRADIATED:
      I_upw = geometry.Ibottom[nspect][mu];
      break;
    case REFLECTIVE:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[BOTTOM]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  } else {
    switch (geometry.vboundary[TOP]) {
    case ZERO:
      I_upw = 0.0;
      break;
    case THERMALIZED:
      Planck(2, &atmos.T[0], spectrum.lambda[nspect], Bnu, -1);
      I_upw = Bnu[0] - (Bnu[1] - Bnu[0]) / dtau_uw;
      break;
    case IRRADIATED:
      I_upw = geometry.Itop[nspect][mu];
      break;
    case REFLECTIVE:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[BOTTOM]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  }
  
  I[k_start] = I_upw;
  if (Psi) Psi[k_start] = 0.0;
  
  /* --- Set variables for first iteration to allow simple 
         shift for all next iterations --              -------------- */

  k = k_start+dk;
  dsup = fabs(geometry.height[k] - geometry.height[k-dk]) * zmu;
  dsdn = fabs(geometry.height[k+dk] - geometry.height[k]) * zmu;
  dchi_up = (chi[k] - chi[k-dk]) / dsup;
  
  /* --- dchi/ds at central point --                   -------------- */

  dchi_c = cent_deriv(dsup, dsdn, chi[k-dk], chi[k], chi[k+dk]);
  
  /* --- upwind path_length (Bezier3 integration) --   -------------- */

  c1 = MAX(chi[k]    - (dsup/3.0) * dchi_c,   0.0);
  c2 = MAX(chi[k-dk] + (dsup/3.0) * dchi_up,  0.0);
  dtau_uw = dsup * (chi[k] + chi[k-dk] + c1 + c2) * 0.25;

  /* dS/dtau at upwind point */

  dS_up = (S[k] - S[k-dk]) / dtau_uw;

  dI_upw = NULL;
  dZup = NULL;
  dZk = NULL;
  if (input.get_atomic_rfs  && to_obs){
    dZup = (double *) malloc(input.n_atomic_pars * sizeof(double));
    dZk = (double *) malloc(input.n_atomic_pars * sizeof(double));
    dI_upw = (double *) malloc(input.n_atomic_pars * sizeof(double));
    for (idp=0; idp<input.n_atomic_pars; idp++){
      dI[k_start][idp] = 0.0;
      dI_upw[idp] = 0.0;

      Zk = -spectrum.dchi_c_lam[nspect][k][idp]/chi[k] * I[k] + spectrum.deta_c_lam[nspect][k][idp]/chi[k];
      Zkm1 = -spectrum.dchi_c_lam[nspect][k-dk][idp]/chi[k-dk] * I[k-dk] + spectrum.deta_c_lam[nspect][k-dk][idp]/chi[k-dk];
      dZup[idp] = (Zk - Zkm1) / dtau_uw;
    }
  }

  
  /* --- Solve transfer along ray --                   -------------- */

  for (k = k_start+dk;  k != k_end+dk;  k += dk) {
    
    if (k != k_end) {

      /* --- Downwind path length --                   -------------- */
      
       dsdn = fabs(geometry.height[k+dk] - geometry.height[k]) * zmu;
       
      /* --- dchi/ds at downwind point --             -------------- */
       
      if (fabs(k - k_end) > 1) {
        dsdn2 = fabs(geometry.height[k+2*dk] - geometry.height[k+dk]) * zmu;
        dchi_dn = cent_deriv(dsdn, dsdn2, chi[k], chi[k+dk], chi[k+2*dk]);       
      } else {
        dchi_dn=(chi[k+dk]-chi[k])/dsdn;
      }
       
       /* --- Make sure that c1 and c2 don't go below zero -- ------- */
    
       c1 = MAX(chi[k]    + (dsdn/3.0) * dchi_c,  0.0);
       c2 = MAX(chi[k+dk] - (dsdn/3.0) * dchi_dn, 0.0);

       /* --- Downwind optical path length --          -------------- */

       dtau_dw = dsdn * (chi[k] + chi[k+dk] + c1 + c2) * 0.25;
       dt03    = dtau_uw / 3.0;
      
      /* --- Compute interpolation parameters --       -------------- */
       
       Bezier3_coeffs(dtau_uw, &alpha, &beta, &gamma, &theta, &eps);
       
       /* --- dS/dt at central point --                -------------- */
       
       dS_c = cent_deriv(dtau_uw, dtau_dw, S[k-dk], S[k], S[k+dk]);

       /* --- Source function control points --        -------------- */
       
       c1 = MAX(S[k]    - dt03 * dS_c , 0.0);
       c2 = MAX(S[k-dk] + dt03 * dS_up, 0.0);       
     
       /* --- Solve integral in this interval --       -------------- */
       
       I[k]= I_upw*eps + alpha*S[k] + beta*S[k-dk] + gamma * c1 + theta * c2; 
  
       if (input.get_atomic_rfs && to_obs){
        for (idp=0; idp<input.n_atomic_pars; idp++){
          Zk = -spectrum.dchi_c_lam[nspect][k][idp] * I[k] + spectrum.deta_c_lam[nspect][k][idp];
          Zk /= chi[k];
          Zkm1 = -spectrum.dchi_c_lam[nspect][k-dk][idp] * I[k-dk] + spectrum.deta_c_lam[nspect][k-dk][idp];
          Zkm1 /= chi[k-dk];
          Zkp1 = -spectrum.dchi_c_lam[nspect][k+dk][idp] * I[k+dk] + spectrum.deta_c_lam[nspect][k+dk][idp];
          Zkp1 /= chi[k+dk];
          dZk[idp] = cent_deriv(dtau_uw, dtau_dw, Zkm1, Zk, Zkp1);
          c1 = MAX(Zk - dt03 * dZk[idp], 0.0);
          c2 = MAX(Zkm1 + dt03 * dZup[idp], 0.0);
          dI[k][idp] = dI_upw[idp]*eps + alpha*Zk + beta*Zkm1 + gamma*c1 + theta*c2;
        }
       }
       if (nspect==50 && to_obs && k>52) {
        //  printf("k = %d | I = %2.8e\n", k, I[k]);
       }

       /* --- Diagonal operator --                     -------------- */

       if (Psi) Psi[k] = alpha + gamma;
       
    } else { 
      
      /* --- Piecewise linear integration at end of ray -- ---------- */
      
      dtau_uw = 0.5 * zmu * (chi[k] + chi[k-dk]) * fabs(geometry.height[k] - geometry.height[k-dk]);
      
      /* --- Defined negative in Han's implementation -- ------------ */
      
      dS_uw = -(S[k] - S[k-dk]) / dtau_uw;
      w3(dtau_uw, w);
      
      I[k] = (1.0 - w[0])*I_upw + w[0]*S[k] + w[1]*dS_uw;

      if (input.get_atomic_rfs && to_obs){
        for (idp=0; idp<input.n_atomic_pars; idp++){
          Zk = spectrum.dchi_c_lam[nspect][k][idp]/chi[k] * I[k] - spectrum.deta_c_lam[nspect][k][idp]/chi[k];
          Zkm1 = spectrum.dchi_c_lam[nspect][k-dk][idp]/chi[k-dk] * I[k-dk] - spectrum.deta_c_lam[nspect][k-dk][idp]/chi[k-dk];
          dZk[idp] = -(Zk - Zkm1) / dtau_uw;
          dI[k][idp] = (1.0 - w[0])*dI_upw[idp] + w[0]*Zk + w[1]*dZk[idp];
        }
      }

      /* --- Diagonal operator --                      -------------- */
      
      if (Psi) Psi[k] = w[0] - w[1] / dtau_uw;
    }
    
    /* --- Re-use downwind quantities for next upwind position -- --- */
    
    I_upw = I[k];
    dsup=dsdn;
    dchi_up=dchi_c;
    dchi_c=dchi_dn;
    dtau_uw=dtau_dw;
    dS_up = dS_c;
    if (input.get_atomic_rfs && to_obs){
      for (idp=0; idp<input.n_atomic_pars; idp++){
        dI_upw[idp] = dI[k][idp];
        dZup[idp] = dZk[idp];
      }
    }
  }
  if (dI_upw!=NULL) free(dI_upw);
  if (dZup!=NULL) free(dZup);
  if (dZk!=NULL) free(dZk);
}
/* ------- end ---------------------------- Piecewise_Bezier3_1D.c -- */

/* ------- end ---------------------------- Piecewise_Bezier3_1D_RFs.c -- */
// Vukadinovic: RF Berzier3 solver for non-polarized RTE
void Piecewise_Bezier3_1D_RFs(int nspect, int mu, bool_t to_obs,
        double *chi, double *Z, double *I, double **dI)
{
  
  /* --- Cubic Bezier solver for unpolarized atomic RFs
         Coded by D. Vukadinovic (copy-pasted from non-RF function)
         --                                            -------------- */
  
  register int k;
  const char routineName[] = "Piecewise_Bezier3_1D_RFs";

  int    k_start, k_end, dk, Ndep = geometry.Ndep;
  double dtau_uw, dtau_dw, dS_uw, I_upw, c1, c2, w[3], zmu, Bnu[2];
  double dsup, dsdn, dt03, eps=0, alpha=0, beta=0, gamma=0, theta=0;
  double dchi_up, dchi_c, dchi_dn, dsdn2;
  int idp;
  double Zk, Zkm1, Zkp1, *dZk, dZkm1, *dZup, *dI_upw;

  zmu = 1.0 / geometry.muz[mu];

  /* --- Distinguish between rays going from BOTTOM to TOP
         (to_obs == TRUE), and vice versa --           -------------- */

  if (to_obs) {
    dk      = -1;
    k_start = Ndep-1;
    k_end   = 0;
  } else {
    dk      = 1;
    k_start = 0;
    k_end   = Ndep-1;
  }
  
  dtau_uw = 0.5 * zmu * (chi[k_start] + chi[k_start+dk]) *
    fabs(geometry.height[k_start] - geometry.height[k_start+dk]);
  
  /* --- Set variables for first iteration to allow simple 
         shift for all next iterations --              -------------- */

  k = k_start+dk;
  dsup = fabs(geometry.height[k] - geometry.height[k-dk]) * zmu;
  dsdn = fabs(geometry.height[k+dk] - geometry.height[k]) * zmu;
  dchi_up = (chi[k] - chi[k-dk]) / dsup;
  
  /* --- dchi/ds at central point --                   -------------- */

  dchi_c = cent_deriv(dsup, dsdn, chi[k-dk], chi[k], chi[k+dk]);
  
  /* --- upwind path_length (Bezier3 integration) --   -------------- */

  c1 = MAX(chi[k]    - (dsup/3.0) * dchi_c,   0.0);
  c2 = MAX(chi[k-dk] + (dsup/3.0) * dchi_up,  0.0);
  dtau_uw = dsup * (chi[k] + chi[k-dk] + c1 + c2) * 0.25;

  /* dS/dtau at upwind point */

  dZup = (double *) malloc(input.n_atomic_pars * sizeof(double));
  dZk = (double *) malloc(input.n_atomic_pars * sizeof(double));
  dI_upw = (double *) malloc(input.n_atomic_pars * sizeof(double));
  for (idp=0; idp<input.n_atomic_pars; idp++){
    dI[k_start][idp] = 0.0;
    dI_upw[idp] = 0.0;

    Zk = -spectrum.dchi_c_lam[nspect][k][idp]/chi[k] * I[k] + spectrum.deta_c_lam[nspect][k][idp]/chi[k];
    Zkm1 = -spectrum.dchi_c_lam[nspect][k-dk][idp]/chi[k-dk] * I[k-dk] + spectrum.deta_c_lam[nspect][k-dk][idp]/chi[k-dk];
    dZup[idp] = (Zk - Zkm1) / dtau_uw;
  }
  
  /* --- Solve transfer along ray --                   -------------- */

  for (k = k_start+dk;  k != k_end+dk;  k += dk) {
    
    if (k != k_end) {

      /* --- Downwind path length --                   -------------- */
      
       dsdn = fabs(geometry.height[k+dk] - geometry.height[k]) * zmu;
       
      /* --- dchi/ds at downwind point --             -------------- */
       
      if (fabs(k - k_end) > 1) {
        dsdn2 = fabs(geometry.height[k+2*dk] - geometry.height[k+dk]) * zmu;
        dchi_dn = cent_deriv(dsdn, dsdn2, chi[k], chi[k+dk], chi[k+2*dk]);       
      } else {
        dchi_dn=(chi[k+dk]-chi[k])/dsdn;
      }
       
       /* --- Make sure that c1 and c2 don't go below zero -- ------- */
    
       c1 = MAX(chi[k]    + (dsdn/3.0) * dchi_c,  0.0);
       c2 = MAX(chi[k+dk] - (dsdn/3.0) * dchi_dn, 0.0);

       /* --- Downwind optical path length --          -------------- */

       dtau_dw = dsdn * (chi[k] + chi[k+dk] + c1 + c2) * 0.25;
       dt03    = dtau_uw / 3.0;
      
      /* --- Compute interpolation parameters --       -------------- */
       
      Bezier3_coeffs(dtau_uw, &alpha, &beta, &gamma, &theta, &eps);

      for (idp=0; idp<input.n_atomic_pars; idp++){
        Zk = -spectrum.dchi_c_lam[nspect][k][idp] * I[k] + spectrum.deta_c_lam[nspect][k][idp];
        Zk /= chi[k];
        Zkm1 = -spectrum.dchi_c_lam[nspect][k-dk][idp] * I[k-dk] + spectrum.deta_c_lam[nspect][k-dk][idp];
        Zkm1 /= chi[k-dk];
        Zkp1 = -spectrum.dchi_c_lam[nspect][k+dk][idp] * I[k+dk] + spectrum.deta_c_lam[nspect][k+dk][idp];
        Zkp1 /= chi[k+dk];
        dZk[idp] = cent_deriv(dtau_uw, dtau_dw, Zkm1, Zk, Zkp1);
        c1 = Zk - dt03 * dZk[idp];
        c2 = Zkm1 + dt03 * dZup[idp];
        if (nspect==50 && to_obs && k>52) {
          // printf("k = %d | Zk = %2.8e | Zkm1 = %2.8e | dZk = %2.8e | dZkm1 = %2.8e\n", k, Zk, Zkm1, dZk[idp], dZup[idp]);
          printf("k = %d | dchi = %2.8e | deta = %2.8e | I = %2.8e | chi = %2.8e | Zk = %2.8e\n", k, spectrum.dchi_c_lam[nspect][k][idp], spectrum.deta_c_lam[nspect][k][idp], I[k], chi[k], Zk);
        }
        dI[k][idp] = dI_upw[idp]*eps + alpha*Zk + beta*Zkm1 + gamma*c1 + theta*c2;
      }
      // if (nspect==50) printf("k = %2d | Zu = %2.5e | Zk = %2.5e | Iu = %2.5e | I0 = %2.5e\n", k, Zkm1, Zk, I[k-dk], I[k]);
       
    } else { 
      
      /* --- Piecewise linear integration at end of ray -- ---------- */
      
      dtau_uw = 0.5 * zmu * (chi[k] + chi[k-dk]) * fabs(geometry.height[k] - geometry.height[k-dk]);
      
      /* --- Defined negative in Han's implementation -- ------------ */
      
      w3(dtau_uw, w);
      
      for (idp=0; idp<input.n_atomic_pars; idp++){
        Zk = -spectrum.dchi_c_lam[nspect][k][idp]/chi[k] * I[k] + spectrum.deta_c_lam[nspect][k][idp]/chi[k];
        Zkm1 = -spectrum.dchi_c_lam[nspect][k-dk][idp]/chi[k-dk] * I[k-dk] + spectrum.deta_c_lam[nspect][k-dk][idp]/chi[k-dk];
        dZk[idp] = -(Zk - Zkm1) / dtau_uw;
        dI[k][idp] = (1.0 - w[0])*dI_upw[idp] + w[0]*Zk + w[1]*dZk[idp];
      }
    }
    
    /* --- Re-use downwind quantities for next upwind position -- --- */
    dsup=dsdn;
    dchi_up=dchi_c;
    dchi_c=dchi_dn;
    dtau_uw=dtau_dw;
    for (idp=0; idp<input.n_atomic_pars; idp++){
      dI_upw[idp] = dI[k][idp];
      dZup[idp] = dZk[idp];
    }
  }

  free(dI_upw);
  free(dZup);
  free(dZk);
}
/* ------- end ---------------------------- Piecewise_Bezier3_1D_RFs.c -- */
