/* ------- file: -------------------------- bezier.h ----------------

   Cubic DELO-Bezier (polarized) and cubic short-char Bezier solvers.
   
   References: de la Cruz Rodriguez & Piskunov (2013), Auer (2003)
               (Derivatives) Fritsch & Butland (1984),
	       
   Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

       --------------------------                      ----------RH-- */


#ifndef BEZIER_H
#define BEZIER_H


/* --- Prototypes auxiliary functions -                -------------- */

/* --- Matrix inversion Cramer method SSE instructions -- ----------- */

void MatInv(float *mat);

/* --- Matrix inversion Shipley-Coleman (not used) --  -------------- */

void m4inv(double MI[4][4]);


void m4v(float a[4][4], double b[4], double c[4]);
void m4m(double a[4][4], double b[4][4], double c[4][4]);

void cent_deriv_vec(double wprime[4], double dsup, double dsdn,
		    double chiup[4], double chic[4], double chidn[4]);

void cent_deriv_mat(double wprime[4][4], double dsup, double dsdn,
		    double chiup[4][4], double chic[4][4],
		    double chidn[4][4]);

double cent_deriv(double dsup,double dsdn, double chiup,
		  double chic, double chidn);
void Svec(int k, double **S, double *Sf);

void Bezier3_coeffs(double dt, double *alpha, double *beta,
		    double *gamma, double *theta, double *eps);

#endif /* !__BEZIER_H__ */

/* ---------------------------------------- bezier.h ---------------- */
