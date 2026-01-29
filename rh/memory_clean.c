#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "inputs.h"
#include "spectrum.h"
#include "rhf1d/geometry.h"

extern Atmosphere atmos;
extern Geometry geometry;
extern InputData input;
extern Spectrum spectrum;

void clean_from_memory()
{
    int nspect;
    ActiveSet *as;

    if (spectrum.I!=NULL) freeMatrix((void **) spectrum.I);
    if (spectrum.Stokes_Q!=NULL) freeMatrix((void **) spectrum.Stokes_Q);
    if (spectrum.Stokes_U!=NULL) freeMatrix((void **) spectrum.Stokes_U);
    if (spectrum.Stokes_V!=NULL) freeMatrix((void **) spectrum.Stokes_V);

    if (spectrum.J!=NULL) freeMatrix((void **) spectrum.J);
    if (input.backgr_pol){
        if (spectrum.J20!=NULL) freeMatrix((void **) spectrum.J20);
    }

    for (int idl=0; idl<atmos.Nrlk; idl++){
        if (atmos.rlk_lines[idl].zm != NULL){
            free(atmos.rlk_lines[idl].zm->q);
            free(atmos.rlk_lines[idl].zm->strength);
            free(atmos.rlk_lines[idl].zm->shift);
            free(atmos.rlk_lines[idl].zm);
        }
    }
    free(atmos.rlk_lines);

    free(atmos.backgrrecno);
    free(atmos.backgrflags);

    if (geometry.Itop!=NULL) freeMatrix((void **) geometry.Itop);
    if (geometry.Ibottom!=NULL) freeMatrix((void **) geometry.Ibottom);

    if (input.get_atomic_rfs) freeMatrix4d(atmos.atomic_rfs, 4, spectrum.Nspect, atmos.Nrays);

    if (spectrum.wave_inds!=NULL) free(spectrum.wave_inds);
    // might have to free as at each wavelength first (for NLTE)
    for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) {
        as = &spectrum.as[nspect];
    
        free(as->Nlower);
        free(as->Nupper);

        free(as->upper_levels);
        free(as->lower_levels);
    }
    if (spectrum.as!=NULL) free(spectrum.as); 

    free(spectrum.lambda);

    freeOpacityEmissivity();
    if (input.get_atomic_rfs) freeOpacityEmissivityDer();

    if (atmos.Nactiveatom > 0) free(atmos.activeatoms);
    if (atmos.Nactivemol > 0) free(atmos.activemols);

    if (atmos.Stokes){
        freeMatrix((void **) atmos.cos_gamma);
        freeMatrix((void **) atmos.cos_2chi);
        freeMatrix((void **) atmos.sin_2chi);
    }

    freeAtoms();
    freeMolecules();
    // have to clear from memory as its a separate pointer
    free(atmos.nHmin);

    // freed inside freeAtoms() when we free Hydrogen (I believe)
    // if (atmos.nH!=NULL) freeMatrix((void **) atmos.nH);

    // free geometry related data
    free(geometry.mux);
    free(geometry.muy);
    free(geometry.muz);
    free(geometry.wmu);
    
    if (geometry.tau_ref!=NULL) free(geometry.tau_ref); geometry.tau_ref = NULL;
    if (geometry.cmass!=NULL) free(geometry.cmass); geometry.cmass = NULL;
    if (geometry.height!=NULL) free(geometry.height); geometry.height = NULL;
    
    free(atmos.N);    
    freeElements();
}