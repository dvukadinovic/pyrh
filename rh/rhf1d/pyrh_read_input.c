#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../rh.h"
#include "../inputs.h"
#include "../atom.h"
#include "../atmos.h"
// #include "geometry.h"
// #include "spectrum.h"
// #include "background.h"
// #include "statistics.h"
// #include "error.h"
#include "../../headers/xdr.h"
#include "../constant.h"

#include "../atomweights.h"

#include "pyrh_read_input.h"

extern Atmosphere atmos;
extern InputData input;
extern CommandLine commandline;

void test_InputData(InputData pyrh_input, Atmosphere pyrh_atmos){
  printf("%f %f\n", pyrh_atmos.totalAbund, pyrh_input.abundances[0]);
  printf("%f\n", pyrh_atmos.elements[0].abund);
  printf("%e\n", pyrh_atmos.elements[0].ionpot[0]);
}

void set_elements(InputData pyrh_input, Atmosphere *pyrh_atmos){
  Element *element;

  pyrh_atmos->elements = (Element *) malloc(pyrh_atmos->Nelem * sizeof(Element));
  for (int n = 0;  n < pyrh_atmos->Nelem;  n++) {
    element = &pyrh_atmos->elements[n];
    strcpy(element->ID, atomweight[n].ID);
    element->abundance_set = TRUE;;
    element->abund = pyrh_input.abundances[n];
    element->mol_index = NULL;
    element->Nmolecule = 0;
    element->weight = atomweight[n].weight;
    element->n = NULL;
    element->model = NULL;
    // set from Cython env
    element->Nstage = 0;
    element->pf = NULL;
    element->ionpot = NULL; 
  }
}

void read_atom_model(char *cwd, char *filename){
  char *tmp = malloc(300);
  char *absolute_path = malloc(300);
  bool_t active;

  concatenate(tmp, cwd, "/keyword.input");
  strcpy(commandline.keyword_input, tmp);

  int argc = 1;
  char* argv[] = {};
  setOptions(argc, argv);
 
  readInput();
  input.pyrhHSE = FALSE;
  input.read_atom_model = TRUE;

  readAbundance(&atmos, 0, NULL, NULL);
  
  concatenate(tmp, "/rh/Atoms/", filename);
  concatenate(absolute_path, input.pyrh_path, tmp);

  Atom *atom;
  readAtom(&atom, absolute_path, active=TRUE);

  // printf("Done!\n");

  // free(tmp);
  // free(absolute_path);
}