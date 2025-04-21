#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "inputs.h"
#include "atom.h"
#include "atmos.h"
// #include "geometry.h"
// #include "spectrum.h"
// #include "background.h"
// #include "statistics.h"
// #include "error.h"
// #include "xdr.h"
// #include "constant.h"

#include "pyrh_read_input.h"

extern Atmosphere atmos;
extern InputData input;
extern CommandLine commandline;

void read_keyword(){
	
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