// -----------------------------------------------------------------
// SU(N) pure-gauge polarization setup
#include "puregauge_includes.h"

#define IF_OK if (status == 0)

// Each node has a params structure for passing simulation parameters
#include "params.h"
params par_buf;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// On node zero, read lattice size and send to others
int initial_set() {
  int prompt = 0, status = 0;
  if (mynode() == 0) {
    // Print banner
    printf("Pure-gauge polarization, Nc = %d\n", NCOL);
    printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());
    time_stamp("start");
    status = get_prompt(stdin, &prompt);

    IF_OK status += get_i(stdin, prompt, "nx", &par_buf.nx);
    IF_OK status += get_i(stdin, prompt, "ny", &par_buf.ny);
    IF_OK status += get_i(stdin, prompt, "nz", &par_buf.nz);
    IF_OK status += get_i(stdin, prompt, "nt", &par_buf.nt);

    if (status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;
  }

  // Broadcast parameter buffer from node 0 to all other nodes
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
  if (par_buf.stopflag != 0)
    normal_exit(0);

  nx = par_buf.nx;
  ny = par_buf.ny;
  nz = par_buf.nz;
  nt = par_buf.nt;

  this_node = mynode();
  number_of_nodes = numnodes();
  volume = nx * ny * nz * nt;
  return prompt;
}
// -----------------------------------------------------------------




// -----------------------------------------------------------------
int setup() {
  int prompt;

  // Print banner, get volume
  prompt = initial_set();
  // Initialize the layout functions, which decide where sites live
  setup_layout();
  // Allocate space for lattice, set up coordinate fields
  make_lattice();
  // Set up neighbor pointers and comlink structures
  make_nn_gathers();
  // Allocate space for temporary su3_matrix field
  FIELD_ALLOC(tempmat1, su3_matrix);

  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read in parameters for SU(N) pure-gauge polarization
// prompt=1 indicates prompts are to be given for input
int readin(int prompt) {
  int status;

  // On node zero, read parameters and send to all other nodes
  if (this_node == 0) {
    printf("\n\n");
    status = 0;

    // Wilson flow parameters
    IF_OK status += get_f(stdin, prompt, "epsilon", &par_buf.epsilon);
    if (par_buf.epsilon == 0) {
      node0_printf("ERROR: epsilon=%g won't get you anywhere\n",
                   par_buf.epsilon);
      status++;
    }
    IF_OK status += get_f(stdin, prompt, "tmax", &par_buf.tmax);
    if (par_buf.epsilon * par_buf.tmax < 0)
      node0_printf("WARNING: epsilon and tmax have different signs\n");

    // Find out what kind of starting lattice to use
    IF_OK status += ask_starting_lattice(stdin, prompt, &par_buf.startflag,
                                         par_buf.startfile);

    if (status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;
  }

  // Broadcast parameter buffer from node0 to all other nodes
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
  if (par_buf.stopflag != 0)
    normal_exit(0);

  epsilon = par_buf.epsilon;
  tmax    = par_buf.tmax;

  startflag = par_buf.startflag;
  strcpy(startfile, par_buf.startfile);

  // Do whatever is needed to get lattice
  startlat_p = reload_lattice(startflag, startfile);
  return 0;
}
// -----------------------------------------------------------------
