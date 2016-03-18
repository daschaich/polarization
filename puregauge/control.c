// -----------------------------------------------------------------
// Main procedure for pure-gauge polarization computation
#define CONTROL
#include "puregauge_includes.h"

int main(int argc, char *argv[]) {
  int prompt, s, x, y, z, t, mu, i, j;
  Real re, im;
  double dssplaq, dstplaq;

  // Setup
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  // Remap standard I/O
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  g_sync();
  prompt = setup();

  // Load input and run (loop removed)
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }

  // Serial code!
  if (this_node != 0) {
    node0_printf("ERROR: run this thing in serial!\n");
    terminate(1);
  }

  // Compute field strength tensor at each site
  make_field_strength(F_OFFSET(link), F_OFFSET(FS));
  FORALLSITES(i, s) {
    // TODO...
  }

  fflush(stdout);
  return 0;
}
// -----------------------------------------------------------------
