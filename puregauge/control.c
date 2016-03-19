// -----------------------------------------------------------------
// Main procedure for pure-gauge polarization computation
#define CONTROL
#include "puregauge_includes.h"

int main(int argc, char *argv[]) {
  register int i;
  register site *s;
  int prompt;
  double ss_plaq, st_plaq, td;

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

  // Check plaquette
  plaquette(&ss_plaq, &st_plaq);
  td = 0.5 * (ss_plaq + st_plaq);
  node0_printf("START %.8g %.8g %.8g\n", ss_plaq, st_plaq, td);

  // Compute field strength tensor at each site
  make_field_strength(F_OFFSET(link), F_OFFSET(FS));
  FORALLSITES(i, s) {
    // TODO...
  }

  fflush(stdout);
  return 0;
}
// -----------------------------------------------------------------
