# This file contains hooks to override default MCell4 model
# code behavior for models generated from CellBlender
# WARNING: do not import file parameters.py at the top level of this file,
# only from individual functions if needed, parameters such as SEED
# will not be set correctly if parameters are imported too early
import sys
import os
import shared
import mcell as m

def custom_argparse_and_parameters():
    # Called before the auto-generated parameter values are set, so values
    # placed in shared.parameter_overrides win over the defaults in
    # parameters.py.
    #
    # Behavior:
    #   * preserves the default optional "-seed N" command-line argument
    #   * additionally allows ITERATIONS / TIME_STEP / SEED to be overridden
    #     via environment variables MCELL_ITERATIONS / MCELL_TIME_STEP /
    #     MCELL_SEED. This lets us run a short validation (e.g.
    #     MCELL_ITERATIONS=350000) without editing parameters.py, leaving the
    #     1 s production config (1e-7 step, 10,000,000 iterations) untouched.
    import sys
    import os

    argv = sys.argv
    if len(argv) == 1:
        pass
    elif len(argv) == 3 and argv[1] == '-seed':
        shared.parameter_overrides['SEED'] = int(argv[2])
    else:
        print("Error: invalid command line arguments")
        print("  usage: " + argv[0] + " [-seed N]")
        sys.exit(1)

    for name, cast in (('SEED', int), ('ITERATIONS', int), ('TIME_STEP', float)):
        raw = os.environ.get('MCELL_' + name)
        if raw is not None:
            shared.parameter_overrides[name] = cast(raw)

"""
def custom_config(model):
    # When uncommented, this function is called to set custom
    # model configuration.
    # It is executed after basic parameter setup is done and 
    # before any components are added to the model. 
    pass
"""

"""
def custom_init_and_run(model):
    # When uncommented, this function is called after all the model
    # components defined in CellBlender were added to the model.
    # It allows to add additional model components before initialization 
    # is done and then to customize how simulation is ran.
    # The module parameters must be imported locally otherwise
    # changes to shared.parameter_overrides done elsewhere won't be applied.
    import parameters as parameters
    model.initialize()
    model.run_iterations(parameters.ITERATIONS)
    model.end_simulation()
"""

