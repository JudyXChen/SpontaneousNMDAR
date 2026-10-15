"""Runtime overrides for shared models; no changes to historical exports."""
import os
import sys
import shared
from condition import CONFIG


def custom_argparse_and_parameters():
    if len(sys.argv) == 3 and sys.argv[1] == '-seed':
        seed = int(sys.argv[2])
    elif len(sys.argv) == 1:
        seed = 1
    else:
        raise ValueError('usage: model.py [-seed N]')
    shared.parameter_overrides.update(
        SEED=seed, ITERATIONS=round(CONFIG['duration_s'] / CONFIG['time_step_s']),
        TIME_STEP=CONFIG['time_step_s'],
        n_NMDAR=sum(c['receptors'] for c in CONFIG['clusters']),
        n_glu=CONFIG['clusters'][0]['glutamate'],
    )
    for name, cast in [('SEED', int), ('ITERATIONS', int), ('TIME_STEP', float)]:
        if 'MCELL_' + name in os.environ:
            shared.parameter_overrides[name] = cast(os.environ['MCELL_' + name])
