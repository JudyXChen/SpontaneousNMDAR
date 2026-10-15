"""Runtime overrides for shared models; no changes to historical exports."""
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


def custom_init_and_run(model):
    """Export a directly importable pink/green CellBlender data model."""
    import json
    from pathlib import Path
    from parameters import ITERATIONS, SEED
    from calcium_labels import CARRIERS, labeled
    model.initialize()
    model.export_data_model()
    folder = Path(f'viz_data/seed_{SEED:05d}')
    for path in folder.glob('*.json'):
        data = json.loads(path.read_text())
        for molecule in data['mcell']['define_molecules']['molecule_list']:
            name = molecule['mol_name']
            for tag, color in [('A', [1.0, 0.12, 0.55]), ('B', [0.1, 0.9, 0.2])]:
                if name in [labeled(c, tag) for c in CARRIERS] or name.startswith('NMDAR_' + tag):
                    molecule['display']['color'] = color
                    molecule['display']['glyph'] = 'Sphere_1'
                    molecule['display']['scale'] = 1.0
        path.write_text(json.dumps(data, indent=1) + '\n')
    model.run_iterations(ITERATIONS)
    model.end_simulation()
