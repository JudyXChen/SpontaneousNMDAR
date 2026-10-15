"""Make a Ca_A/Ca_B-only Blender recording without rerunning or altering MCell.

Preserve original ASCII rows byte-for-byte, including particle IDs/coordinates.
Keep empty frames: an empty early frame means no free source calcium yet.
The manual CellBlender loader enumerates every file, so its frames directory
must contain only particle files (the import JSON and manifest live above it).
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from paths import inputs_for, verify_inputs, repository

SPECIES = {'Ca_A', 'Ca_B'}


def export(run, seed=1, destination=None):
    run = Path(run).resolve()
    config = json.loads((inputs_for(run) / 'condition.json').read_text())
    validation = json.loads((run / 'validation.json').read_text())
    if validation['validation'] != 'passed':
        raise ValueError('Use a completed, validated run.')
    source = run / f'viz_data/seed_{seed:05d}'
    paths = sorted(source.glob('Scene.ascii.*.dat'))
    iterations = list(range(0, round(config['duration_s'] / config['time_step_s']) + 1,
                            config['viz_every_n_timesteps']))
    if [int(p.name.split('.')[-2]) for p in paths] != iterations:
        raise ValueError('Incomplete position frame sequence')
    counts = {}
    for name in SPECIES:
        trace = run / f'react_data/seed_{seed:05d}/{name}_World.dat'
        counts[name] = {round(float(t) / config['time_step_s']): int(float(n))
                        for t, n in (line.split() for line in trace.read_text().splitlines())}
    destination = Path(destination).resolve() if destination else run / 'blender_ca_only'
    destination.mkdir(parents=True, exist_ok=False)
    frames = destination / 'frames'
    frames.mkdir()
    manifest = []
    for path, iteration in zip(paths, iterations):
        original = path.read_bytes()
        rows = [line for line in original.splitlines(keepends=True)
                if line.split() and line.split()[0].decode() in SPECIES]
        kept = b''.join(rows)
        tally = Counter(line.split()[0].decode() for line in rows)
        for name in SPECIES:
            if tally[name] != counts[name][iteration]:
                raise ValueError(f'Count mismatch: {path.name} {name}')
        target = frames / path.name
        target.write_bytes(kept)
        if target.read_bytes() != kept:
            raise IOError('Filtered frame verification failed')
        manifest.append(dict(frame=path.name, iteration=iteration,
            time_s=iteration * config['time_step_s'],
            counts={name: tally[name] for name in sorted(SPECIES)},
            source_sha256=hashlib.sha256(original).hexdigest(),
            output_sha256=hashlib.sha256(kept).hexdigest()))
    model_path = next(source.glob('Scene.data_model.*.json'))
    data = json.loads(model_path.read_text())
    model = data['mcell']
    for molecule in model['define_molecules']['molecule_list']:
        molecule['export_viz'] = molecule['mol_name'] in SPECIES
    model['viz_output']['export_all'] = False
    model['mol_viz'].update(manual_select_viz_dir=True, file_dir=str(frames),
        file_index=0, file_start_index=0, file_stop_index=len(paths)-1,
        file_num=len(paths), file_step_index=1, file_name=paths[0].name, viz_list=[])
    (destination / model_path.name).write_text(json.dumps(data, indent=1) + '\n')
    report = dict(species=sorted(SPECIES), frames=len(paths),
        original_run=str(run), seed=seed, validation='passed',
        method='Exact ASCII row filtering; IDs, positions and timing unchanged. Empty frames retained.',
        model_source_sha256=hashlib.sha256(model_path.read_bytes()).hexdigest(),
        frame_manifest=manifest)
    (destination / 'manifest.json').write_text(json.dumps(report, indent=2) + '\n')
    print(f'Validated {len(paths)} Ca_A/Ca_B-only frames: {destination}')
    return destination


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory', type=Path)
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--destination', type=Path)
    args = parser.parse_args()
    export(args.run_directory, args.seed, args.destination)
