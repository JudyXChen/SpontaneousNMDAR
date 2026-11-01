"""Audit actual exported MCell chemistry, particle placement and calcium balance.

Reads completed or initialized snapshots; --initialized-only checks chemistry and
initial placement without requiring a completed response. No simulation edits.
"""
import sys as _sys
from pathlib import Path as _Path
_here = _Path(__file__).resolve().parent
_sys.path[:0] = [str(_here.parent), str(_here.parent / 'chemistry')]  # source layout; snapshots are flat
import argparse
import hashlib
import json
from pathlib import Path
from paths import inputs_for, verify_inputs, repository
import re
import numpy as np
from calcium_labels import CARRIERS, labeled, relabel
from geometry_checks import check_placement


def on_triangle(point, triangle, tolerance=1e-6):
    a, b, c = triangle
    u, v, w = b-a, c-a, point-a
    normal = np.cross(u, v)
    if abs(np.dot(w, normal)) / np.linalg.norm(normal) > tolerance:
        return False
    coef = np.linalg.lstsq(np.column_stack((u,v)), w, rcond=None)[0]
    return bool(np.all(coef >= -tolerance) and coef.sum() <= 1+tolerance)


def audit(run, initialized_only=False):
    run = Path(run).resolve()
    config = json.loads((inputs_for(run)/'condition.json').read_text())
    placement = check_placement(config,inputs_for(run)/'geometry')
    seed = 1
    status = run/'status.json'
    if status.exists():
        seed = json.loads(status.read_text())['seed']
    verify_inputs(inputs_for(run))
    viz = run/f'viz_data/seed_{seed:05d}'
    dm = json.loads(next(viz.glob('Scene.data_model.*.json')).read_text())['mcell']
    molecules = {m['mol_name']:m for m in dm['define_molecules']['molecule_list']}
    reactions = {r['rxn_name']:r for r in dm['define_reactions']['reaction_list']}
    releases = {r['name']:r for r in dm['release_sites']['release_site_list']}
    for index,site in enumerate(placement['sites'],1):
        if config['clusters'][index-1]['glutamate'] == 0:
            assert 'glutamate_cluster_'+str(index) not in releases
            continue
        release = releases['glutamate_cluster_'+str(index)]
        xyz = [float(release['location_'+axis]) for axis in 'xyz']
        assert np.allclose(xyz,site['glutamate_location_um'],atol=1e-12,rtol=0)
        assert int(release['quantity']) == config['clusters'][index-1]['glutamate']
    assert dm['use_bng_units'] is False
    for carrier in CARRIERS:
        for tag in 'AB':
            for key in ('mol_type','diffusion_constant','custom_time_step'):
                assert molecules[labeled(carrier,tag)][key] == molecules[carrier][key], (carrier,tag,key)
    for prefix in ('Ca2_ER', 'SERCA', 'RyR', 'IP3R'):
        assert any(name.startswith(prefix) for name in molecules), prefix
    assert 'SA' in [o['name'] for o in dm['geometrical_objects']['object_list']]
    for tag in 'AB':
        for name, rule in reactions.items():
            if name.endswith('_'+tag) and not name.startswith('rxn_NMDAR'):
                baseline = reactions[name[:-2]]
                for key in ('reactants','products'):
                    assert rule[key] == relabel(baseline[key],tag), name
                assert rule['fwd_rate'] == baseline['fwd_rate'], name
        entry = reactions['rxn_NMDAR_'+tag]
        assert entry['reactants'] == f"NMDAR_{tag}'"
        assert entry['products'] == f"NMDAR_{tag}' + Ca_{tag}@IN"
        for name, rule in reactions.items():
            if name.startswith('NMDAR_A_gating') or name.startswith('rxn_NMDAR_A'):
                other = reactions[name.replace('NMDAR_A','NMDAR_B')]
                for key in ('reactants','products','fwd_rate','variable_rate_text'):
                    assert rule[key].replace('NMDAR_A','NMDAR_B').replace('Ca_A','Ca_B') == other[key], name
    assert reactions['rxn_VGCC']['products'] == "VGCC' + Ca2@IN"
    for name in ('rxn_PMCA_leak','rxn_NCX_leak'):
        assert 'Ca2@IN' in reactions[name]['products']
        assert not any(name+'_'+tag in reactions for tag in 'AB')
    geometry = next(o for o in dm['geometrical_objects']['object_list'] if o['name']=='Mushroom_dendrite')
    vertices = np.asarray(geometry['vertex_list'])
    walls = np.asarray(geometry['element_connections'])
    first = min(viz.glob('Scene.ascii.*.dat'))
    initial = {}
    for line in first.read_text().splitlines():
        name, identifier, *rest = line.split()
        if name.startswith(('NMDAR_A','NMDAR_B')):
            initial[int(identifier)] = (name, np.asarray(rest[:3],dtype=float))
    centers, placed_centers, positions = [], [], {}
    for tag, cluster in zip('AB',config['clusters']):
        population = 'NMDAR_'+tag
        assert float(molecules[population]['diffusion_constant']) == 0
        triangles = vertices[walls[cluster['walls']]]
        pts = [pos for name,pos in initial.values() if name == population+'_0']
        assert len(pts) == cluster['receptors']
        assert all(any(on_triangle(p,t) for t in triangles) for p in pts), population+' placement'
        centers.append(triangles.reshape(-1,3).mean(axis=0))
        placed_centers.append(np.mean(pts,axis=0))
        positions[tag] = [p.tolist() for p in pts]
    report = dict(validation='passed', run=str(run),
        release_placement=placement,
        checks=['Original snapshot hashes intact','Exported Ca_A/B diffusion and occupied-state kinetics equal Ca2',
                'NMDAR populations use identical gating and voltage tables; each emits only its own label',
                'VGCC/pump leak and ER release calcium remain Ca2; pumps, buffers and ER channels shared','ER species and SA object present',
                'Exported glutamate release coordinates match configuration and lie in the local cleft',
                'All initial NMDAR particles lie on their assigned cluster triangles'],
        region_center_distance_nm=float(np.linalg.norm(centers[1]-centers[0])*1000),
        receptor_centroid_distance_nm=float(np.linalg.norm(placed_centers[1]-placed_centers[0])*1000),
        region_centers_um=[p.tolist() for p in centers], receptor_positions_um=positions,
        free_calcium_diffusion_um2_per_ms=float(molecules['Ca_A']['diffusion_constant'])*1e5)
    D = report['free_calcium_diffusion_um2_per_ms']
    report['unbounded_diffusion_scale'] = dict(
        rms_over_output_interval_nm=float(np.sqrt(6*D*config['viz_every_n_timesteps']*config['time_step_s']*1000)*1000),
        time_for_100nm_rms_us=0.1**2/(6*D)*1000,
        time_for_200nm_rms_us=0.2**2/(6*D)*1000,
        caveat='Free unbounded 3D RMS scale, not a first-passage or measured spine mixing time; excludes binding and boundaries.')
    if not initialized_only:
        data = run/f'react_data/seed_{seed:05d}'
        def trace(name):
            return np.loadtxt(data/(name+'_World.dat'))[:,1]
        balance = {}
        for tag in 'AB':
            entry = trace('NMDAR_'+tag+'_Ca_flux')
            free_inside = np.loadtxt(data/f'Ca_{tag}_Mushroom_dendrite.dat')[:,1]
            assert np.array_equal(trace('Ca_'+tag),free_inside),tag+' free calcium outside cytoplasm'
            inventory = sum(trace(labeled(c,tag)) for c in CARRIERS)
            removed = sum(trace(s+'_'+tag) for s in ('Ca2Null','rxn_PMCA','rxn_NCX'))
            removed = removed + sum(trace(f'Ca_{tag}_to_{p}') for p in ('SERCA','RyR','IP3R'))
            residual = entry-inventory-removed
            assert np.all(residual==0), tag+' balance'
            assert np.all(trace('NMDAR_'+tag+'_total')==config['clusters']['AB'.index(tag)]['receptors'])
            balance[tag] = dict(entry=int(entry[-1]), final_inventory=int(inventory[-1]),
                                final_removed=int(removed[-1]), max_abs_residual=int(abs(residual).max()))
        for partner in ('Bm','Bf','PMCA','NCX','PLC','RyR','IP3R','SERCA'):
            total = trace(partner+'_total')
            assert np.all(total==total[0]), partner
        # Surface populations must stay fixed through ALL gating states/frames.
        # MCell may assign new particle IDs on a state-changing reaction.
        # Compare spatial multisets by population, not IDs across reactions.
        for path in viz.glob('Scene.ascii.*.dat'):
            seen = {'A':[], 'B':[]}
            for line in path.read_text().splitlines():
                fields = line.split()
                if fields[0].startswith(('NMDAR_A','NMDAR_B')):
                    tag = fields[0][len('NMDAR_')]
                    seen[tag].append([float(x) for x in fields[2:5]])
            for tag in 'AB':
                assert len(seen[tag])==len(positions[tag]),path.name
                assert np.allclose(sorted(seen[tag]),sorted(positions[tag]),atol=1e-8,rtol=0),path.name
        report['balance'] = balance
        report['checks'].extend(['Exact calcium balance and shared-partner conservation at every count sample',
                                 'All free Ca_A/B remains inside the spine plus dendrite cytoplasmic mesh',
                                 'Receptor populations and positions conserved across all saved frames; IDs may change on reactions'])
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory',type=Path)
    parser.add_argument('--initialized-only',action='store_true')
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    report = audit(args.run_directory,args.initialized_only)
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='receptor_positions_um'},indent=2))
