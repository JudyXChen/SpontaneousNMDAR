"""Compare source-side occupancy using 3D positions at matched 0.1 ms times.

The plane bisecting the two cluster centers defines side A. The difference
P(Ca_A on side A)-P(Ca_B on side A) is 1 for perfect source-side separation and
0 for equal occupancy of this half-space. Zero is not proof of identical 3D
distributions. This descriptive statistic is not an ensemble test of mixing.
"""
import sys as _sys
from pathlib import Path as _Path
_here = _Path(__file__).resolve().parent
_sys.path[:0] = [str(_here.parent), str(_here.parent / 'chemistry')]  # source layout; snapshots are flat
import argparse
import csv
import json
from pathlib import Path
from paths import inputs_for, verify_inputs, repository
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from geometry_checks import check_placement, read_mesh

sys.path.insert(0,str(repository()/'scripts'))
import plot_style as style


def positions(path):
    result = {'A':[], 'B':[]}
    for line in path.read_text().splitlines():
        fields = line.split()
        if fields[0] in ('Ca_A','Ca_B'):
            result[fields[0][-1]].append([float(s) for s in fields[2:5]])
    return {tag:np.asarray(rows).reshape(-1,3) for tag,rows in result.items()}


def summarize(run):
    config = json.loads((inputs_for(run)/'condition.json').read_text())
    assert json.loads((run/'validation.json').read_text())['validation']=='passed'
    layout = check_placement(config,inputs_for(run)/'geometry')
    centers = np.asarray([s['region_center_um'] for s in layout['sites']])
    mid = centers.mean(axis=0)
    axis = centers[1]-centers[0]; axis /= np.linalg.norm(axis)
    rows = []
    seed = json.loads((run/'status.json').read_text())['seed']
    frames = sorted((run/f'viz_data/seed_{seed:05d}').glob('Scene.ascii.*.dat'))
    for path in frames:
        iteration = int(path.name.split('.')[-2])
        if iteration % 1000:
            continue
        xyz = positions(path)
        fractions = {tag:float(np.mean((points-mid)@axis<0)) if len(points) else np.nan
                     for tag,points in xyz.items()}
        rows.append(dict(time_ms=iteration*config['time_step_s']*1000,
                         n_A=len(xyz['A']),n_B=len(xyz['B']),
                         fraction_A_on_A_side=fractions['A'], fraction_B_on_A_side=fractions['B'],
                         side_fraction_difference=fractions['A']-fractions['B']))
    return rows, centers, frames, config


def preview(run, centers, frames, config, destination):
    # A matched time avoids comparing one condition's early frame to another's
    # late peak. Existing per-run peak previews remain separate artifacts.
    iteration = round((config['waveform_start_s']+0.010)/config['time_step_s'])  # 10 ms after release
    path = next(p for p in frames if int(p.name.split('.')[-2])==iteration)
    xyz = positions(path)
    mid = centers.mean(axis=0)
    horizontal = (centers[1]-centers[0])[:2]
    horizontal /= np.linalg.norm(horizontal)
    vertices,walls,_ = read_mesh(inputs_for(run)/'geometry','Mushroom_transparent')
    vertices,walls = np.asarray(vertices),np.asarray(walls)
    projected = np.column_stack(((vertices[:,:2]-mid[:2])@horizontal,vertices[:,2]))
    edges = np.concatenate([walls[:,[0,1]],walls[:,[1,2]],walls[:,[2,0]]])
    fig,ax = style.figure()
    ax.add_collection(LineCollection(projected[edges],color='#999999',linewidths=.2,alpha=.15))
    colors={'A':'#E52B88','B':'#169B46'}
    for tag,points in xyz.items():
        ax.scatter((points[:,:2]-mid[:2])@horizontal,points[:,2],s=7,color=colors[tag],alpha=.65,linewidths=0,label='Ca_'+tag)
    for tag,center in zip('AB',centers):
        ax.scatter([(center[:2]-mid[:2])@horizontal],[center[2]],marker='x',s=35,color=colors[tag],linewidths=1.3)
    ax.set(xlabel='Along cluster separation (µm)',ylabel='z (µm)',xlim=(-.75,.75),ylim=(.35,1.75))
    ax.set_aspect('equal')
    time_ms=iteration*config['time_step_s']*1000
    style.save(fig,destination/f'source_calcium_{config["distance_nm"]}nm_{time_ms:.3f}ms.png')
    return dict(time_ms=time_ms,iteration=iteration,
                selection='Matched 10 ms snapshot for both conditions.',
                projection='Horizontal axis follows the xy cluster-separation direction; orthogonal horizontal coordinate collapsed.',
                cluster_centers_um=centers.tolist())


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_100nm',type=Path)
    parser.add_argument('run_200nm',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    style.apply_style()
    fig,ax=style.figure()
    summaries={}
    for run,color in [(args.run_100nm,'#2E6F9D'),(args.run_200nm,'#CC762B')]:
        rows,centers,frames,config=summarize(run)
        distance=config['distance_nm']
        with (args.output/f'spatial_{distance}nm.csv').open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        # Sparse frames are retained in the CSV but masked in the plot.
        values=[r['side_fraction_difference'] if min(r['n_A'],r['n_B'])>=20 else np.nan for r in rows]
        ax.plot([r['time_ms'] for r in rows],values,color=color,label=f'{distance} nm')
        summaries[str(distance)]=dict(final=rows[-1],preview=preview(run,centers,frames,config,args.output))
    ax.axhline(0,color='#777777',linewidth=.7,linestyle='--')
    ax.set(xlabel='Time (ms)',ylabel='Source-side fraction difference',xlim=(0,35),ylim=(-1,1))
    style.save(fig,args.output/'source_side_comparison.png')
    summaries['interpretation']=__doc__
    summaries['sampling']='Both runs sampled every 0.1 ms for comparison; full 200 nm movie is 10 us. Plot masks frames with fewer than 20 free ions of either label.'
    (args.output/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
