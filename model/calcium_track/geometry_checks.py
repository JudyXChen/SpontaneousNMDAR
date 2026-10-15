"""Verify cluster triangles and glutamate locations against the actual cleft."""
import ast
import math
from pathlib import Path


def read_mesh(directory, name):
    tree = ast.parse((Path(directory)/(name+'.py')).read_text())
    values = {n.targets[0].id: ast.literal_eval(n.value) for n in tree.body
              if isinstance(n, ast.Assign) and isinstance(n.value, ast.List)}
    return values[name+'_vertex_list'], values[name+'_wall_list'], values


def heights(vertices, walls, x, y):
    result = []
    for indices in walls:
        a,b,c = [vertices[i] for i in indices]
        ux,uy = b[0]-a[0],b[1]-a[1]
        vx,vy = c[0]-a[0],c[1]-a[1]
        determinant = ux*vy-uy*vx
        if abs(determinant)<1e-14:
            continue
        wx,wy = x-a[0],y-a[1]
        u = (wx*vy-wy*vx)/determinant
        v = (ux*wy-uy*wx)/determinant
        if u>=-1e-8 and v>=-1e-8 and u+v<=1+1e-8:
            result.append(a[2]+u*(b[2]-a[2])+v*(c[2]-a[2]))
    return sorted(result)


def check_placement(config, geometry):
    vertices,walls,values = read_mesh(geometry,'Mushroom_dendrite')
    pre_vertices,pre_walls,_ = read_mesh(geometry,'presynapse')
    psd = set(values['Mushroom_dendrite_PSD_wall_indices'])
    sites = []
    for cluster in config['clusters']:
        if not set(cluster['walls'])<=psd:
            raise ValueError('Cluster wall is outside the PSD')
        points = [vertices[i] for w in cluster['walls'] for i in walls[w]]
        center = [sum(p[k] for p in points)/len(points) for k in range(3)]
        lower = max(heights(vertices,walls,*center[:2]))
        upper = min(z for z in heights(pre_vertices,pre_walls,*center[:2]) if z>lower)
        release_z = cluster.get('glutamate_z_um',config['glutamate_z_um'])
        if not lower+1e-8 < release_z < upper-1e-8:
            raise ValueError(f"{cluster['population']} glutamate z={release_z:.9g} is outside the local cleft ({lower:.9g}, {upper:.9g}) um")
        sites.append(dict(population=cluster['population'],walls=cluster['walls'],
            region_center_um=center,glutamate_location_um=[*center[:2],release_z],
            spine_surface_z_um=lower,presynaptic_surface_z_um=upper,
            release_above_spine_nm=(release_z-lower)*1000,
            release_below_presynapse_nm=(upper-release_z)*1000))
    a,b = [s['region_center_um'] for s in sites]
    return dict(sites=sites,center_distance_3d_nm=math.dist(a,b)*1000,
                center_distance_xy_nm=math.dist(a[:2],b[:2])*1000)
