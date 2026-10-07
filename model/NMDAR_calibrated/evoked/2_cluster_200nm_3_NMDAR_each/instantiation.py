# WARNING: This is an automatically generated file and will be overwritten
#          by CellBlender on the next model export.

import os
import shared
import mcell as m

from parameters import *
from subsystem import *
from geometry import *
MODEL_PATH = os.path.dirname(os.path.abspath(__file__))

# Two single-triangle PSD regions are 200.40 nm apart and each supports three
# NMDARs. The first region is shared by all 3-NMDAR distance cases.
CLUSTER_WALL_INDICES = ((1681,), (1335,))
GLU_Z = 1.655
GLU_SITES = os.environ.get('MCELL_GLU_SITES', 'both')
if GLU_SITES not in ('1', '2', 'both'):
    raise ValueError('MCELL_GLU_SITES must be 1, 2, or both')
if n_NMDAR != 6:
    raise ValueError('The two-cluster layout requires n_NMDAR = 6')


def region_center(indices):
    vertices = [Mushroom_dendrite_vertex_list[i]
                for index in indices for i in Mushroom_dendrite_wall_list[index]]
    return tuple(sum(point[axis] for point in vertices) / len(vertices) for axis in range(3))


if not set(CLUSTER_WALL_INDICES[0] + CLUSTER_WALL_INDICES[1]).issubset(Mushroom_dendrite_PSD_wall_indices):
    raise ValueError('NMDAR cluster walls must belong to the PSD')
CLUSTER_CENTERS = tuple(region_center(indices) for indices in CLUSTER_WALL_INDICES)
NMDAR_cluster_regions = [
    m.SurfaceRegion(name='NMDAR_cluster_' + str(i), wall_indices=list(indices))
    for i, indices in enumerate(CLUSTER_WALL_INDICES, 1)
]
Mushroom_dendrite.surface_regions = (
    list(Mushroom_dendrite.surface_regions) + NMDAR_cluster_regions
)


# ---- instantiation ----

# ---- release sites ----

# ---- surface classes assignment ----

Mushroom_dendrite_PSD.surface_class = PSD
# OuterCube
Cube.surface_class = Surface_cube
# Mushroom_dendrite
Spine_base.surface_class = Ions_mushroom
# Cylinder_side1_surface
cylinder_side1.surface_class = Cylinder_side
# Cylinder_side2_surface
cylinder_side2.surface_class = Cylinder_side
head_neck.surface_class = Ions_mushroom
# SA_reflective
SA.surface_class = SA_membrane
# Mushroom_transparent
Mushroom_transparent.surface_class = Ions_mushroom
# Cylinder_reflective
Mushroom_dendrite_cylinder_reflect.surface_class = Cylinder_reflective
# ---- compartments assignment ----

Release_Site_1 = m.ReleaseSite(
    name = 'Release_Site_1',
    complex = m.Complex('Ca2'),
    region = Mushroom_transparent - SA,
    concentration = Ca_cyto_IC
)

Release_Site_2 = m.ReleaseSite(
    name = 'Release_Site_2',
    complex = m.Complex('PMCA', orientation = m.Orientation.UP),
    region = Mushroom_dendrite_head,
    density = PMCA_IC_HEAD
)

Release_Site_3 = m.ReleaseSite(
    name = 'Release_Site_3',
    complex = m.Complex('PMCA', orientation = m.Orientation.UP),
    region = Mushroom_dendrite_neck,
    density = PMCA_IC_NECK
)

Release_Site_4 = m.ReleaseSite(
    name = 'Release_Site_4',
    complex = m.Complex('NCX', orientation = m.Orientation.UP),
    region = Mushroom_dendrite_head,
    density = NCX_IC_HEAD
)

Release_Site_5 = m.ReleaseSite(
    name = 'Release_Site_5',
    complex = m.Complex('NCX', orientation = m.Orientation.UP),
    region = Mushroom_dendrite_neck,
    density = NCX_IC_NECK
)

Release_Site_6 = m.ReleaseSite(
    name = 'Release_Site_6',
    complex = m.Complex('Bm'),
    region = Mushroom_transparent - SA,
    concentration = Bm_IC
)

Release_Site_7 = m.ReleaseSite(
    name = 'Release_Site_7',
    complex = m.Complex('Bf', orientation = m.Orientation.UP),
    region = Mushroom_dendrite_spine,
    density = Bf_IC
)

Release_Site_8 = m.ReleaseSite(
    name = 'Release_Site_8',
    complex = m.Complex('VSCC_0', orientation = m.Orientation.UP),
    region = Mushroom_dendrite_spine,
    density = VSCC_IC
)

NMDAR_cluster_sites = [
    m.ReleaseSite(
        name = 'NMDAR_cluster_' + str(cluster_number),
        complex = m.Complex('NMDAR_0', orientation = m.Orientation.UP),
        region = region,
        number_to_release = n_NMDAR // 2
    )
    for cluster_number, region in enumerate(NMDAR_cluster_regions, 1)
]

Release_Site_10 = m.ReleaseSite(
    name = 'Release_Site_10',
    complex = m.Complex('AMPAR_0', orientation = m.Orientation.UP),
    region = Mushroom_dendrite_PSD,
    number_to_release = n_AMPAR
)

Release_Site_11 = m.ReleaseSite(
    name = 'Release_Site_11',
    complex = m.Complex('glu'),
    shape = m.Shape.SPHERICAL,
    location = (CLUSTER_CENTERS[0][0], CLUSTER_CENTERS[0][1], GLU_Z),
    site_diameter = 0,
    number_to_release = n_glu
)

Release_Site_23 = m.ReleaseSite(
    name = 'Release_Site_23',
    complex = m.Complex('glu'),
    shape = m.Shape.SPHERICAL,
    location = (CLUSTER_CENTERS[1][0], CLUSTER_CENTERS[1][1], GLU_Z),
    site_diameter = 0,
    number_to_release = n_glu
)

Release_Site_12 = m.ReleaseSite(
    name = 'Release_Site_12',
    complex = m.Complex('Ca2_ER'),
    region = SA,
    concentration = Ca_ER_IC
)

Release_Site_13 = m.ReleaseSite(
    name = 'Release_Site_13',
    complex = m.Complex('SERCAx', orientation = m.Orientation.UP),
    region = SA,
    density = SERCA_IC
)

Release_Site_14 = m.ReleaseSite(
    name = 'Release_Site_14',
    complex = m.Complex('RyR_C(a~0, i~0)', orientation = m.Orientation.UP),
    region = SA,
    number_to_release = 60
)

Release_Site_17 = m.ReleaseSite(
    name = 'Release_Site_17',
    complex = m.Complex('PMCA', orientation = m.Orientation.UP),
    region = Mushroom_dendrite_cylinder_receptors,
    density = PMCA_IC_HEAD
)

Release_Site_18 = m.ReleaseSite(
    name = 'Release_Site_18',
    complex = m.Complex('NCX', orientation = m.Orientation.UP),
    region = Mushroom_dendrite_cylinder_receptors,
    density = NCX_IC_HEAD
)

Release_Site_19 = m.ReleaseSite(
    name = 'Release_Site_19',
    complex = m.Complex('IP3'),
    region = Mushroom_dendrite- SA,
    concentration = 50e-9
)

Release_Site_20 = m.ReleaseSite(
    name = 'Release_Site_20',
    complex = m.Complex('mGlur', orientation = m.Orientation.UP),
    region = Mushroom_dendrite- Mushroom_dendrite_PSD ,
    density = 602.2
)

Release_Site_21 = m.ReleaseSite(
    name = 'Release_Site_21',
    complex = m.Complex('PLC', orientation = m.Orientation.UP),
    region = Mushroom_dendrite - Mushroom_dendrite_PSD ,
    density = 120.44
)

# ORIGINAL 8-subunit IP3R (single subunit conducts) — replaced by 24-state channel model, see README
#Release_Site_22 = m.ReleaseSite(
#    name = 'Release_Site_22',
#    complex = m.Complex('IP3R_000', orientation = m.Orientation.UP),
#    region = SA,
#    density = 26.389
#)
# density = 26.389 now seeds CHANNELS, not subunits -> 3x the IP3 binding sites of
# the published model. This is deliberate; see README ("Channel density").
Release_Site_22 = m.ReleaseSite(
    name = 'Release_Site_22',
    complex = m.Complex('IP3R(s1~000,s2~000,s3~000)', orientation = m.Orientation.UP),
    region = SA,
    density = 26.389
)

# ---- create instantiation object and add components ----

instantiation = m.Instantiation()
instantiation.add_geometry_object(presynapse)
instantiation.add_geometry_object(Cube)
instantiation.add_geometry_object(cylinder_side1)
instantiation.add_geometry_object(cylinder_side2)
instantiation.add_geometry_object(head_neck)
instantiation.add_geometry_object(Spine_base)
instantiation.add_geometry_object(Mushroom_dendrite)
instantiation.add_geometry_object(Mushroom_transparent)
instantiation.add_geometry_object(SA)
# Place the spatially restricted NMDARs before other surface molecules can
# occupy their cluster tiles; preserve all counts, regions, and release times.
for site in NMDAR_cluster_sites:
    instantiation.add_release_site(site)
instantiation.add_release_site(Release_Site_1)
instantiation.add_release_site(Release_Site_2)
instantiation.add_release_site(Release_Site_3)
instantiation.add_release_site(Release_Site_4)
instantiation.add_release_site(Release_Site_5)
instantiation.add_release_site(Release_Site_6)
instantiation.add_release_site(Release_Site_7)
instantiation.add_release_site(Release_Site_8)
instantiation.add_release_site(Release_Site_10)
if GLU_SITES in ('1', 'both'):
    instantiation.add_release_site(Release_Site_11)
if GLU_SITES in ('2', 'both'):
    instantiation.add_release_site(Release_Site_23)
instantiation.add_release_site(Release_Site_12)
instantiation.add_release_site(Release_Site_13)
instantiation.add_release_site(Release_Site_14)
instantiation.add_release_site(Release_Site_17)
instantiation.add_release_site(Release_Site_18)
instantiation.add_release_site(Release_Site_19)
instantiation.add_release_site(Release_Site_20)
instantiation.add_release_site(Release_Site_21)
instantiation.add_release_site(Release_Site_22)

# load seed species information from bngl file
instantiation.load_bngl_compartments_and_seed_species(os.path.join(MODEL_PATH, 'model.bngl'), None, shared.parameter_overrides)
