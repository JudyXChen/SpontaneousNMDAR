# NMDAR cluster geometry

## How cluster placement is encoded

The spine and PSD meshes were exported from CellBlender, but the NMDAR
clusters are defined in each model leaf's `instantiation.py`. They are not
stored as manually positioned individual receptors in the Blender scene.

Each cluster is a MCell `SurfaceRegion` constructed from hard-coded PSD mesh
triangle indices. For example, the calibrated 100 nm two-cluster model uses:

```python
CLUSTER_WALL_INDICES = (
    (1681, 1686, 1688),
    (1687, 1691, 1694),
)
```

The selected wall indices are the actual geometry definition. A directory
name such as `2_cluster_100nm` records the intended condition, but does not
itself control the separation.

MCell places an exact number of receptors on each selected surface region. In
a two-cluster model, the current implementation uses:

```python
number_to_release = n_NMDAR // 2
```

Thus, `n_NMDAR = 20` means 20 receptors in total and 10 receptors in each
cluster. It does not mean 20 receptors per cluster. The individual receptor
positions are sampled over the selected triangles, so receptor-to-receptor
distances are not prescribed.

## Receptor count and effective density

NMDAR placement is count-based rather than density-based. The effective
surface density is derived after choosing the receptor count and mesh patch:

$$
\rho_{\mathrm{NMDAR}} = \frac{N_{\mathrm{NMDAR}}}{A_{\mathrm{patch}}}.
$$

The matched one-cluster controls use approximately the same density by
increasing the patch area with receptor count:

| Condition | PSD wall indices | Patch area (nm²) | NMDAR | Effective density (µm⁻²) |
| --- | --- | ---: | ---: | ---: |
| One cluster, 10 NMDAR | 1681, 1686, 1688 | 2,546.8 | 10 | about 3,926 |
| One cluster, 20 NMDAR | 1681, 1686, 1688, 1690, 1696, 1719 | 5,044.5 | 20 | about 3,965 |

The 20-receptor patch therefore has only about 0.9% higher density than the
10-receptor patch; its density is not doubled.

For the 100 nm two-cluster model, the two patches contain 10 receptors each:

| Cluster | PSD wall indices | Patch area (nm²) | NMDAR | Effective density (µm⁻²) |
| --- | --- | ---: | ---: | ---: |
| 1 | 1681, 1686, 1688 | 2,546.8 | 10 | about 3,926 |
| 2 | 1687, 1691, 1694 | 2,313.2 | 10 | about 4,323 |

The second cluster is therefore about 10% denser because the discrete mesh
does not provide an exactly equal-area patch at the requested separation.

The two-triangle patch sometimes recalled in this context belongs to the
spontaneous 5-NMDAR condition. In that comparison, 5 NMDARs use two triangles
(`1688, 1696`; 1,328.5 nm²), whereas 10 NMDARs use three triangles
(`1681, 1686, 1688`; 2,546.8 nm²). Their effective densities differ by about
4.2%.

## Definition of the cluster center

The center is not the midpoint of a shared vertex or shared triangle edge.
For every selected triangle, the code collects its three vertex coordinates
and averages all coordinate entries:

```python
vertices = [
    mesh_vertices[vertex_index]
    for wall_index in cluster_wall_indices
    for vertex_index in mesh_walls[wall_index]
]
center = tuple(
    sum(point[axis] for point in vertices) / len(vertices)
    for axis in range(3)
)
```

Consequently:

- the center of a one-triangle cluster is that triangle's centroid;
- for multiple triangles, it is the equally weighted average of their
  triangle centroids;
- a vertex shared by multiple selected triangles appears once for each
  triangle containing it;
- this is not an area-weighted centroid of the union of the triangles when
  their areas differ.

The calculated center is used to measure cluster separation and to set the
horizontal coordinates of the corresponding glutamate release site. The
glutamate release height is fixed separately at `z = 1.655 µm`. NMDARs are
distributed across the entire surface region rather than placed at its
calculated center.

## How target separation is converted into a model

The distance is an intended center-to-center separation. It is converted to a
mesh layout before the simulation:

1. Anchor the first cluster to a selected PSD patch.
2. Search or inspect candidate PSD triangles for the second patch.
3. Calculate the candidate patch center, center-to-center separation, and
   patch area.
4. Choose a patch that most closely matches the target separation while
   keeping patch area comparable when possible.
5. Hard-code the selected wall indices in `instantiation.py` and record the
   realized separation.

There is no run-time back-calculation from a requested distance or density.
The discretized mesh limits the available patch locations and areas, so the
realized values differ slightly from the nominal conditions:

| Nominal separation | Receptors per cluster | Realized separation |
| ---: | ---: | ---: |
| 50 nm | 10 + 10 | 49.469 nm |
| 100 nm | 10 + 10 | 99.422 nm |
| 200 nm | 10 + 10 | 200.083 nm |

The separate 3 + 3-NMDAR series uses one PSD triangle per cluster and realizes
46.672, 99.978, and 200.400 nm for its nominal 50, 100, and 200 nm conditions,
respectively. A small 20 nm condition realizes 21.008 nm with 3 + 3 NMDARs.

Because `instantiation.py` is derived from a CellBlender export, a subsequent
CellBlender re-export may overwrite these code-level cluster definitions. The
wall-index definitions and their validation should therefore be preserved or
reapplied after an export.
