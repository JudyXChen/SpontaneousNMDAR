# WARNING: This is an automatically generated file and will be overwritten
#          by CellBlender on the next model export.

import os
import shared
import mcell as m

from parameters import *

# ---- subsystem ----

MODEL_PATH = os.path.dirname(os.path.abspath(__file__))

PSD = m.SurfaceClass(
    name = 'PSD',
    type = m.SurfacePropertyType.REACTIVE
)

Surface_base = m.SurfaceClass(
    name = 'Surface_base',
    type = m.SurfacePropertyType.TRANSPARENT,
    affected_complex_pattern = m.Complex('Ca2', orientation = m.Orientation.DOWN)
)

Surface_cube = m.SurfaceClass(
    name = 'Surface_cube',
    type = m.SurfacePropertyType.ABSORPTIVE,
    affected_complex_pattern = m.AllMolecules
)

Ions_mushroom = m.SurfaceClass(
    name = 'Ions_mushroom',
    type = m.SurfacePropertyType.TRANSPARENT,
    affected_complex_pattern = m.AllMolecules
)

Cylinder_reflective = m.SurfaceClass(
    name = 'Cylinder_reflective',
    type = m.SurfacePropertyType.REFLECTIVE,
    affected_complex_pattern = m.AllMolecules
)

Cylinder_side = m.SurfaceClass(
    name = 'Cylinder_side',
    type = m.SurfacePropertyType.TRANSPARENT,
    affected_complex_pattern = m.AllMolecules
)

SA_membrane = m.SurfaceClass(
    name = 'SA_membrane',
    type = m.SurfacePropertyType.REFLECTIVE,
    affected_complex_pattern = m.AllMolecules
)

rxn_SERCA_leak = m.ReactionRule(
    name = 'rxn_SERCA_leak',
    reactants = [ m.Complex('Ca2_ER@IN', compartment_name = 'IN'), m.Complex('SA_membrane') ],
    products = [ m.Complex('Ca2@OUT', compartment_name = 'OUT') ],
    fwd_rate = k_SERCA_leak
)

# ---- create subsystem object and add components ----

subsystem = m.Subsystem()
subsystem.add_surface_class(PSD)
subsystem.add_surface_class(Surface_base)
subsystem.add_surface_class(Surface_cube)
subsystem.add_surface_class(Ions_mushroom)
subsystem.add_surface_class(Cylinder_reflective)
subsystem.add_surface_class(Cylinder_side)
subsystem.add_surface_class(SA_membrane)
subsystem.add_reaction_rule(rxn_SERCA_leak)

# load subsystem information from bngl file
subsystem.load_bngl_molecule_types_and_reaction_rules(os.path.join(MODEL_PATH, 'model.bngl'), shared.parameter_overrides)


from channels import add_channel_rules, channel_species

BACKGROUND_DIFFUSION = {'Ca2': ['diffusion_constant_3d', 2.2e-06], 'Bm': ['diffusion_constant_3d', 2e-07], 'CaBm': ['diffusion_constant_3d', 2e-07], 'Bf': ['diffusion_constant_2d', 0], 'CaBf': ['diffusion_constant_2d', 0], 'PMCA': ['diffusion_constant_2d', 0], 'PMCA1': ['diffusion_constant_2d', 0], 'NCX': ['diffusion_constant_2d', 0], 'NCX1': ['diffusion_constant_2d', 0], 'SERCAx': ['diffusion_constant_2d', 0], 'SERCAx1': ['diffusion_constant_2d', 0], 'SERCAx2': ['diffusion_constant_2d', 0], 'SERCAy': ['diffusion_constant_2d', 0], 'SERCAy1': ['diffusion_constant_2d', 0], 'SERCAy2': ['diffusion_constant_2d', 0], 'AMPAR_0': ['diffusion_constant_2d', 0], 'AMPAR_1': ['diffusion_constant_2d', 0], 'AMPAR_2': ['diffusion_constant_2d', 0], 'AMPAR_3': ['diffusion_constant_2d', 0], 'AMPAR_4': ['diffusion_constant_2d', 0], 'AMPAR_5': ['diffusion_constant_2d', 0], 'AMPAR': ['diffusion_constant_2d', 0], 'glu': ['diffusion_constant_3d', 2.2e-06], 'Ca2_ER': ['diffusion_constant_3d', 2.2e-06], 'RyR_C': ['diffusion_constant_2d', 0], 'RyR_O': ['diffusion_constant_2d', 0], 'RyR_CI': ['diffusion_constant_2d', 0], 'RyR_OI': ['diffusion_constant_2d', 0], 'IP3R': ['diffusion_constant_2d', 0], 'IP3': ['diffusion_constant_3d', 2.8e-06], 'mGlur': ['diffusion_constant_2d', 0], 'mGlur_Glu': ['diffusion_constant_2d', 0], 'mGlur_Glu_Glu': ['diffusion_constant_2d', 0], 'PLC': ['diffusion_constant_2d', 0], 'PLC_Ca': ['diffusion_constant_2d', 0], 'PLC_Ca_PIP2': ['diffusion_constant_2d', 0]}
for name, (attribute, value) in BACKGROUND_DIFFUSION.items():
    molecule = subsystem.find_elementary_molecule_type(name)
    if molecule is None:
        raise ValueError('Missing background molecule: ' + name)
    setattr(molecule, attribute, value)
for name in channel_species():
    subsystem.find_elementary_molecule_type(name).diffusion_constant_2d = 0
NMDAR_FLUX_RULES, VGCC_FLUX_RULES = add_channel_rules(subsystem)
