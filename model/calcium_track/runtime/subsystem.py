"""Calcium chemistry with ER; labeled calcium has exactly the inherited diffusion rates."""
import json
from pathlib import Path
import mcell as m
import shared
from calcium_labels import CARRIERS, labeled
from channels import add_channel_rules, channel_species
from condition import CONFIG
from parameters import k_SERCA_leak

ROOT = Path(__file__).resolve().parent
PSD = m.SurfaceClass(name='PSD', type=m.SurfacePropertyType.REACTIVE)
Surface_cube = m.SurfaceClass(name='Surface_cube', type=m.SurfacePropertyType.ABSORPTIVE,
                             affected_complex_pattern=m.AllMolecules)
Ions_mushroom = m.SurfaceClass(name='Ions_mushroom', type=m.SurfacePropertyType.TRANSPARENT,
                             affected_complex_pattern=m.AllMolecules)
Cylinder_reflective = m.SurfaceClass(name='Cylinder_reflective', type=m.SurfacePropertyType.REFLECTIVE,
                                    affected_complex_pattern=m.AllMolecules)
Cylinder_side = m.SurfaceClass(name='Cylinder_side', type=m.SurfacePropertyType.TRANSPARENT,
                               affected_complex_pattern=m.AllMolecules)
SA_membrane = m.SurfaceClass(name='SA_membrane', type=m.SurfacePropertyType.REFLECTIVE,
                             affected_complex_pattern=m.AllMolecules)
# BNGL cannot express this surface-class reaction. It releases unlabeled calcium.
rxn_SERCA_leak = m.ReactionRule(name='rxn_SERCA_leak',
    reactants=[m.Complex('Ca2_ER@IN', compartment_name='IN'), m.Complex('SA_membrane')],
    products=[m.Complex('Ca2@OUT', compartment_name='OUT')], fwd_rate=k_SERCA_leak)
subsystem = m.Subsystem()
for surface in (PSD, Surface_cube, Ions_mushroom, Cylinder_reflective, Cylinder_side, SA_membrane):
    subsystem.add_surface_class(surface)
subsystem.add_reaction_rule(rxn_SERCA_leak)
subsystem.load_bngl_molecule_types_and_reaction_rules(str(ROOT / 'model.bngl'), shared.parameter_overrides)
diffusion = json.loads((ROOT / 'diffusion.json').read_text())
if CONFIG['immobile_buffer'] == 'cytosol':
    diffusion['Bf'] = diffusion['CaBf'] = ['diffusion_constant_3d', 0]
for carrier in CARRIERS:
    for tag in 'AB':
        diffusion[labeled(carrier, tag)] = diffusion[carrier]
for name, (attribute, value) in diffusion.items():
    setattr(subsystem.find_elementary_molecule_type(name), attribute, value)
for name in channel_species():
    subsystem.find_elementary_molecule_type(name).diffusion_constant_2d = 0
NMDAR_FLUX_RULES, VGCC_FLUX_RULES = add_channel_rules(subsystem)
