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

# set additional information about species and molecule types that cannot be stored in BNGL,
# elementary molecule types are already in the subsystem after they were loaded from BNGL
def set_bngl_molecule_types_info(subsystem):
    Ca2 = subsystem.find_elementary_molecule_type('Ca2')
    assert Ca2, "Elementary molecule type 'Ca2' was not found"
    Ca2.diffusion_constant_3d = 2.2e-6

    Bm = subsystem.find_elementary_molecule_type('Bm')
    assert Bm, "Elementary molecule type 'Bm' was not found"
    Bm.diffusion_constant_3d = 2e-7

    CaBm = subsystem.find_elementary_molecule_type('CaBm')
    assert CaBm, "Elementary molecule type 'CaBm' was not found"
    CaBm.diffusion_constant_3d = 2e-7

    Bf = subsystem.find_elementary_molecule_type('Bf')
    assert Bf, "Elementary molecule type 'Bf' was not found"
    Bf.diffusion_constant_2d = 0

    CaBf = subsystem.find_elementary_molecule_type('CaBf')
    assert CaBf, "Elementary molecule type 'CaBf' was not found"
    CaBf.diffusion_constant_2d = 0

    VSCC = subsystem.find_elementary_molecule_type('VSCC')
    assert VSCC, "Elementary molecule type 'VSCC' was not found"
    VSCC.diffusion_constant_2d = 0

    PMCA = subsystem.find_elementary_molecule_type('PMCA')
    assert PMCA, "Elementary molecule type 'PMCA' was not found"
    PMCA.diffusion_constant_2d = 0

    PMCA1 = subsystem.find_elementary_molecule_type('PMCA1')
    assert PMCA1, "Elementary molecule type 'PMCA1' was not found"
    PMCA1.diffusion_constant_2d = 0

    NCX = subsystem.find_elementary_molecule_type('NCX')
    assert NCX, "Elementary molecule type 'NCX' was not found"
    NCX.diffusion_constant_2d = 0

    NCX1 = subsystem.find_elementary_molecule_type('NCX1')
    assert NCX1, "Elementary molecule type 'NCX1' was not found"
    NCX1.diffusion_constant_2d = 0

    SERCAx = subsystem.find_elementary_molecule_type('SERCAx')
    assert SERCAx, "Elementary molecule type 'SERCAx' was not found"
    SERCAx.diffusion_constant_2d = 0

    SERCAx1 = subsystem.find_elementary_molecule_type('SERCAx1')
    assert SERCAx1, "Elementary molecule type 'SERCAx1' was not found"
    SERCAx1.diffusion_constant_2d = 0

    SERCAx2 = subsystem.find_elementary_molecule_type('SERCAx2')
    assert SERCAx2, "Elementary molecule type 'SERCAx2' was not found"
    SERCAx2.diffusion_constant_2d = 0

    SERCAy = subsystem.find_elementary_molecule_type('SERCAy')
    assert SERCAy, "Elementary molecule type 'SERCAy' was not found"
    SERCAy.diffusion_constant_2d = 0

    SERCAy1 = subsystem.find_elementary_molecule_type('SERCAy1')
    assert SERCAy1, "Elementary molecule type 'SERCAy1' was not found"
    SERCAy1.diffusion_constant_2d = 0

    SERCAy2 = subsystem.find_elementary_molecule_type('SERCAy2')
    assert SERCAy2, "Elementary molecule type 'SERCAy2' was not found"
    SERCAy2.diffusion_constant_2d = 0

    VSCC_0 = subsystem.find_elementary_molecule_type('VSCC_0')
    assert VSCC_0, "Elementary molecule type 'VSCC_0' was not found"
    VSCC_0.diffusion_constant_2d = 0

    VSCC_1 = subsystem.find_elementary_molecule_type('VSCC_1')
    assert VSCC_1, "Elementary molecule type 'VSCC_1' was not found"
    VSCC_1.diffusion_constant_2d = 0

    VSCC_2 = subsystem.find_elementary_molecule_type('VSCC_2')
    assert VSCC_2, "Elementary molecule type 'VSCC_2' was not found"
    VSCC_2.diffusion_constant_2d = 0

    VSCC_3 = subsystem.find_elementary_molecule_type('VSCC_3')
    assert VSCC_3, "Elementary molecule type 'VSCC_3' was not found"
    VSCC_3.diffusion_constant_2d = 0
    VSCC_3.custom_time_step = 1e-7

    NMDAR_0 = subsystem.find_elementary_molecule_type('NMDAR_0')
    assert NMDAR_0, "Elementary molecule type 'NMDAR_0' was not found"
    NMDAR_0.diffusion_constant_2d = 0
    NMDAR_0.custom_time_step = 1e-7

    NMDAR_1 = subsystem.find_elementary_molecule_type('NMDAR_1')
    assert NMDAR_1, "Elementary molecule type 'NMDAR_1' was not found"
    NMDAR_1.diffusion_constant_2d = 0

    NMDAR_2 = subsystem.find_elementary_molecule_type('NMDAR_2')
    assert NMDAR_2, "Elementary molecule type 'NMDAR_2' was not found"
    NMDAR_2.diffusion_constant_2d = 0

    NMDAR_D = subsystem.find_elementary_molecule_type('NMDAR_D')
    assert NMDAR_D, "Elementary molecule type 'NMDAR_D' was not found"
    NMDAR_D.diffusion_constant_2d = 0

    NMDAR_0B = subsystem.find_elementary_molecule_type('NMDAR_0B')
    assert NMDAR_0B, "Elementary molecule type 'NMDAR_0B' was not found"
    NMDAR_0B.diffusion_constant_2d = 0

    NMDAR_1B = subsystem.find_elementary_molecule_type('NMDAR_1B')
    assert NMDAR_1B, "Elementary molecule type 'NMDAR_1B' was not found"
    NMDAR_1B.diffusion_constant_2d = 0

    NMDAR_2B = subsystem.find_elementary_molecule_type('NMDAR_2B')
    assert NMDAR_2B, "Elementary molecule type 'NMDAR_2B' was not found"
    NMDAR_2B.diffusion_constant_2d = 0

    NMDAR_DB = subsystem.find_elementary_molecule_type('NMDAR_DB')
    assert NMDAR_DB, "Elementary molecule type 'NMDAR_DB' was not found"
    NMDAR_DB.diffusion_constant_2d = 0

    NMDAR = subsystem.find_elementary_molecule_type('NMDAR')
    assert NMDAR, "Elementary molecule type 'NMDAR' was not found"
    NMDAR.diffusion_constant_2d = 0

    NMDAR_B = subsystem.find_elementary_molecule_type('NMDAR_B')
    assert NMDAR_B, "Elementary molecule type 'NMDAR_B' was not found"
    NMDAR_B.diffusion_constant_2d = 0

    AMPAR_0 = subsystem.find_elementary_molecule_type('AMPAR_0')
    assert AMPAR_0, "Elementary molecule type 'AMPAR_0' was not found"
    AMPAR_0.diffusion_constant_2d = 0
    AMPAR_0.custom_time_step = 1e-7

    AMPAR_1 = subsystem.find_elementary_molecule_type('AMPAR_1')
    assert AMPAR_1, "Elementary molecule type 'AMPAR_1' was not found"
    AMPAR_1.diffusion_constant_2d = 0

    AMPAR_2 = subsystem.find_elementary_molecule_type('AMPAR_2')
    assert AMPAR_2, "Elementary molecule type 'AMPAR_2' was not found"
    AMPAR_2.diffusion_constant_2d = 0

    AMPAR_3 = subsystem.find_elementary_molecule_type('AMPAR_3')
    assert AMPAR_3, "Elementary molecule type 'AMPAR_3' was not found"
    AMPAR_3.diffusion_constant_2d = 0

    AMPAR_4 = subsystem.find_elementary_molecule_type('AMPAR_4')
    assert AMPAR_4, "Elementary molecule type 'AMPAR_4' was not found"
    AMPAR_4.diffusion_constant_2d = 0

    AMPAR_5 = subsystem.find_elementary_molecule_type('AMPAR_5')
    assert AMPAR_5, "Elementary molecule type 'AMPAR_5' was not found"
    AMPAR_5.diffusion_constant_2d = 0

    AMPAR = subsystem.find_elementary_molecule_type('AMPAR')
    assert AMPAR, "Elementary molecule type 'AMPAR' was not found"
    AMPAR.diffusion_constant_2d = 0

    glu = subsystem.find_elementary_molecule_type('glu')
    assert glu, "Elementary molecule type 'glu' was not found"
    glu.diffusion_constant_3d = 2.2e-6
    glu.custom_time_step = 1e-7

    Ca2_ER = subsystem.find_elementary_molecule_type('Ca2_ER')
    assert Ca2_ER, "Elementary molecule type 'Ca2_ER' was not found"
    Ca2_ER.diffusion_constant_3d = 2.2e-6

    RyR_C = subsystem.find_elementary_molecule_type('RyR_C')
    assert RyR_C, "Elementary molecule type 'RyR_C' was not found"
    RyR_C.diffusion_constant_2d = 0

    RyR_O = subsystem.find_elementary_molecule_type('RyR_O')
    assert RyR_O, "Elementary molecule type 'RyR_O' was not found"
    RyR_O.diffusion_constant_2d = 0

    RyR_CI = subsystem.find_elementary_molecule_type('RyR_CI')
    assert RyR_CI, "Elementary molecule type 'RyR_CI' was not found"
    RyR_CI.diffusion_constant_2d = 0

    RyR_OI = subsystem.find_elementary_molecule_type('RyR_OI')
    assert RyR_OI, "Elementary molecule type 'RyR_OI' was not found"
    RyR_OI.diffusion_constant_2d = 0

    # ORIGINAL 8-subunit IP3R (single subunit conducts) — replaced by 24-state channel model, see README
    #    IP3R_000 = subsystem.find_elementary_molecule_type('IP3R_000')
    #    assert IP3R_000, "Elementary molecule type 'IP3R_000' was not found"
    #    IP3R_000.diffusion_constant_2d = 0
    #
    #    IP3R_001 = subsystem.find_elementary_molecule_type('IP3R_001')
    #    assert IP3R_001, "Elementary molecule type 'IP3R_001' was not found"
    #    IP3R_001.diffusion_constant_2d = 0
    #
    #    IP3R_010 = subsystem.find_elementary_molecule_type('IP3R_010')
    #    assert IP3R_010, "Elementary molecule type 'IP3R_010' was not found"
    #    IP3R_010.diffusion_constant_2d = 0
    #
    #    IP3R_011 = subsystem.find_elementary_molecule_type('IP3R_011')
    #    assert IP3R_011, "Elementary molecule type 'IP3R_011' was not found"
    #    IP3R_011.diffusion_constant_2d = 0
    #
    #    IP3R_100 = subsystem.find_elementary_molecule_type('IP3R_100')
    #    assert IP3R_100, "Elementary molecule type 'IP3R_100' was not found"
    #    IP3R_100.diffusion_constant_2d = 0
    #
    #    IP3R_101 = subsystem.find_elementary_molecule_type('IP3R_101')
    #    assert IP3R_101, "Elementary molecule type 'IP3R_101' was not found"
    #    IP3R_101.diffusion_constant_2d = 0
    #
    #    IP3R_110 = subsystem.find_elementary_molecule_type('IP3R_110')
    #    assert IP3R_110, "Elementary molecule type 'IP3R_110' was not found"
    #    IP3R_110.diffusion_constant_2d = 0
    #
    #    IP3R_111 = subsystem.find_elementary_molecule_type('IP3R_111')
    #    assert IP3R_111, "Elementary molecule type 'IP3R_111' was not found"
    #    IP3R_111.diffusion_constant_2d = 0

    IP3R = subsystem.find_elementary_molecule_type('IP3R')
    assert IP3R, "Elementary molecule type 'IP3R' was not found"
    IP3R.diffusion_constant_2d = 0

    IP3 = subsystem.find_elementary_molecule_type('IP3')
    assert IP3, "Elementary molecule type 'IP3' was not found"
    IP3.diffusion_constant_3d = 2.8e-6

    mGlur = subsystem.find_elementary_molecule_type('mGlur')
    assert mGlur, "Elementary molecule type 'mGlur' was not found"
    mGlur.diffusion_constant_2d = 0

    mGlur_Glu = subsystem.find_elementary_molecule_type('mGlur_Glu')
    assert mGlur_Glu, "Elementary molecule type 'mGlur_Glu' was not found"
    mGlur_Glu.diffusion_constant_2d = 0

    mGlur_Glu_Glu = subsystem.find_elementary_molecule_type('mGlur_Glu_Glu')
    assert mGlur_Glu_Glu, "Elementary molecule type 'mGlur_Glu_Glu' was not found"
    mGlur_Glu_Glu.diffusion_constant_2d = 0

    PLC = subsystem.find_elementary_molecule_type('PLC')
    assert PLC, "Elementary molecule type 'PLC' was not found"
    PLC.diffusion_constant_2d = 0

    PLC_Ca = subsystem.find_elementary_molecule_type('PLC_Ca')
    assert PLC_Ca, "Elementary molecule type 'PLC_Ca' was not found"
    PLC_Ca.diffusion_constant_2d = 0

    PLC_Ca_PIP2 = subsystem.find_elementary_molecule_type('PLC_Ca_PIP2')
    assert PLC_Ca_PIP2, "Elementary molecule type 'PLC_Ca_PIP2' was not found"
    PLC_Ca_PIP2.diffusion_constant_2d = 0

set_bngl_molecule_types_info(subsystem)
