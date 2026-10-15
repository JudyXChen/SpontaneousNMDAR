# WARNING: This is an automatically generated file and will be overwritten
#          by CellBlender on the next model export.

import os
import shared
import mcell as m

from parameters import *
from subsystem import *
from geometry import *
MODEL_PATH = os.path.dirname(os.path.abspath(__file__))


# ---- observables ----

viz_output = m.VizOutput(
    mode = m.VizMode.ASCII,
    output_files_prefix = './viz_data/seed_' + str(SEED).zfill(5) + '/Scene',
    species_list = [m.Species('Ca2')] + [m.Species(p + '_0') for p in NMDAR_FLUX_RULES],
    every_n_timesteps = 1000  # 0.1 ms at the default timestep
)

# ---- declaration of rxn rules defined in BNGL and used in counts ----

cterm_count_Ca2_species_Mushroom_transparent = m.CountTerm(
    species_pattern = m.Complex('Ca2'),
    region = Mushroom_transparent
)

count_Ca2_Mushroom_transparent = m.Count(
    name = 'Ca2_Mushroom_transparent',
    expression = cterm_count_Ca2_species_Mushroom_transparent,
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/Ca2_Mushroom_transparent.dat',
    every_n_timesteps = 1e-6/1e-8
)

cterm_count_Ca2_species = m.CountTerm(
    species_pattern = m.Complex('Ca2')
)

count_Ca2_World = m.Count(
    name = 'Ca2_World',
    expression = cterm_count_Ca2_species,
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/Ca2_World.dat',
    every_n_timesteps = 1e-6/1e-8
)

cterm_count_CaBm_species_Mushroom_transparent = m.CountTerm(
    species_pattern = m.Complex('CaBm'),
    region = Mushroom_transparent
)

count_CaBm_Mushroom_transparent = m.Count(
    name = 'CaBm_Mushroom_transparent',
    expression = cterm_count_CaBm_species_Mushroom_transparent,
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/CaBm_Mushroom_transparent.dat',
    every_n_timesteps = 1e-6/1e-8
)

cterm_count_CaBf_species = m.CountTerm(
    species_pattern = m.Complex('CaBf')
)

count_CaBf_World = m.Count(
    name = 'CaBf_World',
    expression = cterm_count_CaBf_species,
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/CaBf_World.dat',
    every_n_timesteps = 1e-6/1e-8
)

cterm_count_Bf_species_Mushroom_transparent = m.CountTerm(
    species_pattern = m.Complex('Bf'),
    region = Mushroom_transparent
)

count_Bf_Mushroom_transparent = m.Count(
    name = 'Bf_Mushroom_transparent',
    expression = cterm_count_Bf_species_Mushroom_transparent,
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/Bf_Mushroom_transparent.dat',
    every_n_timesteps = 1e-6/1e-8
)

cterm_count_Bm_species_Mushroom_transparent = m.CountTerm(
    species_pattern = m.Complex('Bm'),
    region = Mushroom_transparent
)

count_Bm_Mushroom_transparent = m.Count(
    name = 'Bm_Mushroom_transparent',
    expression = cterm_count_Bm_species_Mushroom_transparent,
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/Bm_Mushroom_transparent.dat',
    every_n_timesteps = 1e-6/1e-8
)

from channel_counts import make_channel_counts
channel_counts = make_channel_counts(NMDAR_FLUX_RULES, VGCC_FLUX_RULES, SEED)
count_VSCC_World = channel_counts['VSCC_World']
count_NMDAR_World = channel_counts['NMDAR_World']
count_NMDAR_Ca_flux_World = channel_counts['NMDAR_Ca_flux_World']

cterm_count_PMCA1_species = m.CountTerm(
    species_pattern = m.Complex('PMCA1')
)

count_PMCA1_World = m.Count(
    name = 'PMCA1_World',
    expression = cterm_count_PMCA1_species,
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/PMCA1_World.dat',
    every_n_timesteps = 1e-6/1e-8
)

cterm_count_NCX1_species = m.CountTerm(
    species_pattern = m.Complex('NCX1')
)

count_NCX1_World = m.Count(
    name = 'NCX1_World',
    expression = cterm_count_NCX1_species,
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/NCX1_World.dat',
    every_n_timesteps = 1e-6/1e-8
)

cterm_count_Ca2_species_Spine_base = m.CountTerm(
    species_pattern = m.Complex('Ca2'),
    region = Spine_base
)

count_Ca2_Spine_base = m.Count(
    name = 'Ca2_Spine_base',
    expression = cterm_count_Ca2_species_Spine_base,
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/Ca2_Spine_base.dat',
    every_n_timesteps = 1e-6/1e-8
)

cterm_count_Ca2_species_cylinder_side1 = m.CountTerm(
    species_pattern = m.Complex('Ca2'),
    region = cylinder_side1
)

count_Ca2_cylinder_side1 = m.Count(
    name = 'Ca2_cylinder_side1',
    expression = cterm_count_Ca2_species_cylinder_side1,
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/Ca2_cylinder_side1.dat',
    every_n_timesteps = 1e-6/1e-8
)

cterm_count_Ca2_species_cylinder_side2 = m.CountTerm(
    species_pattern = m.Complex('Ca2'),
    region = cylinder_side2
)

count_Ca2_cylinder_side2 = m.Count(
    name = 'Ca2_cylinder_side2',
    expression = cterm_count_Ca2_species_cylinder_side2,
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/Ca2_cylinder_side2.dat',
    every_n_timesteps = 1e-6/1e-8
)

cterm_count_Ca2_species_head_neck = m.CountTerm(
    species_pattern = m.Complex('Ca2'),
    region = head_neck
)

count_Ca2_head_neck = m.Count(
    name = 'Ca2_head_neck',
    expression = cterm_count_Ca2_species_head_neck,
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/Ca2_head_neck.dat',
    every_n_timesteps = 1e-6/1e-8
)

cterm_count_Ca2_species_SA = m.CountTerm(
    species_pattern = m.Complex('Ca2'),
    region = SA
)

count_Ca2_SA = m.Count(
    name = 'Ca2_SA',
    expression = cterm_count_Ca2_species_SA,
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/Ca2_SA.dat',
    every_n_timesteps = 1e-6/1e-8
)

cterm_count_Ca2_ER_species_SA = m.CountTerm(
    species_pattern = m.Complex('Ca2_ER'),
    region = SA
)

count_Ca2_ER_SA = m.Count(
    name = 'Ca2_ER_SA',
    expression = cterm_count_Ca2_ER_species_SA,
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/Ca2_ER_SA.dat',
    every_n_timesteps = 1e-6/1e-8
)

cterm_count_RyR_O_species = m.CountTerm(
    species_pattern = m.Complex('RyR_O')
)

count_RyR_O_World = m.Count(
    name = 'RyR_O_World',
    expression = cterm_count_RyR_O_species,
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/RyR_O_World.dat',
    every_n_timesteps = 1e-6/1e-8
)

cterm_count_Ca2_species_Mushroom_dendrite = m.CountTerm(
    species_pattern = m.Complex('Ca2'),
    region = Mushroom_dendrite
)

count_Ca2_Mushroom_dendrite = m.Count(
    name = 'Ca2_Mushroom_dendrite',
    expression = cterm_count_Ca2_species_Mushroom_dendrite,
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/Ca2_Mushroom_dendrite.dat',
    every_n_timesteps = 1e-6/1e-8
)

# ORIGINAL 8-subunit IP3R (single subunit conducts) — replaced by 24-state channel model, see README
#cterm_count_IP3R_110_species = m.CountTerm(
#    species_pattern = m.Complex('IP3R_110')
#)
#
#count_IP3R_110_World = m.Count(
#    name = 'IP3R_110_World',
#    expression = cterm_count_IP3R_110_species,
#    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/IP3R_110_World.dat',
#    every_n_timesteps = 1e-6/1e-8
#)

# ---- IP3R: 24-state channel model ----

# Open channels: all three subunits in 110. This is the DYK (x110)^3 gate, and the
# quantity that actually conducts.
count_IP3R_open_World = m.Count(
    name = 'IP3R_open_World',
    expression = m.CountTerm(
        species_pattern = m.Complex('IP3R(s1~110,s2~110,s3~110)')
    ),
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/IP3R_open_World.dat',
    every_n_timesteps = 1e-6/1e-8
)

# Total channels: must be constant at density * A_SA for the whole run.
count_IP3R_total_World = m.Count(
    name = 'IP3R_total_World',
    expression = m.CountTerm(
        species_pattern = m.Complex('IP3R')
    ),
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/IP3R_total_World.dat',
    every_n_timesteps = 1e-6/1e-8
)

# Per-subunit-in-110 populations. IP3R_s1_110 is the direct analogue of the old
# IP3R_110_World observable (one subunit in 110, regardless of the other two), so it
# makes the old-vs-new comparison possible. s2/s3 are recorded to check that the
# three subunits are statistically equivalent.
for _s in ('s1', 's2', 's3'):
    observables_count = m.Count(
        name = 'IP3R_' + _s + '_110_World',
        expression = m.CountTerm(
            species_pattern = m.Complex('IP3R(' + _s + '~110)')
        ),
        file_name = './react_data/seed_' + str(SEED).zfill(5) + '/IP3R_' + _s + '_110_World.dat',
        every_n_timesteps = 1e-6/1e-8
    )
    globals()['count_IP3R_' + _s + '_110_World'] = observables_count

# ---- directly measured Ca fluxes (reaction-firing counts) ----
#
# m.CountTerm(reaction_rule=...) yields a cumulative count of firings since t=0, so
# these traces are monotonically increasing; a rate is np.gradient(count, t).
# Replaces plot_flux.py's reconstruction of these fluxes from rate laws, and (for
# SERCA) its ER mass-balance residual.

# IP3R Ca flux: ER -> cytosol, one firing = one Ca ion.
count_IP3R_flux_World = m.Count(
    name = 'IP3R_flux_World',
    expression = m.CountTerm(
        reaction_rule = subsystem.find_reaction_rule('rxn_IP3R_flux')
    ),
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/IP3R_flux_World.dat',
    every_n_timesteps = 1e-6/1e-8
)

# RyR Ca flux: the reversible rule was split into named f/r rules, so gross flux in
# each direction is separable and the net is f - r.
count_RyR_flux_f_World = m.Count(
    name = 'RyR_flux_f_World',
    expression = m.CountTerm(
        reaction_rule = subsystem.find_reaction_rule('rxn_RyR_flux_f')
    ),
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/RyR_flux_f_World.dat',
    every_n_timesteps = 1e-6/1e-8
)

count_RyR_flux_r_World = m.Count(
    name = 'RyR_flux_r_World',
    expression = m.CountTerm(
        reaction_rule = subsystem.find_reaction_rule('rxn_RyR_flux_r')
    ),
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/RyR_flux_r_World.dat',
    every_n_timesteps = 1e-6/1e-8
)

count_RyR_flux_net_World = m.Count(
    name = 'RyR_flux_net_World',
    expression = (m.CountTerm(reaction_rule = subsystem.find_reaction_rule('rxn_RyR_flux_f'))
                - m.CountTerm(reaction_rule = subsystem.find_reaction_rule('rxn_RyR_flux_r'))),
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/RyR_flux_net_World.dat',
    every_n_timesteps = 1e-6/1e-8
)

# SERCA leak: ER -> cytosol. Declared in Python (subsystem.py), object already in
# scope via 'from subsystem import *', so no lookup is needed.
count_SERCA_leak_World = m.Count(
    name = 'SERCA_leak_World',
    expression = m.CountTerm(reaction_rule = rxn_SERCA_leak),
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/SERCA_leak_World.dat',
    every_n_timesteps = 1e-6/1e-8
)

# SERCA uptake: the two ER-side release steps deposit Ca into the lumen; their
# reverses remove it. Net uptake = (y2y1f + y1y0f) - (y2y1r + y1y0r).
# The four components are also emitted separately so gross uptake and gross
# backflux are recoverable, not just the net.
for _rxn in ('rxn_SERCA_y2y1f', 'rxn_SERCA_y2y1r',
             'rxn_SERCA_y1y0f', 'rxn_SERCA_y1y0r',
             'rxn_SERCA_x2y2f', 'rxn_SERCA_x2y2r'):
    _c = m.Count(
        name = _rxn + '_World',
        expression = m.CountTerm(reaction_rule = subsystem.find_reaction_rule(_rxn)),
        file_name = './react_data/seed_' + str(SEED).zfill(5) + '/' + _rxn + '_World.dat',
        every_n_timesteps = 1e-6/1e-8
    )
    globals()['count_' + _rxn + '_World'] = _c

count_SERCA_flux_World = m.Count(
    name = 'SERCA_flux_World',
    expression = (m.CountTerm(reaction_rule = subsystem.find_reaction_rule('rxn_SERCA_y2y1f'))
                + m.CountTerm(reaction_rule = subsystem.find_reaction_rule('rxn_SERCA_y1y0f'))
                - m.CountTerm(reaction_rule = subsystem.find_reaction_rule('rxn_SERCA_y2y1r'))
                - m.CountTerm(reaction_rule = subsystem.find_reaction_rule('rxn_SERCA_y1y0r'))),
    file_name = './react_data/seed_' + str(SEED).zfill(5) + '/SERCA_flux_World.dat',
    every_n_timesteps = 1e-6/1e-8
)

# ---- create observables object and add components ----

observables = m.Observables()
observables.add_viz_output(viz_output)
observables.add_count(count_Ca2_Mushroom_transparent)
observables.add_count(count_Ca2_World)
observables.add_count(count_CaBm_Mushroom_transparent)
observables.add_count(count_CaBf_World)
observables.add_count(count_Bf_Mushroom_transparent)
observables.add_count(count_Bm_Mushroom_transparent)
observables.add_count(count_VSCC_World)
observables.add_count(count_NMDAR_World)
observables.add_count(count_NMDAR_Ca_flux_World)
observables.add_count(count_PMCA1_World)
observables.add_count(count_NCX1_World)
observables.add_count(count_Ca2_Spine_base)
observables.add_count(count_Ca2_cylinder_side1)
observables.add_count(count_Ca2_cylinder_side2)
observables.add_count(count_Ca2_head_neck)
observables.add_count(count_Ca2_SA)
observables.add_count(count_Ca2_ER_SA)
observables.add_count(count_RyR_O_World)
observables.add_count(count_Ca2_Mushroom_dendrite)
# ORIGINAL 8-subunit IP3R (single subunit conducts) — replaced by 24-state channel model, see README
#observables.add_count(count_IP3R_110_World)
observables.add_count(count_IP3R_open_World)
observables.add_count(count_IP3R_total_World)
observables.add_count(count_IP3R_s1_110_World)
observables.add_count(count_IP3R_s2_110_World)
observables.add_count(count_IP3R_s3_110_World)
observables.add_count(count_IP3R_flux_World)
observables.add_count(count_RyR_flux_f_World)
observables.add_count(count_RyR_flux_r_World)
observables.add_count(count_RyR_flux_net_World)
observables.add_count(count_SERCA_leak_World)
observables.add_count(count_rxn_SERCA_y2y1f_World)
observables.add_count(count_rxn_SERCA_y2y1r_World)
observables.add_count(count_rxn_SERCA_y1y0f_World)
observables.add_count(count_rxn_SERCA_y1y0r_World)
observables.add_count(count_rxn_SERCA_x2y2f_World)
observables.add_count(count_rxn_SERCA_x2y2r_World)
observables.add_count(count_SERCA_flux_World)

# Population-specific counts, receptor totals, and aggregate channel fluxes.
for name, count in channel_counts.items():
    if name not in ('VSCC_World', 'NMDAR_World', 'NMDAR_Ca_flux_World'):
        observables.add_count(count)
