"""Free and bound source calcium, entry/removal counts and particle positions."""
from functools import reduce
import operator
import mcell as m
from parameters import SEED
from subsystem import subsystem, NMDAR_FLUX_RULES, VGCC_FLUX_RULES
from geometry import Mushroom_transparent, Mushroom_dendrite
from calcium_labels import CARRIERS, labeled
from channel_counts import make_channel_counts
from channels import NMDAR_SUFFIXES
from condition import CONFIG

observables = m.Observables()
all_carriers = list(CARRIERS) + [labeled(c, tag) for tag in 'AB' for c in CARRIERS]
viz_species = all_carriers + [p + s for p in NMDAR_FLUX_RULES for s in NMDAR_SUFFIXES]
observables.add_viz_output(m.VizOutput(
    mode=m.VizMode.ASCII, output_files_prefix=f'./viz_data/seed_{SEED:05d}/Scene',
    species_list=[m.Species(s) for s in viz_species],
    every_n_timesteps=CONFIG['viz_every_n_timesteps']))


def summed(expressions):
    return reduce(operator.add, expressions)


def species(name, region=None):
    kwargs = dict(species_pattern=m.Complex(name))
    if region is not None:
        kwargs['region'] = region
    return m.CountTerm(**kwargs)


def count(name, expression):
    observables.add_count(m.Count(name=name, expression=expression,
        file_name=f'./react_data/seed_{SEED:05d}/{name}.dat', every_n_timesteps=100))


for name in all_carriers + ['glu', 'Bm', 'Bf', 'PMCA', 'NCX', 'PLC']:
    count(name + '_World', species(name))
for tag, free in [('background', 'Ca2'), ('A', 'Ca_A'), ('B', 'Ca_B')]:
    carriers = list(CARRIERS) if tag == 'background' else [labeled(c, tag) for c in CARRIERS]
    count('Ca_' + tag + '_inventory_World', summed([species(c) for c in carriers]))
    for region in (Mushroom_transparent, Mushroom_dendrite):
        count(free + '_' + region.name, species(free, region))
    for sink in ('Ca2Null', 'rxn_PMCA', 'rxn_NCX'):
        name = sink if tag == 'background' else sink + '_' + tag
        count(name + '_World', m.CountTerm(reaction_rule=subsystem.find_reaction_rule(name)))
count('Ca_free_total_World', summed([species(c) for c in ('Ca2', 'Ca_A', 'Ca_B')]))
for name, expression in [('Bm_total', [species(s) for s in ('Bm','CaBm','CaBm_A','CaBm_B')]),
                         ('Bf_total', [species(s) for s in ('Bf','CaBf','CaBf_A','CaBf_B')]),
                         ('PMCA_total', [species(s) for s in ('PMCA','PMCA1','PMCA1_A','PMCA1_B')]),
                         ('NCX_total', [species(s) for s in ('NCX','NCX1','NCX1_A','NCX1_B')]),
                         ('PLC_total', [species(s) for s in ('PLC','PLC_Ca','PLC_Ca_PIP2','PLC_Ca_A','PLC_Ca_PIP2_A','PLC_Ca_B','PLC_Ca_PIP2_B')])]:
    count(name + '_World', summed(expression))
for c in make_channel_counts(NMDAR_FLUX_RULES, VGCC_FLUX_RULES, SEED).values():
    observables.add_count(c)
