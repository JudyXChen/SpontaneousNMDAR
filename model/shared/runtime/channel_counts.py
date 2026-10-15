"""Separate and combined channel counts, sampled on exactly the same time grid."""
from functools import reduce
import operator
import mcell as m
from channels import NMDAR_SUFFIXES, VGCC_SUFFIXES


def make_channel_counts(nmda_flux, vgcc_flux, seed):
    counts = {}

    def summed(expressions):
        return reduce(operator.add, expressions)

    def species(name):
        return m.CountTerm(species_pattern=m.Complex(name))

    def count(name, expression):
        counts[name] = m.Count(name=name, expression=expression,
                               file_name=f'./react_data/seed_{seed:05d}/{name}.dat',
                               every_n_timesteps=100)

    for fluxes, suffixes, prefix in [(nmda_flux, NMDAR_SUFFIXES, 'NMDAR'),
                                    (vgcc_flux, VGCC_SUFFIXES, 'VGCC')]:
        for population, reaction in fluxes.items():
            for suffix in suffixes:
                count(population + suffix + '_World', species(population + suffix))
            count(population + '_total_World', summed([species(population + s) for s in suffixes]))
            count(population + '_Ca_flux_World', m.CountTerm(reaction_rule=reaction))
        count(prefix + '_World', summed([species(p) for p in fluxes]))
        count(prefix + '_total_World', summed([species(p + s) for p in fluxes for s in suffixes]))
        count(prefix + '_Ca_flux_World', summed([m.CountTerm(reaction_rule=r) for r in fluxes.values()]))
    # Existing plotters use the historical VSCC spelling.
    count('VSCC_World', summed([species(p) for p in vgcc_flux]))
    count('VSCC_Ca_flux_World', summed([m.CountTerm(reaction_rule=r) for r in vgcc_flux.values()]))
    return counts
