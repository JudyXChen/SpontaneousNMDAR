"""Conserved NMDAR/VGCC populations with separate imposed voltage drives.

Glutamate is shared; NMDAR populations emit their corresponding calcium label.
Only channel names and the voltage-dependent rates differ between populations.
"""
import json
from pathlib import Path
import mcell as m
from condition import CONFIG, PROTOCOL

ROOT = Path(__file__).resolve().parent
REST = json.loads((ROOT / 'resting_rates.json').read_text())
EVOKED = json.loads((ROOT / 'evoked_rates.json').read_text())
KINETICS = json.loads((ROOT / 'nmda_kinetics.json').read_text())
NMDAR_SUFFIXES = [name[len('NMDAR'):] for name in KINETICS['states']]
VGCC_SUFFIXES = ['', '_0', '_1', '_2', '_3']


def channel_species():
    return ([p + s for p in PROTOCOL['nmda'] for s in NMDAR_SUFFIXES]
            + [p + s for p in PROTOCOL['vgcc'] for s in VGCC_SUFFIXES])


def shifted_rates(table_name, resting_rate):
    """Hold at the resting rate until the prescribed waveform begins."""
    start = CONFIG['waveform_start_s']
    result = [[0.0, resting_rate]]
    for time, value in EVOKED[table_name]:
        point = [time + start, value]
        if point[0] == result[-1][0]:
            result[-1] = point
        else:
            result.append(point)
    return result


def add_channel_rules(subsystem):
    def rule(name, reactant, product, mode, resting_key, table, calcium=False):
        products = [m.Complex(product)]
        if calcium:
            source = CONFIG['calcium_sources'].get(product, 'Ca2')
            products.append(m.Complex(source + '@IN', compartment_name='IN'))
        kwargs = dict(name=name, reactants=[m.Complex(reactant)], products=products)
        if mode == 'rest':
            kwargs['fwd_rate'] = REST[resting_key]
        elif mode == 'evoked':
            kwargs['variable_rate'] = shifted_rates(table, REST[resting_key])
        else:
            raise ValueError('Unknown channel drive: ' + mode)
        reaction = m.ReactionRule(**kwargs)
        subsystem.add_reaction_rule(reaction)
        return reaction

    nmda_flux = {}
    for population, mode in PROTOCOL['nmda'].items():
        rule('rxn_' + population + '_block', population, population + '_B', mode,
             'k_NMDAR_block', 'k_B_rate')
        rule('rxn_' + population + '_unblock', population + '_B', population, mode,
             'k_NMDAR_unblock', 'k_U_rate')
        nmda_flux[population] = rule('rxn_' + population, population, population, mode,
                                     'k_NMDAR_Ca', 'nmdar_rate_scaled', calcium=True)
    vgcc_flux = {}
    for population, mode in PROTOCOL['vgcc'].items():
        states = [population + '_' + str(i) for i in range(4)] + [population]
        for i in range(1, 5):
            for direction, source, target in [('forward', states[i-1], states[i]),
                                              ('reverse', states[i], states[i-1])]:
                rule(f'rxn_{population}_{i}_{direction}', source, target, mode,
                     f'k_VSCC_{i}_{direction}', f'vscc_{i}_{direction}_scaled')
        vgcc_flux[population] = rule('rxn_' + population, population, population, mode,
                                     'k_VSCC_Ca', 'vscc_rate_scaled', calcium=True)
    return nmda_flux, vgcc_flux
