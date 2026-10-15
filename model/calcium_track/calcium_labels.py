"""Passive source labels: free Ca plus every occupied non-ER calcium state.

Each listed carrier holds ONE calcium in the inherited chemistry. Unoccupied
buffers/pumps/PLC stay shared, so labels compete for the original capacity.
"""
import re

CARRIERS = ('Ca2', 'CaBm', 'CaBf', 'PMCA1', 'NCX1', 'PLC_Ca', 'PLC_Ca_PIP2')


def labeled(name, tag):
    return 'Ca_' + tag if name == 'Ca2' else name + '_' + tag


def relabel(text, tag):
    return re.sub(r'\b(' + '|'.join(CARRIERS) + r')\b',
                  lambda match: labeled(match[0], tag), text)


def expand_background(text):
    """Duplicate calcium-involving rules, except unlabelled external sources.

Pump leaks inject new background Ca2; duplicating them would triple entry.
Names of all directed rules are explicit for validation and sink counters.
"""
    output, section, index = [], None, 0
    for line in text.splitlines():
        s = line.strip()
        if s.startswith('begin '):
            section = s[6:]
        if section == 'molecule types' and s in CARRIERS:
            output.extend('    ' + labeled(s, tag) for tag in 'AB')
        if section == 'reaction rules' and s and not s.startswith(('#', 'begin', 'end')):
            index += 1
            if ':' in s:
                name, s = s.split(':', 1)
            else:
                name = f'background_{index}'
            name, s = name.strip(), s.strip()
            if '<->' in s:
                left, right = s.split('<->')
                match = re.fullmatch(r'\s*(.*?)\s+(\S+),\s*(\S+)\s*', right)
                product, forward, reverse = match.groups()
                directed = [(name + '_f', left.strip(), product, forward),
                            (name + '_r', product, left.strip(), reverse)]
            else:
                left, right = s.split('->')
                product, rate = right.strip().rsplit(None, 1)
                directed = [(name, left.strip(), product, rate)]
            # Calcium must already be present to propagate a source label.
            # Inspect the whole reversible reaction so dissociation is included.
            has_carrier = any(re.search(r'\b' + c + r'\b', s) for c in CARRIERS)
            is_source = name in ('rxn_PMCA_leak', 'rxn_NCX_leak')
            for rule_name, left, right, rate in directed:
                output.append(f'    {rule_name}: {left} -> {right} {rate}')
                if has_carrier and not is_source:
                    for tag in 'AB':
                        output.append(f'    {rule_name}_{tag}: {relabel(left, tag)} -> {relabel(right, tag)} {rate}')
            continue
        output.append(line)
    return '\n'.join(output) + '\n'
