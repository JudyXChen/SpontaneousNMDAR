"""Passive source labels: free Ca plus every occupied non-ER calcium state.

Each listed carrier holds ONE calcium in the inherited chemistry. Unoccupied
buffers/pumps/PLC stay shared, so labels compete for the original capacity.

Labels stop at the ER. Labeled calcium binds SERCA, RyR and IP3R and passes through
open RyRs exactly as Ca2 does, but those rules end the label: anything that later
returns to the cytosol from the ER or its channels is unlabeled Ca2.
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
                # Release toward the cytosol from the ER side stays unlabeled.
                from_er = re.search(r'\bCa2@OUT\b', right)
                if has_carrier and not is_source and not from_er:
                    for tag in 'AB':
                        output.append(f'    {rule_name}_{tag}: {relabel(left, tag)} -> {relabel(right, tag)} {rate}')
            continue
        output.append(line)
    return '\n'.join(output) + '\n'


def er_handoffs(bngl_text, tag):
    """Names of the rules that end a label at the ER, grouped by partner."""
    groups = {'SERCA': [], 'RyR': [], 'IP3R': []}
    pattern = re.compile(r'^\s*(\S+_' + tag + r'):\s*(.*?)\s*->')
    for line in bngl_text.splitlines():
        match = pattern.match(line)
        if match and re.search(r'\bCa_' + tag + r'@OUT\b', match[2]):
            groups[next(g for g in groups if g in match[2])].append(match[1])
    return groups
