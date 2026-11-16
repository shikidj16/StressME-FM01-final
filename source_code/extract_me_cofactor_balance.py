"""Extract steady-state cofactor reaction balances from two StressME solutions.

Run inside queensysbio/stressme:v1.1 with the result parent mounted at /results.
"""
from __future__ import print_function
import csv, json, os, pickle
from collections import defaultdict
from sympy import Symbol

ROOT = os.environ.get('ME_RESULTS', '/results')
IJO_PATH = os.environ.get(
    'IJO_REACTION_METADATA',
    '/source/ecolime/ecolime/building_data/iJO1366.json')
COFACTORS = ['atp_c', 'nadh_c', 'nadph_c', 'fadh2_c', 'q8h2_c']
LABELS = {'sim_01': 'fmdA + ammonium', 'sim_02': 'fmdA + formamide'}
mu_symbol = Symbol('mu')

state_map_file = os.environ.get('ME_STATE_MAP_FILE')
state_map_json = os.environ.get('ME_STATE_MAP_JSON')
if state_map_file:
    STATE_MAP = json.load(open(state_map_file))
elif state_map_json:
    STATE_MAP = json.loads(state_map_json)
else:
    STATE_MAP = {
        'sim_01': {'path': 'sim_01', 'label': LABELS['sim_01']},
        'sim_02': {'path': 'sim_02', 'label': LABELS['sim_02']},
    }


def fallback_process_name(reaction_id):
    prefixes = {
        'translation_': 'Translation of ',
        'transcription_': 'Transcription of ',
        'charging_': 'tRNA charging: ',
        'formation_': 'Complex formation: ',
    }
    for prefix, label in prefixes.items():
        if reaction_id.startswith(prefix):
            return label + reaction_id[len(prefix):]
    return 'ME process: ' + reaction_id


def evaluated_reaction_equation(reaction, mu):
    reactants, products = [], []
    for metabolite, coefficient in reaction._metabolites.items():
        try:
            if hasattr(coefficient, 'subs'):
                coefficient = float(coefficient.subs(mu_symbol, mu))
            else:
                coefficient = float(coefficient)
        except Exception:
            continue
        if abs(coefficient) < 1e-15:
            continue
        magnitude = abs(coefficient)
        token = ('' if abs(magnitude - 1.0) < 1e-12 else '{:.8g} '.format(magnitude)) + metabolite.id
        (reactants if coefficient < 0 else products).append(token)
    return ' + '.join(reactants) + ' -> ' + ' + '.join(products)


def load_ijo_metadata():
    names, equations = {}, {}
    if not os.path.exists(IJO_PATH):
        return names, equations
    data = json.load(open(IJO_PATH))
    for reaction in data.get('reactions', []):
        reaction_id = reaction['id']
        names[reaction_id] = reaction.get('name', '') or ''
        reactants, products = [], []
        for metabolite, coefficient in reaction.get('metabolites', {}).items():
            magnitude = abs(float(coefficient))
            token = ('' if magnitude == 1 else '{} '.format(magnitude)) + metabolite
            (reactants if coefficient < 0 else products).append(token)
        equations[reaction_id] = ' + '.join(reactants) + ' -> ' + ' + '.join(products)
    return names, equations


ijo_names, ijo_equations = load_ijo_metadata()

rows = []
summaries = []
for sim_id, state_spec in STATE_MAP.items():
    state_path = state_spec.get('path', sim_id)
    condition_label = state_spec.get('label', sim_id)
    path = os.path.join(ROOT, state_path, 'solution.pickle')
    if not os.path.exists(path):
        continue
    solution = pickle.load(open(path, 'rb'))
    fluxes = getattr(solution, 'x_dict', None) or solution
    mu = float(fluxes.get('biomass_dilution', 0))
    model_path = os.path.join(ROOT, state_path, 'configured_model.pickle')
    if not os.path.exists(model_path):
        raise IOError('Missing exact configured model: {}'.format(model_path))
    me = pickle.load(open(model_path, 'rb'))
    reaction_to_root = {}
    root_equations = {}
    root_names = {}
    for reaction in me.reactions:
        root = reaction.id
        if type(reaction).__name__ == 'MetabolicReaction':
            data = getattr(reaction, 'stoichiometric_data', None)
            if data is not None:
                root = data.id
        reaction_to_root[reaction.id] = root
        root_names.setdefault(root, getattr(reaction, 'name', '') or '')
    for data in me.process_data:
        if type(data).__name__ != 'StoichiometricData':
            continue
        reactants, products = [], []
        for metabolite, coefficient in (getattr(data, '_stoichiometry', {}) or {}).items():
            if not isinstance(coefficient, (int, float)):
                continue
            token = ('{} '.format(abs(coefficient)) if abs(coefficient) != 1 else '') + metabolite
            (reactants if coefficient < 0 else products).append(token)
        root_equations[data.id] = ' + '.join(reactants) + ' -> ' + ' + '.join(products)
    for root in set(reaction_to_root.values()):
        if not root_names.get(root) and ijo_names.get(root):
            root_names[root] = ijo_names[root]
        if not root_equations.get(root) and ijo_equations.get(root):
            root_equations[root] = ijo_equations[root]
    for reaction in me.reactions:
        root = reaction_to_root[reaction.id]
        if not root_names.get(root):
            root_names[root] = fallback_process_name(root)
        if not root_equations.get(root):
            equation = evaluated_reaction_equation(reaction, mu)
            if equation != ' -> ':
                root_equations[root] = equation
    coefficient_index = {}
    for reaction in me.reactions:
        values = {metabolite.id: coefficient for metabolite, coefficient in reaction._metabolites.items()
                  if metabolite.id in COFACTORS}
        if values:
            coefficient_index[reaction.id] = values
    for cofactor in COFACTORS:
        aggregate = defaultdict(lambda: {'flux': 0.0, 'contribution': 0.0, 'variants': 0})
        variant_contributions = []
        for reaction_id, coefficients in coefficient_index.items():
            if cofactor not in coefficients:
                continue
            flux = float(fluxes.get(reaction_id, 0.0))
            if abs(flux) < 1e-15:
                continue
            coefficient = coefficients[cofactor]
            if hasattr(coefficient, 'subs'):
                try:
                    coefficient = float(coefficient.subs(mu_symbol, mu))
                except Exception:
                    continue
            else:
                coefficient = float(coefficient)
            contribution = coefficient * flux
            root = reaction_to_root.get(reaction_id, reaction_id)
            aggregate[root]['flux'] += flux
            aggregate[root]['contribution'] += contribution
            aggregate[root]['variants'] += 1
            variant_contributions.append(contribution)

        production = sum(value for value in variant_contributions if value > 0)
        consumption = sum(-value for value in variant_contributions if value < 0)
        turnover = 0.5 * (production + consumption)
        imbalance = production - consumption
        summaries.append({
            'simulation': sim_id, 'condition': condition_label, 'cofactor': cofactor,
            'growth_rate_h-1': mu, 'steady_state_turnover': turnover,
            'production_total_raw': production, 'consumption_total_raw': consumption,
            'numerical_imbalance': imbalance,
        })
        for root, info in aggregate.items():
            contribution = info['contribution']
            if abs(contribution) < 1e-15:
                continue
            role = 'production' if contribution > 0 else 'consumption'
            rows.append({
                'simulation': sim_id, 'condition': condition_label, 'cofactor': cofactor,
                'reaction_id': root, 'reaction_name': root_names.get(root, ''),
                'equation': root_equations.get(root, ''),
                'flux_sum_across_variants': info['flux'],
                'cofactor_contribution': contribution, 'role': role,
                'percent_of_turnover_side': (100.0 * abs(contribution) / turnover
                                             if turnover > 1e-15 else 0.0),
                'n_variants': info['variants'],
            })

analysis = os.environ.get('ME_ANALYSIS_OUT', os.path.join(ROOT, 'cofactor_analysis'))
if not os.path.exists(analysis):
    os.makedirs(analysis)
with open(os.path.join(analysis, 'ME_cofactor_reaction_balance_long.csv'), 'w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=[
        'simulation', 'condition', 'cofactor', 'reaction_id', 'reaction_name',
        'equation', 'flux_sum_across_variants', 'cofactor_contribution', 'role',
        'percent_of_turnover_side', 'n_variants'])
    writer.writeheader(); writer.writerows(rows)
with open(os.path.join(analysis, 'ME_cofactor_turnover_summary.csv'), 'w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=[
        'simulation', 'condition', 'cofactor', 'growth_rate_h-1',
        'steady_state_turnover', 'production_total_raw', 'consumption_total_raw',
        'numerical_imbalance'])
    writer.writeheader(); writer.writerows(summaries)
print('Wrote {} reaction rows and {} turnover rows'.format(len(rows), len(summaries)))
