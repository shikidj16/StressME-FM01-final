"""
Parameterized StressME simulation runner.

Usage (inside Docker container):
  python run_stressme.py --ros 1.0 [--no-mfa] --out /mount_point/

Always WT, T=37C, pH=7.0. Output goes to --out directory as:
  solution.pickle, phenotype.csv, fluxes.csv, proteome.csv,
  mfa_comparison.csv (only if --no-mfa not set).

MFA constraints come from the wt-mid sheet, scaled by /10 so glucose
uptake target = 10 mmol/gDW/h. PGM-direction sign fix: R9 (G3P->PEP)
maps to ENO (whose BiGG forward direction is glycolytic), not PGM
(whose BiGG forward is gluconeogenic).
"""
from __future__ import print_function
import argparse
import csv
import math
import os
import pickle
import re
import sys

import numpy as np
import pandas as pd
import openpyxl

import cobrame
import ecolime
from cobrame.core.reaction import SummaryVariable
from cobrame.core.component import Constraint

from AcidifyME.periplasmic_proteome import add_protein_stability_constraint
from AcidifyME.membrane_lipid_constraint import add_membrane_constraint
from AcidifyME.membrane_protein_activity import modify_membrane_protein_activity
from AcidifyME.proton_influx import add_proton_leak_rxn
from oxidizeme.model import StressME
from ecolime import chaperones
from qminospy.me1 import ME_NLP1


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--ros', type=float, required=True,
                    help='ROS multiplier (e.g. 0.01, 1.0, 10.0)')
parser.add_argument('--no-mfa', action='store_true',
                    help='If set, skip MFA constraints (baseline run)')
parser.add_argument('--out', required=True,
                    help='Output directory (mounted from host)')
parser.add_argument('--mfa-sheet', default='1.# wt-mid')
parser.add_argument('--mfa-workbook', default='/mfa/mfa_copy.xlsx',
                    help='MFA workbook path inside the container. The default preserves '
                         'historical runs; updated analyses can select another mounted workbook.')
parser.add_argument('--mfa-scale', type=float, default=10.0)
parser.add_argument('--mfa-tol', type=float, default=0.10)
parser.add_argument('--n-source', default='ammonium',
                    choices=['ammonium', 'formamide'],
                    help='Nitrogen source for the cell')
parser.add_argument('--add-formamide-direct', action='store_true',
                    help='Add formamide -> formate via FmdA (formamidase)')
parser.add_argument('--add-formamide-twostep', action='store_true',
                    help='Add formamide -> formaldehyde -> formate (two-step)')
parser.add_argument('--fmd-uptake', type=float, default=10.0,
                    help='Max formamide uptake rate (mmol/gDW/h)')
parser.add_argument('--add-betaxanthin', action='store_true',
                    help='Add betaxanthin synthesis pathway (R75/R76/R77) + EX_tybx_e')
parser.add_argument('--constrain-acetate', action='store_true',
                    help='Constrain EX_ac_e to MFA R70 from selected --mfa-sheet (overflow constraint)')
parser.add_argument('--sim-number', type=int, default=None,
                    help='Optional simulation number for config.json traceability')
parser.add_argument('--strain-label', default='wt',
                    help='Strain label written to config.json (informational)')
parser.add_argument('--medium-label', default='M9-MOPS-glc-NH4',
                    help='Medium label written to config.json (informational)')
parser.add_argument('--temperature', type=float, default=37.0,
                    help='Temperature in Celsius (default 37)')
parser.add_argument('--ph', type=float, default=7.0,
                    help='Extracellular pH (default 7.0)')
parser.add_argument('--hard-glucose-uptake', action='store_true',
                    help='Pin glucose uptake to exactly -scaled (lb=ub), ignoring mfa_tol')
parser.add_argument('--mfa-biomass', action='store_true',
                    help='Constrain biomass_dilution to MFA R74 (scaled to glucose=10)')
parser.add_argument('--biomass-tol', type=float, default=0.0,
                    help='Tolerance around the MFA biomass target (e.g. 0.01 = +/-1%). '
                         '0 (default) = single-point pin (lb=ub=target).')
parser.add_argument('--bisect-precision', type=float, default=1e-4,
                    help='qMINOS bisection precision (default 1e-4; tighter = slower)')
parser.add_argument('--basis-file', default=None,
                    help='Optional qMINOS basis .npy file used to warm-start bisection.')
parser.add_argument('--save-basis', action='store_true',
                    help='Save the final qMINOS basis as qminos_basis.npy in --out.')
parser.add_argument('--fixed-growth-feasibility', action='store_true',
                    help='With --pin-growth, solve one direct qMINOS LP at that same '
                         'symbolic mu instead of independently bisecting mu.')
parser.add_argument('--maximize-product', default=None,
                    help='Reaction ID to maximize (e.g. EX_tybx_e for betaxanthin). '
                         'When set, run LP at fixed mu (from --mfa-biomass) with this '
                         'reaction as objective instead of bisectmu(max growth).')
parser.add_argument('--betax-mfa-lower-only', action='store_true',
                    help='For the R75/R76/R77 betaxanthin MFA constraints, keep the '
                         'lower bound (must produce at least the MFA value) but leave '
                         'the upper bound unbounded. Lets max-product go higher than the '
                         'measured MFA flux.')
parser.add_argument('--skip-betax-mfa', action='store_true',
                    help='Add the betaxanthin pathway but do not constrain R75/R76/R77 '
                         'from MFA. Use when EX_tybx_e is pinned independently from '
                         'experimental exchange data.')
parser.add_argument('--mfa-expanded-tiered', action='store_true',
                    help='Apply the direction-audited MFA37 core: 22 additional '
                         'direct intracellular constraints plus R39 and R66-net. '
                         'R72 CO2 exchange is deliberately excluded. Near-zero '
                         'estimates use absolute bounds and R44 uses +/-25%%.')
parser.add_argument('--acetate-lower-only', action='store_true',
                    help='For the acetate secretion constraint (R70 -> EX_ac_e), keep '
                         'the lower bound (must secrete at least MFA measured value) but '
                         'leave the upper bound unbounded.')
parser.add_argument('--only-biomass-mfa', action='store_true',
                    help='Keep ONLY the biomass (R74) MFA constraint. Skip all '
                         'per-reaction MFA flux bounds (R2-R37, R75-R77) and skip the '
                         'acetate overflow constraint (R70). --hard-glucose-uptake still '
                         'works via R1 special-case. Use to compare against fully MFA-'
                         'constrained sims.')
parser.add_argument('--target-growth', type=float, default=None,
                    help='Anchor MFA scaling to a measured growth rate mu (h^-1). '
                         'Reads R74 raw value from the selected MFA sheet, then sets '
                         'mfa_scale = R74_raw / target_growth so that biomass (R74/mfa_scale) '
                         '== target_growth, and every other MFA flux is rescaled by the '
                         'same factor. Implies --mfa-biomass. Overrides any --mfa-scale value.')
parser.add_argument('--free-glucose', action='store_true',
                    help='Do NOT constrain glucose uptake (R1) from MFA. Leaves '
                         'EX_glc__D_e at model default [-1000, 1000] so the model '
                         'chooses GUR. Overrides --hard-glucose-uptake.')
parser.add_argument('--exchange-json', default=None,
                    help='Path to a JSON file specifying data-derived exchange bounds. '
                         'Format: {"exchanges": {"EX_glc__D_e": {"q": -10.9, "sigma": 0.65}, '
                         '"EX_ac_e": {"q": +2.66, "sigma": 0.59}, ...}, ...}. '
                         'Applied as [q - N*sigma, q + N*sigma] with N controlled by '
                         '--exchange-ci-mult. Sign convention is BiGG (negative = uptake, '
                         'positive = secretion). Overrides any conflicting MFA-derived bound.')
parser.add_argument('--exchange-ci-mult', type=float, default=1.0,
                    help='Confidence-interval multiplier for exchange bounds (default 1.0 '
                         '= +/-1 sigma; use 1.96 for 95%% CI).')
parser.add_argument('--pin-growth', type=float, default=None,
                    help='Pin biomass_dilution to this mu (h^-1), independent of any MFA '
                         'logic. If both --pin-growth and --target-growth are set, '
                         '--pin-growth wins.')
parser.add_argument('--skip-exchange', action='append', default=[],
                    help='Reaction ID to SKIP when applying --exchange-json bounds '
                         '(can be repeated: --skip-exchange EX_nh4_e --skip-exchange EX_tybx_e). '
                         'Used to build "no-NH4" or "tybx-free" variants from the same JSON.')
parser.add_argument('--close-nh4-exchange', action='store_true',
                    help='Fix EX_nh4_e to exactly zero after all nitrogen-source and '
                         'exchange-JSON bounds are applied. Required for formamide-only '
                         'nitrogen simulations when ammonium secretion is also prohibited.')
parser.add_argument('--close-exchange', action='append', default=[],
                    help='Fix a named exchange reaction to exactly [0, 0] after all '
                         'other bounds are applied. May be repeated for sensitivity tests.')
parser.add_argument('--reaction-bound', action='append', default=[],
                    metavar='REACTION,LB,UB',
                    help='Override bounds for any named reaction after all other '
                         'constraints are applied. Repeatable; intended for isolated '
                         'reaction-direction and cycle sensitivity tests.')
parser.add_argument('--block-reaction-direction', action='append', default=[],
                    metavar='ROOT,DIRECTION',
                    help='Block every enzyme-coupled ME reaction variant for a metabolic '
                         'root in the selected BiGG direction: forward, reverse, or all. '
                         'Repeatable; e.g. LDH_D,reverse.')
parser.add_argument('--close-nonacetate-carbon-overflow', action='store_true',
                    help='Block conventional non-acetate overflow secretion through '
                         'ethanol, D/L-lactate, formate, pyruvate, acetaldehyde, and '
                         'both 1,2-propanediol exchanges. Succinate remains open because '
                         'its closure makes the MFA37 formamide state infeasible.')
parser.add_argument('--overflow-percentage-limits', action='store_true',
                    help='Constrain selected secreted organic-carbon shares: acetate '
                         '>=85%%; lactate, ethanol, formate, succinate, and grouped '
                         'pyruvate/acetaldehyde/1,2-propanediol each <=5%%. This is '
                         'separate from the hard overflow-closure option.')
parser.add_argument('--save-configured-model', action='store_true',
                    help='Save configured_model.pickle before rebuilding the qMINOS solver, '
                         'then continue with the simulation. Used for exact post-solution '
                         'stoichiometric and cofactor accounting.')
parser.add_argument('--save-configured-model-only', action='store_true',
                    help='Save the fully configured ME model as configured_model.pickle '
                         'and exit before rebuilding or solving. Intended for exact '
                         'post-solution stoichiometric audits.')
args = parser.parse_args()

TEMPERATURE_C = args.temperature
pH = args.ph
ROS = args.ros
USE_MFA = not args.no_mfa
OUT = args.out
os.makedirs(OUT, exist_ok=True)

print('=' * 70)
print('StressME run')
print('  Sim number:  {}'.format(args.sim_number))
print('  Strain:      {}'.format(args.strain_label))
print('  Medium:      {}'.format(args.medium_label))
print('  Temperature: {} C'.format(TEMPERATURE_C))
print('  pH:          {}'.format(pH))
print('  ROS:         {}x basal'.format(ROS))
print('  N source:    {}'.format(args.n_source))
print('  Fmd direct:  {}'.format(args.add_formamide_direct))
print('  Fmd 2-step:  {}'.format(args.add_formamide_twostep))
print('  Betaxanthin: {}'.format(args.add_betaxanthin))
print('  MFA:         {}'.format('CONSTRAINED ({})'.format(args.mfa_sheet) if USE_MFA else 'unconstrained'))
print('  Acetate ovf: {}'.format(args.constrain_acetate))
print('  Output dir:  {}'.format(OUT))
print('=' * 70)

# ----------------------------------------------------------------------
# MFA-to-ME mapping
# ----------------------------------------------------------------------
# IMPORTANT BiGG sign conventions:
#   PGM BiGG forward: 2pg -> 3pg  (gluconeogenic) -- glycolysis runs in reverse
#   PGK BiGG forward: 3pg + ATP -> 13dpg + ADP  (gluconeogenic) -- glycolysis runs in reverse
#   ENO BiGG forward: 2pg -> pep + h2o  (GLYCOLYTIC, correct sign)
#   GAPD BiGG forward: g3p + nad + pi -> 13dpg + nadh + h  (GLYCOLYTIC, correct sign)
# MFA reports fluxes in glycolytic direction (positive).
# Lumped pairs: constrain ONE of each lumped pair; mass balance handles the other.
MAPPINGS = {
    # Glycolysis (8 reactions)
    'R1':      'EX_glc__D_e',  # exchange, special handling (negative for uptake)
    'R2 net':  'PGI',
    'R3':      'PFK',
    'R6 net':  'FBA',
    'R7 net':  'TPI',
    'R8 net':  'GAPD',   # PGK lumped; GAPD has correct BiGG forward (glycolytic)
    'R9 net':  'ENO',    # FIXED: use ENO (glycolytic forward) not PGM (gluconeogenic forward)
    'R10':     'PYK',
    # PDH + upper TCA (4 reactions)
    'R14':     'PDH',
    'R15':     'CS',
    'R16 net': 'ACONTa', # ACONTb lumped
    'R17 net': 'ICDHyr',
    # Lower TCA -- sign-flipped (see SIGN_FLIPS below)
    'R19 net': 'SUCOAS', # MFA forward = succoa->succ (TCA fwd, ATP+);
                         # BiGG forward = succ->succoa (gluconeogenic, ATP-)
    # PPP -- sign-flipped (see SIGN_FLIPS below)
    'R34 net': 'RPI',    # MFA forward = ru5p->r5p; BiGG forward = r5p->ru5p
    # Still dropped (model topology disagrees with MFA, would re-introduce infeasibility):
    #   R18 AKGDH, R22 MDH, R23 ME1
    #   R30 G6PDH2r, R31 GND, R35 TKT1, R37 TKT2
}  # 14 total

# MFA reactions whose forward direction is opposite to the BiGG-ME forward
# direction. Constraint and parity comparison both apply -1 multiplier so
# MFA "+X (forward)" becomes ME "-X (= BiGG reverse net)".
SIGN_FLIPS = {'R19 net', 'R34 net'}

EXPANDED_TIERED_MAPPINGS = {
    'R18': 'AKGDH', 'R20 net': 'SUCDi', 'R21 net': 'FUM', 'R22 net': 'MDH',
    'R23': 'ME1', 'R25': 'PPC', 'R26': 'PPCK', 'R28': 'ICL', 'R29': 'MALS',
    'R30': 'G6PDH2r', 'R31': 'GND', 'R33 net': 'RPE', 'R35 net': 'TKT1',
    'R36 net': 'TALA', 'R37 net': 'TKT2', 'R38': 'EDD',
    'R40': 'GLUDy', 'R44': 'ASPTA', 'R46': 'ALATA_L',
    'R48 net': 'GHMT2r', 'R49 net': 'GLYCL', 'R50': 'MTHFR2',
    # Tier-2 accepted gross acetate-formation representation.
    'R39': 'ACKr',
}
EXPANDED_SIGN_FLIPS = {'R39', 'R40', 'R44', 'R46'}
EXPANDED_ABSOLUTE_ZERO = {'R26', 'R28', 'R29', 'R38'}

# Betaxanthin reactions - added to MAPPINGS only when --add-betaxanthin is set
# Reaction names match the ones added by add_betaxanthin_pathway()
BETAXANTHIN_MAPPINGS = {
    'R75': 'BETAX_TYR_DOPA',
    'R76': 'BETAX_DOPA_BTLA',
    'R77': 'BETAX_BTLA_TYBX',
}


def add_formamide_metabolites_and_transport(me):
    """Add formamide metabolite + transports + exchange. Idempotent."""
    from cobra import Metabolite
    from cobrame.core.reaction import MEReaction

    # Metabolites (formula CH3NO for formamide HCONH2)
    for compt in ['c', 'p', 'e']:
        mid = 'fmd_' + compt
        if mid not in [m.id for m in me.metabolites]:
            m = Metabolite(mid)
            m.name = 'Formamide'
            m.formula = 'CH3NO'
            m.compartment = compt
            me.add_metabolites([m])
            print('  added metabolite', mid)

    # Outer-membrane transport (extracellular <-> periplasm), passive diffusion
    if 'FMDtex' not in [r.id for r in me.reactions]:
        r = MEReaction('FMDtex')
        me.add_reaction(r)
        r.add_metabolites({me.metabolites.fmd_e: -1, me.metabolites.fmd_p: 1},
                          combine=False)
        r.lower_bound = -1000; r.upper_bound = 1000  # reversible
        print('  added reaction FMDtex')

    # Inner-membrane transport (periplasm <-> cytoplasm), passive diffusion
    if 'FMDtpp' not in [r.id for r in me.reactions]:
        r = MEReaction('FMDtpp')
        me.add_reaction(r)
        r.add_metabolites({me.metabolites.fmd_p: -1, me.metabolites.fmd_c: 1},
                          combine=False)
        r.lower_bound = -1000; r.upper_bound = 1000
        print('  added reaction FMDtpp')

    # Boundary exchange
    if 'EX_fmd_e' not in [r.id for r in me.reactions]:
        r = MEReaction('EX_fmd_e')
        me.add_reaction(r)
        r.add_metabolites({me.metabolites.fmd_e: -1}, combine=False)
        r.lower_bound = 0; r.upper_bound = 1000  # bounds set by --n-source later
        print('  added reaction EX_fmd_e')


def add_formamide_pathway_direct(me):
    """Pathway A: formamide -> formate via formamidase (FmdA, EC 3.5.1.49)."""
    from cobrame.core.reaction import MEReaction
    add_formamide_metabolites_and_transport(me)
    if 'FMDA' not in [r.id for r in me.reactions]:
        r = MEReaction('FMDA')
        me.add_reaction(r)
        # HCONH2 + H2O -> HCOOH + NH3 (released as NH4+ at physiol pH)
        r.add_metabolites({
            me.metabolites.fmd_c:  -1,
            me.metabolites.h2o_c:  -1,
            me.metabolites.for_c:   1,
            me.metabolites.nh4_c:   1,
        }, combine=False)
        r.lower_bound = 0; r.upper_bound = 1000
        print('  added reaction FMDA (formamidase): fmd + h2o -> for + nh4')


def add_formamide_pathway_twostep(me):
    """Pathway B: formamide -> formaldehyde (FMDR) then formaldehyde -> formate
    via existing GSH-dependent machinery (FALDH2 chain) already in iJL1678.
    Plus an explicit FALDH_simple as backup."""
    from cobrame.core.reaction import MEReaction
    add_formamide_metabolites_and_transport(me)

    # Step 1: formamide -> formaldehyde + NH4 (hypothetical reductase)
    if 'FMDR' not in [r.id for r in me.reactions]:
        r = MEReaction('FMDR')
        me.add_reaction(r)
        # HCONH2 + NADH + H+ -> HCHO + NH3 + NAD+
        # Atom balance: LHS C1H3N1O1 + 2H(NADH) + 1H(h+) = C1H6N1O1
        #               RHS C1H2O1 + N1H4 = C1H6N1O1 ✓
        r.add_metabolites({
            me.metabolites.fmd_c:   -1,
            me.metabolites.nadh_c:  -1,
            me.metabolites.h_c:     -1,
            me.metabolites.fald_c:   1,
            me.metabolites.nh4_c:    1,
            me.metabolites.nad_c:    1,
        }, combine=False)
        r.lower_bound = 0; r.upper_bound = 1000
        print('  added reaction FMDR: fmd + nadh + h -> fald + nh4 + nad')

    # Step 2: explicit formaldehyde dehydrogenase (backup to GSH cycle)
    if 'FALDH_simple' not in [r.id for r in me.reactions]:
        r = MEReaction('FALDH_simple')
        me.add_reaction(r)
        # HCHO + NAD+ + H2O -> HCOOH + NADH + H+
        r.add_metabolites({
            me.metabolites.fald_c:  -1,
            me.metabolites.nad_c:   -1,
            me.metabolites.h2o_c:   -1,
            me.metabolites.for_c:    1,
            me.metabolites.nadh_c:   1,
            me.metabolites.h_c:      1,
        }, combine=False)
        r.lower_bound = 0; r.upper_bound = 1000
        print('  added reaction FALDH_simple: fald + nad + h2o -> for + nadh + h')


def set_nitrogen_source(me, source, fmd_uptake_max):
    """Configure N-source: ammonium (default) or formamide.

    'formamide': close EX_nh4_e uptake (lb=0), open EX_fmd_e uptake (lb=-X).
    'ammonium' : default behavior, no changes.
    """
    if source == 'formamide':
        nh4 = me.reactions.EX_nh4_e
        nh4.lower_bound = 0       # no NH4 uptake
        nh4.upper_bound = 1000    # secretion still allowed
        if 'EX_fmd_e' in [r.id for r in me.reactions]:
            fmd = me.reactions.EX_fmd_e
            fmd.lower_bound = -fmd_uptake_max
            fmd.upper_bound = 1000
            print('  N-source: formamide  (EX_nh4_e lb=0, EX_fmd_e lb={})'.format(
                -fmd_uptake_max))
        else:
            print('  WARN: --n-source formamide but EX_fmd_e not in model;'
                  ' did you forget --add-formamide-direct/--add-formamide-twostep?')


def add_betaxanthin_pathway(me):
    """Add R75/R76/R77 betaxanthin biosynthesis with real biochemistry, plus
    transports + EX_tybx_e secretion. No enzyme coupling (plain MEReactions).

    Pathway:
      R75 (HpaBC, tyrosine hydroxylase):
          tyr__L_c + o2_c + nadph_c + h_c -> 34dhphe_c + nadp_c + h2o_c
      R76 (DOD, DOPA dioxygenase):
          34dhphe_c + o2_c -> btla_c + h2o_c
      R77 (spontaneous Schiff base):
          btla_c + tyr__L_c -> tybx_c + h2o_c
      Transport: tybx_c -> tybx_p -> tybx_e
      Exchange:  EX_tybx_e (secretion only)
    """
    from cobra import Metabolite
    from cobrame.core.reaction import MEReaction

    # New metabolites
    new_mets = [
        ('34dhphe_c', 'L-DOPA',                  'C9H11NO4', 'c'),
        ('btla_c',    'Betalamic acid',          'C9H9NO5',  'c'),
        ('tybx_c',    'Tyrosine-betaxanthin',    'C18H18N2O7', 'c'),
        ('tybx_p',    'Tyrosine-betaxanthin',    'C18H18N2O7', 'p'),
        ('tybx_e',    'Tyrosine-betaxanthin',    'C18H18N2O7', 'e'),
    ]
    for mid, name, formula, compt in new_mets:
        if mid not in [m.id for m in me.metabolites]:
            m = Metabolite(mid)
            m.name = name; m.formula = formula; m.compartment = compt
            me.add_metabolites([m])
            print('  added metabolite', mid)

    # R75: TYR + O2 + NADPH + H -> DOPA + NADP + H2O
    if 'BETAX_TYR_DOPA' not in [r.id for r in me.reactions]:
        r = MEReaction('BETAX_TYR_DOPA')
        r.name = 'HpaBC tyrosine hydroxylation (R75)'
        me.add_reaction(r)
        r.add_metabolites({
            me.metabolites.tyr__L_c:  -1,
            me.metabolites.o2_c:      -1,
            me.metabolites.nadph_c:   -1,
            me.metabolites.h_c:       -1,
            me.metabolites.get_by_id('34dhphe_c'): 1,
            me.metabolites.nadp_c:     1,
            me.metabolites.h2o_c:      1,
        }, combine=False)
        r.lower_bound = 0; r.upper_bound = 1000
        print('  added reaction BETAX_TYR_DOPA (R75)')

    # R76: DOPA + O2 -> BTLA + H2O
    if 'BETAX_DOPA_BTLA' not in [r.id for r in me.reactions]:
        r = MEReaction('BETAX_DOPA_BTLA')
        r.name = 'DOD DOPA dioxygenase (R76)'
        me.add_reaction(r)
        r.add_metabolites({
            me.metabolites.get_by_id('34dhphe_c'): -1,
            me.metabolites.o2_c:                   -1,
            me.metabolites.btla_c:                  1,
            me.metabolites.h2o_c:                   1,
        }, combine=False)
        r.lower_bound = 0; r.upper_bound = 1000
        print('  added reaction BETAX_DOPA_BTLA (R76)')

    # R77: BTLA + TYR -> TYBX + H2O
    if 'BETAX_BTLA_TYBX' not in [r.id for r in me.reactions]:
        r = MEReaction('BETAX_BTLA_TYBX')
        r.name = 'Tyrosine-betaxanthin condensation (R77, spontaneous)'
        me.add_reaction(r)
        r.add_metabolites({
            me.metabolites.btla_c:   -1,
            me.metabolites.tyr__L_c: -1,
            me.metabolites.tybx_c:    1,
            me.metabolites.h2o_c:     1,
        }, combine=False)
        r.lower_bound = 0; r.upper_bound = 1000
        print('  added reaction BETAX_BTLA_TYBX (R77)')

    # Transport tybx_c -> tybx_p -> tybx_e
    if 'TYBXtipp' not in [r.id for r in me.reactions]:
        r = MEReaction('TYBXtipp')
        r.name = 'Tyrosine-betaxanthin transport (cytoplasm -> periplasm)'
        me.add_reaction(r)
        r.add_metabolites({
            me.metabolites.tybx_c: -1,
            me.metabolites.tybx_p:  1,
        }, combine=False)
        r.lower_bound = 0; r.upper_bound = 1000
        print('  added reaction TYBXtipp')

    if 'TYBXtex' not in [r.id for r in me.reactions]:
        r = MEReaction('TYBXtex')
        r.name = 'Tyrosine-betaxanthin transport (periplasm -> extracellular)'
        me.add_reaction(r)
        r.add_metabolites({
            me.metabolites.tybx_p: -1,
            me.metabolites.tybx_e:  1,
        }, combine=False)
        r.lower_bound = 0; r.upper_bound = 1000
        print('  added reaction TYBXtex')

    if 'EX_tybx_e' not in [r.id for r in me.reactions]:
        r = MEReaction('EX_tybx_e')
        r.name = 'Tyrosine-betaxanthin exchange (secretion only)'
        me.add_reaction(r)
        r.add_metabolites({me.metabolites.tybx_e: -1}, combine=False)
        r.lower_bound = 0; r.upper_bound = 1000
        print('  added reaction EX_tybx_e')


def constrain_acetate_overflow(me, mfa_xlsx, sheet_name, scale, tolerance,
                                lower_only=False):
    """Bound EX_ac_e using R70 from chosen MFA sheet. EX_ac_e is positive for
    secretion in BiGG, so both bounds are positive.

    If lower_only=True, the upper bound is left unbounded (secretion must be
    at least the MFA value; no cap)."""
    import openpyxl
    wb = openpyxl.load_workbook(mfa_xlsx, read_only=True, data_only=True)
    if sheet_name not in wb.sheetnames:
        print('  WARN: sheet {} not in workbook; skipping acetate constraint'.format(sheet_name))
        return
    mfa_values, _ = read_mfa_fitted_table(wb[sheet_name])
    r70_raw = mfa_values.get('R70')
    if r70_raw is None:
        print('  WARN: R70 not found in sheet {}; skipping acetate constraint'.format(sheet_name))
        return
    scaled = r70_raw / scale
    lb = scaled * (1 - tolerance)
    ub = 1000.0 if lower_only else scaled * (1 + tolerance)
    ex = me.reactions.EX_ac_e
    ex.lower_bound = lb
    ex.upper_bound = ub
    tag = ' [LOWER-BOUND ONLY]' if lower_only else ''
    print('  acetate overflow: EX_ac_e in [{:.3f}, {:.3f}]  (R70 raw={:.3f}, scaled={:.3f}){}'.format(
        lb, ub, r70_raw, scaled, tag))


def read_mfa_fitted_table(worksheet):
    """Read fitted net fluxes from legacy A:C or updated E:J INCA sheets."""
    rows = list(worksheet.iter_rows(values_only=True))
    for header_index, row in enumerate(rows[:10]):
        normalized = [str(value or '').strip().lower() for value in row]
        for id_index, value in enumerate(normalized):
            if value != 'id':
                continue
            equation_index = id_index + 1
            value_index = id_index + 2
            if (value_index >= len(normalized) or
                    normalized[equation_index] != 'equation' or
                    normalized[value_index] != 'value'):
                continue
            values, equations = {}, {}
            for data_row in rows[header_index + 1:]:
                if len(data_row) <= value_index:
                    continue
                reaction_id = data_row[id_index]
                raw_value = data_row[value_index]
                if reaction_id is None or raw_value is None:
                    continue
                try:
                    reaction_id = str(reaction_id).strip()
                    values[reaction_id] = float(raw_value)
                    equations[reaction_id] = str(data_row[equation_index] or '').strip()
                except (TypeError, ValueError):
                    continue
            if values:
                return values, equations
    raise ValueError('Could not locate fitted MFA ID/Equation/Value table in sheet {}'.format(
        worksheet.title))


def add_mfa_constraint(me, root_id, target_value, tolerance):
    """Force (sum_FWD - sum_REV) = target_value +/- tolerance via Constraint + drain.

    For native split reactions (e.g. PGI), finds _FWD_/_REV_ variants.
    For plain MEReactions (e.g. BETAX_TYR_DOPA, FMDA), uses the single reaction.
    """
    fwd = [r for r in me.reactions if r.id.startswith(root_id + '_FWD_')]
    rev = [r for r in me.reactions if r.id.startswith(root_id + '_REV_')]
    if not fwd and not rev:
        # Try the reaction itself (plain MEReaction case)
        try:
            r = me.reactions.get_by_id(root_id)
            fwd = [r]  # treat as a single forward variant
        except KeyError:
            print('  WARNING: no variants found for {}'.format(root_id))
            return None

    c_id = 'mfa_C_' + root_id
    try:
        constraint = me.metabolites.get_by_id(c_id)
    except KeyError:
        constraint = Constraint(c_id)
        me.add_metabolites([constraint])

    for r in fwd:
        r.add_metabolites({constraint: 1}, combine=False)
    for r in rev:
        r.add_metabolites({constraint: -1}, combine=False)

    drain_id = 'mfa_drain_' + root_id
    try:
        drain = me.reactions.get_by_id(drain_id)
    except KeyError:
        drain = SummaryVariable(drain_id)
        me.add_reactions([drain])
    drain.add_metabolites({constraint: -1}, combine=False)
    if target_value >= 0:
        drain.lower_bound = target_value * (1.0 - tolerance)
        drain.upper_bound = target_value * (1.0 + tolerance)
    else:
        drain.lower_bound = target_value * (1.0 + tolerance)
        drain.upper_bound = target_value * (1.0 - tolerance)
    return drain


def add_mfa_linear_constraint(me, constraint_id, weighted_roots,
                              target_value, tolerance):
    """Constrain a weighted sum of net ME metabolic reaction roots."""
    constraint = Constraint('mfa_C_' + constraint_id)
    me.add_metabolites([constraint])
    for root_id, weight in weighted_roots.items():
        fwd = [r for r in me.reactions if r.id.startswith(root_id + '_FWD_')]
        rev = [r for r in me.reactions if r.id.startswith(root_id + '_REV_')]
        if not fwd and not rev:
            try:
                fwd = [me.reactions.get_by_id(root_id)]
            except KeyError:
                print('  WARNING: no variants found for linear term {}'.format(root_id))
                return None
        for reaction in fwd:
            reaction.add_metabolites({constraint: weight}, combine=False)
        for reaction in rev:
            reaction.add_metabolites({constraint: -weight}, combine=False)
    drain = SummaryVariable('mfa_drain_' + constraint_id)
    me.add_reactions([drain])
    drain.add_metabolites({constraint: -1}, combine=False)
    lower = target_value * (1 + tolerance) if target_value < 0 else target_value * (1 - tolerance)
    upper = target_value * (1 - tolerance) if target_value < 0 else target_value * (1 + tolerance)
    drain.lower_bound, drain.upper_bound = lower, upper
    return drain


# ----------------------------------------------------------------------
# 1. Load model + WT setup
# ----------------------------------------------------------------------
print('Loading FoldME pickle...')
me = pickle.load(open('/home/meuser/FoldME_Ali_keff.pickle', 'rb'))

rxn = me.reactions.get_by_id('DXPRIi_FWD_DXPREDISOM-CPLX_mod_cobalt2')
rxn.keff = 0.01
rxn.update()

# ----------------------------------------------------------------------
# 1b. Add heterologous pathways (before stress modules so they're in the
#     model when make_stressme() and substitute_ros() run)
# ----------------------------------------------------------------------
if args.add_formamide_direct:
    print('Adding formamide -> formate (direct, FmdA):')
    add_formamide_pathway_direct(me)
if args.add_formamide_twostep:
    print('Adding formamide -> formaldehyde -> formate (two-step):')
    add_formamide_pathway_twostep(me)
if args.add_betaxanthin:
    print('Adding betaxanthin biosynthesis (R75/R76/R77 + EX_tybx_e):')
    add_betaxanthin_pathway(me)

# Switch N source AFTER pathways are added so EX_fmd_e exists
set_nitrogen_source(me, args.n_source, args.fmd_uptake)

# ----------------------------------------------------------------------
# 2. Apply stress modules
# ----------------------------------------------------------------------
print('Applying stress (T={}C, pH={}, ROS={}x)...'.format(TEMPERATURE_C, pH, ROS))
chaperones.change_temperature(me, TEMPERATURE_C)

temperatureK = 273.15 + TEMPERATURE_C
fold_change = (5.6666666666 * pH**4 - 139.6666666666 * pH**3
               + 1286.5833333 * pH**2 - 5254.0833333 * pH + 8037.0)
add_membrane_constraint(me, 'MJR_habituated')
add_protein_stability_constraint(me, pH, temperatureK, constrain_periplasm_biomass=False)
add_proton_leak_rxn(me, pH)
modify_membrane_protein_activity(me, fold_change)
me.unmodeled_protein_fraction = 0.1

solver_tmp = ME_NLP1(me)
stress = StressME(solver_tmp)
stress.make_stressme()

subs_dict = {'h2o2_c': 5e-08 * ROS, 'o2s_c': 2e-10 * ROS,
             'h2o2_e': 1.5e-06, 'o2s_e': 0, 'pq2_c': 0.0}

for pid in ['SPODM', 'SPODMpp', 'CAT']:
    for r in me.process_data.get_by_id(pid).parent_reactions:
        r.upper_bound = 1000
me.reactions.MOX_REV_CPLX_dummy.upper_bound = 0

# ----------------------------------------------------------------------
# 3. (optional) Apply MFA constraints
# ----------------------------------------------------------------------
if USE_MFA:
    print('\nLoading MFA data from sheet "{}"...'.format(args.mfa_sheet))
    wb = openpyxl.load_workbook(args.mfa_workbook, read_only=True, data_only=True)
    ws = wb[args.mfa_sheet]
    mfa, mfa_equations = read_mfa_fitted_table(ws)
    print('  Loaded {} MFA reactions'.format(len(mfa)))

    # --target-growth: rescale MFA data so R74/mfa_scale == target_growth.
    # This overrides --mfa-scale and forces --mfa-biomass on.
    if args.target_growth is not None:
        if 'R74' not in mfa:
            print('*** ERROR: --target-growth {:.4f} but R74 (biomass) not in sheet "{}"'.format(
                args.target_growth, args.mfa_sheet))
            sys.exit(1)
        if args.target_growth <= 0:
            print('*** ERROR: --target-growth must be > 0, got {}'.format(args.target_growth))
            sys.exit(1)
        r74_raw = mfa['R74']
        old_scale = args.mfa_scale
        args.mfa_scale = r74_raw / args.target_growth
        args.mfa_biomass = True
        print('  --target-growth {:.4f} h^-1: R74_raw={:.4f}  =>  mfa_scale = {:.4f} '
              '(was {:.4f}); --mfa-biomass forced ON.'.format(
                  args.target_growth, r74_raw, args.mfa_scale, old_scale))
        expected_glc = mfa.get('R1', 100.0) / args.mfa_scale
        print('  (implied glucose uptake if MFA-scaled: {:.3f} mmol gDW^-1 h^-1)'.format(
            expected_glc))
        if args.free_glucose:
            print('  --free-glucose is set: glucose bound will NOT be applied from MFA '
                  '(model chooses GUR).')

    # Build the active mapping from the reactions present in this workbook.
    active_mappings = dict(MAPPINGS)
    updated_pts_r1 = ('R1' in mfa and 'PEP' in mfa_equations.get('R1', '').upper()
                      and 'PYR' in mfa_equations.get('R1', '').upper())
    updated_glk_r1 = ('R1' in mfa and 'ATP' in mfa_equations.get('R1', '').upper()
                      and 'G6P' in mfa_equations.get('R1', '').upper()
                      and not updated_pts_r1)
    if updated_pts_r1:
        active_mappings['R1'] = 'GLCptspp'
        active_mappings.pop('R10', None)
        print('  Updated R1 stoichiometry detected: constraining summed GLCptspp '
              '(glucose + PEP -> G6P + PYR); obsolete R10/PYK removed.')
    elif updated_glk_r1:
        # BiGG/iJL1678b-ME calls ATP-dependent glucokinase HEX1.
        active_mappings['R1'] = 'HEX1'
        print('  Updated R1 stoichiometry detected: constraining summed HEX1 '
              '(glucose + ATP -> G6P); R10/PYK retained as an independent flux.')
    if args.mfa_expanded_tiered:
        active_mappings.update(EXPANDED_TIERED_MAPPINGS)
        print('  Expanded tiered profile enabled: {} mapped MFA rows before R66-net'.format(
            len(active_mappings)))
    if args.add_betaxanthin and not args.skip_betax_mfa:
        active_mappings.update(BETAXANTHIN_MAPPINGS)
        print('  Mapping size: {} ({} metabolic + 3 betaxanthin)'.format(
            len(active_mappings), len(active_mappings) - 3))
    elif args.add_betaxanthin:
        print('  Mapping size: {} (central carbon only; R75/R76/R77 skipped by '
              '--skip-betax-mfa)'.format(len(active_mappings)))
    else:
        print('  Mapping size: {} (central carbon only)'.format(len(active_mappings)))

    print('\nApplying MFA constraints (scale=/{}, tol=+/-{:.0f}%):'.format(
        int(args.mfa_scale), args.mfa_tol * 100))
    if args.only_biomass_mfa:
        print('  --only-biomass-mfa is set: skipping ALL per-reaction MFA bounds '
              '(will apply R1 glucose + R74 biomass only).')
    for mfa_id, me_root in active_mappings.items():
        if mfa_id not in mfa:
            continue
        scaled = mfa[mfa_id] / args.mfa_scale
        flipped = mfa_id in (SIGN_FLIPS | EXPANDED_SIGN_FLIPS)
        if flipped:
            scaled = -scaled
        flip_tag = ' [sign-flipped: MFA fwd = BiGG rev]' if flipped else ''
        # R1 may be represented either as glucose exchange (legacy workbook)
        # or as the PEP-dependent PTS reaction (revised workbook). CLI glucose
        # controls must have identical semantics in both representations.
        if mfa_id == 'R1' and args.free_glucose:
            print('  {:<10s} -> {:<22s} target {:7.3f}  [SKIPPED: --free-glucose]'.format(
                mfa_id, me_root, scaled))
            continue
        if me_root == 'EX_glc__D_e':
            r = me.reactions.EX_glc__D_e
            if args.hard_glucose_uptake:
                r.lower_bound = -scaled
                r.upper_bound = -scaled
                hard_tag = ' [HARD: no deviation]'
            else:
                r.lower_bound = -scaled * (1 + args.mfa_tol)
                r.upper_bound = -scaled * (1 - args.mfa_tol)
                hard_tag = ''
            print('  {:<10s} -> {:<22s} target {:7.3f}  bounds [{:.3f}, {:.3f}]{}{}'.format(
                mfa_id, me_root, scaled, r.lower_bound, r.upper_bound, flip_tag, hard_tag))
        else:
            if args.only_biomass_mfa and mfa_id != 'R1':
                continue   # skip per-reaction bounds (R2-R37, R75-R77)
            effective_tolerance = (0.0 if (mfa_id == 'R1' and
                                            args.hard_glucose_uptake) else
                                   0.25 if (args.mfa_expanded_tiered and
                                             mfa_id in ('R44', 'R48 net'))
                                   else args.mfa_tol)
            drain = add_mfa_constraint(me, me_root, scaled, effective_tolerance)
            if drain is not None:
                if (args.mfa_expanded_tiered and mfa_id in EXPANDED_ABSOLUTE_ZERO
                        and abs(scaled) < 1e-5):
                    drain.lower_bound = 0.0
                    drain.upper_bound = 0.01
                    print('  {:<10s} -> {:<22s} target near zero; absolute bounds [0, 0.01]'.format(
                        mfa_id, me_root))
                print('  {:<10s} -> {:<22s} target {:7.3f}  drain bounds [{:.3f}, {:.3f}]{}'.format(
                    mfa_id, me_root, scaled, drain.lower_bound, drain.upper_bound, flip_tag))

    if args.mfa_expanded_tiered and 'R66 net' in mfa:
        r66_target = mfa['R66 net'] / args.mfa_scale
        r66_drain = add_mfa_linear_constraint(
            me, 'R66_net', {'THD2pp': 1.0, 'NADTRHD': -1.0},
            r66_target, args.mfa_tol)
        if r66_drain is not None:
            print('  R66 net    -> THD2pp - NADTRHD       target {:7.3f}  '
                  'drain bounds [{:.3f}, {:.3f}]'.format(
                      r66_target, r66_drain.lower_bound, r66_drain.upper_bound))

    # Optionally relax UPPER bound on R75/R76/R77 betaxanthin drains so those
    # reactions can carry more flux than the measured MFA value. Lower bound
    # stays as-is (must produce at least (1 - tol) x measured target).
    if args.betax_mfa_lower_only and args.add_betaxanthin:
        print('\nRelaxing UPPER bounds on R75/R76/R77 drains (lower bound only):')
        for me_root in BETAXANTHIN_MAPPINGS.values():
            drain_id = 'mfa_drain_' + me_root
            try:
                drain = me.reactions.get_by_id(drain_id)
                old_ub = drain.upper_bound
                drain.upper_bound = 1000.0
                print('  {:<28s} lb={:.4f}  ub {:.4f} -> {:.1f}'.format(
                    drain_id, drain.lower_bound, old_ub, drain.upper_bound))
            except KeyError:
                print('  WARNING: {} not found (betaxanthin drain missing)'.format(drain_id))

    # Optional MFA biomass constraint (R74). If --biomass-tol > 0, allow a
    # tolerance band [target*(1-tol), target*(1+tol)]; else pin (lb=ub).
    if args.mfa_biomass:
        if 'R74' in mfa:
            target_mu = mfa['R74'] / args.mfa_scale
            try:
                bm = me.reactions.get_by_id('biomass_dilution')
                btol = args.biomass_tol
                if btol > 0:
                    bm.lower_bound = target_mu * (1 - btol)
                    bm.upper_bound = target_mu * (1 + btol)
                    tag = '[+/-{:.1%} band]'.format(btol)
                else:
                    bm.lower_bound = target_mu
                    bm.upper_bound = target_mu
                    tag = '[HARD: pinned]'
                print('  R74        -> biomass_dilution        target {:7.4f}  bounds [{:.4f}, {:.4f}]  {}'.format(
                    target_mu, bm.lower_bound, bm.upper_bound, tag))
            except KeyError:
                print('  WARNING: biomass_dilution reaction not found in model')
        else:
            print('  WARNING: R74 not in MFA sheet, biomass constraint skipped')

# Acetate overflow constraint (independent of MFA reaction constraints).
# Suppressed when --only-biomass-mfa is set, per that flag's contract.
if args.constrain_acetate and not args.only_biomass_mfa:
    print('\nApplying acetate overflow constraint from {}:'.format(args.mfa_sheet))
    constrain_acetate_overflow(me, args.mfa_workbook, args.mfa_sheet,
                                args.mfa_scale, args.mfa_tol,
                                lower_only=args.acetate_lower_only)
elif args.only_biomass_mfa and args.constrain_acetate:
    print('\n--only-biomass-mfa is set: SKIPPING acetate overflow constraint even '
          'though --constrain-acetate was requested.')

# ----------------------------------------------------------------------
# 3b. Apply data-derived exchange bounds from --exchange-json
#     (overrides any MFA-derived bound on the same reaction)
# ----------------------------------------------------------------------
if args.exchange_json:
    import json as _json
    with open(args.exchange_json) as _f:
        _payload = _json.load(_f)
    _exchanges = _payload.get('exchanges', {})
    _skip = set(args.skip_exchange or [])
    _ci = args.exchange_ci_mult
    print('\nApplying --exchange-json bounds (mult={:g} sigma):'.format(_ci))
    print('  file: {}'.format(args.exchange_json))
    if _skip:
        print('  --skip-exchange (dropped): {}'.format(sorted(_skip)))
    for _rid, _spec in _exchanges.items():
        if _rid in _skip:
            print('  {:<15s} SKIPPED (per --skip-exchange)'.format(_rid))
            continue
        try:
            _r = me.reactions.get_by_id(_rid)
        except KeyError:
            print('  {:<15s} MISSING from model - skipped'.format(_rid))
            continue
        _q = float(_spec['q']); _sig = float(_spec.get('sigma', 0.0))
        _mode = _spec.get('mode', 'both')     # 'both' (default) | 'lb_only' | 'ub_only'
        _lb = _q - _ci * _sig
        _ub = _q + _ci * _sig
        if _mode == 'lb_only':
            _r.lower_bound = _lb
            print('  {:<15s} q={:+.4f}  sigma={:.4f}  ->  lb-only [{:+.4f}, model-default]'.format(
                _rid, _q, _sig, _lb))
        elif _mode == 'ub_only':
            _r.upper_bound = _ub
            print('  {:<15s} q={:+.4f}  sigma={:.4f}  ->  ub-only [model-default, {:+.4f}]'.format(
                _rid, _q, _sig, _ub))
        else:
            _r.lower_bound = _lb
            _r.upper_bound = _ub
            print('  {:<15s} q={:+.4f}  sigma={:.4f}  ->  bounds [{:+.4f}, {:+.4f}]'.format(
                _rid, _q, _sig, _lb, _ub))

# Revised MFA R1 is the PEP-dependent PTS step and is scaled to the measured
# glucose uptake. Couple the two fluxes so alternate glucose-entry routes
# cannot separate the MFA anchor from the experimental exchange measurement.
if USE_MFA and (updated_pts_r1 or updated_glk_r1) and not args.free_glucose:
    _r1_entry = 'GLCptspp' if updated_pts_r1 else 'HEX1'
    _r1_label = 'PTS' if updated_pts_r1 else 'HEX1'
    _r1_exchange_coupling = add_mfa_linear_constraint(
        me, 'updated_R1_{}_equals_glucose_uptake'.format(_r1_label),
        {_r1_entry: 1.0, 'EX_glc__D_e': 1.0}, 0.0, 0.0)
    if _r1_exchange_coupling is None:
        raise RuntimeError('Could not create revised R1 {}/glucose-exchange equality'.format(
            _r1_label))
    print('  Revised R1 coupling: {} + EX_glc__D_e = 0'.format(_r1_entry))

if args.close_nh4_exchange:
    _nh4 = me.reactions.get_by_id('EX_nh4_e')
    _nh4.lower_bound = 0.0
    _nh4.upper_bound = 0.0
    print('\n--close-nh4-exchange: EX_nh4_e fixed to [0, 0]')

OVERFLOW_SECRETION_CLOSURES = {
    'EX_etoh_e', 'EX_lac__D_e', 'EX_lac__L_e', 'EX_for_e', 'EX_pyr_e',
    'EX_acald_e', 'EX_12ppd__R_e', 'EX_12ppd__S_e',
}
OVERFLOW_FEASIBILITY_EXCEPTION = {
    'EX_succ_e': ('left open because closing succinate makes the MFA37 '
                  'formamide growth state infeasible')
}
closed_carbon_overflow = {}
if args.close_nonacetate_carbon_overflow:
    for _reaction_id in sorted(OVERFLOW_SECRETION_CLOSURES):
        try:
            _reaction = me.reactions.get_by_id(_reaction_id)
        except KeyError:
            continue
        _old_bounds = (float(_reaction.lower_bound), float(_reaction.upper_bound))
        _reaction.lower_bound = min(_old_bounds[0], 0.0)
        _reaction.upper_bound = 0.0
        closed_carbon_overflow[_reaction.id] = {
            'old_bounds': _old_bounds,
            'new_bounds': (float(_reaction.lower_bound),
                           float(_reaction.upper_bound)),
        }
    print('\n--close-nonacetate-carbon-overflow: closed {} carbon exchanges'.format(
        len(closed_carbon_overflow)))
    print('  closed reactions: {}'.format(sorted(closed_carbon_overflow)))
    print('  feasibility exception: {}'.format(OVERFLOW_FEASIBILITY_EXCEPTION))

OVERFLOW_CARBON_GROUPS = {
    'acetate': {'EX_ac_e': 2.0},
    'lactate': {'EX_lac__D_e': 3.0, 'EX_lac__L_e': 3.0},
    'ethanol': {'EX_etoh_e': 2.0},
    'formate': {'EX_for_e': 1.0},
    'succinate': {'EX_succ_e': 4.0},
    'minor': {
        'EX_pyr_e': 3.0, 'EX_acald_e': 2.0,
        'EX_12ppd__R_e': 3.0, 'EX_12ppd__S_e': 3.0,
    },
}
OVERFLOW_SHARE_LIMITS = {
    'acetate_min': 0.85, 'lactate_max': 0.05, 'ethanol_max': 0.05,
    'formate_max': 0.05, 'succinate_max': 0.05, 'minor_max': 0.05,
}
overflow_percentage_drains = []
if args.overflow_percentage_limits:
    _group_weights = {}
    for _group, _members in OVERFLOW_CARBON_GROUPS.items():
        for _reaction_id, _carbon_count in _members.items():
            try:
                _reaction = me.reactions.get_by_id(_reaction_id)
            except KeyError:
                continue
            _reaction.lower_bound = max(float(_reaction.lower_bound), 0.0)
            _group_weights.setdefault(_group, {})[_reaction_id] = _carbon_count

    _total_weights = {}
    for _members in _group_weights.values():
        for _reaction_id, _weight in _members.items():
            _total_weights[_reaction_id] = _total_weights.get(_reaction_id, 0.0) + _weight

    def _add_overflow_inequality(_name, _weights, _lb, _ub):
        _constraint = Constraint('overflow_C_' + _name)
        me.add_metabolites([_constraint])
        for _reaction_id, _weight in _weights.items():
            me.reactions.get_by_id(_reaction_id).add_metabolites(
                {_constraint: _weight}, combine=False)
        _drain = SummaryVariable('overflow_drain_' + _name)
        me.add_reactions([_drain])
        _drain.add_metabolites({_constraint: -1.0}, combine=False)
        _drain.lower_bound, _drain.upper_bound = _lb, _ub
        overflow_percentage_drains.append(_drain.id)

    _acetate_weights = dict(_group_weights.get('acetate', {}))
    for _reaction_id, _weight in _total_weights.items():
        _acetate_weights[_reaction_id] = (_acetate_weights.get(_reaction_id, 0.0)
                                           - 0.85 * _weight)
    _add_overflow_inequality('acetate_share_min_85pct', _acetate_weights, 0.0, 1e6)

    for _group in ('lactate', 'ethanol', 'formate', 'succinate', 'minor'):
        _weights = dict(_group_weights.get(_group, {}))
        for _reaction_id, _weight in _total_weights.items():
            _weights[_reaction_id] = _weights.get(_reaction_id, 0.0) - 0.05 * _weight
        _add_overflow_inequality(_group + '_share_max_5pct', _weights, -1e6, 0.0)
    print('\n--overflow-percentage-limits: acetate >=85%; each non-acetate group <=5%')
    print('  constrained drains: {}'.format(overflow_percentage_drains))

explicitly_closed_exchanges = {}
for _reaction_id in args.close_exchange:
    try:
        _reaction = me.reactions.get_by_id(_reaction_id)
    except KeyError:
        raise ValueError('--close-exchange reaction not found: {}'.format(_reaction_id))
    _old_bounds = (float(_reaction.lower_bound), float(_reaction.upper_bound))
    _reaction.lower_bound = 0.0
    _reaction.upper_bound = 0.0
    explicitly_closed_exchanges[_reaction_id] = {
        'old_bounds': _old_bounds, 'new_bounds': (0.0, 0.0)}
if explicitly_closed_exchanges:
    print('\n--close-exchange: {}'.format(sorted(explicitly_closed_exchanges)))

# Apply general reaction overrides last so a sensitivity test changes only the
# nominated reaction while preserving the configured MFA and exchange design.
reaction_bound_overrides = {}
for _bound_spec in args.reaction_bound:
    _parts = [part.strip() for part in _bound_spec.split(',')]
    if len(_parts) != 3:
        raise ValueError('--reaction-bound must be REACTION,LB,UB: {}'.format(
            _bound_spec))
    _reaction_id, _lb_text, _ub_text = _parts
    if _reaction_id not in me.reactions:
        raise ValueError('--reaction-bound reaction not found: {}'.format(
            _reaction_id))
    _reaction = me.reactions.get_by_id(_reaction_id)
    _old_bounds = (float(_reaction.lower_bound), float(_reaction.upper_bound))
    _new_bounds = (float(_lb_text), float(_ub_text))
    if _new_bounds[0] > _new_bounds[1]:
        raise ValueError('--reaction-bound lower bound exceeds upper bound: {}'.format(
            _bound_spec))
    _reaction.lower_bound, _reaction.upper_bound = _new_bounds
    reaction_bound_overrides[_reaction_id] = {
        'old_bounds': _old_bounds, 'new_bounds': _new_bounds}
if reaction_bound_overrides:
    print('\n--reaction-bound overrides: {}'.format(reaction_bound_overrides))

reaction_direction_blocks = {}
for _block_spec in args.block_reaction_direction:
    _parts = [part.strip() for part in _block_spec.split(',')]
    if len(_parts) != 2 or _parts[1].lower() not in ('forward', 'reverse', 'all'):
        raise ValueError('--block-reaction-direction must be ROOT,forward|reverse|all: '
                         '{}'.format(_block_spec))
    _root, _direction = _parts[0], _parts[1].lower()
    _prefixes = []
    if _direction in ('forward', 'all'):
        _prefixes.append(_root + '_FWD_')
    if _direction in ('reverse', 'all'):
        _prefixes.append(_root + '_REV_')
    _variants = [reaction for reaction in me.reactions
                 if any(reaction.id.startswith(prefix) for prefix in _prefixes)]
    # Some custom/plain reactions do not have enzyme-coupled variants.
    if not _variants and _direction == 'all' and _root in me.reactions:
        _variants = [me.reactions.get_by_id(_root)]
    if not _variants:
        raise ValueError('--block-reaction-direction found no variants for {} ({})'.format(
            _root, _direction))
    reaction_direction_blocks[_block_spec] = []
    for _reaction in _variants:
        _old_bounds = (float(_reaction.lower_bound), float(_reaction.upper_bound))
        _reaction.lower_bound = 0.0
        _reaction.upper_bound = 0.0
        reaction_direction_blocks[_block_spec].append({
            'reaction_id': _reaction.id,
            'old_bounds': _old_bounds,
            'new_bounds': (0.0, 0.0),
        })
if reaction_direction_blocks:
    print('\n--block-reaction-direction:')
    for _spec, _variants in sorted(reaction_direction_blocks.items()):
        print('  {} -> blocked {} variants'.format(_spec, len(_variants)))

# ----------------------------------------------------------------------
# 3c. Pin biomass_dilution to a specific mu (independent of MFA logic).
#     --pin-growth takes precedence over --target-growth if both are set.
# ----------------------------------------------------------------------
if args.pin_growth is not None:
    _pmu = float(args.pin_growth)
    _bm = me.reactions.get_by_id('biomass_dilution')
    _bm.lower_bound = _pmu
    _bm.upper_bound = _pmu
    print('\n--pin-growth {:.4f} h^-1: biomass_dilution bounds pinned to [{:.4f}, {:.4f}]'.format(
        _pmu, _pmu, _pmu))

if args.save_configured_model or args.save_configured_model_only:
    _configured_path = os.path.join(OUT, 'configured_model.pickle')
    with open(_configured_path, 'wb') as _configured_handle:
        pickle.dump(me, _configured_handle, protocol=pickle.HIGHEST_PROTOCOL)
    print('\nSaved fully configured model to {}'.format(_configured_path))

if args.save_configured_model_only:
    print('Exiting before solver rebuild (--save-configured-model-only).')
    sys.exit(0)

# ----------------------------------------------------------------------
# 4. Rebuild solver, re-apply ROS+metal substitutions, solve
# ----------------------------------------------------------------------
print('\nRebuilding solver after model modifications...')
solver = ME_NLP1(me)
# Point existing stress object at the new solver so force_o2s_from_pq()
# reads the right substitution_dict
stress.me_nlp = solver
stress.substitute_ros(solver, subs_dict)
stress.substitute_metal(solver)
stress.force_o2s_from_pq()

if args.fixed_growth_feasibility:
    if args.pin_growth is None:
        print('*** ERROR: --fixed-growth-feasibility requires --pin-growth')
        sys.exit(1)
    for r in me.reactions:
        r.objective_coefficient = 0
    target_mu = float(args.pin_growth)
    _basis = np.load(args.basis_file) if args.basis_file else None
    print('\nSolving direct fixed-growth feasibility LP at mu = {:.6f} '
          '(quad-precision qMINOS)...'.format(target_mu))
    xopt, status, hs = solver.solvelp(
        muf=target_mu, basis=_basis, quad=True, verbosity=1)
    if status != 'optimal':
        print('*** FIXED-GROWTH LP STATUS: {} ***'.format(status))
        with open(os.path.join(OUT, 'INFEASIBLE.txt'), 'w') as f:
            f.write('Fixed-growth qMINOS solve failed at mu={} with status={}\n'.format(
                target_mu, status))
        sys.exit(1)
    from types import SimpleNamespace
    rxn_ids = [r.id for r in me.reactions]
    x_dict = dict(zip(rxn_ids, xopt))
    me.solution = SimpleNamespace(
        x=list(xopt), x_dict=x_dict, status=str(status), f=0.0)
    if args.save_basis:
        np.save(os.path.join(OUT, 'qminos_basis.npy'), hs)
    print('  fixed-growth feasibility status: {}'.format(status))
    print('  symbolic mu: {:.6f}'.format(target_mu))
    print('  biomass_dilution: {:.6f}'.format(
        x_dict.get('biomass_dilution', 0.0)))
elif args.maximize_product:
    # LP mode: pin mu via --mfa-biomass or --pin-growth, set product reaction as
    # objective, solve LP with quad-precision qMINOS at that fixed mu.
    if not (args.mfa_biomass or args.pin_growth is not None):
        print('*** ERROR: --maximize-product requires biomass to be pinned via '
              '--mfa-biomass or --pin-growth')
        sys.exit(1)
    prod_id = args.maximize_product
    if prod_id not in me.reactions:
        print('*** ERROR: --maximize-product reaction "{}" not in model'.format(prod_id))
        sys.exit(1)
    # Clear all objective coefficients, set this one to 1 (maximize)
    for r in me.reactions:
        r.objective_coefficient = 0
    me.reactions.get_by_id(prod_id).objective_coefficient = 1

    # Pull the pinned mu from biomass_dilution bounds
    bm = me.reactions.get_by_id('biomass_dilution')
    target_mu = float(bm.lower_bound)  # lower==upper when biomass-tol==0
    print('\nSolving LP: max {} at fixed mu = {:.4f} (quad-precision qMINOS)...'.format(
        prod_id, target_mu))

    xopt, status, hs = solver.solvelp(muf=target_mu, quad=True, verbosity=1)
    if status != 'optimal':
        print('*** LP status: {} (not optimal)'.format(status))
        # still keep going -- LP may report suboptimal but yield usable values

    # Attach solution to me.solution using SimpleNamespace (picklable, provides
    # the .x / .x_dict / .status / .f attributes the rest of this script reads).
    from types import SimpleNamespace
    rxn_ids = [r.id for r in me.reactions]
    x_dict = dict(zip(rxn_ids, xopt))
    me.solution = SimpleNamespace(
        x=list(xopt),
        x_dict=x_dict,
        status=str(status),
        f=x_dict.get(prod_id, 0.0),
    )
    print('  product flux ({}): {:.6f}'.format(prod_id, x_dict.get(prod_id, 0)))
    print('  mu (biomass_dilution) : {:.6f}'.format(x_dict.get('biomass_dilution', 0)))
else:
    _basis = None
    if args.basis_file:
        _basis = np.load(args.basis_file)
        print('Solving via bisection (precision {:.0e}, warm-start basis {})...'.format(
            args.bisect_precision, args.basis_file))
    else:
        print('Solving via bisection (precision {:.0e}, no warm-start)...'.format(
            args.bisect_precision))
    sol, hs, xopt, cache = solver.bisectmu(
        precision=args.bisect_precision, basis=_basis)
    if args.save_basis:
        np.save(os.path.join(OUT, 'qminos_basis.npy'), hs)
        print('Saved qMINOS basis to {}'.format(
            os.path.join(OUT, 'qminos_basis.npy')))

if me.solution is None or me.solution.x_dict is None:
    print('\n*** NO FEASIBLE SOLUTION ***')
    with open(os.path.join(OUT, 'INFEASIBLE.txt'), 'w') as f:
        f.write('No feasible solution found for ROS={}x, MFA={}\n'.format(
            ROS, 'on' if USE_MFA else 'off'))
    sys.exit(1)

mu = me.solution.x_dict['biomass_dilution']
gur = -1 * me.solution.x_dict.get('EX_glc__D_e', 0)
our = -1 * me.solution.x_dict.get('EX_o2_e', 0)
apr = me.solution.x_dict.get('EX_ac_e', 0)
co2 = me.solution.x_dict.get('EX_co2_e', 0)
lac = me.solution.x_dict.get('EX_lac__D_e', 0)

print('\n=== SOLUTION ===')
print('  mu:   {:.4f} h^-1'.format(mu))
print('  GUR:  {:.4f}'.format(gur))
print('  OUR:  {:.4f}'.format(our))
print('  APR:  {:.4f}'.format(apr))

# ----------------------------------------------------------------------
# 5. Save outputs (standardized filenames)
# ----------------------------------------------------------------------
print('\nWriting outputs to {}'.format(OUT))

# solution.pickle
with open(os.path.join(OUT, 'solution.pickle'), 'wb') as f:
    pickle.dump(me.solution, f)

# phenotype.csv
with open(os.path.join(OUT, 'phenotype.csv'), 'w') as f:
    w = csv.writer(f)
    w.writerow(['metric', 'value', 'units'])
    w.writerow(['ros_multiplier', ROS, 'x basal'])
    w.writerow(['mfa_constrained', 'yes' if USE_MFA else 'no', ''])
    w.writerow(['growth_rate', mu, 'h^-1'])
    w.writerow(['GUR', gur, 'mmol/gDW/h'])
    w.writerow(['OUR', our, 'mmol/gDW/h'])
    w.writerow(['APR_acetate', apr, 'mmol/gDW/h'])
    w.writerow(['LPR_lactate', lac, 'mmol/gDW/h'])
    w.writerow(['CO2_secretion', co2, 'mmol/gDW/h'])
    if gur > 1e-12:
        w.writerow(['biomass_yield', mu / gur / 0.1800634, 'gDCW/g-glucose'])

# fluxes.csv (M-model-equivalent net fluxes, non-zero only)
# Columns: reaction, name, equation, flux

def _format_equation_from_stoich(stoich_dict):
    reactants, products = [], []
    for met_id, coeff in stoich_dict.items():
        if not isinstance(coeff, (int, float)):
            continue
        c = abs(coeff)
        token = '{} {}'.format(c, met_id) if c != 1 else met_id
        if coeff < 0:
            reactants.append(token)
        elif coeff > 0:
            products.append(token)
    arrow = ' <=> ' if (reactants and products) else ' -> '
    return ' + '.join(reactants) + arrow + ' + '.join(products)


def _get_name_and_equation(me, rxn_id):
    from cobrame import StoichiometricData
    name = ''
    equation = ''
    # Prefer clean M-model stoichiometry from process_data
    try:
        pd_obj = me.process_data.get_by_id(rxn_id)
        if isinstance(pd_obj, StoichiometricData):
            equation = _format_equation_from_stoich(pd_obj._stoichiometry)
    except KeyError:
        pass
    # Get name from reaction instance if available
    try:
        r = me.reactions.get_by_id(rxn_id)
        name = getattr(r, 'name', None) or ''
        if not equation:
            stoich = {m.id: c for m, c in r.metabolites.items()}
            equation = _format_equation_from_stoich(stoich)
    except KeyError:
        # Fall back to a representative variant for name lookup
        variants = [r for r in me.reactions if r.id.startswith(rxn_id + '_FWD_')]
        if variants and not name:
            name = getattr(variants[0], 'name', None) or ''
    return name, equation


try:
    met_flux = me.get_metabolic_flux(me.solution)
    # COBRAme omits custom MEReaction objects without process_data from
    # get_metabolic_flux(). Export their exact solved values explicitly so
    # whole-network analyses do not silently lose the heterologous pathway.
    for custom_id in ('BETAX_TYR_DOPA', 'BETAX_DOPA_BTLA',
                      'BETAX_BTLA_TYBX', 'TYBXtipp', 'TYBXtex',
                      'FMDA', 'FMDtex', 'FMDtpp'):
        if custom_id in me.solution.x_dict:
            met_flux[custom_id] = float(me.solution.x_dict[custom_id])
    with open(os.path.join(OUT, 'fluxes.csv'), 'w') as f:
        w = csv.writer(f)
        w.writerow(['reaction', 'name', 'equation', 'flux'])
        for rxn, flux in sorted(met_flux.items()):
            if abs(flux) > 1e-20 and 'DM_RNA' not in rxn:
                name, equation = _get_name_and_equation(me, rxn)
                w.writerow([rxn, name, equation, flux])
except Exception as e:
    print('  WARN: could not compute metabolic fluxes:', e)

# proteome.csv (protein mass fractions)
MW_dict = {}
for gene in me.metabolites.query(re.compile('RNA_b[0-9]')):
    prot_id = gene.id.replace('RNA', 'protein')
    try:
        m = me.metabolites.get_by_id(prot_id)
        MW_dict[prot_id.replace('protein_', '')] = m.formula_weight
    except KeyError:
        continue
MW_dict['dummy'] = 33210.96
total_protein = 0.0
mf = {}
for rid, flux in me.solution.x_dict.items():
    if rid.startswith('translation_') and flux > 1e-20:
        gene = rid.replace('translation_', '')
        if gene in MW_dict:
            mass = flux * MW_dict[gene]
            total_protein += mass
            mf[gene] = mass
with open(os.path.join(OUT, 'proteome.csv'), 'w') as f:
    w = csv.writer(f)
    w.writerow(['gene', 'mass_fraction'])
    if total_protein > 0:
        for g in sorted(mf.keys(), key=lambda g: -mf[g]):
            w.writerow([g, mf[g] / total_protein])

# mfa_comparison.csv (only for constrained runs)
# Reports central-carbon mappings + acetate (R70) and biomass (R74) when
# those constraints are active. R70/R74 are applied via separate code paths
# (constrain_acetate_overflow / --mfa-biomass) so they aren't in
# active_mappings; we add them here for reporting parity.
if USE_MFA:
    # Build a reporting-only mapping by adding R70 / R74 (and ensuring R75-R77
    # already in active_mappings stay) for the comparison file.
    reporting_mappings = dict(active_mappings)
    if args.constrain_acetate:
        reporting_mappings['R70'] = 'EX_ac_e'
    if args.mfa_biomass:
        reporting_mappings['R74'] = 'biomass_dilution'
    if args.mfa_expanded_tiered and 'R66 net' in mfa:
        reporting_mappings['R66 net'] = 'THD2pp - NADTRHD'

    with open(os.path.join(OUT, 'mfa_comparison.csv'), 'w') as f:
        w = csv.DictWriter(f, fieldnames=['mfa_id', 'me_root', 'mfa_raw',
                                          'mfa_scaled', 'me_predicted', 'diff_pct'])
        w.writeheader()
        for mfa_id, me_root in reporting_mappings.items():
            if mfa_id not in mfa:
                continue
            raw = mfa[mfa_id]
            target = raw / args.mfa_scale
            if mfa_id in (SIGN_FLIPS | EXPANDED_SIGN_FLIPS):
                target = -target  # compare in BiGG-forward convention
            if me_root == 'EX_glc__D_e':
                predicted = -1 * me.solution.x_dict.get('EX_glc__D_e', 0)
            elif me_root == 'THD2pp - NADTRHD':
                thd2pp = sum(me.solution.x_dict.get(r.id, 0)
                             for r in me.reactions
                             if r.id == 'THD2pp' or r.id.startswith('THD2pp_'))
                nadtrhd = sum(me.solution.x_dict.get(r.id, 0)
                              for r in me.reactions
                              if r.id == 'NADTRHD' or r.id.startswith('NADTRHD_'))
                predicted = thd2pp - nadtrhd
            elif me_root in ('EX_ac_e', 'biomass_dilution'):
                # Plain reactions: no FWD/REV split, single id lookup
                predicted = me.solution.x_dict.get(me_root, 0)
            else:
                fwd_flux = sum(me.solution.x_dict.get(r.id, 0)
                               for r in me.reactions if r.id.startswith(me_root + '_FWD_'))
                rev_flux = sum(me.solution.x_dict.get(r.id, 0)
                               for r in me.reactions if r.id.startswith(me_root + '_REV_'))
                predicted = fwd_flux - rev_flux
                # Plain MEReaction case (no FWD/REV variants)
                if predicted == 0 and me_root in me.reactions:
                    predicted = me.solution.x_dict.get(me_root, 0)
            diff = (predicted - target) / target * 100 if abs(target) > 1e-12 else 0
            w.writerow({'mfa_id': mfa_id, 'me_root': me_root,
                        'mfa_raw': raw, 'mfa_scaled': target,
                        'me_predicted': predicted, 'diff_pct': diff})

# config.json — full traceability of run parameters
import json
config = {
    'sim_number': args.sim_number,
    'strain_label': args.strain_label,
    'medium_label': args.medium_label,
    'temperature_C': TEMPERATURE_C,
    'pH': pH,
    'ros_multiplier': ROS,
    'n_source': args.n_source,
    'fmd_uptake': args.fmd_uptake,
    'add_formamide_direct': args.add_formamide_direct,
    'add_formamide_twostep': args.add_formamide_twostep,
    'add_betaxanthin': args.add_betaxanthin,
    'mfa_constrained': USE_MFA,
    'mfa_sheet': args.mfa_sheet if USE_MFA else None,
    'mfa_workbook': args.mfa_workbook if USE_MFA else None,
    'mfa_scale': args.mfa_scale,
    'mfa_tol': args.mfa_tol,
    'hard_glucose_uptake': args.hard_glucose_uptake,
    'free_glucose': args.free_glucose,
    'target_growth': args.target_growth,
    'exchange_json': args.exchange_json,
    'exchange_ci_mult': args.exchange_ci_mult,
    'pin_growth': args.pin_growth,
    'skip_exchange': args.skip_exchange,
    'close_nh4_exchange': args.close_nh4_exchange,
    'close_exchange': args.close_exchange,
    'explicitly_closed_exchanges': explicitly_closed_exchanges,
    'reaction_bound': args.reaction_bound,
    'reaction_bound_overrides': reaction_bound_overrides,
    'block_reaction_direction': args.block_reaction_direction,
    'reaction_direction_blocks': reaction_direction_blocks,
    'close_nonacetate_carbon_overflow': args.close_nonacetate_carbon_overflow,
    'overflow_percentage_limits': args.overflow_percentage_limits,
    'overflow_carbon_groups': OVERFLOW_CARBON_GROUPS,
    'overflow_share_limits': OVERFLOW_SHARE_LIMITS,
    'overflow_percentage_constraint_drains': overflow_percentage_drains,
    'overflow_secretion_closure_candidates': sorted(OVERFLOW_SECRETION_CLOSURES),
    'overflow_secretion_closures': sorted(closed_carbon_overflow),
    'overflow_feasibility_exception': OVERFLOW_FEASIBILITY_EXCEPTION,
    'closed_carbon_overflow_exchanges': closed_carbon_overflow,
    'save_configured_model': args.save_configured_model,
    'mfa_biomass': args.mfa_biomass,
    'biomass_tol': args.biomass_tol,
    'bisect_precision': args.bisect_precision,
    'constrain_acetate': args.constrain_acetate,
    'solution': {
        'growth_rate': mu,
        'GUR': gur,
        'OUR': our,
        'APR_acetate': apr,
        'CO2_secretion': co2,
        'biomass_yield': mu / gur / 0.1800634 if gur > 1e-12 else None,
    },
}
with open(os.path.join(OUT, 'config.json'), 'w') as f:
    json.dump(config, f, indent=2)

print('\nDone.')
