# ME-model publication figure captions

## Combined carbon and cofactor figure (Figure 5)

**Figure 5. ME-model carbon allocation and cofactor use under different
nitrogen sources.** **a,** Carbon allocation predicted by StressME for FM01 +
NH4Cl and FM01 + Formamide after applying condition-specific 13C-MFA,
experimental growth, and exchange constraints. Each exchange flux was weighted
by the number of carbon atoms in the exchanged metabolite and expressed as a
percentage of total quantified carbon input. Carbon import is plotted below
zero, whereas quantified carbon export and retained carbon are plotted above
zero. Retained carbon is the arithmetic remainder after subtracting all
quantified carbon-containing exchange exports from total quantified carbon
input. It represents carbon not resolved among the plotted extracellular
products and must not be interpreted as a direct measurement of biomass or as
a separately predicted intracellular carbon pool. **b-c,** Reaction-level ATP
(**b**) and NADPH (**c**) production and consumption calculated from the exact
StressME flux solutions. Contributions are expressed as percentages of the
corresponding production or consumption total, with each side normalized
separately to 100%. The accompanying proteome panels show the summed mass
fraction of unique model genes associated with each reaction. Reaction labels
are model reaction identifiers. ATPS4rpp, ATP synthase (four protons per ATP);
PGK, phosphoglycerate kinase; PYK, pyruvate kinase; ACKr, acetate kinase;
ATPM, ATP maintenance requirement; ADNK1, adenosine kinase; ADK1, adenylate
kinase; HEX1, ATP-dependent hexokinase; NDPK1, nucleoside-diphosphate kinase
(ATP:GDP); ICDHyr, NADP-dependent isocitrate dehydrogenase; THD2pp,
proton-translocating NAD(P) transhydrogenase; ALDD2y, NADP-dependent aldehyde
dehydrogenase; G6PDH2r, glucose-6-phosphate dehydrogenase; GND,
6-phosphogluconate dehydrogenase; MTHFD, NADP-dependent
methylenetetrahydrofolate dehydrogenase; 3OAR1401, 3-oxoacyl-[acyl-carrier
protein] reductase for the saturated C14 intermediate; 3OAR1201,
3-oxoacyl-[acyl-carrier protein] reductase for the saturated C12 intermediate;
3OAR1001, 3-oxoacyl-[acyl-carrier protein] reductase for the saturated C10
intermediate; G3PD2, NADP-dependent glycerol-3-phosphate dehydrogenase;
NADTRHD, soluble NAD transhydrogenase; FADRx2, FAD reductase; and GLUDy,
NADP-dependent glutamate dehydrogenase.

## Carbon balance

**ME-model carbon import, export, and retention.** Carbon-equivalent exchange
fluxes are shown for FM01 + NH4Cl and FM01 + Formamide after replacing the
ammonium state with the specified LDH directional constraints. Fluxes were multiplied by
metabolite carbon number and are presented either in absolute units or as a
percentage of total carbon input. Import is plotted below zero and export or
retained carbon above zero. `Retained carbon` is the accounting remainder after
quantified export and must not be interpreted as a direct biomass measurement.
The complete exchange-reaction breakdown is supplied with the figure data.

## Cofactor turnover and reaction distributions

**ME-model cofactor turnover and reaction-level contributions.** ATP, NADH,
NADPH, and ubiquinol turnover was recomputed from the exact configured models
and qMINOS solution vectors. Turnover is the mean of total production and total
consumption and is normalized to absolute glucose uptake where indicated.
Reaction-distribution panels report each reaction's percentage contribution to
one production or consumption side; complete unfiltered ledgers are provided
in `publication_data/raw_me_cofactor/`.

## Flux and proteome investment

**Cofactor-linked flux distribution and associated proteome investment.** For
each cofactor, the left panel shows signed reaction contributions, with
consumption and production on opposite sides of zero. The right panel reports
the summed mass fraction of unique model genes associated with each reaction.
Separate production-only and consumption-only panels retain reactions meeting
the stated 10% contribution threshold in at least one condition; the CSV files
also retain excluded rows and the inclusion flag.

## CSD-to-MFA distance

**ME-model distance between MFA-constrained growth and non-MFA maximum-product
states.** Panel a reports the root-mean-square difference between paired
glucose-normalized flux vectors over reactions active in both states. Numbers
above bars give the number of shared-active reactions. Panel b shows metric MDS
of whole-network Euclidean distances calculated in one fixed reaction space.
The FM01 + NH4Cl MFA-growth point uses the LDH-constrained solution; WT,
FM01 + Formamide, and historical maximum-product comparator states retain their
previous definitions.
