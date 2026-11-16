# Methods and provenance

The two ME-model states retained the MFA0831-updated intracellular and
experimental exchange constraints, fixed experimental growth rates, and the
default-exchange configuration used in the preceding publication analysis.
For FM01 + NH4Cl, reverse flux through `LDH_D` and forward flux through
`LDH_D2` were set to zero to eliminate the diagnosed internal lactate-
dehydrogenase cycle. The resulting model was solved from the configured ME
model with qMINOS and satisfied the fixed-growth problem with `inform = 0`,
`ninf = 0`, and `sinf = 0`. FM01 + Formamide was retained unchanged.

Cofactor production and consumption were recomputed from the exact configured
ME models and solution vectors. For each reaction, flux was multiplied by the
stoichiometric coefficient of ATP, NADH, NADPH, or ubiquinol. Positive and
negative contributions were summed separately; turnover was their arithmetic
mean and was divided by absolute glucose uptake where shown. Reaction-linked
proteome investment was calculated by summing the solved mass fractions of
the unique genes associated with each reaction. Carbon import and export were
calculated from exchange flux multiplied by metabolite carbon number, with a
complete reaction-level breakdown retained in the publication data.

All figures in this package were regenerated using the FM01 + NH4Cl state with
the stated LDH directional constraints. No values were copied from the earlier
ammonium state.
