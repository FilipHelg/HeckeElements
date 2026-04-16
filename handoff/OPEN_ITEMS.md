# Open Items and Suggested Next Actions

## 1) Resolve output convention explicitly
Problem:
- Current emitted labels show asymmetric unit behavior (`e*y=y`, `x*e=x^{-1}`).
- This causes confusion for inversion and associativity checks on raw CSV labels.

Action options:
1. Patch generator output mapping to standard labeling.
2. Patch validator to validate current convention exactly.
3. Offer both mappings and add a `--convention` argument.

## 2) Add strict regression tests
Suggested tests to add:
- Unit tests for identity behavior for all elements in S4/S5.
- Associativity exhaustive test on S4 for chosen convention.
- Spot-check known products with multi-term expansions.

## 3) Document data schema for CSV consumers
Add a short schema note near outputs or README:
- meaning of `x`, `y`, `z`
- polynomial format in `h`
- convention used for multiplication labels

## 4) Optional publication prep
Before sharing externally:
- regenerate S4/S5/S6 with final chosen convention
- re-run validator with strict associativity mode under matching semantics
- freeze tool versions/flags in a release note
