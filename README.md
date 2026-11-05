# 🧪 ChemML — Your little molecule lab

A friendly Streamlit workbench for exploring molecular structures with RDKit.

## What you can do

- Explore a molecule from SMILES or choose aspirin, caffeine, ethanol, or ibuprofen.
- See its 2D structure, formula, canonical SMILES, and eight molecular descriptors.
- Review a transparent Lipinski rule-of-five violation count.
- Compare two molecules with side-by-side structures, properties, and Morgan/Tanimoto similarity.
- Analyze a UTF-8 CSV of up to 500 molecules (2 MB), preserving invalid rows and input row numbers.
- Download structures as PNG and results as CSV.
- Use native keyboard-accessible controls, visible focus indicators, responsive layouts, and an in-app quick guide.

## Open the app online

[Launch ChemML in Binder](https://mybinder.org/v2/gh/tishhuy99/ChemML-Molecular-Property-Analyzer/main?urlpath=proxy%2F8501%2F) — a public, no-login launch of this Python app. The first start can take several minutes while Binder builds the environment. Sessions are temporary and may be unavailable when the free service is busy; work and uploads are not saved between sessions. For a reliable always-on link, a dedicated app host is needed.

## Run locally

Use Python 3.11 or newer. From this project directory:

```bash
python -m venv .venv
source .venv/bin/activate
# Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

Open the localhost URL printed by Streamlit. The `.streamlit/config.toml` theme is loaded when starting from this directory.

## Tests

```bash
python -m pip install -r requirements-dev.txt
python -m pytest -q
```

Tests cover known descriptors, invalid/empty inputs, equivalent SMILES similarity, mixed-validity batches, row limits, safe CSV exports, and the Streamlit example/validation flow.

## Files

- `app.py`: interactive UI and downloads
- `chemistry.py`: shared RDKit calculations, validation, and export
- `style.css` and `.streamlit/config.toml`: lavender theme
- `tests/`: calculation and app regression checks
- `molecular_analysis.ipynb`: original exploratory notebook

## Scientific scope

Descriptors are calculated from the supplied structure; they are not experimental measurements or trained ML predictions. No automatic salt removal, protonation adjustment, or tautomer standardization is performed. The rule-of-five count checks molecular weight >500, Crippen LogP >5, H-bond donors >5, and H-bond acceptors >10. It does not establish safety, efficacy, or absorption. Similarity uses radius-2, 2,048-bit Morgan fingerprints and the Tanimoto coefficient; matching fingerprints do not prove identical biological behavior.

Calculations run on the hosting server. This app does not call an external chemistry API or persist uploaded files. On shared deployments, users must be allowed to send their data to that host. CSV exports neutralize leading spreadsheet-formula characters in text values.

## Author

**Tish Huynh** — M.S. Chemistry, University of Connecticut; B.A. Chemistry, Drexel University.

Future work: curated training data, evaluated property-prediction models, and applicability-domain analysis.

## Proton map and acidity references

Explore now lists every explicit or implicit hydrogen, with one-based atom IDs matching a labeled, downloadable structure. Choose a hydrogen to highlight its bond site and download all proton rows as CSV. For structures above 120 atoms including hydrogens, the diagram is omitted for readability; the full table remains available.

`pka.py` matches conservative functional-group rules to the [Indiana University teaching chart](https://www.chem.indiana.edu/wp-content/uploads/2018/03/pka-chart.pdf). **These are family references, not molecule-specific pKa predictions, experimentally measured values for the submitted compound, or uncertainty intervals.** The chart does not specify a uniform solvent/temperature scale. Do not calculate equilibrium constants from these mixed teaching references.

Supported families include carboxylic acid, phenol, alcohol, aliphatic thiol, amide, selected carbonyl-alpha C–H, terminal alkyne, simple neutral amine N–H, and limited hydrocarbon environments. Unsupported hydrogens remain visible as “Not estimated.” Charged/radical structures, isotope-specific acidity, and multiply activated sites deliberately abstain. Substituent and long-range electronic effects are not modeled; a matching functional group does not validate a compound-specific estimate. No tautomer enumeration, microstate coupling, stereotopic hydrogen resolution, or sequential deprotonation is performed. Hydrogens on one atom share a first-deprotonation family reference.

Simple aliphatic amine proton acceptors have a separate conjugate-acid **pKaH** table; this is different from the acidity of an existing N–H bond. Other basic sites are outside the method. Atom IDs refer to the supplied SMILES after adding hydrogens and can change with alternative equivalent SMILES. CSV rows include method and source provenance.

Acidity alone does not predict nucleophilicity, electrophilicity, radical reactivity, reaction rates or products. Assessing a reaction requires a partner, solvent, temperature and intended transformation. Use measured solvent-matched data or a validated pKa predictor for quantitative decisions.
