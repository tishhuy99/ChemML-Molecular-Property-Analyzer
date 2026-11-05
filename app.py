"""ChemML: a friendly, local molecular workbench."""
import io
from pathlib import Path
import pandas as pd
import streamlit as st
from rdkit.Chem import Draw
from pka import site_analysis, site_image, SOURCE
from chemistry import EXAMPLES, MAX_ROWS, analyze_batch, describe, export_csv, molecule, similarity

st.set_page_config(page_title="ChemML · Your little molecule lab", page_icon="🧪", layout="wide")
st.markdown(f"<style>{Path(__file__).with_name('style.css').read_text()}</style>", unsafe_allow_html=True)
st.markdown('<div class="lab-eyebrow">🧪 CHEMML / MOLECULAR PLAYGROUND</div>', unsafe_allow_html=True)
st.title("Small molecules. Big discoveries.")
st.caption("Your little molecule lab — explore structures, compare properties, and follow your curiosity.")
explore, compare, batch, guide = st.tabs(["🌷 Explore a molecule", "🔎 Compare two", "📂 Batch analysis", "📖 Quick guide"])

def png(smiles):
    buffer = io.BytesIO()
    Draw.MolToImage(molecule(smiles), size=(600, 380)).save(buffer, format="PNG")
    return buffer.getvalue()

def set_example(name):
    st.session_state["smiles"] = EXAMPLES[name]

def details(values):
    keys = ["Molecular weight (g/mol)", "LogP (Crippen)", "H-bond donors", "H-bond acceptors"]
    for pair in (keys[:2], keys[2:]):
        for col, key in zip(st.columns(2), pair):
            col.metric(key, f"{values[key]:.2f}" if isinstance(values[key], float) else values[key])
    st.caption("LogP estimates preference for oil over water; larger values indicate greater lipophilicity.")
    with st.expander("More descriptors & rule of five"):
        st.write({key: values[key] for key in ["TPSA (Å²)", "Rotatable bonds", "Rings", "Rule-of-five violations"]})
        st.caption("Counts thresholds exceeded: molecular weight >500, LogP >5, donors >5, acceptors >10. This is a screening heuristic, not a prediction of safety, efficacy, or absorption.")

with explore:
    st.subheader("Meet your next molecule")
    st.write("Start with a familiar favorite, or bring your own SMILES.")
    for col, name in zip(st.columns(4), EXAMPLES):
        col.button(name, on_click=set_example, args=(name,), use_container_width=True)
    if "smiles" not in st.session_state:
        st.session_state.smiles = EXAMPLES["Aspirin"]
    smiles = st.text_input("Molecule SMILES", key="smiles", help="SMILES describes a molecular structure as text. For example, CCO is ethanol.", max_chars=2000)
    try:
        values = describe(smiles)
    except ValueError as exc:
        st.info(str(exc)) if not smiles.strip() else st.error(str(exc))
    else:
        structure, properties = st.columns([1, 1])
        with structure, st.container(border=True):
            st.subheader("A closer look")
            picture = png(smiles)
            st.image(picture, caption=f"2D structure · {values['Formula']}", use_container_width=True)
            st.download_button("↓ Save structure", picture, "molecule.png", "image/png")
        with properties:
            st.subheader("The molecule at a glance")
            details(values)
            st.download_button("↓ Download properties", export_csv([values]), "molecule-properties.csv", "text/csv")
        with st.expander("Canonical SMILES · copy or reuse"):
            st.code(values["Canonical SMILES"], language=None)

        st.subheader("🧭 Proton map & acidity")
        st.info("Educational pKa references — not molecule-specific predictions. Each hydrogen is listed; a number is shown only where a supported functional-group rule applies.")
        st.caption("These teaching-table values do not specify a consistent solvent or temperature. Do not use them to calculate reaction equilibria. Substitution, tautomerism and coupled sites can change acidity substantially.")
        hydrogenated, proton_rows, basic_rows = site_analysis(smiles)
        if proton_rows:
            supported = sum(row["Reference pKa"] != "Not estimated" for row in proton_rows)
            st.caption(f"{len(proton_rows)} hydrogens · {supported} have a family reference · {len(proton_rows) - supported} not estimated")
            options = [row["Hydrogen"] for row in proton_rows]
            labels = {row["Hydrogen"]: f"{row['Hydrogen']} on {row['Attached atom']} · {row['Environment']} · {row['Reference pKa']}" for row in proton_rows}
            selected = st.selectbox("Highlight a hydrogen", options, format_func=lambda value: labels[value])
            diagram = site_image(smiles, selected)
            if diagram is not None:
                st.image(diagram, caption="Atom labels match the table for this input SMILES. Selected hydrogen and its attached atom are highlighted.", use_container_width=True)
                image_buffer = io.BytesIO()
                diagram.save(image_buffer, format="PNG")
                st.download_button("↓ Save proton map", image_buffer.getvalue(), "proton-map.png", "image/png")
            else:
                st.info("The labeled structure is limited to 120 atoms for readability. Every hydrogen is still included in the table and CSV.")
            chosen = next(row for row in proton_rows if row["Hydrogen"] == selected)
            st.write(f"**{selected}: {chosen['Environment']} · {chosen['Reference pKa']}**")
            st.caption(chosen["Interpretation"])
            st.dataframe(pd.DataFrame(proton_rows)[["Hydrogen", "Attached atom", "Environment", "Reference pKa", "Interpretation"]], hide_index=True, use_container_width=True)
            st.download_button("↓ Download proton references", export_csv(proton_rows), "proton-pka-references.csv", "text/csv")
        else:
            st.info("This supplied structure has no hydrogens to map.")
        with st.expander("Proton-accepting sites · pKaH is different"):
            st.write("A base accepts a proton; pKaH describes loss of that proton from its conjugate acid. It is not the acidity of an existing N–H bond.")
            if basic_rows:
                st.dataframe(pd.DataFrame(basic_rows), hide_index=True, use_container_width=True)
            else:
                st.caption("No supported simple aliphatic amine sites found. This does not mean the molecule has no basic sites; aromatic nitrogens and other basic groups are outside this method.")
        with st.expander("How to interpret the proton map"):
            st.markdown("""Lower pKa means easier **heterolytic proton loss**, when values refer to the same solvent and conditions. Hydrogens attached to one atom share a family reference here; these are not independent or sequential dissociation constants. The map does not resolve stereotopic hydrogens or calculate microscopic pKa values.

Acidity alone does **not** identify a reaction product or rank nucleophilicity, electrophilicity, radical reactivity, or reaction rates. To assess a partner molecule, specify the reagent, solvent, temperature and intended transformation.

Numbering follows the supplied SMILES after adding hydrogens. It can change when you enter a different SMILES, even for the same molecule. Charged/radical structures, isotope effects and uncovered environments are deliberately left unestimated.
""")
            st.markdown(f"Reference: [Indiana University functional-group pKa chart]({SOURCE}). Values are illustrative family references, not confidence intervals or validated estimates for your compound.")

with compare:
    st.subheader("A little side-by-side science")
    st.caption("Compare calculated properties and Morgan fingerprint similarity.")
    left, right = st.columns(2)
    a = left.text_input("Molecule A · SMILES", EXAMPLES["Aspirin"], max_chars=2000)
    b = right.text_input("Molecule B · SMILES", EXAMPLES["Ibuprofen"], max_chars=2000)
    try:
        va, vb = describe(a), describe(b)
        score = similarity(a, b)
    except ValueError as exc:
        st.error(str(exc))
    else:
        left.image(png(a), caption=f"Molecule A · {va['Formula']}", use_container_width=True)
        right.image(png(b), caption=f"Molecule B · {vb['Formula']}", use_container_width=True)
        st.metric("Fingerprint similarity", f"{score:.3f}")
        st.caption("Tanimoto coefficient · Morgan radius 2, 2,048 bits. 0 = no shared bits; 1 = identical fingerprints. This does not establish biological equivalence.")
        comparison = pd.DataFrame({"Molecule A": va, "Molecule B": vb}).astype(str)
        st.dataframe(comparison, use_container_width=True)
        st.download_button("↓ Download comparison", export_csv([{"Molecule": "A", **va}, {"Molecule": "B", **vb}]), "comparison.csv", "text/csv")

with batch:
    st.subheader("Bring the whole molecule family")
    st.write(f"Upload a UTF-8 CSV with a SMILES column. Up to {MAX_ROWS} rows and 2 MB; invalid rows stay in your results with an explanation.")
    st.download_button("↓ Example CSV", "name,smiles\nAspirin,CC(=O)OC1=CC=CC=C1C(=O)O\nEthanol,CCO\n", "example-molecules.csv", "text/csv")
    uploaded = st.file_uploader("Choose a molecule CSV", type=["csv"])
    if uploaded is not None:
        if uploaded.size > 2_000_000:
            st.error("This file is too large. Please upload a CSV smaller than 2 MB.")
        else:
            try:
                frame = pd.read_csv(uploaded, dtype=str, keep_default_na=False)
                if frame.empty or len(frame) > MAX_ROWS:
                    raise ValueError(f"Please provide between 1 and {MAX_ROWS} data rows.")
                columns = list(frame.columns)
                preferred = next((i for i, name in enumerate(columns) if name.strip().lower() == "smiles"), 0)
                column = st.selectbox("Which column contains SMILES?", columns, index=preferred)
                if st.button("Analyze molecules", type="primary"):
                    with st.spinner("Getting to know your molecules…"):
                        rows = analyze_batch(frame[column].tolist())
                    success = sum(row["Status"] == "OK" for row in rows)
                    st.success(f"{success} analyzed · {len(rows) - success} need attention")
                    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
                    st.download_button("↓ Download all results", export_csv(rows), "chemml-batch.csv", "text/csv", on_click="ignore")
            except (ValueError, UnicodeDecodeError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
                st.error(f"We couldn’t read this CSV. {exc}")

with guide:
    st.subheader("A friendly place to start")
    st.markdown("""**1. Pick a molecule.** Choose an example in Explore or paste a SMILES string. `CCO` is ethanol.

**2. Explore its properties.** Molecular weight is in g/mol. Donors and acceptors count hydrogen-bonding features. TPSA is topological polar surface area in Å².

**3. Take your work with you.** Save a structure as PNG, download descriptors as CSV, or compare two molecules.

**Working with a dataset?** Open Batch analysis, upload a CSV, select the SMILES column, and analyze. The result retains input row numbers, including invalid rows.

**What is calculated?** RDKit calculates descriptors from the structure you provide. These are not experimental measurements or trained ML predictions. Salts, protonation, and stereochemistry can affect results; this app does not standardize them automatically.

**Where does my data go?** Calculations run on the server hosting this app. The app does not send molecules to an external chemistry API or save uploads to disk. On a shared deployment, use only data you are allowed to send to that host.
""")
st.divider()
st.caption("Made for curious minds · Powered by RDKit · Research & learning")
