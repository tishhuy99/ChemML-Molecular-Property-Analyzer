import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

st.set_page_config(
    page_title="Molecule Analyzer",
    page_icon="🧪"
)

st.title("🧪 Molecule Analyzer")

smiles = st.text_input(
    "Enter a SMILES string:",
    "CC(=O)OC1=CC=CC=C1C(=O)O"
)

mol = Chem.MolFromSmiles(smiles)

if mol:
    st.image(Draw.MolToImage(mol))

    st.subheader("Molecular Properties")

    st.metric("Molecular Weight", f"{Descriptors.MolWt(mol):.2f}")
    st.metric("LogP", f"{Descriptors.MolLogP(mol):.2f}")

    st.write("H-Bond Donors:", Descriptors.NumHDonors(mol))
    st.write("H-Bond Acceptors:", Descriptors.NumHAcceptors(mol))

else:
    st.error("Invalid SMILES string")
    