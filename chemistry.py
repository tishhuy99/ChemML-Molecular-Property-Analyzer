"""Local RDKit calculations shared by the interactive and batch workflows."""
import csv
import io
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import Descriptors, rdMolDescriptors, rdFingerprintGenerator

EXAMPLES = {
    "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Caffeine": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
    "Ethanol": "CCO",
    "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
}
MAX_ROWS = 500

def molecule(smiles):
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("Enter a SMILES string, or choose an example to get started.")
    if len(smiles) > 2000:
        raise ValueError("Please use a SMILES string of 2,000 characters or fewer.")
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError("We couldn’t read this molecule. Check brackets, ring numbers, and element symbols, or try an example.")
    return mol

def describe(smiles):
    mol = molecule(smiles)
    values = {
        "Canonical SMILES": Chem.MolToSmiles(mol),
        "Formula": rdMolDescriptors.CalcMolFormula(mol),
        "Molecular weight (g/mol)": Descriptors.MolWt(mol),
        "LogP (Crippen)": Descriptors.MolLogP(mol),
        "H-bond donors": Descriptors.NumHDonors(mol),
        "H-bond acceptors": Descriptors.NumHAcceptors(mol),
        "TPSA (Å²)": Descriptors.TPSA(mol),
        "Rotatable bonds": Descriptors.NumRotatableBonds(mol),
        "Rings": rdMolDescriptors.CalcNumRings(mol),
    }
    values["Rule-of-five violations"] = sum([
        values["Molecular weight (g/mol)"] > 500,
        values["LogP (Crippen)"] > 5,
        values["H-bond donors"] > 5,
        values["H-bond acceptors"] > 10,
    ])
    return values

def similarity(left, right):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    return DataStructs.TanimotoSimilarity(
        generator.GetFingerprint(molecule(left)), generator.GetFingerprint(molecule(right)))

def analyze_batch(values):
    if len(values) > MAX_ROWS:
        raise ValueError(f"Please upload {MAX_ROWS} molecules or fewer at a time.")
    rows = []
    for number, value in enumerate(values, start=1):
        row = {"Row": number, "Input SMILES": value if isinstance(value, str) else ""}
        try:
            row.update(describe(value))
            row["Status"] = "OK"
        except ValueError as exc:
            row["Status"] = str(exc)
        rows.append(row)
    return rows

def export_csv(rows):
    """Neutralize spreadsheet formula cells in user-controlled text."""
    fields = list(dict.fromkeys(key for row in rows for key in row))
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=fields)
    writer.writeheader()
    for row in rows:
        writer.writerow({key: "'" + value if isinstance(value, str) and value.lstrip().startswith(("=", "+", "-", "@")) else value for key, value in row.items()})
    return output.getvalue()
