"""Educational functional-group acidity references, not a pKa prediction model.

Atom IDs refer to the parsed input after AddHs (one-based), never canonical SMILES.
Each H is listed, but H atoms on one parent share a first-deprotonation reference.
No independent microscopic or sequential dissociation constants are calculated.
"""
from rdkit import Chem
from rdkit.Chem import Draw
from chemistry import molecule

SOURCE = "https://www.chem.indiana.edu/wp-content/uploads/2018/03/pka-chart.pdf"
METHOD = "Functional-group reference v1; solvent unspecified; not a molecular prediction"
# Priority is most specific first. Atom map 1 identifies the proton-bearing atom.
RULES = [
    ("Carboxylic acid O–H", "[O;H1;+0:1]-[C](=[O])-[#6]", "4–5"),
    ("Phenol O–H", "[O;H1;+0:1]-[c]", "≈10"),
    ("Alcohol O–H", "[O;H1;+0:1]-[C;X4]", "16–18"),
    ("Thiol S–H", "[S;H1;X2;+0:1]-[C;X4]", "≈10"),
    ("Amide N–H", "[N;!a;+0;!H0:1]-[C](=[O])", "≈18"),
    ("Alpha C–H of ester", "[C;X4;!H0;+0:1]-[C](=[O])-[O;H0]-[#6]", "≈25"),
    ("Alpha C–H of ketone", "[C;X4;!H0;+0:1]-[C](=[O])-[#6]", "≈20"),
    ("Alpha C–H of aldehyde", "[C;X4;!H0;+0:1]-[C;H1](=[O])", "≈20"),
    ("Terminal alkyne C–H", "[C;H1;+0:1]#[C]", "≈25"),
]
PATTERNS = [(name, Chem.MolFromSmarts(pattern), value) for name, pattern, value in RULES]


def _simple_amine(atom):
    """Exclude conjugated, aromatic, heteroatom-substituted and acylated N."""
    neighbors = [a for a in atom.GetNeighbors() if a.GetAtomicNum() != 1]
    return (atom.GetAtomicNum() == 7 and not atom.GetIsAromatic()
            and atom.GetFormalCharge() in (0, 1)
            and bool(neighbors)
            and all(a.GetAtomicNum() == 6 and a.GetHybridization() == Chem.HybridizationType.SP3 for a in neighbors)
            and all(b.GetBondType() == Chem.BondType.SINGLE for b in atom.GetBonds()))


def site_analysis(smiles):
    mol = molecule(smiles)
    hydrogenated = Chem.AddHs(mol)
    matches = {}
    for name, pattern, value in PATTERNS:
        anchor = next(a.GetIdx() for a in pattern.GetAtoms() if a.GetAtomMapNum() == 1)
        for match in mol.GetSubstructMatches(pattern, uniquify=False, maxMatches=100000):
            matches.setdefault(match[anchor], (name, value))
    rows = []
    complex_state = any(a.GetFormalCharge() or a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    for hydrogen in hydrogenated.GetAtoms():
        if hydrogen.GetAtomicNum() != 1:
            continue
        neighbors = list(hydrogen.GetNeighbors())
        parent = neighbors[0] if len(neighbors) == 1 else None
        kind, value = "Uncovered environment", "Not estimated"
        note = "No supported functional-group rule. Use a validated predictor or measured data."
        if parent is not None and parent.GetIdx() < mol.GetNumAtoms():
            atom = mol.GetAtomWithIdx(parent.GetIdx())
            # Multiple electron-withdrawing groups invalidate simple alpha-C–H rules.
            carbonyls = sum(any(b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(n).GetAtomicNum() == 8 for b in n.GetBonds()) for n in atom.GetNeighbors())
            if not complex_state and carbonyls <= 1 and hydrogen.GetIsotope() == 0:
                kind, value = matches.get(atom.GetIdx(), (kind, value))
                if value == "Not estimated" and _simple_amine(atom):
                    kind, value = "Amine N–H (neutral acid)", "38–40"
                # Restrict hydrocarbon references to C/H-only molecules.
                if value == "Not estimated" and all(a.GetAtomicNum() in (1, 6) for a in mol.GetAtoms()):
                    if atom.GetIsAromatic():
                        kind, value = "Aryl C–H", "≈43"
                    elif atom.GetHybridization() == Chem.HybridizationType.SP3:
                        if any(n.GetIsAromatic() for n in atom.GetNeighbors()):
                            kind, value = "Benzylic C–H", "≈41"
                        elif all(n.GetHybridization() == Chem.HybridizationType.SP3 for n in atom.GetNeighbors()):
                            kind, value = "Alkane C–H", ">50"
            if value != "Not estimated":
                note = "Typical family value only; substituents, solvent, temperature and other sites can shift pKa substantially."
            elif complex_state:
                note = "Charged or radical structure: this reference method abstains; specify a relevant microstate and validated method."
            elif carbonyls > 1:
                note = "Multiple adjacent carbonyls: simple single-group acidity rule does not apply."
            elif hydrogen.GetIsotope():
                note = "Isotope-specific acidity is not modeled."
        rows.append({
            "Hydrogen": f"H{hydrogen.GetIdx()+1}",
            "Attached atom": f"{parent.GetSymbol()}{parent.GetIdx()+1}" if parent is not None else "Unbound",
            "Environment": kind, "Reference pKa": value,
            "Interpretation": note, "Method": METHOD,
            "Source": SOURCE if value != "Not estimated" else "",
        })
    basic = []
    if not complex_state:
        for atom in mol.GetAtoms():
            if _simple_amine(atom) and atom.GetDegree() < 4:
                basic.append({"Atom": f"N{atom.GetIdx()+1}", "Site": "Aliphatic amine · proton acceptor",
                              "Conjugate-acid reference pKaH": "9–10", "Source": SOURCE})
    return hydrogenated, rows, basic


def site_image(smiles, hydrogen_id=None):
    mol, rows, _ = site_analysis(smiles)
    if mol.GetNumAtoms() > 120:
        return None
    highlighted = []
    for atom in mol.GetAtoms():
        atom.SetProp("atomNote", f"{atom.GetSymbol()}{atom.GetIdx()+1}")
        if f"H{atom.GetIdx()+1}" == hydrogen_id:
            highlighted = [atom.GetIdx()] + [n.GetIdx() for n in atom.GetNeighbors()]
    return Draw.MolToImage(mol, size=(1000, 650), highlightAtoms=highlighted)
