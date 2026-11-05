from rdkit import Chem
from pka import site_analysis, site_image

def rows(smiles):
    return site_analysis(smiles)[1]

def test_every_hydrogen_is_accounted_for():
    for smiles in ['CCO', 'CC(=O)Oc1ccccc1C(=O)O', 'CN', 'c1ccccc1', '[2H]OC', '[H][H]']:
        mol, data, _ = site_analysis(smiles)
        assert len(data) == sum(a.GetAtomicNum() == 1 for a in mol.GetAtoms())
        assert len({r['Hydrogen'] for r in data}) == len(data)

def test_acid_alcohol_phenol_and_carbonyl_distinction():
    assert any(r['Reference pKa'] == '4–5' for r in rows('CC(=O)O'))
    assert any(r['Environment'] == 'Alcohol O–H' for r in rows('CCO'))
    assert any(r['Environment'] == 'Phenol O–H' for r in rows('Oc1ccccc1'))
    assert sum(r['Environment'] == 'Alpha C–H of ketone' for r in rows('CC(=O)C')) == 6
    assert sum(r['Environment'] == 'Alpha C–H of ester' for r in rows('CC(=O)OC')) == 3
    # Carboxylic acid alpha hydrogens must not be mislabeled ketone.
    assert not any(r['Environment'] == 'Alpha C–H of ketone' for r in rows('CC(=O)O'))

def test_multiple_carbonyls_are_not_treated_as_single_ketone():
    mol, data, _ = site_analysis('CC(=O)CC(C)=O')
    center = [r for r in data if r['Attached atom'] == 'C4']
    assert len(center) == 2
    assert all(r['Reference pKa'] == 'Not estimated' for r in center)

def test_acid_and_base_are_separate():
    _, protons, bases = site_analysis('CN')
    assert any(r['Reference pKa'] == '38–40' for r in protons)
    assert bases[0]['Conjugate-acid reference pKaH'] == '9–10'
    assert site_analysis('CN(C)C')[2]
    for smiles in ['CC(=O)N', 'Nc1ccccc1', 'n1ccccc1', 'C[N+](C)(C)C']:
        assert not site_analysis(smiles)[2]

def test_abstention_and_zero_hydrogens():
    assert rows('O=C=O') == []
    assert all(r['Reference pKa'] == 'Not estimated' for r in rows('C[NH3+]'))
    isotope = next(r for r in rows('[2H]OC') if r['Hydrogen'] == 'H1')
    assert isotope['Reference pKa'] == 'Not estimated'

def test_image_and_size_guard():
    assert site_image('CCO', rows('CCO')[0]['Hydrogen']).size == (1000,650)
    assert site_image('C'*50) is None

def test_proton_map_app_selection_and_no_hydrogen_state():
    from pathlib import Path
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_file(str(Path(__file__).parents[1] / 'app.py')).run(timeout=30)
    assert not app.exception
    app.text_input(key='smiles').set_value('CCO').run()
    picker = next(s for s in app.selectbox if s.label == 'Highlight a hydrogen')
    picker.select('H9').run()
    assert not app.exception
    assert any('Alcohol O–H' in m.value for m in app.markdown)
    app.text_input(key='smiles').set_value('O=C=O').run()
    assert not app.exception
    assert any('no hydrogens' in info.value for info in app.info)
