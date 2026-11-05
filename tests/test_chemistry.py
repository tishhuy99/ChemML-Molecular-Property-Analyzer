import pytest
from chemistry import describe, molecule, similarity, analyze_batch, export_csv

def test_ethanol_descriptors():
    result = describe("CCO")
    assert result["Formula"] == "C2H6O"
    assert result["Molecular weight (g/mol)"] == pytest.approx(46.069, abs=.001)
    assert result["H-bond donors"] == 1
    assert result["Rule-of-five violations"] == 0

@pytest.mark.parametrize("value", ["", " ", "nonsense", "C1CC", None, "C" * 2001])
def test_invalid_inputs(value):
    with pytest.raises(ValueError):
        molecule(value)

def test_similarity_equivalent_smiles():
    assert similarity("CCO", "OCC") == 1
    assert 0 <= similarity("CCO", "c1ccccc1") < 1

def test_batch_preserves_errors_and_row_numbers():
    rows = analyze_batch(["CCO", "bad", "", "O"])
    assert [row["Row"] for row in rows] == [1, 2, 3, 4]
    assert [row["Status"] == "OK" for row in rows] == [True, False, False, True]
    with pytest.raises(ValueError):
        analyze_batch(["O"] * 501)

def test_export_neutralizes_formulas():
    result = export_csv([{"Input": "=1+1", "LogP": -1.2}])
    assert "'=1+1" in result
    assert "-1.2" in result
