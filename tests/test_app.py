from pathlib import Path
from streamlit.testing.v1 import AppTest

APP = str(Path(__file__).parents[1] / "app.py")

def test_explore_and_validation():
    app = AppTest.from_file(APP).run(timeout=30)
    assert not app.exception
    assert any(metric.value == "180.16" for metric in app.metric)
    next(button for button in app.button if button.label == "Ethanol").click().run()
    assert app.text_input(key="smiles").value == "CCO"
    assert any(metric.value == "46.07" for metric in app.metric)
    app.text_input(key="smiles").set_value("invalid").run()
    assert app.error
    assert not app.exception
    app.text_input(key="smiles").set_value("").run()
    assert app.info
    assert not app.exception
