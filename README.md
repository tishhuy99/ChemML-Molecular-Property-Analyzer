# 🧪 ChemML — Molecular Property Analyzer

An interactive cheminformatics application for visualizing molecules and calculating molecular properties from SMILES strings.

Built with **Python, RDKit, and Streamlit**, this project is designed as a foundation for integrating cheminformatics with machine learning for molecular property prediction.

## Features

- Enter molecules using SMILES notation
- Generate 2D molecular structures with RDKit
- Calculate molecular weight
- Calculate LogP
- Calculate hydrogen-bond donors
- Calculate hydrogen-bond acceptors
- Interactive web interface built with Streamlit
- Foundation for future ML-based molecular property prediction

## Example

Enter a SMILES string such as aspirin:

```text
CC(=O)OC1=CC=CC=C1C(=O)O
```

The application generates the molecular structure and calculates properties including:

| Property | Description |
|---|---|
| Molecular Weight | Molecular mass calculated by RDKit |
| LogP | Estimated octanol/water partition coefficient |
| H-Bond Donors | Number of hydrogen-bond donor groups |
| H-Bond Acceptors | Number of hydrogen-bond acceptor groups |

## Tech Stack

- Python
- RDKit
- Streamlit
- NumPy
- Pandas
- scikit-learn
- Jupyter

## Project Structure

```text
chem-app/
├── app.py
├── Test-1.ipynb
├── requirements.txt
├── .gitignore
└── README.md
```

`app.py` contains the Streamlit application.

`Test-1.ipynb` is used for cheminformatics and machine-learning experimentation.

## Installation

Clone the repository:

```bash
git clone <your-repository-url>
cd chem-app
```

Create a Conda environment:

```bash
conda create -n chem-ml python=3.11
conda activate chem-ml
```

Install the dependencies:

```bash
pip install -r requirements.txt
```

## Running the Application

Start the Streamlit application:

```bash
python -m streamlit run app.py
```

Then open the local address displayed by Streamlit, typically port `8501`.

## Machine Learning Roadmap

The project will be expanded from molecular descriptor calculation into a molecular machine-learning platform.

Planned features include:

- Morgan molecular fingerprints
- Lipinski Rule of Five analysis
- Molecular similarity calculations
- Molecular dataset upload and analysis
- Descriptor generation for ML datasets
- Exploratory data analysis
- Molecular property prediction
- Model training with scikit-learn
- Random Forest and gradient-boosting models
- Model performance visualization
- Applicability-domain analysis
- PyTorch-based molecular models
- Graph neural networks

Potential prediction targets include aqueous solubility, LogP, permeability, and other physicochemical or ADME-related molecular properties.

## Scientific Purpose

This project explores the intersection of **chemistry, cheminformatics, data science, and machine learning**.

RDKit converts molecular structures into computational representations such as descriptors and fingerprints. These representations can subsequently be used as features for statistical and machine-learning models.

A future workflow will follow:

```text
SMILES
   ↓
RDKit Molecular Structure
   ↓
Descriptors / Molecular Fingerprints
   ↓
Data Processing
   ↓
Machine Learning Model
   ↓
Molecular Property Prediction
```

## Author

**Tish Huynh**

M.S. Chemistry — University of Connecticut  
B.A. Chemistry — Drexel University

Interests: cheminformatics, pharmaceutical development, analytical chemistry, process development, and machine learning.

## Disclaimer

This application is intended for educational, research, and software-development purposes. Computational predictions should not be interpreted as experimental measurements or medical advice.