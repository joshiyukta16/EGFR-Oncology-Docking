import streamlit as st
from docking_utils import smiles_to_3d_block

st.set_page_config(page_title="EGFR Molecular Docking", layout="wide")
st.title("EGFR Lung Cancer Molecular Docking Tool")

st.sidebar.header("Target Configuration")
target_variant = st.sidebar.selectbox(
    "Select EGFR Variant", 
    ["EGFR Wild-Type (1M17)", "EGFR T790M Mutant", "EGFR L858R Mutant"]
)

smiles_input = st.text_input(
    "Enter Ligand SMILES", 
    value="CN1CCN(CC1)C2=CC3=C(C=C2)N=CN=C3NC4=C(C=C(C=C4)F)Cl"
)

if st.button("Generate 3D Structure"):
    mol_block = smiles_to_3d_block(smiles_input)
    if mol_block:
        st.success("3D Conformation Generated Successfully!")
    else:
        st.error("Invalid SMILES input. Please check the structure.") 