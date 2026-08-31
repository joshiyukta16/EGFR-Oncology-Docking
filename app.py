import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
from docking_utils import smiles_to_3d_molblock, predict_binding_affinity

st.set_page_config(page_title="EGFR Molecular Docking", layout="wide")
st.title("EGFR Lung Cancer Molecular Docking Dashboard")

# Sidebar
st.sidebar.header("Target Configuration")
target_variant = st.sidebar.selectbox(
    "Select EGFR Variant", 
    ["EGFR Wild-Type (1M17)", "EGFR T790M Mutant", "EGFR L858R Mutant"]
)

# Main Interface
col_input, col_view = st.columns([1, 1])

with col_input:
    st.subheader("1. Ligand Input")
    smiles_input = st.text_input(
        "Enter Ligand SMILES (Default: Gefitinib):", 
        value="COC1=C(OCCCN2CCOCC2)C=C3C(=C1)N=CN=C3NC4=CC(=C(C=C4)F)Cl"
    )
    
    if st.button("Generate & Dock Ligand", type="primary"):
        molblock = smiles_to_3d_molblock(smiles_input)
        if molblock:
            st.success("3D Structure Successfully Generated!")
            score = predict_binding_affinity(target_variant, smiles_input)
            
            st.subheader("2. Predicted Results")
            st.metric(label=f"Binding Affinity ({target_variant})", value=f"{score} kcal/mol")
            
            st.session_state["molblock"] = molblock
        else:
            st.error("Invalid SMILES string. Please check the structure.")

with col_view:
    st.subheader("3. Interactive 3D Viewer")
    if "molblock" in st.session_state:
        # Create py3Dmol view
        view = py3Dmol.view(width=500, height=400)
        view.addModel(st.session_state["molblock"], "mol")
        view.setStyle({'stick': {'colorscheme': 'greenCarbon'}})
        view.zoomTo()
        
        # Generate 3Dmol HTML content
        html_content = view._make_html()
        
        # Render HTML component safely without deprecation warnings
        components.html(html_content, height=410, width=510)
    else:
        st.info("Click 'Generate & Dock Ligand' to render the 3D structure.") 
        