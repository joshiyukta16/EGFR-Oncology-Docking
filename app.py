import streamlit as st
import streamlit.components.v1 as components  # Fixes 3Dmol.js rendering
import py3Dmol
import pandas as pd
from docking_utils import (
    smiles_to_3d_molblock, 
    get_all_variant_scores, 
    fetch_smiles_from_pubchem,
    load_receptor_pdb
)

st.set_page_config(page_title="EGFR Molecular Docking Dashboard", layout="wide")
st.title("EGFR Lung Cancer Molecular Docking Platform")

# Sidebar Configuration
st.sidebar.header("Target Configuration")
target_variant = st.sidebar.selectbox(
    "Select Target EGFR Structure", 
    ["EGFR Wild-Type (1M17)", "EGFR T790M Mutant", "EGFR L858R Mutant"]
)

st.sidebar.markdown("---")
st.sidebar.subheader("Presets / PubChem Search")
search_query = st.sidebar.text_input("Search Drug by Name (e.g. Erlotinib, Gefitinib, Osimertinib):")

selected_smiles = "COC1=C(OCCCN2CCOCC2)C=C3C(=C1)N=CN=C3NC4=CC(=C(C=C4)F)Cl"

if search_query:
    fetched = fetch_smiles_from_pubchem(search_query)
    if fetched:
        st.sidebar.success(f"Found structure for '{search_query}'")
        selected_smiles = fetched
    else:
        st.sidebar.error("Drug not found on PubChem. Using default structure.")

# Main Interface Tabs
tab1, tab2, tab3 = st.tabs(["3D Docking Viewer", "Comparative Analysis Table", "Receptor Info"])

with tab1:
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("1. Ligand Input & Processing")
        smiles_input = st.text_area("SMILES String:", value=selected_smiles, height=100)
        
        show_receptor = st.checkbox("Display Processed Receptor Target in 3D Viewer", value=True)
        run_btn = st.button("Run Virtual Screening & Docking", type="primary")
        
        if run_btn:
            molblock = smiles_to_3d_molblock(smiles_input)
            if molblock:
                st.session_state["molblock"] = molblock
                st.session_state["scores"] = get_all_variant_scores(smiles_input)
                st.success("Structure and docking pose calculated!")
            else:
                st.error("Invalid SMILES input string.")

        if "scores" in st.session_state:
            st.subheader("2. Target Binding Affinity")
            selected_score = st.session_state["scores"][target_variant]
            st.metric(label=f"Binding Energy ({target_variant})", value=f"{selected_score} kcal/mol")
            
            # Export Option
            st.download_button(
                label="Download Generated 3D Ligand (.sdf/.mol)",
                data=st.session_state["molblock"],
                file_name="docked_ligand.mol",
                mime="chemical/x-mdl-molfile"
            )

    with col2:
        st.subheader("3. Interactive 3D Complex Viewer")
        if "molblock" in st.session_state:
            view = py3Dmol.view(width=550, height=450)
            
            # Render Ligand (Green Sticks)
            view.addModel(st.session_state["molblock"], "mol")
            view.setStyle({'model': -1}, {'stick': {'colorscheme': 'greenCarbon'}})
            
            # Render Receptor Backbone (Spectrum Ribbon/Cartoon)
            if show_receptor:
                receptor_data, file_ext = load_receptor_pdb(target_variant)
                if receptor_data:
                    fmt = "pdbqt" if file_ext == "pdbqt" else "pdb"
                    view.addModel(receptor_data, fmt)
                    view.setStyle({'model': -1}, {'cartoon': {'color': 'spectrum'}})
                else:
                    st.info(f"Receptor file for '{target_variant}' not found in /receptors directory.")

            view.zoomTo()
            html_content = view._make_html()
            
            # Replaced st.html with components.html to allow 3Dmol.js CDN scripts
            components.html(html_content, height=470, width=560)
        else:
            st.info("Click 'Run Virtual Screening & Docking' to view 3D interaction poses.")

with tab2:
    st.subheader("Multi-Variant Binding Energy Comparison Table")
    if "scores" in st.session_state:
        df_scores = pd.DataFrame(list(st.session_state["scores"].items()), columns=["EGFR Variant Target", "Predicted Affinity (kcal/mol)"])
        
        st.dataframe(df_scores.style.highlight_min(subset=["Predicted Affinity (kcal/mol)"], color="lightgreen"), width="stretch")
        
        st.bar_chart(df_scores.set_index("EGFR Variant Target"))
    else:
        st.info("Run docking on Tab 1 to generate comparison statistics.")

with tab3:
    st.subheader("Target Structure Metadata")
    st.markdown("""
    * **EGFR Wild-Type (1M17)**: Active kinase domain conformation bound to Erlotinib/Gefitinib class[cite: 1].
    * **EGFR T790M Mutant**: Gatekeeper mutation causing drug resistance against 1st-generation TKIs.
    * **EGFR L858R Mutant**: Activating point mutation with heightened sensitivity to TKIs.
    """) 