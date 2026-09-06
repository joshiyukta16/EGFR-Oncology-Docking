import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
import pandas as pd
from docking_utils import (
    smiles_to_3d_molblock, 
    get_all_variant_scores, 
    fetch_smiles_from_pubchem,
    load_receptor_pdb,
    calculate_lipinski_properties,
    get_receptor_sequence
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
tab1, tab2, tab3 = st.tabs(["3D Docking Viewer & Analysis", "Comparative Analysis Table", "Receptor Info"])

with tab1:
    # --- TOP HORIZONTAL SECTION: 3D VISUALIZER & SEQUENCE DISPLAY ---
    st.subheader("1. Interactive 3D Complex & Sequence Inspection")
    
    # Load receptor data first to extract sequence
    receptor_data, file_ext = load_receptor_pdb(target_variant)
    seq_list, seq_str = get_receptor_sequence(receptor_data) if receptor_data else ([], "")

    # Expandable Horizontal Receptor Sequence Viewer
    with st.expander("🔬 View Target Receptor Amino Acid Sequence", expanded=True):
        if seq_str:
            st.text_area("Receptor Primary Sequence (1-Letter Code):", value=seq_str, height=70, disabled=True)
            
            # Interactive Residue Picker
            res_options = [item["label"] for item in seq_list]
            selected_custom_residues = st.multiselect(
                "Select Residues to Highlight in 3D Structure:",
                options=res_options,
                default=[r for r in ["MET793", "LYS745", "THR790", "MET790", "LEU858", "ARG858"] if r in res_options]
            )
        else:
            st.info("Receptor file not found. Load receptor to extract sequence.")
            selected_custom_residues = []

    # 3D Visualizer Box across full horizontal width
    if "molblock" in st.session_state:
        view = py3Dmol.view(width=1100, height=500)
        
        # 1. Render Ligand (Green Carbon Sticks)
        view.addModel(st.session_state["molblock"], "mol")
        view.setStyle({'model': -1}, {'stick': {'colorscheme': 'greenCarbon'}})
        
        # 2. Render Receptor Backbone (Spectrum Cartoon)
        if receptor_data:
            fmt = "pdbqt" if file_ext == "pdbqt" else "pdb"
            view.addModel(receptor_data, fmt)
            view.setStyle({'model': -1}, {'cartoon': {'color': 'spectrum'}})
            
            # 3. Highlight Selected Residues in Red/Yellow Sticks
            if selected_custom_residues:
                # Extract numerical residue indices
                res_nums = [r[3:] for r in selected_custom_residues]
                view.addStyle(
                    {'model': -1, 'resi': res_nums},
                    {'stick': {'colorscheme': 'redCarbon'}, 'cartoon': {}}
                )
                
        view.zoomTo()
        html_content = view._make_html()
        components.html(html_content, height=520, width=1120)
    else:
        st.info("Enter ligand SMILES below and click 'Run Virtual Screening & Docking' to generate 3D complex.")

    st.markdown("---")

    # --- BOTTOM SECTION: INPUTS, BINDING AFFINITY, & DESCRIPTORS ---
    b_col1, b_col2 = st.columns([1, 1])

    with b_col1:
        st.subheader("2. Ligand Input & Processing")
        smiles_input = st.text_area("SMILES String Format:", value=selected_smiles, height=100)
        run_btn = st.button("Run Virtual Screening & Docking", type="primary")
        
        if run_btn:
            molblock = smiles_to_3d_molblock(smiles_input, target_variant)
            if molblock:
                st.session_state["molblock"] = molblock
                st.session_state["scores"] = get_all_variant_scores(smiles_input)
                st.success("Structure and docking pose calculated!")
                st.rerun()
            else:
                st.error("Invalid SMILES input string.")

        if "scores" in st.session_state:
            st.markdown("---")
            st.subheader("3. Predicted Binding Affinity")
            selected_score = st.session_state["scores"][target_variant]
            st.metric(label=f"Binding Energy ({target_variant})", value=f"{selected_score} kcal/mol")
            
            st.download_button(
                label="Download Generated 3D Ligand (.sdf/.mol)",
                data=st.session_state["molblock"],
                file_name="docked_ligand.mol",
                mime="chemical/x-mdl-molfile"
            )

    with b_col2:
        if "scores" in st.session_state:
            st.subheader("4. Comprehensive Molecular Descriptors")
            props = calculate_lipinski_properties(smiles_input)
            
            if props:
                df_descriptors = pd.DataFrame(list(props.items()), columns=["Descriptor Property", "Calculated Value"])
                st.dataframe(df_descriptors, width="stretch", hide_index=True)
                
                csv_data = df_descriptors.to_csv(index=False).encode('utf-8')
                st.download_button(
                    label="Download Descriptors (.csv)",
                    data=csv_data,
                    file_name="ligand_descriptors.csv",
                    mime="text/csv"
                )

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
    * **Interactive Sequence Search**: Select any amino acid from the upper panel to highlight its 3D sidechain in red inside the active site.
    """) 