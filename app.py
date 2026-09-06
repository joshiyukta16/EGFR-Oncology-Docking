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

# 1. Initialize Default Session State Variables
DEFAULT_SMILES = "COC1=C(OCCCN2CCOCC2)C=C3C(=C1)N=CN=C3NC4=CC(=C(C=C4)F)Cl"

if "smiles" not in st.session_state:
    st.session_state["smiles"] = DEFAULT_SMILES

# 2. Sidebar Configuration
st.sidebar.header("Target Configuration")
target_variant = st.sidebar.selectbox(
    "Select Target EGFR Structure", 
    ["EGFR Wild-Type (1M17)", "EGFR T790M Mutant", "EGFR L858R Mutant"]
)

st.sidebar.markdown("---")
st.sidebar.subheader("Presets / PubChem Search")
search_query = st.sidebar.text_input("Search Drug by Name:")

if search_query:
    fetched = fetch_smiles_from_pubchem(search_query)
    if fetched:
        st.sidebar.success(f"Found structure for '{search_query}'")
        st.session_state["smiles"] = fetched
    else:
        st.sidebar.error("Drug not found on PubChem.")

# 3. Main Interface Tabs
tab1, tab2, tab3 = st.tabs(["3D Docking Viewer & Analysis", "Comparative Analysis Table", "Receptor Info"])

with tab1:
    st.subheader("1. Interactive 3D Complex & Sequence Inspection")
    
    # Load Receptor Data
    receptor_data, file_ext = load_receptor_pdb(target_variant)
    
    # Extract sequence list for multi-select options
    seq_list, _ = get_receptor_sequence(receptor_data) if receptor_data else ([], "")
    res_options = [item["label"] for item in seq_list]

    selected_custom_residues = st.multiselect(
        "Select Residues to Highlight in Sequence & 3D Ribbon Structure:",
        options=res_options,
        default=[]
    )
    
    selected_nums = [r[3:] for r in selected_custom_residues]
    _, html_sequence = get_receptor_sequence(receptor_data, selected_nums) if receptor_data else ([], "")

    # Expandable Primary Sequence Viewer
    with st.expander("🔬 View Primary Sequence (Highlighted Residues in Red)", expanded=True):
        if html_sequence:
            st.markdown(
                f'<div style="font-family: monospace; font-size: 14px; line-height: 24px; word-wrap: break-word; background: #f8f9fa; padding: 12px; border-radius: 8px; border: 1px solid #e0e0e0; max-height: 120px; overflow-y: auto;">{html_sequence}</div>',
                unsafe_allow_html=True
            )
        else:
            st.info("Receptor structure file missing from /receptors folder.")

    # Calculate current 3D molblock for the active SMILES & variant
    active_smiles = st.session_state.get("smiles", DEFAULT_SMILES)
    molblock = smiles_to_3d_molblock(active_smiles, target_variant)

    # 4. 3D Visualizer Canvas (Guaranteed to render continuously)
    if receptor_data and molblock:
        view = py3Dmol.view(width=1100, height=520)
        
        # Model 0: Ligand
        view.addModel(molblock, "mol")
        view.setStyle({'model': 0}, {'stick': {'colorscheme': 'greenCarbon'}})
        
        # Model 1: Receptor Target Structure
        fmt = "pdbqt" if file_ext == "pdbqt" else "pdb"
        view.addModel(receptor_data, fmt)
        view.setStyle({'model': 1}, {'cartoon': {'color': 'spectrum'}})
        
        # Highlight selected residues on the cartoon ribbon backbone
        if selected_custom_residues:
            for res_label in selected_custom_residues:
                res_num = res_label[3:]
                
                # Turn ribbon segment solid red
                view.addStyle(
                    {'model': 1, 'resi': str(res_num)},
                    {'cartoon': {'color': 'red'}}
                )
                
                # Floating 3D Text Label
                view.addResidueLabel(
                    {'model': 1, 'resi': str(res_num)},
                    res_label,
                    {'fontColor': 'black', 'backgroundColor': 'yellow', 'backgroundOpacity': 0.85}
                )
        
        # Keep camera focused on docked ligand inside binding pocket
        view.zoomTo({'model': 0})
            
        html_content = view._make_html()
        components.html(html_content, height=540, width=1120)
    else:
        st.error("Unable to load 3D visualizer canvas. Check receptor PDB files or SMILES format.")

    st.markdown("---")

    # 5. Input Controls & Descriptor Sections
    b_col1, b_col2 = st.columns([1, 1])

    with b_col1:
        st.subheader("2. Ligand Input & Processing")
        smiles_input = st.text_area("SMILES String Format:", value=active_smiles, height=100)
        run_btn = st.button("Run Virtual Screening & Docking", type="primary")
        
        if run_btn:
            st.session_state["smiles"] = smiles_input
            st.rerun()

        scores = get_all_variant_scores(active_smiles)
        selected_score = scores.get(target_variant, -7.5)
        
        st.markdown("---")
        st.subheader("3. Predicted Binding Affinity")
        st.metric(label=f"Binding Energy ({target_variant})", value=f"{selected_score} kcal/mol")
        
        st.download_button(
            label="Download Generated 3D Ligand (.sdf/.mol)",
            data=molblock if molblock else "",
            file_name="docked_ligand.mol",
            mime="chemical/x-mdl-molfile"
        )

    with b_col2:
        st.subheader("4. Comprehensive Molecular Descriptors")
        props = calculate_lipinski_properties(active_smiles)
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
    active_smiles = st.session_state.get("smiles", DEFAULT_SMILES)
    scores = get_all_variant_scores(active_smiles)
    df_scores = pd.DataFrame(list(scores.items()), columns=["EGFR Variant Target", "Predicted Affinity (kcal/mol)"])
    st.dataframe(df_scores.style.highlight_min(subset=["Predicted Affinity (kcal/mol)"], color="lightgreen"), width="stretch")
    st.bar_chart(df_scores.set_index("EGFR Variant Target"))

with tab3:
    st.subheader("Target Structure Metadata")
    st.markdown("""
    * **EGFR Wild-Type (1M17)**: Active kinase domain conformation bound to Erlotinib/Gefitinib class.
    * **EGFR T790M Mutant**: Gatekeeper mutation causing drug resistance against 1st-generation TKIs.
    * **EGFR L858R Mutant**: Activating point mutation with heightened sensitivity to TKIs.
    """) 