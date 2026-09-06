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

# Sidebar
st.sidebar.header("Target Configuration")
target_variant = st.sidebar.selectbox(
    "Select Target EGFR Structure", 
    ["EGFR Wild-Type (1M17)", "EGFR T790M Mutant", "EGFR L858R Mutant"]
)

st.sidebar.markdown("---")
st.sidebar.subheader("Presets / PubChem Search")
search_query = st.sidebar.text_input("Search Drug by Name:")
selected_smiles = "COC1=C(OCCCN2CCOCC2)C=C3C(=C1)N=CN=C3NC4=CC(=C(C=C4)F)Cl"

if search_query:
    fetched = fetch_smiles_from_pubchem(search_query)
    if fetched:
        st.sidebar.success(f"Found structure for '{search_query}'")
        selected_smiles = fetched

tab1, tab2, tab3 = st.tabs(["3D Docking Viewer & Analysis", "Comparative Analysis Table", "Receptor Info"])

with tab1:
    st.subheader("1. Interactive 3D Complex & Sequence Inspection")
    
    receptor_data, file_ext = load_receptor_pdb(target_variant)
    
    dummy_seq, _ = get_receptor_sequence(receptor_data) if receptor_data else ([], "")
    res_options = [item["label"] for item in dummy_seq]

    selected_custom_residues = st.multiselect(
        "Select Residues to Highlight in Sequence & 3D Ribbon Structure:",
        options=res_options,
        default=[]
    )
    
    selected_nums = [r[3:] for r in selected_custom_residues]
    seq_list, html_sequence = get_receptor_sequence(receptor_data, selected_nums) if receptor_data else ([], "")

    with st.expander("🔬 View Primary Sequence (Highlighted Residues in Red)", expanded=True):
        if html_sequence:
            st.markdown(
                f'<div style="font-family: monospace; font-size: 14px; line-height: 24px; word-wrap: break-word; background: #f8f9fa; padding: 12px; border-radius: 8px; border: 1px solid #e0e0e0; max-height: 120px; overflow-y: auto;">{html_sequence}</div>',
                unsafe_allow_html=True
            )

    # Force recalculate molblock whenever SMILES or target variant changes
    molblock = smiles_to_3d_molblock(selected_smiles, target_variant)

    # 3D Visualizer Box (Always renders regardless of residue selection)
    if molblock and receptor_data:
        view = py3Dmol.view(width=1100, height=520)
        
        # Model 0: Ligand
        view.addModel(molblock, "mol")
        view.setStyle({'model': 0}, {'stick': {'colorscheme': 'greenCarbon'}})
        
        # Model 1: Receptor Backbone
        fmt = "pdbqt" if file_ext == "pdbqt" else "pdb"
        view.addModel(receptor_data, fmt)
        view.setStyle({'model': 1}, {'cartoon': {'color': 'spectrum'}})
        
        # Highlight selected residues directly on the ribbon backbone
        if selected_custom_residues:
            for res_label in selected_custom_residues:
                res_num = res_label[3:]
                
                # Color cartoon ribbon segment red
                view.addStyle(
                    {'model': 1, 'resi': str(res_num)},
                    {'cartoon': {'color': 'red'}}
                )
                
                # Add 3D Floating Label
                view.addResidueLabel(
                    {'model': 1, 'resi': str(res_num)},
                    res_label,
                    {'fontColor': 'black', 'backgroundColor': 'yellow', 'backgroundOpacity': 0.85}
                )
        
        # Always center camera on the ligand inside the active site
        view.zoomTo({'model': 0})
            
        html_content = view._make_html()
        components.html(html_content, height=540, width=1120)

    st.markdown("---")

    b_col1, b_col2 = st.columns([1, 1])

    with b_col1:
        st.subheader("2. Ligand Input & Processing")
        smiles_input = st.text_area("SMILES String Format:", value=selected_smiles, height=100)
        run_btn = st.button("Run Virtual Screening & Docking", type="primary")
        
        scores = get_all_variant_scores(smiles_input)
        selected_score = scores[target_variant]
        
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
    scores = get_all_variant_scores(selected_smiles)
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