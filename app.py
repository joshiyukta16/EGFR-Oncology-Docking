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
    get_receptor_sequence,
    compute_dynamic_binding_score
)

st.set_page_config(page_title="EGFR Molecular Docking Dashboard", layout="wide")
st.title("EGFR Lung Cancer Molecular Docking Platform")

DEFAULT_SMILES = "COC1=C(OCCCN2CCOCC2)C=C3C(=C1)N=CN=C3NC4=CC(=C(C=C4)F)Cl"

# Session State Initialization
if "docking_ran" not in st.session_state:
    st.session_state["docking_ran"] = False
if "smiles" not in st.session_state:
    st.session_state["smiles"] = DEFAULT_SMILES
if "selected_residues" not in st.session_state:
    st.session_state["selected_residues"] = []

# Sidebar Controls
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

tab1, tab2, tab3 = st.tabs(["3D Docking Viewer & Analysis", "Comparative Analysis Table", "Receptor Info"])

with tab1:
    st.subheader("1. Interactive 3D Complex & Sequence Inspection")
    
    receptor_data, file_ext = load_receptor_pdb(target_variant)
    seq_list, _ = get_receptor_sequence(receptor_data) if receptor_data else ([], "")
    res_options = [item["label"] for item in seq_list]

    selected_custom_residues = st.multiselect(
        "Select Residues to Highlight in Sequence & 3D Structure:",
        options=res_options,
        key="selected_residues"
    )
    
    selected_nums = [r[3:] for r in selected_custom_residues]
    _, html_sequence = get_receptor_sequence(receptor_data, selected_nums) if receptor_data else ([], "")

    with st.expander("🔬 View Primary Sequence (Highlighted Residues in Red)", expanded=True):
        if html_sequence:
            st.markdown(
                f'<div style="font-family: monospace; font-size: 14px; line-height: 24px; word-wrap: break-word; background: #ffffff; padding: 12px; border-radius: 8px; border: 1px solid #e0e0e0; max-height: 120px; overflow-y: auto;">{html_sequence}</div>',
                unsafe_allow_html=True
            )
        else:
            st.info("Receptor structure file missing from /receptors folder.")

    active_smiles = st.session_state.get("smiles", DEFAULT_SMILES)

    # Visualizer Canvas (Renders docked pose ONLY when Virtual Screening is executed)
    if st.session_state["docking_ran"]:
        molblock = smiles_to_3d_molblock(active_smiles, receptor_data)
        
        if receptor_data and molblock:
            view = py3Dmol.view(width=1100, height=520)
            
            # Model 0: Ligand
            view.addModel(molblock, "mol")
            view.setStyle({'model': 0}, {'stick': {'colorscheme': 'greenCarbon'}})
            
            # Model 1: Receptor
            fmt = "pdbqt" if file_ext == "pdbqt" else "pdb"
            view.addModel(receptor_data, fmt)
            view.setStyle({'model': 1}, {'cartoon': {'color': 'spectrum'}})
            
            # Ribbon Highlights
            if selected_custom_residues:
                for res_label in selected_custom_residues:
                    res_num = res_label[3:]
                    view.addStyle(
                        {'model': 1, 'resi': str(res_num)},
                        {'cartoon': {'color': 'red'}}
                    )
                    view.addResidueLabel(
                        {'model': 1, 'resi': str(res_num)},
                        res_label,
                        {'fontColor': 'black', 'backgroundColor': 'yellow', 'backgroundOpacity': 0.85}
                    )
            
            # Click Interaction Callback
            view.setClickable(
                {'model': 1},
                True,
                """
                function(atom, viewer, event, container) {
                    if(atom && atom.resn && atom.resi) {
                        var tag = atom.resn + atom.resi;
                        viewer.addResidueLabel({'resi': atom.resi}, tag, {fontColor: 'black', backgroundColor: 'yellow'});
                        viewer.render();
                    }
                }
                """
            )
            
            view.zoomTo({'model': 0})
            html_content = view._make_html()
            components.html(html_content, height=540, width=1120)
        else:
            st.error("Error generating 3D complex. Verify structure files and SMILES format.")
    else:
        st.info("Visualizer canvas clean. Click 'Run Virtual Screening & Docking' below to execute docking.")

    st.markdown("---")

    b_col1, b_col2 = st.columns([1, 1])

    with b_col1:
        st.subheader("2. Ligand Input & Processing")
        smiles_input = st.text_area("SMILES String Format:", value=active_smiles, height=100)
        
        run_btn = st.button("Run Virtual Screening & Docking", type="primary")
        
        if run_btn:
            st.session_state["smiles"] = smiles_input
            st.session_state["docking_ran"] = True
            st.rerun()

        selected_score = compute_dynamic_binding_score(active_smiles, target_variant)
        
        st.markdown("---")
        st.subheader("3. Predicted Binding Affinity")
        st.metric(label=f"Binding Energy ({target_variant})", value=f"{selected_score} kcal/mol")

    with b_col2:
        st.subheader("4. Comprehensive Molecular Descriptors")
        props = calculate_lipinski_properties(active_smiles)
        if props:
            df_descriptors = pd.DataFrame(list(props.items()), columns=["Descriptor Property", "Calculated Value"])
            st.dataframe(df_descriptors, width="stretch", hide_index=True)

with tab2:
    st.subheader("Multi-Variant Binding Energy Comparison Table")
    active_smiles = st.session_state.get("smiles", DEFAULT_SMILES)
    scores = get_all_variant_scores(active_smiles)
    df_scores = pd.DataFrame(list(scores.items()), columns=["EGFR Variant Target", "Predicted Affinity (kcal/mol)"])
    st.dataframe(df_scores.style.highlight_min(subset=["Predicted Affinity (kcal/mol)"], color="lightgreen"), width="stretch")

with tab3:
    st.subheader("Target Structure Metadata")
    st.markdown("""
    * **EGFR Wild-Type (1M17)**: Active kinase domain.
    * **EGFR T790M Mutant**: Gatekeeper mutation.
    * **EGFR L858R Mutant**: Activating point mutation.
    """) 