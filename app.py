import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
import pandas as pd
from docking_utils import (
    generate_docked_poses, 
    fetch_smiles_from_pubchem,
    load_receptor_pdb,
    calculate_lipinski_properties,
    get_receptor_sequence,
    extract_ligand_atoms,
    compute_dynamic_binding_score
)

# Dark Navy Blue and Neon Green Theme Setup
st.set_page_config(page_title="EGFR Molecular Docking Dashboard", layout="wide")

st.markdown("""
    <style>
        .stApp {
            background-color: #0b132b;
            color: #ffffff;
        }
        h1, h2, h3, h4 {
            color: #ffffff !important;
            font-weight: 700 !important;
        }
        .neon-box {
            background-color: #001219;
            border: 2px solid #00f5d4;
            box-shadow: 0 0 10px #00f5d4;
            padding: 12px;
            border-radius: 8px;
            max-height: 140px;
            overflow-y: auto;
            margin-bottom: 15px;
        }
        .stButton>button {
            border-radius: 6px;
        }
    </style>
""", unsafe_allow_html=True)

st.title("EGFR Lung Cancer Molecular Docking Platform")

DEFAULT_SMILES = "COC1=C(OCCCN2CCOCC2)C=C3C(=C1)N=CN=C3NC4=CC(=C(C=C4)F)Cl"

# Initialize Session States
if "docking_completed" not in st.session_state:
    st.session_state["docking_completed"] = False
if "smiles" not in st.session_state:
    st.session_state["smiles"] = DEFAULT_SMILES
if "selected_receptor_res" not in st.session_state:
    st.session_state["selected_receptor_res"] = None
if "selected_ligand_atom" not in st.session_state:
    st.session_state["selected_ligand_atom"] = None

# Sidebar
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
    receptor_data, file_ext = load_receptor_pdb(target_variant)
    seq_list = get_receptor_sequence(receptor_data) if receptor_data else []
    
    # --- UPPER PANEL: INTERACTIVE RECEPTOR SEQUENCE (NEON STYLE) ---
    st.subheader("1. Interactive Target Receptor Sequence")
    
    if seq_list:
        # Render Sequence in Neon Box
        seq_html = '<div class="neon-box">'
        for res in seq_list:
            is_selected = (st.session_state["selected_receptor_res"] == res["label"])
            bg_color = "#ff0055" if is_selected else "#00f5d4"
            text_color = "#ffffff" if is_selected else "#000000"
            seq_html += f'<span style="background-color: {bg_color}; color: {text_color}; padding: 3px 6px; margin: 2px; border-radius: 4px; font-weight: bold; display: inline-block; font-size: 13px;" title="{res["label"]}">{res["single"]}</span>'
        seq_html += '</div>'
        st.markdown(seq_html, unsafe_allow_html=True)
        
        # Interactive Clickable Sequence Buttons
        res_cols = st.columns(10)
        selected_res_button = st.selectbox(
            "Click/Select Residue to Focus on 3D Ribbon Pose:",
            options=["None"] + [r["label"] for r in seq_list]
        )
        if selected_res_button != "None":
            st.session_state["selected_receptor_res"] = selected_res_button
        else:
            st.session_state["selected_receptor_res"] = None
    else:
        st.info("Load receptor file to view amino acid sequence.")

    st.markdown("---")

    # --- MAIN VISUALIZER PANEL ---
    st.subheader("2. Docked Pose Visualizer Canvas")
    
    # Check if screening has been executed
    if st.session_state["docking_completed"]:
        active_smiles = st.session_state["smiles"]
        poses_data = generate_docked_poses(active_smiles, target_variant)
        
        # Pose Selector Controls (Top 3 Poses)
        p_col1, p_col2 = st.columns([1, 3])
        with p_col1:
            selected_pose_key = st.radio("Select Docking Conformer:", ["Pose 1 (Top Rank)", "Pose 2", "Pose 3"])
            pose_idx_map = {"Pose 1 (Top Rank)": "Pose 1", "Pose 2": "Pose 2", "Pose 3": "Pose 3"}
            current_pose = poses_data[pose_idx_map[selected_pose_key]]
            
            st.metric("Binding Affinity", f"{current_pose['score']} kcal/mol")
            st.metric("RMSD Value", f"{current_pose['rmsd']} Å")
            
        with p_col2:
            view = py3Dmol.view(width=850, height=500)
            
            # Model 0: Ligand Pose
            view.addModel(current_pose["molblock"], "mol")
            view.setStyle({'model': 0}, {'stick': {'colorscheme': 'greenCarbon'}})
            
            # Highlight Specific Ligand Atom if clicked in Lower Panel
            if st.session_state["selected_ligand_atom"] is not None:
                atom_idx = st.session_state["selected_ligand_atom"]
                view.addStyle(
                    {'model': 0, 'serial': atom_idx},
                    {'sphere': {'color': 'cyan', 'radius': 0.8}}
                )
            
            # Model 1: Receptor Structure
            if receptor_data:
                fmt = "pdbqt" if file_ext == "pdbqt" else "pdb"
                view.addModel(receptor_data, fmt)
                view.setStyle({'model': 1}, {'cartoon': {'color': 'spectrum'}})
                
                # Highlight Selected Residue directly on Cartoon Ribbon
                if st.session_state["selected_receptor_res"]:
                    res_label = st.session_state["selected_receptor_res"]
                    res_num = res_label[3:]
                    
                    # Highlight Cartoon Ribbon in Solid Red
                    view.addStyle(
                        {'model': 1, 'resi': str(res_num)},
                        {'cartoon': {'color': 'red'}}
                    )
                    
                    # Add 3D Residue Text Label directly on the ribbon
                    view.addResidueLabel(
                        {'model': 1, 'resi': str(res_num)},
                        res_label,
                        {'fontColor': 'white', 'backgroundColor': '#ff0055', 'backgroundOpacity': 0.9}
                    )
            
            # Direct 3D Surface Clickability (Shows Residue Tag on Click)
            view.setClickable(
                {'model': 1},
                True,
                """
                function(atom, viewer, event, container) {
                    if(atom && atom.resn && atom.resi) {
                        var tag = atom.resn + atom.resi;
                        viewer.addResidueLabel({'resi': atom.resi}, tag, {fontColor: 'white', backgroundColor: '#ff0055'});
                        viewer.render();
                    }
                }
                """
            )
            
            view.zoomTo({'model': 0})
            html_content = view._make_html()
            components.html(html_content, height=520, width=860)

        st.markdown("---")

        # --- LOWER PANEL: INTERACTIVE LIGAND SEQUENCE / ATOM EXPLORER ---
        st.subheader("3. Interactive Ligand Atom Sequence")
        lig_atoms = extract_ligand_atoms(active_smiles)
        
        if lig_atoms:
            atom_options = ["None"] + [f"Atom {a['index']}: {a['symbol']}" for a in lig_atoms]
            selected_atom_str = st.selectbox("Click Atom to Highlight in 3D Ligand Pose:", options=atom_options)
            
            if selected_atom_str != "None":
                selected_idx = int(selected_atom_str.split(":")[0].replace("Atom ", ""))
                st.session_state["selected_ligand_atom"] = selected_idx
            else:
                st.session_state["selected_ligand_atom"] = None

    else:
        # Initial Clear Visualizer State
        st.info("Visualizer Canvas clear. Enter SMILES string below and click 'Run Virtual Screening & Docking' to display docking poses.")

    st.markdown("---")

    # --- BOTTOM PANEL: INPUT CONTROLS & DESCRIPTORS ---
    b_col1, b_col2 = st.columns([1, 1])

    with b_col1:
        st.subheader("4. Input Ligand Processing")
        smiles_input = st.text_area("SMILES String Format:", value=st.session_state["smiles"], height=100)
        
        # Red Screening Button
        run_btn = st.button("Run Virtual Screening & Docking", type="primary")
        
        if run_btn:
            st.session_state["smiles"] = smiles_input
            st.session_state["docking_completed"] = True
            st.rerun()

    with b_col2:
        st.subheader("5. Lipinski Physicochemical Descriptors")
        props = calculate_lipinski_properties(st.session_state["smiles"])
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
    
    variants = ["EGFR Wild-Type (1M17)", "EGFR T790M Mutant", "EGFR L858R Mutant"]
    comparison_data = []
    for var in variants:
        score, rmsd = compute_dynamic_binding_score(active_smiles, var, 0)
        comparison_data.append({
            "EGFR Variant Target": var,
            "Top Pose Affinity (kcal/mol)": score,
            "Reference RMSD (Å)": rmsd
        })
        
    df_scores = pd.DataFrame(comparison_data)
    st.dataframe(df_scores.style.highlight_min(subset=["Top Pose Affinity (kcal/mol)"], color="#00f5d4"), width="stretch")
    st.bar_chart(df_scores.set_index("EGFR Variant Target")["Top Pose Affinity (kcal/mol)"])

with tab3:
    st.subheader("Target Structure Metadata")
    st.markdown("""
    * **EGFR Wild-Type (1M17)**: Active kinase domain conformation bound to Erlotinib/Gefitinib class.
    * **EGFR T790M Mutant**: Gatekeeper mutation causing drug resistance against 1st-generation TKIs.
    * **EGFR L858R Mutant**: Activating point mutation with heightened sensitivity to TKIs.
    * **Dynamic Scoring Algorithm**: Docking scores change dynamically based on MW, LogP, rotatable bonds, and specific mutation pocket geometry.
    """) 