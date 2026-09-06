import streamlit as st
import py3Dmol
import pandas as pd
from docking_utils import (
    smiles_to_3d_molblock, 
    generate_standalone_ligand_3d,
    get_all_variant_scores, 
    fetch_smiles_from_pubchem,
    load_receptor_pdb,
    calculate_lipinski_properties,
    get_receptor_sequence,
    compute_dynamic_binding_score
)

st.set_page_config(page_title="EGFR Molecular Docking Dashboard", layout="wide")

DEFAULT_SMILES = "COC1=C(OCCCN2CCOCC2)C=C3C(=C1)N=CN=C3NC4=CC(=C(C=C4)F)Cl"
COVER_IMAGE_URL = "https://images.unsplash.com/photo-1532187863486-abf9dbad1b69?auto=format&fit=crop&w=1200&q=80"

if "docking_ran" not in st.session_state:
    st.session_state["docking_ran"] = False
if "smiles" not in st.session_state:
    st.session_state["smiles"] = DEFAULT_SMILES

if "complex_spin" not in st.session_state:
    st.session_state["complex_spin"] = False
if "complex_dark" not in st.session_state:
    st.session_state["complex_dark"] = False
if "complex_zoom" not in st.session_state:
    st.session_state["complex_zoom"] = 1.0

if "ligand_spin" not in st.session_state:
    st.session_state["ligand_spin"] = False
if "ligand_dark" not in st.session_state:
    st.session_state["ligand_dark"] = False
if "ligand_zoom" not in st.session_state:
    st.session_state["ligand_zoom"] = 1.0
if "highlight_element" not in st.session_state:
    st.session_state["highlight_element"] = "ALL"

st.markdown("""
    <style>
        html { scroll-behavior: smooth; }
        .hero-card {
            background-color: #ffffff;
            border-radius: 12px;
            padding: 24px;
            box-shadow: 0 4px 16px rgba(0,0,0,0.08);
            border: 1px solid #e0e0e0;
            margin-bottom: 25px;
        }
        .hero-title {
            font-size: 36px;
            font-weight: 800;
            color: #1e293b;
            margin-bottom: 8px;
        }
        .hero-subtitle {
            font-size: 16px;
            color: #64748b;
            margin-bottom: 20px;
        }
    </style>
""", unsafe_allow_html=True)

st.markdown('<div class="hero-card">', unsafe_allow_html=True)
hero_col1, hero_col2 = st.columns([1.2, 1])

with hero_col1:
    st.markdown('<div class="hero-title">EGFR Molecular Docking Dashboard</div>', unsafe_allow_html=True)
    st.markdown(
        '<div class="hero-subtitle">Interactive structure-based virtual screening platform targeting '
        'EGFR Wild-Type (1M17), T790M Resistance Mutant (2ITY), and L858R Activating Mutant (9VOX).</div>',
        unsafe_allow_html=True
    )
    
    st.markdown("""
        <a href="#dashboard-main" style="text-decoration: none;">
            <div style="background-color: #ff4b4b; color: white; padding: 14px 28px; border-radius: 8px; 
                        font-weight: bold; font-size: 16px; text-align: center; display: inline-block; 
                        box-shadow: 0 4px 10px rgba(255, 75, 75, 0.3);">
                🚀 Start Docking Analysis
            </div>
        </a>
    """, unsafe_allow_html=True)

with hero_col2:
    st.image(COVER_IMAGE_URL, caption="Protein–Ligand Structural Interaction Interface", width="stretch")

st.markdown('</div>', unsafe_allow_html=True)
st.markdown('<div id="dashboard-main"></div>', unsafe_allow_html=True)
st.markdown("---")

st.sidebar.header("Target Configuration")
target_variant = st.sidebar.selectbox(
    "Select Target EGFR Structure", 
    ["EGFR Wild-Type (1M17)", "EGFR T790M Mutant", "EGFR L858R Mutant"]
)

st.sidebar.markdown("---")
st.sidebar.subheader("Presets / PubChem Search")
search_query = st.sidebar.text_input(
    "Search Drug by Name (e.g. Gefitinib, Erlotinib, Osimertinib):",
    placeholder="e.g. Gefitinib"
)

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
                f'<div style="font-family: monospace; font-size: 14px; line-height: 24px; word-wrap: break-word; background: rgba(227, 242, 253, 0.6); padding: 12px; border-radius: 8px; border: 1px solid #bbdefb; max-height: 120px; overflow-y: auto;">{html_sequence}</div>',
                unsafe_allow_html=True
            )
        else:
            st.info("Receptor structure file missing from /receptors folder.")

    active_smiles = st.session_state.get("smiles", DEFAULT_SMILES)

    if st.session_state["docking_ran"]:
        molblock = smiles_to_3d_molblock(active_smiles, receptor_data)
        
        if receptor_data and molblock:
            view = py3Dmol.view(width=1100, height=480)
            bg_color = "#111111" if st.session_state["complex_dark"] else "#ffffff"
            view.setBackgroundColor(bg_color)
            
            view.addModel(molblock, "mol")
            view.setStyle({'model': 0}, {'stick': {'colorscheme': 'greenCarbon'}})
            
            fmt = "pdbqt" if file_ext == "pdbqt" else "pdb"
            view.addModel(receptor_data, fmt)
            view.setStyle({'model': 1}, {'cartoon': {'color': 'spectrum'}})
            
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
            
            if st.session_state["complex_spin"]:
                view.spin(True)
            else:
                view.spin(False)
            
            view.zoomTo({'model': 0})
            if st.session_state["complex_zoom"] != 1.0:
                view.zoom(st.session_state["complex_zoom"])
                
            html_content = view._make_html()
            st.iframe(html_content, height=490, width=1120)
            
            c_btn1, c_btn2, c_btn3, c_btn4, c_btn5 = st.columns([1.2, 1.2, 1, 1, 1.2])
            with c_btn1:
                spin_lbl = "Pause Spin ⏸️" if st.session_state["complex_spin"] else "Toggle Spin 🔄"
                if st.button(spin_lbl, key="c_spin", width="stretch"):
                    st.session_state["complex_spin"] = not st.session_state["complex_spin"]
                    st.rerun()
            with c_btn2:
                thm_lbl = "Light Mode ☀️" if st.session_state["complex_dark"] else "Dark Mode 🌙"
                if st.button(thm_lbl, key="c_dark", width="stretch"):
                    st.session_state["complex_dark"] = not st.session_state["complex_dark"]
                    st.rerun()
            with c_btn3:
                if st.button("Zoom In 🔍+", key="c_zin", width="stretch"):
                    st.session_state["complex_zoom"] *= 1.2
                    st.rerun()
            with c_btn4:
                if st.button("Zoom Out 🔍-", key="c_zout", width="stretch"):
                    st.session_state["complex_zoom"] *= 0.8
                    st.rerun()
            with c_btn5:
                if st.button("Reset View 🎯", key="c_reset", width="stretch"):
                    st.session_state["complex_zoom"] = 1.0
                    st.session_state["complex_spin"] = False
                    st.rerun()
        else:
            st.error("Error generating 3D complex. Verify structure files and SMILES format.")
    else:
        st.info("Visualizer canvas clean. Click 'Run Virtual Screening & Docking' below to execute docking.")

    st.markdown("---")

    st.subheader("2. Ligand Input & Standalone 3D Conformer Inspection")
    
    l_col1, l_col2 = st.columns([1, 1])

    with l_col1:
        st.caption("ℹ️ *The default sample SMILES format provided below represents **Gefitinib**.*")
        smiles_input = st.text_area("SMILES String Format:", value=active_smiles, height=110)
        
        run_btn = st.button("Run Virtual Screening & Docking", type="primary")
        if run_btn:
            st.session_state["smiles"] = smiles_input
            st.session_state["docking_ran"] = True
            st.rerun()

        if st.session_state["docking_ran"]:
            selected_score = compute_dynamic_binding_score(active_smiles, target_variant)
            st.markdown("---")
            st.metric(label=f"Predicted Binding Energy ({target_variant})", value=f"{selected_score} kcal/mol")

    with l_col2:
        st.markdown("##### 🧬 Standalone 3D Ligand Conformer Visualizer")
        lig_molblock = generate_standalone_ligand_3d(active_smiles)
        
        if lig_molblock:
            l_view = py3Dmol.view(width=540, height=300)
            l_bg = "#111111" if st.session_state["ligand_dark"] else "#ffffff"
            l_view.setBackgroundColor(l_bg)
            
            l_view.addModel(lig_molblock, "mol")
            
            if st.session_state["highlight_element"] == "NITROGEN":
                l_view.setStyle({'elem': 'N'}, {'sphere': {'color': 'blue', 'radius': 0.6}})
                l_view.setStyle({'elem': 'C'}, {'stick': {'colorscheme': 'greenCarbon'}})
            elif st.session_state["highlight_element"] == "OXYGEN":
                l_view.setStyle({'elem': 'O'}, {'sphere': {'color': 'red', 'radius': 0.6}})
                l_view.setStyle({'elem': 'C'}, {'stick': {'colorscheme': 'greenCarbon'}})
            elif st.session_state["highlight_element"] == "HALOGEN":
                l_view.setStyle({'elem': 'Cl'}, {'sphere': {'color': 'green', 'radius': 0.7}})
                l_view.setStyle({'elem': 'F'}, {'sphere': {'color': 'cyan', 'radius': 0.7}})
                l_view.setStyle({'elem': 'C'}, {'stick': {'colorscheme': 'greenCarbon'}})
            else:
                l_view.setStyle({}, {'stick': {'colorscheme': 'greenCarbon'}, 'sphere': {'radius': 0.3}})
            
            l_view.setClickable(
                {},
                True,
                """
                function(atom, viewer, event, container) {
                    if(atom) {
                        var lbl = (atom.elem || 'Atom') + ' (Idx: ' + atom.serial + ')';
                        viewer.addLabel(lbl, {position: {x: atom.x, y: atom.y, z: atom.z}, backgroundColor: 'black', fontColor: 'white'});
                        viewer.render();
                    }
                }
                """
            )
            
            if st.session_state["ligand_spin"]:
                l_view.spin(True)
            else:
                l_view.spin(False)
                
            l_view.zoomTo()
            if st.session_state["ligand_zoom"] != 1.0:
                l_view.zoom(st.session_state["ligand_zoom"])
                
            l_html = l_view._make_html()
            st.iframe(l_html, height=310, width=550)
            
            l_btn1, l_btn2, l_btn3, l_btn4, l_btn5 = st.columns([1, 1, 1, 1, 1])
            with l_btn1:
                l_spin_lbl = "Pause ⏸️" if st.session_state["ligand_spin"] else "Spin 🔄"
                if st.button(l_spin_lbl, key="l_spin", width="stretch"):
                    st.session_state["ligand_spin"] = not st.session_state["ligand_spin"]
                    st.rerun()
            with l_btn2:
                l_thm_lbl = "Light ☀️" if st.session_state["ligand_dark"] else "Dark 🌙"
                if st.button(l_thm_lbl, key="l_dark", width="stretch"):
                    st.session_state["ligand_dark"] = not st.session_state["ligand_dark"]
                    st.rerun()
            with l_btn3:
                if st.button("In 🔍+", key="l_zin", width="stretch"):
                    st.session_state["ligand_zoom"] *= 1.2
                    st.rerun()
            with l_btn4:
                if st.button("Out 🔍-", key="l_zout", width="stretch"):
                    st.session_state["ligand_zoom"] *= 0.8
                    st.rerun()
            with l_btn5:
                if st.button("Reset 🎯", key="l_reset", width="stretch"):
                    st.session_state["ligand_zoom"] = 1.0
                    st.session_state["ligand_spin"] = False
                    st.session_state["highlight_element"] = "ALL"
                    st.rerun()

            st.markdown("<small><b>Interactive Element Filter Badges:</b></small>", unsafe_allow_html=True)
            f_col1, f_col2, f_col3, f_col4 = st.columns(4)
            with f_col1:
                if st.button("All Atoms", width="stretch"):
                    st.session_state["highlight_element"] = "ALL"
                    st.rerun()
            with f_col2:
                if st.button("Nitrogens (N)", width="stretch"):
                    st.session_state["highlight_element"] = "NITROGEN"
                    st.rerun()
            with f_col3:
                if st.button("Oxygens (O)", width="stretch"):
                    st.session_state["highlight_element"] = "OXYGEN"
                    st.rerun()
            with f_col4:
                if st.button("Halogens (F/Cl)", width="stretch"):
                    st.session_state["highlight_element"] = "HALOGEN"
                    st.rerun()

    st.markdown("---")

    st.subheader("3. Lipinski Rule & Comprehensive Chemical Descriptors")
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
    st.subheader("Target Structure Metadata & Structural Notes")
    st.markdown("""
    * **EGFR Wild-Type (1M17)**: Represents the active kinase domain conformation bound to 1st-generation TKIs (Erlotinib/Gefitinib).
    * **EGFR T790M Mutant**: Contains the critical gatekeeper point mutation (Threonine 790 to Methionine), creating steric hindrance causing TKI drug resistance.
    * **EGFR L858R Mutant**: Features the activating mutation in exon 21 (Leucine 858 to Arginine), causing constitutive kinase domain activation.
    """)
    
    st.markdown("---")
    st.subheader("📌 Note & Structural Disclaimer")
    st.info(
        "The three receptor target structures provided in this platform are pre-processed structures initially downloaded "
        "from the Protein Data Bank (PDB). These target structures have been cleaned, optimized, and pre-processed "
        "specifically for molecular docking purposes."
    ) 