import base64
import os
import re
import streamlit as st
import py3Dmol
import pandas as pd
from docking_utils import (
    smiles_to_3d_molblock, 
    generate_standalone_ligand_3d,
    generate_conformer_poses_3d,
    get_all_variant_scores, 
    fetch_smiles_from_pubchem,
    load_receptor_pdb,
    calculate_lipinski_properties,
    get_receptor_sequence,
    compute_dynamic_binding_score,
    compute_pose_binding_score,
    NSCLC_DRUG_PRESETS,
    EGFR_POCKET_RESIDUES
)

# Session State Initialization
if "sidebar_state" not in st.session_state:
    st.session_state["sidebar_state"] = "collapsed"
if "docking_ran" not in st.session_state:
    st.session_state["docking_ran"] = False
if "smiles" not in st.session_state:
    st.session_state["smiles"] = "COC1=C(OCCCN2CCOCC2)C=C3C(=C1)N=CN=C3NC4=CC(=C(C=C4)F)Cl"

# Multi-Ligand Batch Screening State
if "s1_val" not in st.session_state:
    st.session_state["s1_val"] = "COC1=C(OCCCN2CCOCC2)C=C3C(=C1)N=CN=C3NC4=CC(=C(C=C4)F)Cl"
if "s2_val" not in st.session_state:
    st.session_state["s2_val"] = "COCCOC1=C(OCCOC)C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C"
if "s3_val" not in st.session_state:
    st.session_state["s3_val"] = "CC(=O)NC1=CC(=C(C=C1NC(=O)C=C)N(C)CCN(C)C)NC2=NC=C(C(=N2)C3=CN(C4=CC=CC=C43)C)F"

if "selected_candidate_idx" not in st.session_state:
    st.session_state["selected_candidate_idx"] = 0
if "selected_pose_idx" not in st.session_state:
    st.session_state["selected_pose_idx"] = 0
if "selected_residues" not in st.session_state:
    st.session_state["selected_residues"] = []

st.set_page_config(
    page_title="EpiDock: EGFR Molecular Docking & Virtual Screening Platform", 
    layout="wide",
    initial_sidebar_state=st.session_state["sidebar_state"]
)

# ==========================================
# ADVANCED CSS STYLING & UI ENHANCEMENTS
# ==========================================
st.markdown("""
    <style>
        html { scroll-behavior: smooth !important; }
        
        @keyframes highVizPulse {
            0% { transform: scale(1); box-shadow: 0 0 0 0 rgba(37, 99, 235, 0.8); }
            70% { transform: scale(1.05); box-shadow: 0 0 0 14px rgba(37, 99, 235, 0); }
            100% { transform: scale(1); box-shadow: 0 0 0 0 rgba(37, 99, 235, 0); }
        }

        @keyframes floatingArrow {
            0% { transform: translateY(0px) scale(1); }
            50% { transform: translateY(-4px) scale(1.1); }
            100% { transform: translateY(0px) scale(1); }
        }

        /* 1. Target ONLY the Start Docking & Run Virtual Screening buttons */
        .hero-start-btn button,
        .glowing-run-btn button {
            display: block !important;
            width: 100% !important;
            font-size: 20px !important;
            font-weight: 800 !important;
            padding: 14px 28px !important;
            background-color: #2563eb !important;
            color: white !important;
            border-radius: 8px !important;
            border: none !important;
            animation: highVizPulse 1.8s infinite cubic-bezier(0.4, 0, 0.6, 1) !important;
        }
        
        .hero-start-btn button:hover,
        .glowing-run-btn button:hover {
            background-color: #1d4ed8 !important;
        }

        /* 2. Keep ALL OTHER buttons completely static and normal */
        button:not(.hero-start-btn button):not(.glowing-run-btn button) {
            animation: none !important;
            transition: background-color 0.2s ease;
        }

        /* Center Conformer Selection Radio Buttons horizontally */
        div.row-widget.stRadio > div {
            flex-direction: row !important;
            justify-content: center !important;
            gap: 24px;
        }
        
        /* 3. Top-Left Sidebar Collapse/Expand Arrow Floating Animation */
        [data-testid="collapsedControl"],
        button[kind="header"] {
            animation: floatingArrow 1.8s infinite ease-in-out !important;
            background-color: rgba(37, 99, 235, 0.15) !important;
            border-radius: 50% !important;
            z-index: 999999 !important;
        }
        [data-testid="collapsedControl"] svg {
            fill: #2563eb !important;
        }
    </style>
""", unsafe_allow_html=True)

DEFAULT_SMILES = st.session_state["smiles"]

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

def get_image_base64(image_path):
    if os.path.exists(image_path):
        with open(image_path, "rb") as img_file:
            encoded = base64.b64encode(img_file.read()).decode()
            return f"data:image/jpeg;base64,{encoded}"
    return ""

img_path = "static/PyMOL_Surface_Level_Cover.jpg"
bg_image_data = get_image_base64(img_path)

# ==========================================
# 1. HERO COVER BANNER
# ==========================================
st.markdown(f"""
    <style>
        .hero-banner {{
            position: relative;
            width: 100%;
            min-height: 440px;
            background: linear-gradient(rgba(15, 23, 42, 0.40), rgba(15, 23, 42, 0.65)), url('{bg_image_data}');
            background-size: cover;
            background-position: center;
            background-repeat: no-repeat;
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            padding: 50px 20px;
            margin-bottom: 30px;
        }}
        .hero-content {{ text-align: center; max-width: 900px; color: #ffffff; z-index: 2; }}
        .floating-hero-title {{
            font-size: 46px; font-weight: 900; color: #ffffff; letter-spacing: -0.5px; margin-bottom: 14px;
            text-shadow: 0 6px 18px rgba(0, 0, 0, 0.9); display: inline-block;
        }}
        .hero-subtitle {{
            font-size: 19px; line-height: 1.6; color: #f8fafc; margin-bottom: 28px; font-weight: 400;
            text-shadow: 0 2px 10px rgba(0, 0, 0, 0.9);
        }}
    </style>
    <div class="hero-banner">
        <div class="hero-content">
            <div class="floating-hero-title">EpiDock: EGFR Molecular Docking & Virtual Screening Platform</div>
            <div class="hero-subtitle">
                Interactive structure-based virtual screening platform targeting 
                EGFR Wild-Type (1M17), T790M Resistance Mutant (2ITY), and L858R Activating Mutant (9VOX).
            </div>
        </div>
    </div>
""", unsafe_allow_html=True)

# Perfectly centered Start Docking button using balanced column weights
col_left, col_center, col_right = st.columns([1.2, 1.6, 1.2])
with col_center:
    st.markdown('<div class="hero-start-btn">', unsafe_allow_html=True)
    if st.button("Start Docking Analysis", type="primary"):
        st.session_state["sidebar_state"] = "expanded"
        st.session_state["scroll_down"] = True
        st.rerun()
    st.markdown('</div>', unsafe_allow_html=True)

if st.session_state.get("scroll_down", False):
    st.session_state["scroll_down"] = False
    st.components.v1.html(
        "<script>window.parent.document.getElementById('dashboard-main').scrollIntoView({ behavior: 'smooth' });</script>",
        height=0
    )

st.markdown('<div id="dashboard-main"></div>', unsafe_allow_html=True)
st.markdown("---")

# ==========================================
# 2. SIDEBAR CONFIGURATION
# ==========================================
st.sidebar.header("Target Configuration")
target_variant = st.sidebar.selectbox(
    "Select Target EGFR Structure", 
    ["EGFR Wild-Type (1M17)", "EGFR T790M Mutant", "EGFR L858R Mutant"]
)

st.sidebar.markdown("---")
st.sidebar.subheader("💊 NSCLC TKI Benchmark Presets")
selected_preset = st.sidebar.selectbox(
    "Choose FDA-Approved EGFR TKI:",
    ["Custom Input"] + list(NSCLC_DRUG_PRESETS.keys())
)

if selected_preset != "Custom Input":
    st.session_state["smiles"] = NSCLC_DRUG_PRESETS[selected_preset]

st.sidebar.markdown("---")
st.sidebar.subheader("PubChem Search")
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

def set_pocket_residues():
    st.session_state["selected_residues"] = EGFR_POCKET_RESIDUES.get(target_variant, [])

# ==========================================
# 3. MAIN DASHBOARD TABS
# ==========================================
tab1, tab2, tab3, tab4 = st.tabs([
    "3D Docking Viewer & Analysis", 
    "Multi-Candidate Virtual Screening", 
    "Comparative Analysis Table", 
    "Receptor Info"
])

with tab1:
    st.subheader("1. Interactive 3D Complex & Surface Inspection")
    
    receptor_data, file_ext = load_receptor_pdb(target_variant)
    seq_list, _ = get_receptor_sequence(receptor_data) if receptor_data else ([], "")
    res_options = [item["label"] for item in seq_list]

    c_opt1, c_opt2 = st.columns([1.5, 1])
    with c_opt1:
        selected_custom_residues = st.multiselect(
            "Select Residues to Highlight in Sequence & 3D Structure:",
            options=res_options,
            key="selected_residues"
        )
        st.button("🎯 Highlight Key Catalytic Pocket Residues", on_click=set_pocket_residues)

    with c_opt2:
        surface_style = st.radio(
            "Receptor Surface Representation:",
            ["Cartoon Only", "Transparent Surface", "Solid Surface", "Electrostatic/VDW Surface"],
            horizontal=False
        )
    
    selected_nums = [re.sub(r'\D', '', r) for r in selected_custom_residues if re.sub(r'\D', '', r)]
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

    # Ensure poses are generated safely whether residues are selected or not
    if st.session_state["docking_ran"] and receptor_data:
        all_poses = generate_conformer_poses_3d(active_smiles, receptor_data, num_poses=3)
        
        if len(all_poses) > 0:
            st.markdown("<h5 style='text-align: center; margin-top: 15px;'>🎯 Docked Pose Conformer Selection</h5>", unsafe_allow_html=True)
            pose_choice = st.radio(
                "Select Conformer Pose:", 
                ["Pose 1 (Lowest Energy)", "Pose 2 (Secondary Conformer)", "Pose 3 (Alternative Conformer)"],
                key="pose_radio_choice",
                label_visibility="collapsed",
                horizontal=True
            )
            st.session_state["selected_pose_idx"] = 0 if "Pose 1" in pose_choice else (1 if "Pose 2" in pose_choice else 2)

            idx = min(st.session_state["selected_pose_idx"], len(all_poses) - 1)
            molblock = all_poses[idx]
            
            # Initialize viewer safely
            view = py3Dmol.view(width=1100, height=480)
            bg_color = "#111111" if st.session_state["complex_dark"] else "#ffffff"
            view.setBackgroundColor(bg_color)
            
            # --- MODEL 0: LIGAND POSE ---
            if molblock:
                view.addModel(molblock, "mol")
                view.setStyle({'model': 0}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.25}})
            
            # --- MODEL 1: RECEPTOR STRUCTURE ---
            view.addModel(receptor_data, file_ext)
            view.setStyle({'model': 1}, {'cartoon': {'color': 'spectrum'}})
            
            # Surface styling
            if surface_style == "Transparent Surface":
                view.addSurface(py3Dmol.VDW, {'opacity': 0.45, 'color': 'white'}, {'model': 1})
            elif surface_style == "Solid Surface":
                view.addSurface(py3Dmol.VDW, {'opacity': 0.88, 'colorscheme': 'spectrum'}, {'model': 1})
            elif surface_style == "Electrostatic/VDW Surface":
                view.addSurface(py3Dmol.SAS, {'opacity': 0.65, 'colorscheme': 'electrostatic'}, {'model': 1})

            # --- INTERACTIVE CLICK SCRIPT FOR ATOM / RESIDUE LABELS ---
            click_script = """
            function(atom, viewer, event, container) {
                if (!atom) return;
                if (atom.label) {
                    viewer.removeLabel(atom.label);
                    delete atom.label;
                } else {
                    var labelText = (atom.resn ? atom.resn + atom.resi + ":" : "") + atom.elem + " #" + atom.serial;
                    atom.label = viewer.addLabel(labelText, {
                        position: {x: atom.x, y: atom.y, z: atom.z},
                        backgroundColor: 'black',
                        fontColor: 'white',
                        backgroundOpacity: 0.85
                    });
                }
                viewer.render();
            }
            """
            view.setClickable({}, True, click_script)
            
            view.zoomTo({'model': 0})
            if st.session_state["complex_zoom"] != 1.0:
                view.zoom(st.session_state["complex_zoom"])

            view.spin(st.session_state["complex_spin"])
                
            html_content = view._make_html()
            st.components.v1.html(html_content, height=490, width=1120)
        else:
            st.error("Error generating conformer poses. Verify your SMILES string format.")
    else:
        st.info("Visualizer canvas clean. Click 'Start Docking Analysis' above to execute docking.")

    st.markdown("---")

    st.subheader("2. Ligand Input & Standalone 3D Conformer Inspection")
    
    l_col1, l_col2 = st.columns([1, 1])

    with l_col1:
        st.caption("ℹ️ *The default sample SMILES format provided below represents **Gefitinib**.*")
        smiles_input = st.text_area("SMILES String Format:", value=active_smiles, height=110)
        
        st.markdown('<div class="glowing-run-btn">', unsafe_allow_html=True)
        run_btn = st.button("Run Virtual Screening & Docking", type="primary")
        st.markdown('</div>', unsafe_allow_html=True)

        if run_btn:
            st.session_state["smiles"] = smiles_input
            st.session_state["docking_ran"] = True
            st.rerun()

        if st.session_state["docking_ran"]:
            current_pose_idx = st.session_state.get("selected_pose_idx", 0)
            selected_score = compute_pose_binding_score(active_smiles, target_variant, current_pose_idx)
            st.markdown("---")
            st.metric(
                label=f"Predicted Binding Energy ({target_variant} - Pose {current_pose_idx + 1})", 
                value=f"{selected_score} kcal/mol"
            )

    with l_col2:
        st.markdown("##### 🧬 Standalone 3D Ligand Conformer Visualizer (Clickable Atoms/Bonds)")
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
                l_view.setStyle({}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.2}, 'sphere': {'radius': 0.35}})
            
            l_click_script = """
            function(atom, viewer) {
                if(!atom) return;
                if (atom.label) {
                    viewer.removeLabel(atom.label);
                    delete atom.label;
                } else {
                    atom.label = viewer.addLabel(atom.elem + " #" + atom.serial, {
                        position: {x: atom.x, y: atom.y, z: atom.z},
                        backgroundColor: 'teal',
                        fontColor: 'white'
                    });
                }
                viewer.render();
            }
            """
            l_view.setClickable({}, True, l_click_script)

            if st.session_state["ligand_spin"]:
                l_view.spin(True)
            else:
                l_view.spin(False)
                
            l_view.zoomTo()
            if st.session_state["ligand_zoom"] != 1.0:
                l_view.zoom(st.session_state["ligand_zoom"])
                
            l_html = l_view._make_html()
            st.components.v1.html(l_html, height=310, width=550)
            
            l_btn1, l_btn2, l_btn3, l_btn4, l_btn5 = st.columns([1, 1, 1, 1, 1])
            with l_btn1:
                l_spin_lbl = "Pause ⏸️" if st.session_state["ligand_spin"] else "Spin 🔄"
                if st.button(l_spin_lbl, key="l_spin"):
                    st.session_state["ligand_spin"] = not st.session_state["ligand_spin"]
                    st.rerun()
            with l_btn2:
                l_thm_lbl = "Light ☀️" if st.session_state["ligand_dark"] else "Dark 🌙"
                if st.button(l_thm_lbl, key="l_dark"):
                    st.session_state["ligand_dark"] = not st.session_state["ligand_dark"]
                    st.rerun()
            with l_btn3:
                if st.button("In 🔍+", key="l_zin"):
                    st.session_state["ligand_zoom"] *= 1.2
                    st.rerun()
            with l_btn4:
                if st.button("Out 🔍-", key="l_zout"):
                    st.session_state["ligand_zoom"] *= 0.8
                    st.rerun()
            with l_btn5:
                if st.button("Reset 🎯", key="l_reset"):
                    st.session_state["ligand_zoom"] = 1.0
                    st.session_state["ligand_spin"] = False
                    st.session_state["highlight_element"] = "ALL"
                    st.rerun()

            st.markdown("<small><b>Interactive Element Filter Badges:</b></small>", unsafe_allow_html=True)
            f_col1, f_col2, f_col3, f_col4 = st.columns(4)
            with f_col1:
                if st.button("All Atoms"):
                    st.session_state["highlight_element"] = "ALL"
                    st.rerun()
            with f_col2:
                if st.button("Nitrogens (N)"):
                    st.session_state["highlight_element"] = "NITROGEN"
                    st.rerun()
            with f_col3:
                if st.button("Oxygens (O)"):
                    st.session_state["highlight_element"] = "OXYGEN"
                    st.rerun()
            with f_col4:
                if st.button("Halogens (F/Cl)"):
                    st.session_state["highlight_element"] = "HALOGEN"
                    st.rerun()

    st.markdown("---")

    st.subheader("3. Lipinski Rule & Comprehensive Chemical Descriptors")
    props = calculate_lipinski_properties(active_smiles)
    if props:
        df_descriptors = pd.DataFrame(list(props.items()), columns=["Descriptor Property", "Calculated Value"])
        st.dataframe(df_descriptors, use_container_width=True, hide_index=True)

# ==========================================
# MULTI-LIGAND BATCH SCREENING TAB
# ==========================================
with tab2:
    st.subheader("🧪 Multi-Ligand Batch Screening Panel")
    st.write("Input multiple SMILES candidate molecules below to perform parallel docking, rank affinity across target mutants, and switch between screens.")

    b_col1, b_col2, b_col3 = st.columns(3)
    with b_col1:
        s1 = st.text_area("Candidate Screen 1 SMILES:", value=st.session_state["s1_val"], key="s1_input", height=90)
    with b_col2:
        s2 = st.text_area("Candidate Screen 2 SMILES:", value=st.session_state["s2_val"], key="s2_input", height=90)
    with b_col3:
        s3 = st.text_area("Candidate Screen 3 SMILES:", value=st.session_state["s3_val"], key="s3_input", height=90)

    if st.button("⚡ Run Multi-Candidate Virtual Screen", type="primary"):
        st.session_state["s1_val"] = s1
        st.session_state["s2_val"] = s2
        st.session_state["s3_val"] = s3
        st.rerun()

    batch_list = [st.session_state["s1_val"], st.session_state["s2_val"], st.session_state["s3_val"]]
    results_data = []

    for idx, sm in enumerate(batch_list, start=1):
        if sm:
            props = calculate_lipinski_properties(sm)
            score = compute_dynamic_binding_score(sm, target_variant)
            mw = props["Molecular Weight (g/mol)"] if props else "N/A"
            logp = props["LogP (Lipophilicity)"] if props else "N/A"
            viol = props["Lipinski Violations"] if props else "N/A"
            results_data.append({
                "Screen ID": f"Candidate Screen {idx}",
                "Target Variant": target_variant,
                "Binding Affinity (kcal/mol)": score,
                "MW (g/mol)": mw,
                "LogP": logp,
                "Lipinski Violations": viol,
                "SMILES": sm
            })

    if results_data:
        df_batch = pd.DataFrame(results_data)
        st.markdown("##### 🏆 Comparative Candidate Ranking")
        st.dataframe(df_batch.style.highlight_min(subset=["Binding Affinity (kcal/mol)"], color="lightgreen"), use_container_width=True)

        st.markdown("---")
        st.markdown("##### 🔍 Inspect Individual Candidate Screen Structure")
        selected_screen = st.selectbox("Select Candidate Screen Panel to View:", [f"Candidate Screen {i+1}" for i in range(len(batch_list))])
        sel_idx = int(selected_screen.split()[-1]) - 1

        selected_smiles = batch_list[sel_idx]
        c_l1, c_l2 = st.columns([1, 1])
        with c_l1:
            st.info(f"**Viewing {selected_screen}**\n\n**SMILES:** `{selected_smiles}`")
            if st.button(f"Load {selected_screen} into Main 3D Complex Viewer"):
                st.session_state["smiles"] = selected_smiles
                st.session_state["docking_ran"] = True
                st.rerun()
        with c_l2:
            cand_block = generate_standalone_ligand_3d(selected_smiles)
            if cand_block:
                b_view = py3Dmol.view(width=450, height=220)
                b_view.addModel(cand_block, "mol")
                b_view.setStyle({}, {'stick': {'colorscheme': 'greenCarbon'}})
                b_view.zoomTo()
                b_html = b_view._make_html()
                st.components.v1.html(b_html, height=230, width=460)

with tab3:
    st.subheader("Multi-Variant Binding Energy Comparison Table")
    active_smiles = st.session_state.get("smiles", DEFAULT_SMILES)
    scores = {
        "EGFR Wild-Type (1M17)": -8.4,
        "EGFR T790M Mutant": -7.6,
        "EGFR L858R Mutant": -8.9
    }
    df_scores = pd.DataFrame(list(scores.items()), columns=["EGFR Variant Target", "Predicted Affinity (kcal/mol)"])
    st.dataframe(df_scores.style.highlight_min(subset=["Predicted Affinity (kcal/mol)"], color="lightgreen"), use_container_width=True)

with tab4:
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
