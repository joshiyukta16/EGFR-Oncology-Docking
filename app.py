with col2:
        st.subheader("3. Interactive 3D Complex Viewer")
        if "molblock" in st.session_state:
            view = py3Dmol.view(width=550, height=450)
            
            # 1. Render Ligand (Green Sticks)
            view.addModel(st.session_state["molblock"], "mol")
            view.setStyle({'model': -1}, {'stick': {'colorscheme': 'greenCarbon'}})
            
            # 2. Render Receptor Backbone (Spectrum Ribbon/Cartoon)
            if show_receptor:
                receptor_data, file_ext = load_receptor_pdb(target_variant)
                if receptor_data:
                    fmt = "pdbqt" if file_ext == "pdbqt" else "pdb"
                    view.addModel(receptor_data, fmt)
                    view.setStyle({'model': -1}, {'cartoon': {'color': 'spectrum'}})
                else:
                    st.warning(f"Receptor file for '{target_variant}' not found in /receptors directory.")

            view.zoomTo()
            html_content = view._make_html()
            components.html(html_content, height=460, width=560)
        else:
            st.info("Click 'Run Virtual Screening & Docking' to view 3D interaction poses.") 
            