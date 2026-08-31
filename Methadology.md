# Methodology

## Overview
This study was conducted to analyze the binding interaction between Epidermal Growth Factor Receptor (EGFR) and the anticancer drug Gefitinib using molecular docking. The workflow involved receptor and ligand preparation, docking using PyRx (AutoDock Vina), and visualization using PyMOL and Chimera.

---

## 1. Retrieval of Structures

### Receptor (Protein)
- The 3D structures of EGFR were obtained from the Protein Data Bank:
  - Wild-type EGFR (PDB ID: 1M17)
  - Mutant EGFR (PDB ID: 2ITY / 2ITZ*)

- Structures were downloaded in `.pdb` format.

### Ligand
- Gefitinib structure was retrieved from PubChem.
- Downloaded in `.sdf` format.

---

## 2. Receptor Preparation

Protein preparation was carried out to make the structure suitable for docking.

Steps:
- Removed water molecules and heteroatoms
- Removed co-crystallized ligands (if present)
- Added polar hydrogen atoms
- Assigned Kollman charges
- Saved the prepared protein in `.pdbqt` format using AutoDock Tools

---

## 3. Ligand Preparation

Ligand preprocessing ensures proper geometry and compatibility with docking software.

Steps:
- Converted `.sdf` file to `.pdb`
- Performed energy minimization
- Added hydrogens
- Assigned Gasteiger charges
- Converted ligand into `.pdbqt` format

---

## 4. Molecular Docking

Docking was performed using PyRx, which integrates AutoDock Vina.

Steps:
- Imported receptor and ligand into PyRx
- Converted all files to `.pdbqt` format (if not already)
- Defined the grid box around the active/binding site (ATP-binding pocket of EGFR)
- Grid parameters were adjusted to fully cover the binding region
- Docking simulation was executed using AutoDock Vina

---

## 5. Docking Analysis

After docking, multiple binding poses (modes) were generated.

Evaluation criteria:
- Binding affinity (kcal/mol)
- RMSD values (lower bound and upper bound)

Steps:
- Selected the best pose (Mode 0) based on lowest binding energy and RMSD = 0.0
- Compared docking results for wild-type and mutant EGFR

---

## 6. Visualization of Interactions

Visualization was performed using PyMOL and Chimera.

Steps:
- Loaded docked complexes into PyMOL
- Examined ligand orientation within the binding pocket
- Identified interacting amino acid residues
- Generated images showing:
  - Binding pocket
  - Ligand-receptor interactions
  - Comparison between wild-type and mutant structures

---

## 7. Data Representation

Results were organized into tables and figures:

- Binding affinity values for each docking mode
- RMSD values
- Key interacting residues
- Screenshots of docking poses and interaction maps

---

## 8. Interpretation

- The docking results were analyzed to determine binding strength and stability
- Lower binding affinity values indicated stronger interactions
- Comparative analysis between wild-type and mutant EGFR was performed to evaluate differences in drug binding

---

## Notes

- RMSD values were obtained directly from PyRx output

---

## Software Versions (Optional)

- PyRx (AutoDock Vina)
- PyMOL
- UCSF Chimera
- AutoDock Tools 
