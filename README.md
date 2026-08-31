# Binding And Beyond Docking Molecular study on EGFR mutations in lung cancer-  

## 📄 License 
This repository has been created for academic and portfolio purposes. Commercial use is not permitted.  
The project report, documentation, figures, and other original content are licensed under the **Creative Commons Attribution–NonCommercial 4.0 International (CC BY-NC 4.0)** License. 

See the [LICENSE](LICENSE) file for the complete license terms.

--- 

## Binding & Beyond: Molecular Docking of Gefitinib with EGFR

## 🧬 Project Overview
This project focuses on studying the interaction between Epidermal Growth Factor Receptor (EGFR) and the anticancer drug Gefitinib using computational approaches, including molecular descriptor analysis and molecular docking. 

EGFR plays a critical role in non-small cell lung cancer (NSCLC). Mutations in EGFR can alter drug-binding behavior and influence therapeutic response, making it important to study ligand–receptor interactions computationally. 

Gefitinib was analyzed using RDKit to evaluate its molecular descriptors and drug-like properties based on Lipinski's Rule of Five. Molecular docking was performed using PyRx (AutoDock Vina) to predict binding affinity between Gefitinib and EGFR structures. 

Protein–ligand interactions were analyzed and visualized using PyMOL and UCSF Chimera. 

--- 

## 🎯 Objectives

- Retrieve EGFR protein structures (wild-type and mutant)
- Prepare receptor and ligand structures
- Perform The Descriptor Calculations using RDKit 
- Perform molecular docking using PyRx
- Analyze binding affinity and RMSD
- Visualize ligand-receptor interactions
- Compare wild-type and mutant EGFR binding

--- 

```mermaid
flowchart TD

A[Target Selection EGFR]
--> B[Protein Retrieval PDB]

B --> C[Ligand Retrieval PubChem]

C --> D[Ligand Preparation]

D --> E[Molecular Descriptor Analysis RDKit<br>
Molecular Weight<br>
LogP<br>
TPSA<br>
H-Bond Donors<br>
H-Bond Acceptors<br>
Rotatable Bonds<br>
Ring Count<br>
Heavy Atom Count<br>
Fraction Csp3]

E --> F[Drug-Likeness Assessment<br>
Lipinski's Rule of Five]

F --> G[Receptor Preparation]

G --> H[Molecular Docking<br>
PyRx / AutoDock Vina]

H --> I[Binding Affinity Analysis]

I --> J[Protein-Ligand Interaction Analysis]

J --> K[Molecular Visualization<br>
PyMOL<br>
UCSF Chimera]

K --> L[Result Interpretation & Conclusion]
```

--- 

## 🛠️ Tools & Software Used

- Protein Data Bank (PDB)
- PubChem
- RDKit
- Python
- VS Codes  
- PyRx (AutoDock Vina)
- AutoDock Tools
- Open Babel 
- PyMOL 

---

## Project Repository Structure 

```text 
📁  Binding_And_Beyond_Docking_Molecular_study_on_EGFR_mutations_in_lung_cancer
│
├── 📄 README.md
├── 📄 LICENSE
│
├── 📘 Project_Report_BI.pdf
├── 📑 Methodology.md
├── 📊 results-summary.md
├── 🧪 Lung_Cancer_Molecular_Docking.pdf
│
├── 📁 RDKit
│   │
│   ├── 📘 Gefitinib_RDKit.ipynb
│   │      └── Complete RDKit workflow:
│   │          - SDF file loading
│   │          - Molecular descriptor calculation
│   │          - Lipinski Rule of Five
│   │          - Morgan fingerprint generation
│   │          - CSV export
│   │          - 2D structure generation
│   │
│   ├── 📁 data
│   │   └── 🧬 Gefitinib.sdf
│   │          └── PubChem ligand input file
│   │
│   ├── 🧾 Gefitinib_Descriptors.csv
│   │          └── RDKit calculated molecular properties
│   │
│   └── 📁 Output
│       └── 🖼️ Gefitinib_2D.png
│              └── RDKit generated 2D molecular structure
│ 
├── 🖼️ Open_Babel.png
├── 🖼️ PyMOL.jpg
└── 🖼️ PyMOL_Surface_Level.jpg 

``` 
--- 

## ⚙️ Methodology

### 1. Receptor Preparation
- Downloaded EGFR structures from PDB
- Removed water molecules and unwanted chains

<img width="599" height="400" alt="image" src="https://github.com/user-attachments/assets/1f70a17d-1770-4788-8d0a-826eaf1fca39" />
  
- Added hydrogen atoms
  
<img width="599" height="400" alt="image" src="https://github.com/user-attachments/assets/f1d9cb07-745f-44c9-b286-e10cc17bdc92" />

- Converted to `.pdbqt` format

   ![Open Babel](Open_Babel.png) 

### 2. Ligand Preparation
- Downloaded Gefitinib from PubChem

<img width="599" height="400" alt="image" src="https://github.com/user-attachments/assets/caa1cd19-7766-4c5a-b115-c3c4891d3e8c" />

- Converted SDF → PDB → PDBQT
- Energy minimization performed 

### 3. Molecular Descriptor Analysis (RDKit)

The Gefitinib ligand structure was retrieved from PubChem in **SDF format** and analyzed using **RDKit**, an open-source cheminformatics library.

The molecule was loaded using RDKit's `Chem.SDMolSupplier()` function, followed by molecular descriptor calculation, drug-likeness evaluation, molecular fingerprint generation, and 2D structure visualization.

📘 [RDKit Analysis Notebook](RDKit/Gefitnib_RDKit.ipynb)  

### RDKit Molecular Descriptor Calculation

The following molecular descriptors were calculated:

- Molecular Weight
- LogP (Lipophilicity)
- Topological Polar Surface Area (TPSA)
- Hydrogen Bond Donors (HBD)
- Hydrogen Bond Acceptors (HBA)
- Rotatable Bonds
- Ring Count
- Heavy Atom Count
- Fraction Csp3


### RDKit Descriptor Results

| Descriptor | Value |
|------------|-------|
| Molecular Weight (Da) | 446.91 |
| LogP | 4.28 |
| TPSA (Å²) | 68.74 |
| Hydrogen Bond Donors | 1 |
| Hydrogen Bond Acceptors | 7 |
| Rotatable Bonds | 8 |
| Ring Count | 4 |
| Heavy Atom Count | 31 |
| Fraction Csp3 | 0.36 |


The complete descriptor output generated using RDKit is available:

📄 [Gefitinib_Descriptors.csv](RDKit/Gefitinib_Descriptors.csv)


### Lipinski's Rule of Five Assessment

The calculated descriptors were evaluated using Lipinski's Rule of Five:

| Parameter | Gefitinib Value | Rule | Result |
|-----------|----------------|------|--------|
| Molecular Weight | 446.91 Da | < 500 Da | Pass |
| LogP | 4.28 | < 5 | Pass |
| H-Bond Donors | 1 | ≤ 5 | Pass |
| H-Bond Acceptors | 7 | ≤ 10 | Pass |


Gefitinib satisfies all Lipinski criteria, indicating favorable drug-like properties.


### Molecular Fingerprint Generation

Morgan fingerprints were generated using RDKit to represent the molecular structure as a binary fingerprint vector.

The fingerprint representation can be used for:
- Molecular similarity analysis
- Chemical clustering
- Machine learning applications


### 2D Molecular Structure Visualization

The 2D structure of Gefitinib was generated using RDKit and saved as an image.

<img src="RDKit/Output/Gefitinib_2D.png" alt="Gefitinib 2D Structure" width="599"/> 

The generated RDKit 2D structure image is available here: [Gefitinib_2D.png](RDKit/Output/Gefitinib_2D.png) 

### 4. Docking (PyRx) 
- Loaded receptor and ligand

<img width="599" height="400" alt="image" src="https://github.com/user-attachments/assets/91cb3d7b-b4ac-4aab-8e45-c0d0e2fc6101" /> 

- Defined grid box around active site
  
<img width="599" height="400" alt="image" src="https://github.com/user-attachments/assets/6726e0bd-9b3b-423c-b819-38105453dcac" /> 
  
- Performed docking using AutoDock Vina 

### 5. Analysis
- Selected best pose (Mode 0)
- Evaluated:
  - Binding affinity (kcal/mol)
  - RMSD values
- Visualized interactions using PyMOL 

---

## 📊 Results

| Receptor | Binding Affinity (kcal/mol) | Key Residues |
|----------|----------------------------|-------------|
| Wild-type EGFR | -5.8 | Met793, Leu718, Lys745 |
| Mutant EGFR | -6.5 | Met793, Leu858, Thr790 |

### Key Observations
- Mutant EGFR shows stronger binding with Gefitinib
- Lower binding energy indicates higher stability
- Critical residues involved in binding were identified

---

## 🧪 Visualization

### Best Docked Pose (Mutant) 

![Mutant](PyMOL.jpg) 

### Best Docked Pose (Wild) 

![Wild](PyMOL_Surface_Level.jpg) 

### Color Coding for Clarity:
To differentiate interaction types:
- Hydrogen bonds → yellow dashed lines
-	Hydrophobic residues → green or grey surface
-	Aromatic/π-interactions → magenta sticks 


---

## 🎥 Project Report 

A detailed Project Report showcasing Methodology and Results is available in this repository.

**📄 Project Report:** [View Project Report](Project_Report_BI%20.pdf) 

--- 

## 📌 Key Learnings

- Molecular docking helps predict drug-protein interactions
- EGFR mutations significantly impact drug binding
- Visualization tools are essential for interpreting docking results

---

## 📚 References

- Protein Data Bank (PDB)
- PubChem Database
- AutoDock Vina Documentation

---

## 👩‍🔬 Author

Yukta Joshi  
Bioinformatics | Computational Biology | AI in Healthcare 

🙏 Acknowledgment

Data provided by PDB And PubChem database. 

⭐ If you like this project, give it a star!!  

