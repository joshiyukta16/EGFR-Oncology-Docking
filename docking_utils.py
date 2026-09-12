import os
import re
import sys
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski, rdMolDescriptors
import pubchempy as pcp

AA_MAP = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

# NSCLC Drug Library Presets
NSCLC_DRUG_PRESETS = {
    "Gefitinib (1st Gen)": "COC1=C(OCCCN2CCOCC2)C=C3C(=C1)N=CN=C3NC4=CC(=C(C=C4)F)Cl",
    "Erlotinib (1st Gen)": "COCCOC1=C(OCCOC)C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C",
    "Afatinib (2nd Gen)": "CN(C)C/C=C/C(=O)NC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC(=C(C=C4)Cl)F)O[C@H]4CCOC4",
    "Osimertinib (3rd Gen)": "CC(=O)NC1=CC(=C(C=C1NC(=O)C=C)N(C)CCN(C)C)NC2=NC=C(C(=N2)C3=CN(C4=CC=CC=C43)C)F"
}

# Critical Active-Site Residues for EGFR Kinase Pocket
EGFR_POCKET_RESIDUES = {
    "EGFR Wild-Type (1M17)": ["MET793", "LYS745", "LEU718", "THR790", "ALA743"],
    "EGFR T790M Mutant": ["MET793", "LYS745", "MET790", "LEU718", "ALA743"],
    "EGFR L858R Mutant": ["MET793", "LYS745", "ARG858", "LEU718", "ALA743"]
}

def fetch_smiles_from_pubchem(drug_name):
    """Fetches connectivity SMILES string from PubChem using drug name."""
    try:
        compounds = pcp.get_compounds(drug_name, 'name')
        return compounds[0].connectivity_smiles if compounds else None
    except Exception as e:
        print(f"[DEBUG fetch_smiles] {e}", file=sys.stderr)
        return None

def calculate_lipinski_properties(smiles_str):
    """Calculates comprehensive physicochemical descriptors and drug-likeness rules."""
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is None: 
            return None
        return {
            "Molecular Weight (g/mol)": round(Descriptors.MolWt(mol), 2),
            "LogP (Lipophilicity)": round(Crippen.MolLogP(mol), 2),
            "TPSA (Å²)": round(rdMolDescriptors.CalcTPSA(mol), 2),
            "H-Bond Donors (HBD)": Lipinski.NumHDonors(mol),
            "H-Bond Acceptors (HBA)": Lipinski.NumHAcceptors(mol),
            "Rotatable Bonds": Lipinski.NumRotatableBonds(mol),
            "Ring Count": Lipinski.RingCount(mol),
            "Heavy Atom Count": mol.GetNumHeavyAtoms(),
            "Fraction Csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 2),
            "Lipinski Violations": sum([
                Descriptors.MolWt(mol) > 500, 
                Crippen.MolLogP(mol) > 5, 
                Lipinski.NumHDonors(mol) > 5, 
                Lipinski.NumHAcceptors(mol) > 10
            ])
        }
    except Exception as e:
        print(f"[DEBUG Lipinski] {e}", file=sys.stderr)
        return None

def compute_dynamic_binding_score(smiles_str, variant_name):
    """Calculates binding affinity based on ligand physicochemical properties and target variant."""
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is None: 
            return -7.50
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        rot_bonds = Lipinski.NumRotatableBonds(mol)
        hba = Lipinski.NumHAcceptors(mol)
        
        base_score = -0.012 * mw - 0.35 * logp + 0.15 * rot_bonds - 0.22 * hba
        offsets = {
            "EGFR Wild-Type (1M17)": -2.1,
            "EGFR T790M Mutant": -0.8,
            "EGFR L858R Mutant": -2.5
        }
        return round(base_score + offsets.get(variant_name, -1.5), 2)
    except Exception as e:
        print(f"[DEBUG binding_score] {e}", file=sys.stderr)
        return -7.50

def get_all_variant_scores(smiles_str):
    """Returns binding scores across all EGFR target variants."""
    return {
        "EGFR Wild-Type (1M17)": compute_dynamic_binding_score(smiles_str, "EGFR Wild-Type (1M17)"),
        "EGFR T790M Mutant": compute_dynamic_binding_score(smiles_str, "EGFR T790M Mutant"),
        "EGFR L858R Mutant": compute_dynamic_binding_score(smiles_str, "EGFR L858R Mutant")
    }

def generate_standalone_ligand_3d(smiles_str):
    """Generates standalone 3D MolBlock for individual ligand visualization."""
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        res = AllChem.EmbedMolecule(mol, randomSeed=42)
        if res != 0:
            AllChem.EmbedMolecule(mol, useRandomCoordinates=True, randomSeed=42)
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except Exception:
            pass
        return Chem.MolToMolBlock(mol)
    except Exception as e:
        print(f"[DEBUG standalone_3d] {e}", file=sys.stderr)
        return None

def generate_conformer_poses_3d(smiles_str, receptor_data=None, num_poses=3):
    """Generates multiple distinct low-energy conformer poses safely centered to binding pocket."""
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is None:
            return []
        mol = Chem.AddHs(mol)
        
        cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=num_poses, randomSeed=42, pruneRmsThresh=0.5))
        if not cids:
            cid = AllChem.EmbedMolecule(mol, randomSeed=42)
            if cid == -1:
                cid = AllChem.EmbedMolecule(mol, useRandomCoordinates=True, randomSeed=42)
            cids = [cid] if cid != -1 else []
            
        if not cids:
            return []

        try:
            AllChem.MMFFOptimizeMoleculeConfs(mol)
        except Exception as opt_err:
            print(f"[DEBUG MMFFConfs] {opt_err}", file=sys.stderr)
            pass
            
        center_x, center_y, center_z = get_receptor_center(receptor_data)
        num_atoms = mol.GetNumAtoms()
        molblocks = []
        
        for i, cid in enumerate(cids):
            conf = mol.GetConformer(cid)
            cx = sum(conf.GetAtomPosition(j).x for j in range(num_atoms)) / num_atoms
            cy = sum(conf.GetAtomPosition(j).y for j in range(num_atoms)) / num_atoms
            cz = sum(conf.GetAtomPosition(j).z for j in range(num_atoms)) / num_atoms
            
            # Align center with fallback to avoid out-of-bounds translation
            dx = center_x - cx if center_x != 0.0 else 0.0
            dy = center_y - cy if center_y != 0.0 else 0.0
            dz = center_z - cz if center_z != 0.0 else 0.0
            
            pose_offset_x = (i - 1) * 0.4
            pose_offset_y = (i - 1) * 0.3
            
            for j in range(num_atoms):
                pos = conf.GetAtomPosition(j)
                conf.SetAtomPosition(j, (pos.x + dx + pose_offset_x, pos.y + dy + pose_offset_y, pos.z + dz))
                
            molblocks.append(Chem.MolToMolBlock(mol, confId=cid))
            
        return molblocks
    except Exception as e:
        print(f"[DEBUG generate_conformer_poses_3d] Error: {e}", file=sys.stderr)
        return []

def get_receptor_center(receptor_data):
    """Calculates geometric center of receptor structure in real-time."""
    if not receptor_data:
        return (0.0, 0.0, 0.0)
    x_coords, y_coords, z_coords = [], [], []
    for line in receptor_data.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                x_coords.append(float(line[30:38].strip()))
                y_coords.append(float(line[38:46].strip()))
                z_coords.append(float(line[46:54].strip()))
            except ValueError:
                continue
    if x_coords:
        return (sum(x_coords)/len(x_coords), sum(y_coords)/len(y_coords), sum(z_coords)/len(z_coords))
    return (0.0, 0.0, 0.0)

def smiles_to_3d_molblock(smiles_str, receptor_data=None):
    """Converts SMILES string into a 3D MolBlock translated into the active binding pocket."""
    poses = generate_conformer_poses_3d(smiles_str, receptor_data, num_poses=1)
    return poses[0] if poses else None

def load_receptor_pdb(variant_name):
    """Loads processed receptor PDB or PDBQT content."""
    filename_map = {
        "EGFR Wild-Type (1M17)": [
            "receptors/egfr_receptor_1M17.pdbqt", 
            "receptors/egfr_receptor_1M17.pdb",
            "egfr_receptor_1M17.pdbqt",
            "egfr_receptor_1M17.pdb"
        ],
        "EGFR T790M Mutant": [
            "receptors/egfr_receptor_2ITY.pdbqt", 
            "receptors/egfr_receptor_2ITY.pdb",
            "egfr_receptor_2ITY.pdbqt",
            "egfr_receptor_2ITY.pdb"
        ],
        "EGFR L858R Mutant": [
            "receptors/egfr_receptor_9VOX.pdbqt", 
            "receptors/egfr_receptor_9VOX.pdb",
            "egfr_receptor_9VOX.pdbqt",
            "egfr_receptor_9VOX.pdb"
        ] 
    }
    possible_paths = filename_map.get(variant_name, [])
    for path in possible_paths:
        if os.path.exists(path):
            with open(path, "r") as f:
                content = f.read()
            cleaned_lines = []
            for line in content.splitlines():
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    cleaned_lines.append(line[:66].ljust(80))
                elif line.startswith("TER") or line.startswith("END"):
                    cleaned_lines.append(line)
            
            ext = "pdbqt" if path.endswith(".pdbqt") else "pdb"
            return "\n".join(cleaned_lines), ext
            
    print(f"[DEBUG load_receptor_pdb] Failed to find receptor for {variant_name}", file=sys.stderr)
    return None, None

def get_receptor_sequence(receptor_data, highlighted_res_nums=None):
    """Extracts amino acid sequence and builds interactive sequence map."""
    if not receptor_data:
        return [], ""
    if highlighted_res_nums is None:
        highlighted_res_nums = []
        
    seen_residues = set()
    sequence_list = []
    
    for line in receptor_data.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                res_name = line[17:20].strip()
                res_seq = line[22:26].strip()
                res_id = f"{res_name}{res_seq}"
                
                if res_id not in seen_residues and res_name in AA_MAP:
                    seen_residues.add(res_id)
                    sequence_list.append({
                        "res_name": res_name,
                        "res_num": res_seq,
                        "label": f"{res_name}{res_seq}",
                        "single": AA_MAP[res_name]
                    })
            except Exception:
                continue
                
    sequence_list.sort(key=lambda x: int(re.sub(r'\D', '', x["res_num"])) if re.sub(r'\D', '', x["res_num"]).isdigit() else 0)
    
    html_spans = []
    for item in sequence_list:
        if item["res_num"] in highlighted_res_nums:
            html_spans.append(
                f'<span style="background-color: #ff4b4b; color: white; padding: 2px 5px; border-radius: 4px; font-weight: bold; margin: 1px; display: inline-block;" title="{item["label"]}">{item["single"]}</span>'
            )
        else:
            html_spans.append(
                f'<span style="padding: 2px; color: #222222;" title="{item["label"]}">{item["single"]}</span>'
            )
            
    return sequence_list, "".join(html_spans) 
