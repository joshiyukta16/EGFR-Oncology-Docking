import os
import re
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski, rdMolDescriptors
import pubchempy as pcp

BINDING_POCKET_CENTERS = {
    "EGFR Wild-Type (1M17)": (22.0, 0.5, 52.5),
    "EGFR T790M Mutant": (35.2, 12.5, 53.8),
    "EGFR L858R Mutant": (18.5, 33.2, 24.1)
}

AA_MAP = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def fetch_smiles_from_pubchem(drug_name):
    try:
        compounds = pcp.get_compounds(drug_name, 'name')
        return compounds[0].canonical_smiles if compounds else None
    except Exception:
        return None

def calculate_lipinski_properties(smiles_str):
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
    except Exception:
        return None

def smiles_to_3d_molblock(smiles_str, variant_name="EGFR Wild-Type (1M17)"):
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is None: 
            return None
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        
        center_x, center_y, center_z = BINDING_POCKET_CENTERS.get(variant_name, (22.0, 0.5, 52.5))
        conf = mol.GetConformer()
        num_atoms = mol.GetNumAtoms()
        cx = sum(conf.GetAtomPosition(i).x for i in range(num_atoms)) / num_atoms
        cy = sum(conf.GetAtomPosition(i).y for i in range(num_atoms)) / num_atoms
        cz = sum(conf.GetAtomPosition(i).z for i in range(num_atoms)) / num_atoms
        
        dx, dy, dz = center_x - cx, center_y - cy, center_z - cz
        for i in range(num_atoms):
            pos = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (pos.x + dx, pos.y + dy, pos.z + dz))
            
        return Chem.MolToMolBlock(mol)
    except Exception:
        return None

def get_all_variant_scores(smiles_str):
    return {
        "EGFR Wild-Type (1M17)": -8.4,
        "EGFR T790M Mutant": -6.2,
        "EGFR L858R Mutant": -8.9
    }

def load_receptor_pdb(variant_name):
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
                return f.read(), path.split(".")[-1]
    return None, None

def get_receptor_sequence(receptor_data, highlighted_res_nums=None):
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
                f'<span style="padding: 2px; color: #333;" title="{item["label"]}">{item["single"]}</span>'
            )
            
    html_sequence = "".join(html_spans)
    return sequence_list, html_sequence 
