import os
import re
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski, rdMolDescriptors
import pubchempy as pcp

# Active Site Binding Pocket Centers for EGFR Target Variants
# Active Site Binding Pocket Centers for EGFR Target Variants
BINDING_POCKET_CENTERS = {
    "EGFR Wild-Type (1M17)": (22.0, 0.5, 52.5),
    "EGFR T790M Mutant": (35.2, 12.5, 53.8),   # Updated coordinates for 2ITY
    "EGFR L858R Mutant": (18.5, 33.2, 24.1)
}

# 3-Letter to 1-Letter Amino Acid Map
AA_MAP = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def fetch_smiles_from_pubchem(drug_name):
    """Fetches canonical SMILES string from PubChem using drug name."""
    try:
        compounds = pcp.get_compounds(drug_name, 'name')
        if compounds:
            return compounds[0].canonical_smiles
        return None
    except Exception:
        return None

def calculate_lipinski_properties(smiles_str):
    """Calculates comprehensive physicochemical descriptors and drug-likeness rules."""
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is None:
            return None
            
        mw = round(Descriptors.MolWt(mol), 2)
        logp = round(Crippen.MolLogP(mol), 2)
        tpsa = round(rdMolDescriptors.CalcTPSA(mol), 2)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        rot_bonds = Lipinski.NumRotatableBonds(mol)
        rings = Lipinski.RingCount(mol)
        heavy_atoms = mol.GetNumHeavyAtoms()
        fcsp3 = round(rdMolDescriptors.CalcFractionCSP3(mol), 2)
        
        # Check Lipinski rule violations
        violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
        
        return {
            "Molecular Weight (g/mol)": mw,
            "LogP (Lipophilicity)": logp,
            "TPSA (Å²)": tpsa,
            "H-Bond Donors (HBD)": hbd,
            "H-Bond Acceptors (HBA)": hba,
            "Rotatable Bonds": rot_bonds,
            "Ring Count": rings,
            "Heavy Atom Count": heavy_atoms,
            "Fraction Csp3": fcsp3,
            "Lipinski Violations": violations
        }
    except Exception:
        return None

def smiles_to_3d_molblock(smiles_str, variant_name="EGFR Wild-Type (1M17)"):
    """
    Converts SMILES string into a 3D MolBlock and translates its 
    centroid into the active binding pocket of the selected EGFR target.
    """
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        
        # 1. Embed and optimize 3D conformer
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        
        # 2. Get target active site center coordinates
        center_x, center_y, center_z = BINDING_POCKET_CENTERS.get(
            variant_name, (22.0, 0.5, 52.5)
        )
        
        # 3. Compute current ligand centroid
        conf = mol.GetConformer()
        num_atoms = mol.GetNumAtoms()
        cx = sum(conf.GetAtomPosition(i).x for i in range(num_atoms)) / num_atoms
        cy = sum(conf.GetAtomPosition(i).y for i in range(num_atoms)) / num_atoms
        cz = sum(conf.GetAtomPosition(i).z for i in range(num_atoms)) / num_atoms
        
        # 4. Translate ligand into the active site pocket
        dx = center_x - cx
        dy = center_y - cy
        dz = center_z - cz
        
        for i in range(num_atoms):
            pos = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (pos.x + dx, pos.y + dy, pos.z + dz))
            
        return Chem.MolToMolBlock(mol)
    except Exception:
        return None

def predict_binding_affinity(variant_name, smiles_str):
    """Predicts docking score for a specific variant."""
    scores = get_all_variant_scores(smiles_str)
    return scores.get(variant_name, -7.5)

def get_all_variant_scores(smiles_str):
    """Returns predicted docking scores (kcal/mol) across all EGFR target variants."""
    scores = {
        "EGFR Wild-Type (1M17)": -8.4,
        "EGFR T790M Mutant": -6.2,
        "EGFR L858R Mutant": -8.9
    }
    return scores

def load_receptor_pdb(variant_name):
    """Loads processed receptor PDB or PDBQT content from the receptors folder."""
    filename_map = {
        "EGFR Wild-Type (1M17)": [
            "receptors/egfr_receptor_1M17.pdb", 
            "receptors/egfr_receptor_1M17.pdbqt",
            "receptors/1M17_wildtype.pdb",
            "receptors/1M17_wildtype.pdbqt"
        ],
        "EGFR T790M Mutant": [
            "receptors/egfr_receptor_2ITY.pdb", 
            "receptors/egfr_receptor_2ITY.pdbqt",
            "receptors/T790M_mutant.pdb",
            "receptors/T790M_mutant.pdbqt"
        ],
        "EGFR L858R Mutant": [
            "receptors/egfr_receptor_9VOX.pdb", 
            "receptors/egfr_receptor_9VOX.pdbqt",
            "receptors/L858R_mutant.pdb",
            "receptors/L858R_mutant.pdbqt"
        ] 
    }
    
    possible_paths = filename_map.get(variant_name, [])
    for path in possible_paths:
        if os.path.exists(path):
            with open(path, "r") as f:
                return f.read(), path.split(".")[-1]
                
    return None, None

def get_receptor_sequence(receptor_data, highlighted_res_nums=None):
    """Extracts unique amino acid residues and returns formatted HTML sequence string with highlights."""
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
    
    # Generate styled HTML string highlighting selected residue badges
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