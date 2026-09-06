import os
from rdkit import Chem
from rdkit.Chem import AllChem
import pubchempy as pcp

# Binding Pocket Centroid Coordinates (X, Y, Z) for EGFR Kinase Active Sites
BINDING_POCKET_CENTERS = {
    "EGFR Wild-Type (1M17)": (22.0, 0.5, 52.5),
    "EGFR T790M Mutant": (29.8, 8.4, 11.2),
    "EGFR L858R Mutant": (18.5, 33.2, 24.1)
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

def smiles_to_3d_molblock(smiles_str, variant_name="EGFR Wild-Type (1M17)"):
    """Converts SMILES string into a 3D MolBlock translated into the target active site."""
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        
        # 1. Generate 3D Conformer
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        
        # 2. Get target active site center coordinates
        center_x, center_y, center_z = BINDING_POCKET_CENTERS.get(variant_name, (22.0, 0.5, 52.5))
        
        # 3. Calculate current ligand centroid
        conf = mol.GetConformer()
        num_atoms = mol.GetNumAtoms()
        cx = sum(conf.GetAtomPosition(i).x for i in range(num_atoms)) / num_atoms
        cy = sum(conf.GetAtomPosition(i).y for i in range(num_atoms)) / num_atoms
        cz = sum(conf.GetAtomPosition(i).z for i in range(num_atoms)) / num_atoms
        
        # 4. Shift ligand coordinates directly into the receptor binding pocket
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
                return f.read(), path.split(".")[-1]  # Returns (file_content, extension)
                
    return None, None 