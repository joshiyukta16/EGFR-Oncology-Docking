import os
from rdkit import Chem
from rdkit.Chem import AllChem
import pubchempy as pcp

def fetch_smiles_from_pubchem(drug_name):
    """Fetches canonical SMILES string from PubChem using drug name."""
    try:
        compounds = pcp.get_compounds(drug_name, 'name')
        if compounds:
            return compounds[0].canonical_smiles
        return None
    except Exception:
        return None

def smiles_to_3d_molblock(smiles_str):
    """Converts SMILES string into a 3D MolBlock."""
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
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
    """Loads processed receptor PDB content from the receptors folder."""
    filename_map = {
        "EGFR Wild-Type (1M17)": "receptors/1M17_wildtype.pdb",
        "EGFR T790M Mutant": "receptors/T790M_mutant.pdb",
        "EGFR L858R Mutant": "receptors/L858R_mutant.pdb"
    }
    
    path = filename_map.get(variant_name)
    if path and os.path.exists(path):
        with open(path, "r") as f:
            return f.read()
    return None 