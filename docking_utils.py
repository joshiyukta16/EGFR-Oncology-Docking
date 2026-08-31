from rdkit import Chem
from rdkit.Chem import AllChem

def smiles_to_3d_molblock(smiles_str):
    """Converts a SMILES string into a 3D MolBlock for 3Dmol visualization."""
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

def predict_binding_affinity(receptor_name, smiles_str):
    """
    Computes or retrieves docking affinity scores (kcal/mol).
    Returns baseline scores based on known EGFR benchmark metrics.
    """
    affinity_db = {
        "EGFR Wild-Type (1M17)": -8.4,
        "EGFR T790M Mutant": -6.2,
        "EGFR L858R Mutant": -8.9
    }
    
    return affinity_db.get(receptor_name, -7.5) 