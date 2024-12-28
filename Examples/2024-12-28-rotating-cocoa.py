"""
Simple file to make rotating theobromine
"""
import sys
from rdkit.Chem import MolFromSmiles
from rdkit.Chem.Draw.rdMolDraw2D import PrepareMolForDrawing
sys.path.append("../")
# pylint: disable=wrong-import-position
import mol_from_text

def run():
    """
    Run the code to generate the rotating cocoa logo
    """
    smiles_theobromine = "CN1C=NC2=C1C(=O)NC(=O)N2C"
    mol_theobromine = MolFromSmiles(smiles_theobromine)
    PrepareMolForDrawing(mol_theobromine)
    mol_from_text.animate_mols(output_file="logo.gif",mols=[mol_theobromine],
                               total_time_s=5., rotation_degrees=360,
                               frames_per_second=30., loop_forever=True,
                               letter_size_pixels=800,start_degrees=0,
                               and_reverse=False,letters_per_row=1,
                               scale=0.09)



if __name__ == "__main__":
    run()
