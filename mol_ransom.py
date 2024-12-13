from rdkit import Chem
from rdkit.Chem import Draw
import warnings
from rdkit.Chem.Draw import MolsToGridImage
import molfiles
import imageio
import math
import tempfile
from rdkit.Chem.Draw import rdMolDraw2D

_letter_to_mol_dict = dict([[l, Chem.MolFromMolBlock(mol_file)]
                            for l, mol_file in molfiles.letter_to_molfile.items()])

def single_word_per_line(string_v):
    words = string_v.split(" ")
    max_len = max(len(w) for w in words)
    to_ret = "".join([e for w in words for e in (w + " " * (max_len - len(w)))])
    return to_ret, max_len

def ascii_to_mols(string):
    _letter_to_mols = _letter_to_mol_dict
    string_upper = [s.upper() for s in string]
    missing_letters = sorted(set(string_upper) - set(_letter_to_mols.keys()) - set(" "))
    if len(missing_letters) > 0:
        warnings.warn(UserWarning(f"Didn't understand the following letters; will appear blank:{missing_letters}"))
    mols = [ _letter_to_mols[s] if s in _letter_to_mols else Chem.MolFromSmiles("")
             for s in string_upper]
    for m in mols:
        rdMolDraw2D.PrepareMolForDrawing(m)
    return mols

def to_image(output_file,mols,same_scale=True,padding=0.00,
             black_and_white=False,rotate_degrees=None,
             letters_per_row=5):
    opt = Draw.MolDrawOptions()
    if black_and_white:
        opt.useBWAtomPalette()
    if rotate_degrees is not None:
        opt.rotate = rotate_degrees
    opt.prepareMolsBeforeDrawing = False
    opt.drawMolsSameScale = same_scale
    opt.padding = padding
    opt.fixedScale = 0.05
    opt.centreMoleculesBeforeDrawing = False
    img = MolsToGridImage(mols=mols, molsPerRow=letters_per_row,
                          drawOptions=opt)
    img.save(output_file)

def acscii_to_mol_image(string,output_file,one_word_per_line=False,**kw):
    string, kw = adjust_kw_as_needed(string, one_word_per_line=one_word_per_line,
                                     kw=kw)
    mols = ascii_to_mols(string)
    to_image(output_file,mols,**kw)

def adjust_kw_as_needed(one_word_per_line,string,kw):
    if one_word_per_line:
        string, length = single_word_per_line(string)
        kw['letters_per_row'] = length
    return string, kw

def ascii_to_mol_movie(string,output_file,total_time_s=5,rotation_degrees=360,
                       s_per_frame=1/5,one_word_per_line=False,loop_forever=True,
                       start_degrees=0,and_reverse=False,**kw):
    string, kw = adjust_kw_as_needed(one_word_per_line,string, kw)
    total_frames = int(math.ceil(total_time_s/s_per_frame))
    degrees_per_frame = rotation_degrees/total_frames
    assert total_frames > 1 , "Must have at least one frame"
    mols = ascii_to_mols(string)
    if and_reverse:
        start_directions = [[0, +1],[rotation_degrees,-1]]
    else:
        start_directions = [[0, +1]]
    with imageio.get_writer(output_file, mode='I',loop=0 if loop_forever else None) as writer:
        with tempfile.NamedTemporaryFile(suffix=".png") as temp:
            for start,directions in start_directions:
                for frame_N in range(0,total_frames,1):
                    rotate_degrees = start_degrees + start + degrees_per_frame * frame_N * directions
                    to_image(temp.name, mols, rotate_degrees=rotate_degrees,**kw)
                    image = imageio.v3.imread(temp.name)
                    writer.append_data(image)
    foo = 1

def run():
    pass


if __name__ == "__main__":
    run()
