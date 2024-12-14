from rdkit import Chem
from rdkit.Chem import Draw
import warnings, string
from rdkit.Chem.Draw import MolsToGridImage
import molfiles
import imageio, click
import math
import tempfile
from rdkit.Chem.Draw import rdMolDraw2D

def process_true_false(value):
    """

    :param value: value
    :return: true if the value is like TRUE (case insensitive"
    """
    return str(value).upper() == "TRUE"

def process_arguments(kw_args):
    """

    :param kw_args: arguments for cli
    :return: processed arguments, where boolean-like arguments (e.g., True/False) are transformed
    """
    set_bool = ['one_word_per_line',
                'and_reverse',
                'comic_mode',
                'black_and_white',
                'loop_forever']
    for k,v in kw_args.items():
        if k in set_bool:
            kw_args[k] = process_true_false(v)
    return kw_args


kw_true_false = dict(type=click.Choice([True,False,"FALSE","False","TRUE","True"]),
                     required=False)


_letter_to_mol_dict = dict([[l, Chem.MolFromMolBlock(mol_file)]
                            for l, mol_file in molfiles.letter_to_molfile.items()])

def _valid_characters():
    """

    :return: list of valid characters which have molefile representation
    Note this is equivalent to ASCII code 32 to 127, but not including anything
    lowercase

    See https://www.ascii-code.com/
    """
    set_ok_to_miss = {'\t', '\n', '\r', '\x0b', '\x0c'}
    return [s for s in string.printable if s not in set_ok_to_miss]

def single_word_per_line(string_v):
    """

    :param string_v: string
    :return: new string, where each word (space-delimited) is right-padded
    such that all words are on their own line
    """
    words = string_v.split(" ")
    max_len = max(len(w) for w in words)
    to_ret = "".join([e for w in words for e in (w + " " * (max_len - len(w)))])
    return to_ret, max_len

def ascii_to_mols(string):
    """

    :param string: string, length N
    :return: length N list of mols, one per character in <string>.
     If string contains unsupported characters
    (i.e., characters not in <_valid_characters>), those will be blank mols
    and a warning raied
    """
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

def mols_to_image(output_file,mols,same_scale=True,padding=0.00,
                  black_and_white=False,rotate_degrees=None,
                  letters_per_row=5,comic_mode=False):
    """

    :param output_file: where to output image
    :param mols: list of mols, see ascii_to_mols
    :param same_scale:  see rdkit.Draw.MolDrawOptions https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions
    :param padding:  see rdkit.Draw.MolDrawOptions https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions
    :param black_and_white: if true, render in black and white
    :param rotate_degrees: how much to rotae molecules (0 is 12 o clock, -90 is 9 o clock, etc)
    :param letters_per_row: how many letters per row
    :param comic_mode:  see rdkit.Draw.MolDrawOptions https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions
    :return: nothing
    """
    opt = Draw.MolDrawOptions()
    if black_and_white:
        opt.useBWAtomPalette()
    if rotate_degrees is not None:
        opt.rotate = rotate_degrees
    opt.prepareMolsBeforeDrawing = False
    opt.drawMolsSameScale = same_scale
    if comic_mode:
        opt.comicMode = True
    opt.padding = padding
    opt.fixedScale = 0.05
    opt.centreMoleculesBeforeDrawing = False
    img = MolsToGridImage(mols=mols, molsPerRow=letters_per_row,
                          drawOptions=opt)
    img.save(output_file)

def _image(string,output_file,one_word_per_line=False,**kw):
    """

    :param string:  se mols_to_image
    :param output_file:  see mols_to_image
    :param one_word_per_line:  if true, output image such that there is one word per line
    :param kw:  see mols_to_image
    :return:  Nothing
    """
    string, kw = adjust_kw_as_needed(string=string,
                                     one_word_per_line=one_word_per_line,
                                     kw=kw)
    mols = ascii_to_mols(string)
    mols_to_image(output_file,mols,**kw)

def adjust_kw_as_needed(string,one_word_per_line,kw):
    """

    :param string: input string
    :param one_word_per_line: see _image
    :param kw: additional keywords
    :return:  tuple of (updated string,updated keywords) necessary to have
    one word per line
    """
    if one_word_per_line:
        string, length = single_word_per_line(string)
        kw['letters_per_row'] = length
    return string, kw

def _animate(string,output_file,total_time_s=5,rotation_degrees=90,
             frames_per_second=10,one_word_per_line=False,loop_forever=True,
             start_degrees=-45,and_reverse=False,**kw):
    """

    :param string:  see <_image>
    :param output_file:  see <_image>, but should be gif
    :param total_time_s:  how long the .gif should last
    :param rotation_degrees: how many degrees to use
    :param frames_per_second:  how many frames per second
    :param one_word_per_line:  see _image
    :param loop_forever: if true, loop animation forever
    :param start_degrees:  offset rotation at start of animatnino
    :param and_reverse:  if true, animation advances forward then reverse
    :param kw: passed to _image
    :return:  nothing
    """
    string, kw = adjust_kw_as_needed(string=string,
                                     one_word_per_line=one_word_per_line,
                                     kw=kw)
    # if going in reverse with N frames total, must have 1/2 the frames in
    # forward direction and 1/2 the frames in the reverse direction
    factor = 0.5 if and_reverse else 1
    total_frames = int(math.ceil(total_time_s*frames_per_second) * factor)
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
                    mols_to_image(temp.name, mols, rotate_degrees=rotate_degrees,**kw)
                    image = imageio.v3.imread(temp.name)
                    writer.append_data(image)


@click.group()
def cli():
    pass

@cli.command()
@click.option('--output_file', required=True,type=click.Path(),
              help="Name of output file (only png supported)")
@click.option('--string', required=True,type=str,
              help="What string to convert to molecules")
@click.option('--black_and_white', default=False,**kw_true_false,
              help="If true, use black and white instead of colors")
@click.option('--comic_mode',  default=False,
              help="If true, use comic-drawing mode (xkcd style)",**kw_true_false)
@click.option('--letters_per_row', required=False,default=8,type=int,
              help="How many letters per row (i.e., line character limit)")
@click.option('--one_word_per_line',  default=False,**kw_true_false,
              help="If true, have one word per row/line (equivalent to setting <letters_per_row> to the longest word length)")
def image(**kw):
    _image(**process_arguments(kw))

@cli.command()
@click.option('--output_file', required=True,type=click.Path(),
              help="Name of output file (only gif supported)")
@click.option('--string', required=True,type=str,
              help="What string to convert to molecules")
@click.option('--black_and_white', default=False,**kw_true_false,
              help="If true, use black and white instead of colors")
@click.option('--comic_mode',  default=False,
              help="If true, use comic-drawing mode (xkcd style)",**kw_true_false)
@click.option('--letters_per_row', required=False,default=8,type=int,
              help="How many letters per row (i.e., line character limit)")
@click.option('--one_word_per_line',  default=False,**kw_true_false,
              help="If true, have one word per row/line (equivalent to setting <letters_per_row> to the longest word length)")
@click.option('--loop_forever', default=True,**kw_true_false,
              help="If true, loop for forever, otherwise just will run once")
@click.option('--and_reverse', default=True,**kw_true_false,
              help="If true, animation will run in forward direction then in the reverse direction")
@click.option('--total_time_s', required=False,
              default=5,type=float,help="How long for total animation")
@click.option('--rotation_degrees', required=False,
              default=90,type=float,help="Total rotation in degrees")
@click.option('--frames_per_second', required=False,
              default=10,type=float,help="How many frames per second")
@click.option('--start_degrees', required=False,
              default=-45,type=float,help="Where to start animation (0 = 12 o clock, -90 would be 9 o clock, +90 would be 3 o clock, etc)")
def animate(**kw):
    _animate(**process_arguments(kw))

if __name__ == '__main__':
    cli()