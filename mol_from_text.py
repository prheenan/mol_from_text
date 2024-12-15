import warnings
import math
import tempfile
import cairosvg
import imageio
import click
from rdkit.Chem.Draw import MolsToGridImage, MolDrawOptions
from rdkit.Chem.Draw.rdMolDraw2D import PrepareMolForDrawing
from rdkit.Chem import SDWriter, MolFromMolBlock, MolFromSmiles
from click import ParamType
import molfiles


class BoolType(ParamType):
    def __init__(self):
        pass

    def get_metavar(self, param):
        return 'Choice([TRUE/True/FALSE/False])'


    def convert(self, value, _, __):
        upper = str(value).upper()
        if upper == "TRUE":
            return True
        elif upper == "FALSE":
            return False
        else:
            self.fail(f"Invalid value: {value}. Expected TRUE/FALSE")
            return False



kw_true_false = dict(type=BoolType(),required=False)

_letter_to_mol_dict = dict([[l, MolFromMolBlock(mol_file)]
                            for l, mol_file in molfiles.letter_to_molfile.items()])



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
        msg = f"Don't support these characters; will be blank:{missing_letters}"
        warnings.warn(UserWarning(msg))
    mols = [ _letter_to_mols[s] if s in _letter_to_mols else MolFromSmiles("")
             for s in string_upper]
    for m in mols:
        PrepareMolForDrawing(m)
    return mols

def mols_to_array(mols,same_scale=True,padding=0.00,
                  black_and_white=False,rotate_degrees=None,
                  letters_per_row=5,comic_mode=False,use_svg=False,
                  dots_per_angstrom=300,font_size=12,bond_width=1.25,
                  scale=0.055,letter_size_pixels=200):
    """

    :param mols: list of mols, see ascii_to_mols
    :param same_scale:  see rdkit.Draw.MolDrawOptions https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions
    :param padding:  see rdkit.Draw.MolDrawOptions https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions
    :param black_and_white: if true, render in black and white
    :param rotate_degrees: how much to rotae molecules (0 is 12 o clock, -90 is 9 o clock, etc)
    :param letters_per_row: how many letters per row
    :param comic_mode:  see rdkit.Draw.MolDrawOptions https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions
    :param use_svg: if true, use svg instead
    :param dots_per_angstrom:  see rdkit.Draw.MolDrawOptions https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions
    :param font_size:  see rdkit.Draw.MolDrawOptions https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions
    :param bond_width: see rdkit.Draw.MolDrawOptions https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions
    :param scale:  see rdkit.Draw.MolDrawOptions https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions
    :param letter_size_pixels: size of each pixel
    :return: nothing
    """
    opt = MolDrawOptions()
    if black_and_white:
        opt.useBWAtomPalette()
    if rotate_degrees is not None:
        opt.rotate = rotate_degrees
    opt.prepareMolsBeforeDrawing = False
    opt.drawMolsSameScale = same_scale
    opt.dotsPerAngstrom = dots_per_angstrom
    opt.minFontSize = font_size
    opt.bondLineWidth = bond_width
    opt.comicMode = comic_mode
    opt.padding = padding
    opt.fixedScale = scale
    opt.centreMoleculesBeforeDrawing = False
    img = MolsToGridImage(mols=mols, molsPerRow=letters_per_row,
                          drawOptions=opt,useSVG=use_svg,
                          subImgSize=(letter_size_pixels,letter_size_pixels))
    return img

def _write_svg(output_file,img):
    """

    :param output_file: where to output
    :param img: svg image as string
    :return:  n/a
    """
    with open(output_file, "w",encoding="utf8") as f:
        f.write(img)

def mols_to_image(output_file,**kw):
    """

    :param output_file: where to output image
    :param kw: see mols_to_array
    :return: nothing, saves the file
    """
    img = mols_to_array(use_svg=True,**kw)
    if output_file.endswith(".png"):
        # convert from svg to png first, this results in better images, see
        # https://iwatobipen.wordpress.com/2017/11/03/draw-high-quality-molecular-image-in-rdkit-rdkit/
        with tempfile.NamedTemporaryFile(suffix=".svg") as f:
            _write_svg(f.name, img)
            # scale = 4 seems to give decent results
            cairosvg.svg2png(url=f.name, write_to=output_file, scale=4.0)
    elif output_file.endswith(".svg"):
        _write_svg(output_file, img)
    else:
        raise ValueError(f"Didn't understand file extension for {output_file}")


def _image(string,output_file,one_word_per_line=False,**kw):
    """

    :param string:  se mols_to_image
    :param output_file:  see mols_to_image
    :param one_word_per_line:  if true, output image such that there is one word per line
    :param kw:  see mols_to_image
    :return:  Nothing
    """
    if output_file is None:
        output_file = f"{string}.svg"
    if not (output_file.endswith(".png") or output_file.endswith(".svg")):
        raise TypeError(f"Only .png and .svg supported, but given output file:{output_file}")
    string, kw = adjust_kw_as_needed(string=string,
                                     one_word_per_line=one_word_per_line,
                                     kw=kw)
    mols = ascii_to_mols(string)
    mols_to_image(output_file=output_file,mols=mols,**kw)

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

def _animate(string,output_file,total_time_s=5.,rotation_degrees=90.,
             frames_per_second=10.,one_word_per_line=False,loop_forever=True,
             start_degrees=-45.,and_reverse=False,**kw):
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
    if output_file is None:
        output_file = f"{string}.gif"
    if not output_file.endswith(".gif"):
        raise TypeError(f"Only .gif supported, but given output file:{output_file}")
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
    all_images = []
    with imageio.get_writer(output_file, mode='I',loop=0 if loop_forever else None) as writer:
        for start,directions in start_directions:
            for frame_n in range(0,total_frames,1):
                rotate_degrees = start_degrees + start + degrees_per_frame * frame_n * directions
                png = mols_to_array(mols, rotate_degrees=rotate_degrees,
                                    use_svg=False,**kw)
                writer.append_data(png)
                all_images.append(png)
    return all_images

def _export_structures(output_file):
    with SDWriter(output_file) as writer:
        for id_v, mol in _letter_to_mol_dict.items():
            mol.SetProp(key='ascii_character', val=id_v)
            writer.write(mol)

@click.group()
def cli():
    pass

@cli.command()
@click.option('--string', required=True,type=str,
              help="What string to convert to molecules")
@click.option('--output_file', required=False,type=click.Path(),
              help="Name of output file (only png and svg supported)")
@click.option('--black_and_white', default=False,**kw_true_false,
              help="If true, use black and white instead of colors")
@click.option('--comic_mode',  default=False,
              help="If true, use comic-drawing mode (xkcd style)",**kw_true_false)
@click.option('--letters_per_row', required=False,default=5,type=int,
              help="How many letters per row (i.e., line character limit)")
@click.option('--one_word_per_line',  default=False,**kw_true_false,
              help="If true, have one word per row/line (equivalent to setting <letters_per_row> to the longest word length)")
@click.option('--dots_per_angstrom', required=False,
              default=300.,type=float,help="See rdkit.Draw.MolDrawOptions")
@click.option('--font_size', required=False,
              default=12,type=int,help="See rdkit.Draw.MolDrawOptions")
@click.option('--bond_width', required=False,
              default=1.25,type=float,help="See rdkit.Draw.MolDrawOptions")
@click.option('--scale', required=False,
              default=0.055,type=float,help="See rdkit.Draw.MolDrawOptions")
@click.option('--letter_size_pixels', required=False,
              default=200,type=int,help="Individual size of letters")
def image(**kw):
    _image(**kw)

@cli.command()
@click.option('--output_file', required=False,type=click.Path(),
              help="Name of output file (only gif supported)")
@click.option('--string', required=True,type=str,
              help="What string to convert to molecules")
@click.option('--black_and_white', default=False,**kw_true_false,
              help="If true, use black and white instead of colors")
@click.option('--comic_mode',  default=False,
              help="If true, use comic-drawing mode (xkcd style)",**kw_true_false)
@click.option('--letters_per_row', required=False,default=5,type=int,
              help="How many letters per row (i.e., line character limit)")
@click.option('--one_word_per_line',  default=False,**kw_true_false,
              help="If true, have one word per row/line (equivalent to setting <letters_per_row> to the longest word length)")
@click.option('--loop_forever', default=True,**kw_true_false,
              help="If true, loop for forever, otherwise just will run once")
@click.option('--and_reverse', default=True,**kw_true_false,
              help="If true, animation will run in forward direction then in the reverse direction")
@click.option('--total_time_s', required=False,
              default=5.,type=float,help="How long for total animation")
@click.option('--rotation_degrees', required=False,
              default=90.,type=float,help="Total rotation in degrees")
@click.option('--frames_per_second', required=False,
              default=10.,type=float,help="How many frames per second")
@click.option('--start_degrees', required=False,
              default=-45.,type=float,help="Where to start animation (0 = 12 o clock, -90 would be 9 o clock, +90 would be 3 o clock, etc)")
@click.option('--dots_per_angstrom', required=False,
              default=300.,type=float,help="See rdkit.Draw.MolDrawOptions")
@click.option('--font_size', required=False,
              default=12,type=int,help="See rdkit.Draw.MolDrawOptions")
@click.option('--bond_width', required=False,
              default=1.25,type=float,help="See rdkit.Draw.MolDrawOptions")
@click.option('--scale', required=False,
              default=0.055,type=float,help="See rdkit.Draw.MolDrawOptions")
@click.option('--letter_size_pixels', required=False,
              default=200,type=int,help="Individual size of letters")
def animate(**kw):
    _animate(**kw)

@cli.command()
@click.option('--output_file', required=True,type=click.Path(),
              help="Name of output file (only .sdf supported)")
def export_structures(**kw):
    _export_structures(**kw)

if __name__ == '__main__':
    cli()
