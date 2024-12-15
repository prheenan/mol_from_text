import unittest
import shlex
import tempfile
import imageio
from click.testing import CliRunner
import numpy as np
from rdkit.Chem import MolToSmiles,SDMolSupplier
import mol_from_text
import molfiles


def assert_at_most_percent_elements_different(pct,one,two):
    """

    :param pct: maximum percent which can be different
    :param one:  array 1
    :param two:  array 2
    :return: true if at most pct elements are different
    """
    assert one.shape == two.shape
    n_different = np.sum(one != two)
    pct_different = 100 * n_different / one.size
    assert pct_different < pct

class MyTestCase(unittest.TestCase):
    """
    For pre-commit, see

    https://pre-commit.com/

    For pylint, see

    https://pylint.pycqa.org/en/stable/user_guide/installation/pre-commit-integration.html

    """
    def __init__(self,*args,**kw):
        super().__init__(*args,**kw)
        self.i_subtest = 0

    def _read_output(self,f):
        return imageio.v3.imread(f)

    def _function_and_cli(self,function,function_cli,suffix=".png",**kw):
        args = [e for k, v in kw.items() for e in [f"--{k}", v]]
        # use the normal function
        with tempfile.NamedTemporaryFile(suffix=suffix) as temp_func:
            with self.subTest(self.i_subtest):
                function(output_file=temp_func.name, **kw)
                data_function = self._read_output(temp_func.name)
            self.i_subtest += 1
        # use the cli
        with tempfile.NamedTemporaryFile(suffix=suffix) as temp_cli:
            runner = CliRunner()
            args += ["--output_file", temp_cli.name]
            args_str = [str(s) for s in args]
            logger_msg = f"Running the following command:\n\npython {function_cli} {' '.join(args_str)}"
            print(logger_msg)
            with self.subTest(self.i_subtest):
                result = runner.invoke(function_cli, args)
            self.i_subtest += 1
            with self.subTest(self.i_subtest):
                # exit code should be zero
                if result.exit_code != 0:
                    raise ValueError(result.stderr)
            self.i_subtest += 1
            data_cli = self._read_output(temp_cli.name)
        # Due to image quantization (which isn't detemrinistic)
        # the images won't be exactly the same, but they should be close
        # see here https://stackoverflow.com/questions/68366920/set-imageio-compression-level-in-python
        with self.subTest(self.i_subtest):
            assert_at_most_percent_elements_different(10,data_cli,data_function)
        self.i_subtest += 1


    def test_characters(self):
        self.i_subtest = 0
        for s in molfiles._valid_characters():
            with self.subTest(i=self.i_subtest,msg=s):
                assert (s in molfiles.letter_to_molfile or s.upper() in molfiles.letter_to_molfile), f"Couldn't find character '{s}'"
            self.i_subtest += 1

    def test_lineify(self):
        self.i_subtest = 0
        with self.subTest(self.i_subtest):
            assert mol_from_text.\
                       single_word_per_line(string_v="why hello there") == ("why  hellothere",5)
        self.i_subtest += 1
        with self.subTest(self.i_subtest):
            assert mol_from_text.\
                       single_word_per_line(string_v="oops a short line") == ("oops a    shortline ", 5)
        self.i_subtest += 1

    def test_image(self):
        self.i_subtest = 0
        self._function_and_cli(function=mol_from_text._image,
                               function_cli=mol_from_text.image,
                               letters_per_row=6,suffix=".png",
                               string="hello world")
        self._function_and_cli(function=mol_from_text._image,
                               function_cli=mol_from_text.image,
                               letters_per_row=10,suffix=".png",
                               string=shlex.quote("".join(molfiles.letter_to_molfile.keys())),)

    def test_animate(self):
        self.i_subtest = 0
        self._function_and_cli(function=mol_from_text._animate,
                               function_cli=mol_from_text.animate,and_reverse=True,
                               black_and_white=False,
                               rotation_degrees=90.,start_degrees=-45.,
                               comic_mode=True,suffix=".gif",total_time_s=1.,
                               letters_per_row=2,string="hi")
        self._function_and_cli(function=mol_from_text._animate,
                               function_cli=mol_from_text.animate,
                               one_word_per_line=True,and_reverse=True,
                               rotation_degrees=90.,start_degrees=-45.,
                               comic_mode=True,suffix=".gif",
                               string="i love my mimi and noni")

    def test_export(self):
        self.i_subtest = 0
        with tempfile.NamedTemporaryFile(suffix=".sdf") as f:
            mol_from_text._export_structures(output_file=f.name)
            suppl = list(SDMolSupplier(f.name))
            for mol in suppl:
                ascii_character = mol.GetProp("ascii_character")
                with self.subTest(self.i_subtest):
                    assert MolToSmiles(mol) == MolToSmiles(mol_from_text._letter_to_mol_dict[ascii_character])
                self.i_subtest += 1



if __name__ == '__main__':
    unittest.main()
