import unittest
import mol_ransom, string, molfiles
from molfiles import letter_to_molfile


class MyTestCase(unittest.TestCase):
    def test_something(self):
        """
        mol_ransom.acscii_to_mol_image(string=[s.lower() for s in string.ascii_letters],
                                       output_file="tmp.png")
        mol_ransom.acscii_to_mol_image(string="hellofelix",rotate_degrees=30,
                                       output_file="tmp.png")
        """
        self.i_subtest = 0
        for s in mol_ransom._valid_characters():
            with self.subTest(i=self.i_subtest,msg=s):
                assert (s in molfiles.letter_to_molfile or s.upper() in letter_to_molfile), f"Couldn't find character '{s}'"
            self.i_subtest += 1
        with self.subTest(self.i_subtest):
            mol_ransom.single_word_per_line(string_v="why hello there") == ("why  hellothere",5)
        self.i_subtest += 1
        with self.subTest(self.i_subtest):
            mol_ransom.single_word_per_line(string_v="oops a short line") == ("oops a    shortline ", 5)
        self.i_subtest += 1
        mol_ransom._acscii_to_mol_image(string="".join(mol_ransom._valid_characters()),
                                        letters_per_row=10,output_file="test.png")
        mol_ransom._ascii_to_mol_movie(string="i love my mimi and noni",
                                       one_word_per_line=True,and_reverse=True,
                                       rotation_degrees=90,
                                       start_degrees=-45,comic_mode=True,
                                       output_file="test.gif")



if __name__ == '__main__':
    unittest.main()
