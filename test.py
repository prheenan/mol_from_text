import unittest
import mol_ransom, string




class MyTestCase(unittest.TestCase):
    def test_something(self):
        """
        mol_ransom.acscii_to_mol_image(string=[s.lower() for s in string.ascii_letters],
                                       output_file="tmp.png")
        mol_ransom.acscii_to_mol_image(string="hellofelix",rotate_degrees=30,
                                       output_file="tmp.png")
        """
        self.i_subtest = 0
        with self.subTest(self.i_subtest):
            mol_ransom.single_word_per_line(string_v="why hello there") == ("why  hellothere",5)
        self.i_subtest += 1
        with self.subTest(self.i_subtest):
            mol_ransom.single_word_per_line(string_v="oops a short line") == ("oops a    shortline ", 5)
        self.i_subtest += 1
        mol_ransom.ascii_to_mol_movie(string="snd boobs or wil blowout yu bn wrnd",
                                      one_word_per_line=True,and_reverse=True,
                                      rotation_degrees=90,
                                      start_degrees=-45,
                                      output_file="tmp.gif")
        foo = 1



if __name__ == '__main__':
    unittest.main()
