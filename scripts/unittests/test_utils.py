import unittest
from src.utils import *
import sys
import os
sys.dont_write_bytecode = True


class TestUtils(unittest.TestCase):

    def test_remove_branch_lengths(self):
        input1 = "((t1:0.12029,t2:0.232309)(t3:23.237323,t4:893.829389));"
        expected_output1 = "((t1,t2)(t3,t4));"
        calculated_output1 = remove_branch_lengths(input1)
        self.assertEqual(expected_output1, calculated_output1)
        input2 = "(t3,t12,(t2,((((t13,(t11,t7)1:0.33894013053783395)1:0.45718108964384196,((t15,t10):0.021230360592972005,#H30:0.04634541896266962::0.3539320314293897):1.0721775005863046)1:0.48235201057096455,((t1,t9)1:0.040234728770749704,((t14,t4)0.985724:0.054416203199337806,(t5,((t6,t8)1:0.0)#H30:0.08039797772894099::0.6460679685706103)0.532078:0.012528005514031309)1:0.3308743841186229)0.997299:1.3089997541445408)1:0.4885065344542067,OUT)1:1.3011012767550771)1:0.8502460402338993)1;"
        expected_output2 = "(t3,t12,(t2,((((t13,(t11,t7)),((t15,t10),#H30)),((t1,t9),((t14,t4),(t5,((t6,t8))#H30)))),OUT)));"
        calculated_output2 = remove_branch_lengths(input2)
        self.assertEqual(expected_output2, calculated_output2)

    def test_get_mpl_score(self):
        loc_path = os.path.dirname(__file__)
        expected_mpl = -487565.7058449446
        calculated_mpl = get_mpl_score(os.path.join(loc_path, "calgtprob_example.txt"))
        self.assertEqual(expected_mpl, calculated_mpl)

    def test_get_fn_fp(self):
        loc_path = os.path.dirname(__file__)
        expected_fn = 0.21052631578947367
        expected_fp = 0.11764705882352941
        calculated_fn, calculated_fp = get_fn_fp(os.path.join(loc_path, "cmpnets_example.txt"))
        self.assertEqual(expected_fn, calculated_fn)
        self.assertEqual(expected_fp, calculated_fp)


