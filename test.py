import itertools
import numpy as np

checker = {
    "[1, 1]": (True, [1, 1]),
    "[1, -1]": (True, [1, -1]),
    "[-1, 1]": (True, [-1, 1]),
    "[-1, -1]": (True, [-1, -1]),
}


checker["[1, 1]"][0] = False
for d_check, d_bool in checker.items():
    print(d_bool)
