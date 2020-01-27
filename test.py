import itertools
import numpy as np
import random
import timeit
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import random


a = 0
b = 10
switcher = False

while b > 5:

    if switcher == False:
        a += 1
        index = a
        switcher = True
    elif switcher == True:
        b -= 1
        index = b
        switcher = False

    print(index)

