import math
import random
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist
import statistics as stat

from classes import aminolattice, atom

if __name__ == "__main__":

    amino = "CHCHCHCHCHCHCHCCCCHHHHCHCHCCCCCHHHP"

    chain = aminolattice.AminoLattice(amino)
    chain.generate_nodes()
