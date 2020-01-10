import math
import random
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from classes import aminolattice


# Just get one calculation of the chain
def calculate_one_chain(random_move, optimal_move):
    chain = aminolattice.AminoLattice(amino, random_move, optimal_move)
    chain.generate_chain()

    if chain.chain_stuck:
        print("Chain got stuck! Try again.")
    else:
        chain.plot_chain()


# (bruteforce) ALGORITHM
def iterate_generations_of_chains(random_move, optimal_move):
    print("\nBrute force generating chains:")
    print(f"Generating {iterations} chains...")
    best_stability = 0
    best_chain = aminolattice.AminoLattice(amino, random_move, optimal_move)

    # try n iterations for best stability
    for i in range(iterations):
        chain = aminolattice.AminoLattice(amino, random_move, optimal_move)
        chain.generate_chain()

        # only count non-stuck chains
        if not chain.chain_stuck:
            stability, moves = chain.get_chain_data(False)

            # if this generation is a new record
            if stability <= best_stability:
                best_chain = chain
                best_stability = stability
                print(f"Generation {i}: Stability {best_stability}.")

    # print this chain
    best_chain.get_chain_data(True)

    # plot this chain
    best_chain.plot_chain()


# ----------------------------------------- Main structure
if __name__ == "__main__":
    amino = input("\nEnter desired amino-chain (C's, H's and P's): ")
    str = input("Random generation or optimal generation? (y = optimal / n = random): ")
    str_brute = input("Brute force generation to find optimal amino fold? (y/n): ")

    if (str != "y" and str != "n") or (str_brute != "y" and str_brute != "n"):
        print("Only answer with y or n please.")
        exit(0)

    if str == "n":
        if str_brute == "y":
            iterations = int(input("How many chain generations for brute forcing?: "))
            iterate_generations_of_chains(True, False)
        else:
            calculate_one_chain(True, False)

    if str == "y":
        if str_brute == "y":
            iterations = int(input("How many chain generations for brute forcing?: "))
            optimalization_tries = int(
                input("How many times should a node try for an optimal move?: ")
            )
            iterate_generations_of_chains(False, True)
        else:
            optimalization_tries = int(
                input("How many times should a node try for an optimal move?: ")
            )
            calculate_one_chain(False, True)
# ----------------------------------------- End main structure
