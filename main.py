from classes.AminoLattice import AminoLattice
from algorithms.chaingenerate import generate_chain
from algorithms.bruteforce import bruteforce_chains
from visualisation.data import get_chain_data, get_plot_data
from visualisation.plot3D import plot_chain
import os

optimalization_tries = 10

def generate_one_chain(amino, use_optimize_algorithm, optimalization_tries):
    lattice = AminoLattice(amino)
    lattice = generate_chain(lattice, use_optimize_algorithm, optimalization_tries)

    if lattice.chain_stuck:
        print("Chain got stuck! Try again.")
        return

    return lattice

def clear_terminal():
    os.system('cls' if os.name == 'nt' else 'clear')

if __name__ == "__main__":
    clear_terminal()
    amino = input("Enter desired amino-chain (C's, H's and P's): ")
    str_optimize = input("Random generation or greedy-move algorithm optimization? (y = greedy / n = random): ")
    str_brute = input("Brute force generate chains and find best one? (y/n): ")

    if (str_optimize != "y" and str_optimize != "n") or (str_brute != "y" and str_brute != "n"):
        print("Only answer with y or n please.")
        exit(0)

    if str_optimize == "n":
        if str_brute == "y":
            iterations = int(input("How many chain generations for brute forcing?: "))
            clear_terminal()
            best_chain = bruteforce_chains(amino, iterations, False, 0)

            # stuck chain check (chain is None when stuck)
            if best_chain:
                plot_chain(best_chain)

        else:
            clear_terminal()
            lattice = generate_one_chain(amino, False, 0)

            # stuck chain check (lattice is None when stuck)
            if lattice:
                plot_chain(lattice)


    if str_optimize == "y":
        if str_brute == "y":
            iterations = int(input("How many chain generations for brute forcing?: "))
            clear_terminal()
            best_chain = bruteforce_chains(amino, iterations, True, optimalization_tries)

            if best_chain:
                plot_chain(best_chain)

        else:
            clear_terminal()
            lattice = generate_one_chain(amino, True, optimalization_tries)

            if lattice:
                plot_chain(lattice)
