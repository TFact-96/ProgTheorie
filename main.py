from classes.AminoLattice import AminoLattice
from algorithms.chaingenerate import generate_chain
from algorithms.bruteforce import bruteforce_chains
from visualisation.data import get_chain_data, get_plot_data, get_chain_from_file, write_chain_to_csv
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

def plotting_and_data_handler(lattice):
    plot = input("Do you want to plot the chain? (y/n): ")
    if (plot == "y"):
        plot_chain(lattice)

    save_data = input("Do you want to save the atom moves of the chain into a .csv file? (y/n): ")
    if (save_data == "y"):
        write_chain_to_csv(get_chain_data(lattice, False))


if __name__ == "__main__":
    clear_terminal()

    ################### Loading existing chain data
    from_csv = input("Plot an existing amino chain from a .csv file, or generate a new one? (y = load / n = generate): ")

    if (from_csv == "y"):
        file = input("Give the filename from the data folder (without .csv extension): ")
        lattice = get_chain_from_file(file)
        if lattice:
            plot_chain(lattice)

        exit(0)

    ################### Generating chain data
    amino = input("Enter desired amino-chain (C's, H's and P's): ")

    for atom in amino:
        if atom != "H" and atom != "C" and atom != "P":
            print("Only C, H and P atoms allowed.")
            exit(0)

    str_optimize = input("Random generation or greedy-move algorithm optimization? (y = greedy / n = random): ")
    str_brute = input("Brute force generate chains and find best one? (y/n): ")

    if (str_optimize != "y" and str_optimize != "n") or (str_brute != "y" and str_brute != "n"):
        print("Only answer with y or n please.")
        exit(0)

    if str_optimize == "n":
        if str_brute == "y":
            iterations = int(input("How many chain generations for brute forcing?: "))
            clear_terminal()
            lattice = bruteforce_chains(amino, iterations, False, 0)

            # stuck chain check (chain is None when stuck)
            if lattice:
                plotting_and_data_handler(lattice)

        else:
            clear_terminal()
            lattice = generate_one_chain(amino, False, 0)

            # stuck chain check (lattice is None when stuck)
            if lattice:
                plotting_and_data_handler(lattice)

    if str_optimize == "y":
        if str_brute == "y":
            iterations = int(input("How many chain generations for brute forcing?: "))
            clear_terminal()
            lattice = bruteforce_chains(amino, iterations, True, optimalization_tries)

            if lattice:
                plotting_and_data_handler(lattice)

        else:
            clear_terminal()
            lattice = generate_one_chain(amino, True, optimalization_tries)

            if lattice:
                plotting_and_data_handler(lattice)
