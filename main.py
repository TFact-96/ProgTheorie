from classes.AminoLattice import AminoLattice
from algorithms.chaingeneration.chaingenerate import generate_chain
from algorithms.chaingeneration.bruteforce import bruteforce_chains
from visualisation.data import get_chain_data, get_plot_data, get_chain_from_file, write_chain_to_csv
from visualisation.plot3D import plot_chain
import os

# Get this shitty optimalization_tries out. It slows down the program by a significant deal.
# Optimalization_tries means the amount of moves a new atom will try until stability change is found.
# It does moves in a random way, and an already used move can be used again, consuming a try for no reason.
# So, even if it uses 20+ random moves, it wouldn't even ensure that it uses the move where the bond is found.
#
# TIP: rewrite the generate_random_valid_node() function in AminoLattice class so that an already used move will not be tried again,
# and if all possible moves are used and no stability change has occured; make it do a random move.
# This would make the "amount of tries" obsolete, and only 6 moves max would be used vs 20+ random moves.
# Which makes the runtime of the greedy-moves algorithm about 5x faster than it now is,
# (but still on average 3x slower than the random-moves, albeit infinitely more efficient)
#
# Spent 8 hours on this shit to no avail, will buy coffee and frikandelbroodje for you if you fix it.
optimalization_tries = 10

def clear_terminal():
    os.system('cls' if os.name == 'nt' else 'clear')

def plotting_and_data_handler(lattice):
    plot = input("Do you want to plot the chain? (y/n): ")
    if (plot == "y"):
        plot_chain(lattice)

    save_data = input("Do you want to save the atom moves of the chain into a .csv file? (y/n): ")
    if (save_data == "y"):
        write_chain_to_csv(get_chain_data(lattice))

def get_user_input_for_generating_chain():
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

    if str_brute == "y":
        iterations = int(input("How many chain generations for brute forcing?: "))
        return iterations, amino, str_optimize, str_brute

    return 0, amino, str_optimize, str_brute

if __name__ == "__main__":
    clear_terminal()

    ################### Loading existing chain data or generating own?
    from_csv = input("Plot an existing amino chain from a .csv file, or generate a new one? (y = load / n = generate): ")

    if (from_csv == "y"):
        file = input("Give the filename from the data folder (without .csv extension): ")
        lattice = get_chain_from_file(file)
        if lattice:
            plot_chain(lattice)

        exit(0)

    ################### Get chain data input for generation
    clear_terminal()
    iterations, amino, str_optimize, str_brute = get_user_input_for_generating_chain()

    ################### Chain generation
    clear_terminal()

    # Random generation
    if str_optimize == "n":
        # bruteforce
        if str_brute == "y":
            lattice = bruteforce_chains(amino, iterations, False, 0)
        # non bruteforce
        else:
            lattice = generate_chain(amino, False, 0)

    # Optimized (greedymoves) generation
    if str_optimize == "y":
        # bruteforce
        if str_brute == "y":
            lattice = bruteforce_chains(amino, iterations, True, optimalization_tries)
        # not bruteforce
        else:
            lattice = generate_chain(amino, True, optimalization_tries)

    # Handling plotting and data of chain
    if lattice:
        plotting_and_data_handler(lattice)

    exit(0)
