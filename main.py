from classes.AminoLattice import AminoLattice
from algorithms.chaingenerate import generate_chain
from algorithms.bruteforce import bruteforce_chains
from visualisation.data import get_chain_data, get_plot_data
from visualisation.plot3D import plot_chain

def generate_one_chain(amino, use_optimize_algorithm, optimalization_tries):
    lattice = AminoLattice(amino)
    generated_chain = generate_chain(lattice, use_optimize_algorithm, optimalization_tries)
    return generated_chain

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
            best_chain = bruteforce_chains(amino, iterations, False, 0)
            # plot this chain
            plot_chain(best_chain)
        else:
            lattice = generate_one_chain(amino, False, 0)
            plot_chain(lattice)


    if str == "y":
        if str_brute == "y":
            iterations = int(input("How many chain generations for brute forcing?: "))
            optimalization_tries = int(input("How many times should a node try for an optimal move?: "))
            best_chain = bruteforce_chains(amino, iterations, True, optimalization_tries)
            # plot this chain
            plot_chain(best_chain)

        else:
            optimalization_tries = int(input("How many times should a node try for an optimal move?: "))
            lattice = generate_one_chain(amino, True, optimalization_tries)
            plot_chain(lattice)
