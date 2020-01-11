from classes.AminoLattice import AminoLattice
from algorithms.chaingenerate import generate_chain
from visualisation.data import get_chain_data
import os

# clearing the terminal
def clear_terminal():
    os.system('cls' if os.name == 'nt' else 'clear')

###################################### (bruteforce) ALGORITHM
def bruteforce_chains(amino, iterations, use_optimize_algorithm, optimalization_tries):
    best_stability = 0

    print(f"Brute force generating {iterations} chains:")

    # try n iterations for best stability
    for i in range(iterations):
        lattice = AminoLattice(amino)
        lattice = generate_chain(lattice, use_optimize_algorithm, optimalization_tries)

        # only count non-stuck chains
        if not lattice.chain_stuck:
            stability, moves = get_chain_data(lattice, False)

            # if this generation is a new record
            if stability <= best_stability:
                best_lattice = lattice
                best_stability = stability
                best_generation = i
                print(f"Generation {i}: Stability {best_stability}.")

    # print this chain
    get_chain_data(best_lattice, True)

    return best_lattice
###################################### END (bruteforce) ALGORITHM
