from classes.AminoLattice import AminoLattice
from algorithms.chaingeneration.chaingenerate import generate_chain
from visualisation.data import get_chain_data

###################################### (bruteforce) ALGORITHM
def bruteforce_chains(amino, iterations, use_optimize_algorithm, optimalization_tries):
    best_stability = 0

    print(f"Brute force generating {iterations} chains:")

    # try n iterations for best stability
    for i in range(iterations):
        lattice = generate_chain(amino, use_optimize_algorithm, optimalization_tries)

        # only count non-stuck chains (aka lattice is not None)
        if lattice:
            stability, moves = get_chain_data(lattice)

            # if this generation is a new record
            if stability <= best_stability:
                best_lattice = lattice
                best_stability = stability
                best_generation = i
                print(f"Generation {i}: Stability {best_stability}.")

    print(f"\nFinished! Best generation was {best_generation} with Stability = {best_stability}")
    return best_lattice
###################################### END (bruteforce) ALGORITHM
