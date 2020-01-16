from classes.AminoLattice import AminoLattice
from algorithms.chaingeneration.chaingenerate import generate_chain
from visualisation.data import get_chain_data

###################################### (bruteforce) ALGORITHM
def multiple_chains(amino, iterations, use_optimize_algorithm, optimalization_tries, ThreeD):
    chain_stability = []
    chain_nr = []
    
    print(f"\nBrute force generating {iterations} chains:")

    # try n iterations for best stability
    for i in range(iterations):
        lattice = generate_chain(amino, use_optimize_algorithm, optimalization_tries, ThreeD)

        # only count non-stuck chains (aka lattice is not None)
        if lattice:
            stability, moves = get_chain_data(lattice)
            chain_nr.append(i)
            chain_stability.append(stability)
            

    return chain_nr, chain_stability
###################################### END (bruteforce) ALGORITHM
