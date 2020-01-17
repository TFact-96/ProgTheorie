from classes.ChainLattice import ChainLattice
from algorithms.chaingeneration.chaingenerate import generate_chain
from visualisation.data import get_chain_data

###################################### (bruteforce) ALGORITHM
def multiple_chains(protein, iterations, use_greedy, greedy_tries, ThreeD):
    chain_stability = []
    chain_nr = []

    print(f"\n Generating {iterations} chains:")

    # try n iterations for best stability
    for i in range(iterations):
        Chain = generate_chain(protein, use_greedy, greedy_tries, ThreeD)

        # only count non-stuck chains (aka lattice is not None)
        if lattice:
            stability, moves = get_chain_data(Chain)
            chain_nr.append(i)
            chain_stability.append(stability)


    return chain_nr, chain_stability
###################################### END (bruteforce) ALGORITHM
