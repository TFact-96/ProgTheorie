from classes.ChainLattice import ChainLattice
from algorithms.chaingeneration.chaingenerate import generate_chain
from visualisation.data import get_chain_data

###################################### Generates multiple random/greedy chains and puts their stability vs chain nr in a list
def multiple_chains(protein, iterations, use_greedy, greedy_tries):
    chain_stability = []
    chain_nr = []

    print(f"\n Generating {iterations} chains:")

    # try n iterations for best stability
    for i in range(iterations):
        Chain = generate_chain(protein, use_greedy, greedy_tries)

        # calculate stability and set bonds in bond list
        Chain.set_stability_and_bonds()
        stability, moves = get_chain_data(Chain)
        chain_nr.append(i)
        chain_stability.append(stability)


    return chain_nr, chain_stability
######################################
