from algorithms.optimizingalgorithms.greedymoves import generate_amino_greedy_move
from classes.ChainLattice import ChainLattice

###################################### Generating all aminos one by one onto the lattice
def generate_chain(protein, use_greedy, greedy_tries):
    Chain = ChainLattice(protein)

    # while amount of aminos already generated is smaller than the whole protein
    while len(Chain.state) < len(Chain.protein):

        if use_greedy:
            # generate a optimalized move (new amino always makes bond-making moves if possible)
            new_amino = generate_amino_greedy_move(Chain, greedy_tries)
        else:
            # get a random move per new amino added to chain
            new_amino = Chain.generate_amino_random_move()
        
        # reset chain if chain stuck with this new move
        if Chain.chain_stuck(new_amino.x, new_amino.y, new_amino.z):
            Chain.state = {0: Chain.first_amino}
            
        # append the new amino to the existing chain
        index = len(Chain.state)
        Chain.state[index] = new_amino
    
    # calculate and set the protein stability, and put the bonds in their respective coord list
    Chain.set_stability_and_bonds()

    return Chain
