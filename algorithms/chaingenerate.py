from algorithms.optimize import generate_optimal_move

###################################### Generating all atoms one by one onto the lattice
def generate_chain(lattice, use_optimize_algorithm, optimalization_tries):
    # while amount of nodes already generated is smaller than the whole amino string length
    while len(lattice.chain) < len(lattice.amino) and not lattice.chain_stuck:

        if use_optimize_algorithm:
            # generate a optimalized move (new node always makes bond-making moves if possible)
            new_node = generate_optimal_move(lattice, optimalization_tries)
        else:
            # get a random move per new atom added to chain
            new_node = lattice.generate_random_valid_node()

        # append the new node to the existing chain
        lattice.chain.append(new_node)

    if lattice.chain_stuck:
        print("Chain got stuck! Try again.")

    return lattice
