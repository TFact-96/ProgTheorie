###################################### Greedy ALGORITHM (The more C's and H's, the better it works)
def generate_atom_greedy_move(lattice, optimalization_tries):
    stability = lattice.get_stability_and_bonds(True)
    new_stability = stability
    i = 0
    index = len(lattice.chain)

    # if a move makes stability change; keep that move
    while i < optimalization_tries and stability == new_stability:
        # get a random new atom
        new_atom = lattice.generate_atom_random_move()

        if new_atom:
            # append to chain for testing
            lattice.chain[index] = new_atom

            # calculate new stability with this move
            new_stability = lattice.get_stability_and_bonds(True)

            # pop atom (otherwise it keeps adding)
            lattice.chain.pop(index)

            # break if stability changed
            if new_stability < stability:
                break

        i += 1

    # return the new atom
    return new_atom
###################################### END OPTIMALIZATION ALGORITHM
