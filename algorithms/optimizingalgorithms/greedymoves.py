###################################### Greedy ALGORITHM (The more C's and H's, the better it works)
def generate_greedy_move(lattice, optimalization_tries):
    stability = lattice.get_stability_and_bonds(True)
    new_stability = stability
    i = 0

    # if a move makes stability change; keep that move
    while i < optimalization_tries and stability == new_stability:
        # get a random new atom
        new_atom = lattice.generate_random_valid_node()

        if new_atom:
            # append to chain for testing
            lattice.chain.append(new_atom)

            # calculate new stability with this move
            new_stability = lattice.get_stability_and_bonds(True)

            # pop atom (otherwise it keeps adding)
            lattice.chain.pop(-1)

            # break if stability changed
            if new_stability < stability:
                break

        i += 1

    # return the new atom
    return new_atom
###################################### END OPTIMALIZATION ALGORITHM
