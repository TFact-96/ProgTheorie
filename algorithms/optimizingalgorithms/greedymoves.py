###################################### Greedy ALGORITHM (The more C's and H's, the better it works)
def generate_amino_greedy_move(Chain, greedy_tries):
    stability = Chain.get_stability_and_bonds(True)
    new_stability = stability
    i = 0
    index = len(Chain.state)

    # if a move makes stability change; keep that move
    while i < greedy_tries and stability == new_stability:
        # get a random new amino
        new_amino = Chain.generate_amino_random_move()

        if new_amino:
            # append to chain for testing
            Chain.state[index] = new_amino

            # calculate new stability with this move
            new_stability = Chain.get_stability_and_bonds(True)

            # pop amino (otherwise it keeps adding)
            Chain.state.pop(index)

            # break if stability changed
            if new_stability < stability:
                break

        i += 1

    # return the new amino
    return new_amino
###################################### END OPTIMALIZATION ALGORITHM
