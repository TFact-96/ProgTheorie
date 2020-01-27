import copy
from algorithms.CalcUpperbound import calc_upperbound
from algorithms.RandomChain import random_chain
from algorithms.PullMove import pull_move

def hill_climber(protein, whole_chain_pull_amount, amount_of_reset_checks):
    local_minimum_chains = {}
    stability_over_time = []
    chain_nr = 1

    # set minimal stability the chain has to get
    naive_upperbound = calc_upperbound(protein)

    # create initial grid object with random chain in it
    temp_chain = random_chain(protein)

    # Save as initial local minimum
    best_stability = copy.copy(temp_chain.stability)

    # how many checks if the chain is a local minimum
    for iteration in range(amount_of_reset_checks):
        better_stab_found = False

        # amount of pulling the whole chain
        for it in range(whole_chain_pull_amount):
            stability_over_time.append(best_stability)

            # pulling the whole chain
            for index in range(1, len(temp_chain.grid_chain) - 1):
                # get node object
                node_coords = temp_chain.grid_chain[index][0]
                node = temp_chain.grid[node_coords].nodes[0]

                # perform a pullmove on this node and update stability and bonds
                pull_move(temp_chain, node)

            # update new bonds and stability of this chain
            temp_chain.update_all_bonds()

            # update best current chain if the stability is better
            if (temp_chain.stability < best_stability):

                # only copy the filled gridpoints (the chain) for lower computing time
                temp_chain.set_filled_gridpoints_from_grid()
                best_current_grid = copy.deepcopy(temp_chain.filled_gridpoints)
                best_current_chain = copy.deepcopy(temp_chain.grid_chain)
                best_stability = copy.copy(temp_chain.stability)
                better_stab_found = True

        print(f"Random chain {chain_nr}: Pulled {whole_chain_pull_amount} times. Stability: {best_stability}")

        # if it didnt find any upgrades after the amount of pulls of the whole chain, must be local minimum
        if not better_stab_found:
            # terminal info
            print(f"No stability change after {whole_chain_pull_amount} pulls. Saving chain {chain_nr} as local minima...\n")

            # set the best filled gridpoints configuration into the chain
            temp_chain.filled_gridpoints = best_current_grid
            temp_chain.grid_chain = best_current_chain

            # merge into its grid
            temp_chain.set_grid_from_filled_gridpoints()

            # set bonds and stability
            temp_chain.update_all_bonds()

            # append to local minima list
            local_minimum_chains[temp_chain.stability] = temp_chain

            # reset to new chain
            temp_chain = random_chain(protein)
            chain_nr += 1

            # set this as best
            best_stability = temp_chain.stability

    return local_minimum_chains, stability_over_time
