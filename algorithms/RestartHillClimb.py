import copy
from algorithms.CalcUpperbound import calc_upperbound
from algorithms.RandomChain import random_chain
from algorithms.PullMove import pull_move
from classes.Grid import Grid

def hill_climber(protein, whole_chain_pull_amount, amount_of_reset_checks):
    local_minimum_chains = {}
    stability_over_time = []
    chain_nr = 1

    # set minimal stability the chain has to get
    naive_upperbound = calc_upperbound(protein)

    # create initial grid object with random chain in it
    grid_object = Grid(protein)
    grid_object = random_chain(grid_object)

    # Save as initial local minimum
    best_stability = copy.copy(grid_object.stability)
    
    # only get filled gridpoints and deepcopy that as best chain
    grid_object.set_filled_gridpoints_from_grid()
    best_current_grid = copy.deepcopy(grid_object.filled_gridpoints)
    best_current_chain = copy.deepcopy(grid_object.grid_chain)

    # how many checks if the chain is a local minimum
    for iteration in range(amount_of_reset_checks):
        better_stab_found = False

        # amount of pulling the whole chain
        for it in range(whole_chain_pull_amount):
            stability_over_time.append(best_stability)

            # pulling the whole chain
            for index in range(1, len(grid_object.grid_chain) - 1):
                # get node object
                node_coords = grid_object.grid_chain[index][0]
                node = grid_object.grid[node_coords].nodes[0]

                # perform a pullmove on this node and update stability and bonds
                pull_move(grid_object, node)

            # update new bonds and stability of this chain
            grid_object.update_all_bonds()

            # update best current chain if the stability is better
            if (grid_object.stability < best_stability):

                # only copy the filled gridpoints (the chain) for lower computing time
                grid_object.set_filled_gridpoints_from_grid()
                best_current_grid = copy.deepcopy(grid_object.filled_gridpoints)
                best_current_chain = copy.deepcopy(grid_object.grid_chain)
                best_stability = copy.copy(grid_object.stability)
                better_stab_found = True

        print(f"Random chain {chain_nr}: Pulled {whole_chain_pull_amount} times. Stability: {best_stability}")

        # if it didnt find any upgrades after the amount of pulls of the whole chain, could be local minimum
        if not better_stab_found:
            # terminal info
            print(f"No stability change after {whole_chain_pull_amount} pulls. Saving chain {chain_nr} as local minima...\n")

            # append to local minima list
            local_minimum_chains[copy.copy(grid_object.stability)] = [copy.deepcopy(best_current_grid), copy.deepcopy(best_current_chain)]

            # reset to new chain
            grid_object = random_chain(grid_object)
            chain_nr += 1

            # set this as best
            best_stability = copy.deepcopy(grid_object.stability)

    return local_minimum_chains, stability_over_time
