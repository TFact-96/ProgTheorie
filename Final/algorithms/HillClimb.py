import copy

def hill_climber(whole_chain_pull_amount, amount_of_reset_checks, grid_class):
    
    # Make a protein chain
    current_hilltop, grid = grid_class.create_chain()        
    current_stability = grid_class.update_neighbours(grid, current_hilltop)[0]
    
    # Save as best chain (initialize)
    best_c = current_hilltop
    best_stab_c = current_stability
    best_grid = grid
    grid_class.best_chain[best_stab_c] = [best_c, best_grid]
    
    # for statistic plotting
    stability_over_time = []
    
    # try iteration amount of reset checks
    for iteration in range(amount_of_reset_checks):
        best_c_found = False
        print(f"try: {iteration}")

        # amount the >whole< chain should be pulled
        for it in range(whole_chain_pull_amount):
            # keep track of all stability improvements for plotting statistics
            stability_over_time.append(best_stab_c)
                        
            # try to pull each node
            for index in range(1, len(current_hilltop) - 1):

                temp_grid, temp_chain = grid_class.pull_move(
                    grid[current_hilltop[index][0]].nodes[0], grid, current_hilltop,
                )

                temp_stability = grid_class.update_neighbours(temp_grid, temp_chain)[0]
                
                # if stability is better after pull, save this as best chain
                if temp_stability < best_stab_c:
                    best_c = copy.deepcopy(temp_chain)
                    best_grid = copy.deepcopy(temp_grid)
                    best_stab_c = temp_stability
                    best_c_found = True
                    
        # if better chain is found with these pulls, go on with this one
        if best_c_found:
            current_hilltop = copy.deepcopy(best_c)
            current_stability = best_stab_c
            grid = copy.deepcopy(best_grid)
            best_c_found = False

        # if no stability change after this amount of chainlength pulls, must be local maximum. Save in best_chain list
        # and reset the chain.
        else:
            grid_class.best_chain[best_stab_c] = [best_c, best_grid]

            current_hilltop, grid = grid_class.create_chain()
            current_stability = grid_class.update_neighbours(grid, current_hilltop)[0]

            best_c = current_hilltop
            best_stab_c = current_stability
            best_grid = grid
    
    return grid_class, stability_over_time
