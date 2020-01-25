import copy
from classes.Grid import Grid

def hill_climber(amino, whole_chain_pull_amount, amount_of_reset_checks):
    best_chains = {}
    
    # Make a protein chain
    grid_class = Grid(amino)
    grid_class.create_chain()        
    grid_class.update_neighbours()
    best_stability = grid_class.stability
    
    # Save as best chain (initialize)
    best_chains[best_stability] = grid_class
    
    # for statistic plotting
    stability_over_time = []
    
    # try iteration amount of reset checks
    for iteration in range(amount_of_reset_checks):
        better_stability_found = False
        print(f"try: {iteration}")

        # amount the >whole< chain should be pulled
        for it in range(whole_chain_pull_amount):
            # keep track of all stability improvements for plotting statistics
            stability_over_time.append(best_stability)
                        
            # try to pull each node
            for index in range(1, len(grid_class.grid_chain) - 1):

                grid_class.pull_move(
                    grid_class.grid[grid_class.grid_chain[index][0]].nodes[0]
                )

                grid_class.update_neighbours()
                temp_stability = grid_class.stability

                # if stability is better after pull, save this as best chain
                if temp_stability < best_stability:
                    best_stability = temp_stability
                    better_stability_found = True
                    
        # if better chain is found with these pulls, go on with this one
        if better_stability_found:
            better_stability_found = False

        # if no stability change after this amount of chainlength pulls, must be local maximum. Save in best_chain list
        # and reset the chain.
        else:
            best_chains[best_stability] = grid_class
            
            # reset chain and set best stability
            grid_class = Grid(amino)
            grid_class.create_chain()
            grid_class.update_neighbours()
            best_stability = grid_class.stability
            
    return best_chains, stability_over_time
