import copy
from classes.Grid import Grid
    
def hill_climber(amino, whole_chain_pull_amount, amount_of_reset_checks):
    best_chains = {}
    stability_over_time = []
    
    grid_class = Grid(amino)
    # Current protein chain
    grid_class.create_chain()
    grid_class.update_neighbours()
    current_stability = grid_class.stability

    # Save as best chain
    best_stab_c = current_stability
    best_chains[best_stab_c] = grid_class

    for iteration in range(amount_of_reset_checks):
        best_c_found = False
        print(f"try: {iteration}")
        
        for it in range(whole_chain_pull_amount):
            stability_over_time.append(current_stability)

            for index in range(1, len(grid_class.grid_chain) - 1):
                grid_class.pull_move(
                    grid_class.grid[grid_class.grid_chain[index][0]].nodes[0]
                )
                grid_class.update_neighbours()
                temp_stability = grid_class.stability
                
                if temp_stability < best_stab_c:
                    best_chain_class = copy.deepcopy(grid_class)
                    best_stab_c = best_chain_class.stability
                    best_c_found = True

        if best_c_found:
            current_hilltop = copy.deepcopy(best_chain_class)
            current_stability = best_stab_c
            best_c_found = False

        else:
            best_chains[best_stab_c] = best_chain_class
            
            current_hilltop.create_chain()
            current_hilltop.update_neighbours()
            current_stability = current_hilltop.stability
            best_stab_c = current_stability

    return best_chains, stability_over_time