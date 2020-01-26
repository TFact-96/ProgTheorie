import copy
from classes.Grid import Grid
from algorithms.CalcUpperbound import calc_upperbound
from algorithms.Random import random_chain

def hill_climber(amino, whole_chain_pull_amount, amount_of_reset_checks):
    local_minimum_chains = {}
    stability_over_time = []
    chain_nr = 1
    
    # set minimal stability the chain has to get
    naive_upperbound = calc_upperbound(amino)
    
    # create initial random chain
    temp_chain = random_chain(amino)
    
    # Save as initial local minimum
    best_stability = copy.copy(temp_chain.stability)
    best_current_chain = temp_chain
    
    for iteration in range(amount_of_reset_checks):
        better_stab_found = False
        
        for it in range(whole_chain_pull_amount):
            stability_over_time.append(best_stability)

            for index in range(1, len(temp_chain.grid_chain) - 1):
                temp_chain.pull_move(
                    temp_chain.grid[temp_chain.grid_chain[index][0]].nodes[0]
                )
                temp_chain.update_neighbours()
                
                if (temp_chain.stability < best_stability) and (temp_chain.stability < naive_upperbound):
                    best_current_chain = copy.deepcopy(temp_chain)
                    best_stability = copy.copy(temp_chain.stability)
                    better_stab_found = True

        print(f"Chain {chain_nr}: Pulled {whole_chain_pull_amount} times. Stability: {best_stability}")

        # if it didnt find any upgrades after the amount of pulls of the whole chain, must be local minimum
        if not better_stab_found:
            # terminal info
            print(f"No stability change after {whole_chain_pull_amount} pulls. Saving chain {chain_nr} as local minima...\n")
            
            
            # append to local minima list
            local_minimum_chains[best_current_chain.stability] = best_current_chain
            
            # reset to new chain
            temp_chain = random_chain(amino)
            chain_nr += 1
            
            # set this as best
            best_stability = temp_chain.stability
            best_current_chain = temp_chain
    
    return local_minimum_chains, stability_over_time