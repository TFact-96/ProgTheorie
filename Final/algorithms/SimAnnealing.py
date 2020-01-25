import copy
import random
import numpy as np
from classes.Grid import Grid

def simulated_annealing(
        amino, iterations, start_temperature, use_linear_temp,
        use_exp_temp, linear_temp_coeff, exp_temp_coeff
    ):
    
    # make a random protein chain
    grid_class = Grid(amino)
    grid_class.create_chain()        
    grid_class.update_neighbours()
    
    # current stability
    current_stability = grid_class.stability

    # for statistic plotting
    stability_over_time = []
        
    # set temperature
    temperature = start_temperature
    
    # commence the simulated annealing
    for iteration in range(iterations):
        # keep track of stability for every iteration
        stability_over_time.append(current_stability)
        
        # keep old state
        old_state = copy.deepcopy(grid_class)
        
        # make pullmove on random node
        random_node_index = np.random.randint(1, len(grid_class.grid_chain) - 1)
        grid_class.pull_move(
            grid_class.grid[grid_class.grid_chain[random_node_index][0]].nodes[0],
        )

        # calculate new stability
        grid_class.update_neighbours()
        new_stability = grid_class.stability
        
        # dart shot
        accept_line = 2**((current_stability - new_stability) / temperature)    
        random_shot = random.random()
        
        # undo pullmove
        if random_shot > accept_line:
            grid_class = copy.deepcopy(old_state)
        else:
            current_stability = copy.deepcopy(new_stability)

        # lower temperature
        if use_linear_temp:
            temperature = start_temperature - (linear_temp_coeff * iteration)
            
            if temperature <= 0.001:
                temperature = 0.001
                
            print(temperature)
        if use_exp_temp:
            temperature = start_temperature * (exp_temp_coeff**iteration)
            print(temperature)
    
    print(current_stability)
    # put the chain into the best_chain dict
    best_chains = {}
    best_chains[best_stability] = grid_class
    
    return best_chains, stability_over_time