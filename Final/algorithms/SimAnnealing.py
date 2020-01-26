import copy
import random
import numpy as np
from algorithms.Random import random_chain

def simulated_annealing(
        amino, iterations, start_temperature, use_linear_temp,
        use_exp_temp, linear_temp_coeff, exp_temp_coeff
    ):
    
    # make a random protein chain
    current_state = random_chain(amino)
    
    # current stability
    current_stability = copy.copy(current_state.stability)

    # for statistic plotting
    stability_over_time = []
        
    # set temperature
    temperature = start_temperature
    
    # commence the simulated annealing
    for iteration in range(iterations):
        
        # cmd output
        print(f"Iteration: {iteration}")
        
        # keep track of stability for every iteration
        stability_over_time.append(current_stability)
        
        # keep old state
        old_state = copy.copy(current_state)
        
        # make pullmove on random node
        random_node_index = np.random.randint(1, len(current_state.grid_chain) - 1)
        
        # get node object
        node_coords = current_state.grid_chain[random_node_index][0]
        node = current_state.grid[node_coords].nodes[0]
        
        # perform a pullmove on this node
        current_state.pull_move(node)

        # calculate new stability
        current_state.update_neighbours()
        new_stability = copy.copy(current_state.stability)
        
        # dart shot
        accept_value = 2**((current_stability - new_stability) / temperature)    
        random_shot = random.random()
        
        # undo move if random shot is above accept value
        if random_shot > accept_value:
            current_state = copy.copy(old_state)
            current_stability = copy.copy(old_state.stability)            

        # lower temperature
        if use_linear_temp:
            temperature = start_temperature - (linear_temp_coeff * iteration)

        if use_exp_temp:
            temperature = start_temperature * (exp_temp_coeff**iteration)

        # overflow fix
        if temperature <= 0.001:
            temperature = 0.001
            
    current_state.update_neighbours()
    # put the chain into the best_chain dict for plotting compatibility
    best_chain = {}
    best_chain[current_stability] = current_state
    
    return best_chain, stability_over_time