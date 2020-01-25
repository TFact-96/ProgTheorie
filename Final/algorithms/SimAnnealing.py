import copy
import random
import numpy as np

def simulated_annealing(grid_class, iterations, temperature):
    # make a random protein chain
    grid_chain, grid = grid_class.create_chain()        
    stability = grid_class.update_neighbours(grid, grid_chain)[0]
    
    # set as best chain as initialization
    current_chain = grid_chain
    current_grid = grid
    
    # make stability positive for acceptation chance formula to work
    current_stability = abs(stability)
    
    # for statistic plotting
    stability_over_time = []
    
    # commence the simulated annealing
    for iteration in range(iterations):
        
        # keep track of stability for every iteration
        stability_over_time.append(current_stability)

        # make pullmove on random node
        random_node_index = np.random.randint(len(grid_class.amino))
        temp_grid, temp_chain = grid_class.pull_move(
            grid[grid_chain[random_node_index][0]].nodes[0], grid, current_chain,
        )

        # calculate new stability
        new_stability = abs(grid_class.update_neighbours(temp_grid, temp_chain)[0])
        
        # dart shot
        accept_chance = 2**((current_stability - new_stability) / temperature)
        random_shot = random.random()
        
        # accept this pullmove if random number is below the accept chance 
        if random_shot < accept_chance:
            current_chain = temp_chain
            current_grid = grid
            current_stability = new_stability

        # lower temperature
        temperature = 1000 * (0.977**iteration)
        #temperature = 1000 - (0.5 * iteration)
    
    return grid_class, stability_over_time