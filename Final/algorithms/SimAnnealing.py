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

    # for statistic plotting
    stability_over_time = []

    # set temperature
    temperature = start_temperature

    # commence the simulated annealing
    for iteration in range(iterations):

        # keep track of stability for every iteration
        stability_over_time.append(current_state.stability)

        # remember old state (only filled gridpoints and the chain list)
        current_state.make_filled_gridpoints()
        old_filled_gridpoints = copy.deepcopy(current_state.filled_gridpoints)
        old_chain = copy.deepcopy(current_state.grid_chain)
        old_stability = copy.deepcopy(current_state.stability)

        # perform pullmoves on the whole chain
        random_node_index = np.random.randint(1, len(current_state.grid_chain) - 1)

        # get node object
        node_coords = current_state.grid_chain[random_node_index][0]
        node = current_state.grid[node_coords].nodes[0]

        # perform a pullmove on this node and update stability and bonds
        current_state.pull_move(node)

        # calculate new stability
        current_state.update_neighbours()
        new_stability = copy.deepcopy(current_state.stability)

        # dart shot
        accept_value = 2**((old_stability - new_stability) / temperature)
        random_shot = random.random()

        # undo move if random shot is above accept value
        if random_shot > accept_value:
            # get the old filled gridpoints and put them back into the grid
            current_state.filled_gridpoints = copy.deepcopy(old_filled_gridpoints)
            current_state.grid_chain = copy.deepcopy(old_chain)
            current_state.merge_filled_gridpoints_back()
            current_state.update_neighbours()

        # lower temperature
        if use_linear_temp:
            temperature = start_temperature - (linear_temp_coeff * iteration)

        if use_exp_temp:
            temperature = start_temperature * (exp_temp_coeff**iteration)

        # overflow fix
        if temperature <= 0.001:
            temperature = 0.001

    # return the last iteration
    return current_state, stability_over_time
