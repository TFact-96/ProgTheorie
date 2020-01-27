import copy
import random
from algorithms.RandomChain import random_chain

def simulated_annealing(
        protein, iterations, start_temperature, use_linear_temp,
        use_exp_temp, linear_temp_coeff, exp_temp_coeff
    ):

    # make a grid object with a random chain configuration with this protein
    current_state = random_chain(protein)

    # for statistic plotting
    stability_over_time = []

    # set temperature
    temperature = start_temperature

    # commence the simulated annealing
    for iteration in range(iterations):

        # keep track of stability for every iteration
        stability_over_time.append(current_state.stability)

        # remember old state (only filled gridpoints and the chain list)
        current_state.set_filled_gridpoints_from_grid()
        old_filled_gridpoints = copy.deepcopy(current_state.filled_gridpoints)
        old_chain = copy.deepcopy(current_state.grid_chain)
        old_stability = copy.deepcopy(current_state.stability)

        # choose a random node index
        random_node_index = np.random.randint(1, len(current_state.grid_chain) - 1)

        # get node object
        node_coords = current_state.grid_chain[random_node_index][0]
        node = current_state.grid[node_coords].nodes[0]

        # perform a pullmove on this node and update stability and bonds
        pull_move(current_state, node)

        # calculate new stability
        current_state.update_all_bonds()
        new_stability = copy.deepcopy(current_state.stability)

        # dart shot
        accept_value = 2**((old_stability - new_stability) / temperature)
        random_shot = random.random()

        # undo move if random shot is above accept value
        if random_shot > accept_value:
            # get the old filled gridpoints and put them back into the grid
            current_state.filled_gridpoints = copy.deepcopy(old_filled_gridpoints)
            current_state.grid_chain = copy.deepcopy(old_chain)
            current_state.set_grid_from_filled_gridpoints()
            current_state.update_all_bonds()

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
