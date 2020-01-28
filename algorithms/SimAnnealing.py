import copy
import random
import numpy as np
from algorithms.RandomChain import random_chain
from algorithms.PullMove import pull_move
from classes.Grid import Grid


def initialize(protein, start_temperature):
    """
    Create Grid object of protein and initialize temperature.
    :param protein: protein string
    :param start_temperature: chosen start temperature
    :return: current_state, stability_over_time, temperature
    """
    # make a grid object with a random chain configuration with this protein
    current_state = Grid(protein)
    current_state = random_chain(current_state)

    # for statistic plotting
    stability_over_time = []

    # set initial temperature
    temperature = start_temperature

    return current_state, stability_over_time, temperature
    
def save_current_chain(current_state):
    """
    Deepcopies the current chain.
    :param current_state: the current state
    :return: old_grid, old_chain, old_stability
    """
    # remember old state (only grid and the chain list)
    old_grid = copy.deepcopy(current_state.grid)
    old_chain = copy.deepcopy(current_state.grid_chain)
    old_stability = copy.deepcopy(current_state.stability)

    return old_grid, old_chain, old_stability

def perform_pullmoves(current_state, amount_of_pulls_per_iteration):
    """
    Performs n random pullmoves on the chain and updates all bonds.
    :param current_state: Current state of chain
    :param amount_of_pulls_per_iteration: Amount of pullmoves
    :return: current_state
    """
    # perform pullmoves
    for _ in range(amount_of_pulls_per_iteration):
        # choose a random node index
        random_node_index = np.random.randint(1, len(current_state.grid_chain) - 1)

        # get node object
        node_coords = current_state.grid_chain[random_node_index][0]
        node = current_state.grid[node_coords].nodes[0]

        # perform a pullmove on this node and update stability and bonds
        pull_move(current_state, node)

    # calculate new stability
    current_state.update_all_bonds()
    
    return current_state

def dart_shot(
    current_state, 
    old_stability, 
    new_stability, 
    temperature, 
    old_grid, 
    old_chain, 
    iteration, 
    amount_of_pulls_per_iteration
):
    """
    Accepts or rejects moves based on new score.
    :param current_state: Current state of chain after pullmoves
    :param amount_of_pulls_per_iteration: Amount of pullmoves
    :param old_grid, old_chain: The old state before pullmoves
    :param temperature: chosen temperature
    :param iteration: N'th iteration we're on.
    :return: Current state after reverted back or accepted moves
    """
    accept_value = 2 ** ((old_stability - new_stability) / temperature)
    random_shot = random.random()

    # undo move if random shot is above accept value
    if random_shot > accept_value:
        # get the old grid and put them back into this grid
        current_state.grid = copy.deepcopy(old_grid)
        current_state.grid_chain = copy.deepcopy(old_chain)
        current_state.update_all_bonds()

        print(
            f"Iteration {iteration}: Undo {amount_of_pulls_per_iteration} pulls: Stability = {current_state.stability}"
        )
    else:
        print(
            f"Iteration {iteration}: Accepted {amount_of_pulls_per_iteration} pulls: Stability = {current_state.stability}"
        )
        
    return current_state
    
def lower_temperature(iteration, use_linear_temp, use_exp_temp, linear_temp_coeff, exp_temp_coeff, start_temperature):
    """
    Lowers temperature each iteration either linearily or exponentially.
    :return: new temperature
    """
    # lower temperature either linearily or exponentially
    if use_linear_temp:
        temperature = start_temperature - (linear_temp_coeff * iteration)

    if use_exp_temp:
        temperature = start_temperature * (exp_temp_coeff ** iteration)

    # overflow fix, 0.01 is good enough for good moves to accept almost 100% of the time.
    if temperature <= 0.01:
        temperature = 0.01
        
    return temperature

def simulated_annealing(
    protein,
    iterations,
    amount_of_pulls_per_iteration,
    start_temperature,
    use_linear_temp,
    use_exp_temp,
    linear_temp_coeff,
    exp_temp_coeff,
):
    """
    Simulated annealing.
    :param protein: protein chain
    :param iterations: amount of iterations
    :param amount_of_pulls_per_iteration: amount of pulls per iteration
    :param start_temperature: the start temperature
    :param use_linear_temp: use the linear temperature
    :param use_exp_temp: use exponential temperature
    :param linear_temp_coeff: linear temperature coefficient
    :param exp_temp_coeff: exponential temperature coefficient
    :return: stability, grid, gridchain
    """
    
    # initialize with random chain
    current_state, stability_over_time, temperature = initialize(protein, start_temperature)
    
    # commence the simulated annealing
    for iteration in range(iterations):

        # keep track of stability for every iteration
        stability_over_time.append(current_state.stability)

        # remember state before pulling
        old_grid, old_chain, old_stability = save_current_chain(current_state)
        
        # Perform given amount of pullmoves on the chain
        current_chain = perform_pullmoves(current_state, amount_of_pulls_per_iteration)
        
        # get new stability
        new_stability = copy.deepcopy(current_state.stability)

        # Accepts or rejects the moves based on new score and temperature.
        current_state = dart_shot(
                            current_state, old_stability, 
                            new_stability, temperature, old_grid, 
                            old_chain, iteration, 
                            amount_of_pulls_per_iteration
                        )

        # lower temperature for next iteration
        temperature = lower_temperature(
                            iteration, linear_temp_coeff, 
                            exp_temp_coeff, start_temperature,
                            use_linear_temp, use_exp_temp
                        )

    # return the last iteration
    return (
        current_state.stability,
        current_state.grid,
        current_state.grid_chain,
        stability_over_time,
    )


def annealing_bruteforce(
    protein,
    repeat_amount,
    iteration_amount,
    amount_of_pulls_per_iteration,
    start_temp,
    coeff,
    exponential,
):
    """
    Annealing bruteforce, a simple bruteforce repeating for the best simulated annealing run. 
    :param protein: protein chain
    :param repeat_amount: repeat amount
    :param iteration_amount: amount of iterations
    :param amount_of_pulls_per_iteration: amount of pulls per iteration
    :param start_temp: starting temperature
    :param coeff: given coefficient
    :param exponential: given exponential
    :return: best chain and stability
    """
    best_stability = 0

    # amount of annealing runs
    for iteration in range(repeat_amount):

        # exponential temp decrease or linear temp decrease, then run the sim annealing.
        if exponential == "n":
            stability, new_grid, new_chain, stability_over_time = simulated_annealing(
                protein,
                iteration_amount,
                amount_of_pulls_per_iteration,
                start_temp,
                True,
                False,
                coeff,
                0,
            )
        else:
            stability, new_grid, new_chain, stability_over_time = simulated_annealing(
                protein,
                iteration_amount,
                amount_of_pulls_per_iteration,
                start_temp,
                False,
                True,
                0,
                coeff,
            )

        print(f"Iteration {iteration}: Stability = {stability}")

        # save new best run if found
        if stability < best_stability:
            best_grid, best_chain, best_stability_over_time = (
                new_grid,
                new_chain,
                stability_over_time,
            )
            best_stability = stability

    print(f"Best chain: {best_stability}")

    # for plotting compatibility
    best_chains = {}
    best_chains[best_stability] = [best_grid, best_chain]

    return best_chains, best_stability_over_time
