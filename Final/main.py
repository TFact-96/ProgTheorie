from classes.Grid import Grid
from algorithms.HillClimb import hill_climber
from algorithms.SimAnnealing import simulated_annealing
from algorithms.Random import random_chain
from visualisation.DataPlots import data_plot_hillclimb, data_plot_annealing
from visualisation.PlotBestChain import plot_best_chain

def get_initial_user_input():
    amino = input("Enter desired protein chain: ")
    make_random_chain = input("Do you want to just make a random chain? (y/n): ")
    restart_hill_climb = input("Do you want to use the Restart Hill Climbing algorithm? (y/n): ")
    sim_annealing = input("Do you want to use the Simulated Annealing Hill Climb algorithm? (y/n): ")

    return amino, make_random_chain, restart_hill_climb, sim_annealing

def plot_chain_request(chains):
    plot_request = input("Do you want to plot the chain? (y/n): ")
    if plot_request == "y":
        plot_best_chain(chains)

# Simulated Annealing
def annealing_flow(amino):
    # best start temp around 2 and decrease rate for exponential 0.995
    # for linear: decrease rate start_temp / iterations
    repeat_amount = int(input("Sim Annealing: Amount of annealings: "))
    iteration_amount = int(input("Sim Annealing: How many random pullmove iterations per annealing run?: "))
    start_temp = float(input("Sim Annealing: Enter the start temperature: "))
    exponential = input("Sim Annealing: Linear or exponential temperature decrease over iterations? (y = exponential / n = linear): ")
    coeff = float(input("Sim Annealing: Enter desired temperature decrease coefficient: "))

    # run the annealing
    best_chain, stability_over_time = annealing_repeat(amino, repeat_amount, iteration_amount, start_temp, coeff, exponential)

    # plot stability over time for hillclimb statistics
    data_plot_request = input("Sim Annealing: Do you want to plot the stability over time? (y/n): ")

    if data_plot_request == "y":
        data_plot_annealing(stability_over_time, amino)

    plot_chain_request(best_chain)

# simple bruteforce repeating for the best simulated annealing run
def annealing_repeat(amino, repeat_amount, iteration_amount, start_temp, coeff, exponential):
    best_stability = 0

    for iteration in range(repeat_amount):
        if exponential == "n":
            new_chain, stability_over_time = simulated_annealing(amino, iteration_amount, start_temp, True,
                    False, coeff, 0)
        else:
            new_chain, stability_over_time = simulated_annealing(amino, iteration_amount, start_temp, False,
                    True, 0, coeff)

        print(f"Iteration {iteration}: Stability = {new_chain.stability}")

        # save new best run if found
        if new_chain.stability < best_stability:
            best_chain, best_stability_over_time = new_chain, stability_over_time
            best_stability = new_chain.stability

    print(f"Best chain: {best_chain.stability}")

    # for plotting compatibility
    best_chains = {}
    best_chains[best_chain.stability] = best_chain

    return best_chains, best_stability_over_time


def restart_hill_climb_flow(amino):
    reset_checks = int(input("Restart Hillclimb: Enter the amount of chain restart checks: "))
    chain_pull_amt = int(input("Restart Hillclimb: Enter the amount of times the whole chain should be pulled per reset check: "))

    local_minima_chains, stability_over_time = hill_climber(amino, chain_pull_amt, reset_checks)

    # plot stability over time for hillclimb statistics
    data_plot_request = input("Restart Hillclimb: Do you want to plot the stability over time? (y/n): ")

    if data_plot_request == "y":
        data_plot_hillclimb(stability_over_time, amino)

    plot_chain_request(local_minima_chains)


def main():
    amino, make_random_chain, restart_hill_climb, sim_annealing = get_initial_user_input()

    if make_random_chain == "y":
        # for plotting compatibility
        chains = {}

        # make random chain
        grid = random_chain(amino)

        # for plotting compatibility
        chain[grid.stability] = grid

        plot_chain_request(chain)

    if sim_annealing == "y":
        annealing_flow(amino)

    if restart_hill_climb == "y":
        restart_hill_climb_flow(amino)

    return

if __name__ == "__main__":
    main()
