from classes.Grid import Grid
from algorithms.RestartHillClimb import hill_climber
from algorithms.SimAnnealing import annealing_bruteforce
from algorithms.RandomChain import random_chain
from algorithms.SelfLearning import Agents
from algorithms.SelfLearning import Environment
from visualisation.DataPlots import data_plot_hillclimb, data_plot_annealing
from visualisation.PlotBestChain import plot_best_chain
from visualisation.SelfLearningShow import showProtein

import sys


# Request user at the end if they want to 3D plot the best chain
def plot_chain_request(chains):
    plot_request = input("Do you want to plot the chain? (y/n): ")
    if plot_request == "y":
        plot_best_chain(chains)


# Simulated Annealing Algorithm
def annealing_flow(protein):
    # best start temp around 2 and decrease rate for exponential 0.995
    # for linear: decrease rate start_temp / iterations
    repeat_amount = int(input("Sim Annealing: Amount of bruteforce annealing runs: "))
    iteration_amount = int(input("Sim Annealing: How many iterations per annealing run?: "))
    amount_of_pulls_per_iteration = int(input("Sim Annealing: How many random pullmoves per iteration?: "))
    start_temp = float(input("Sim Annealing: Enter the start temperature: "))
    exponential = input(
        "Sim Annealing: Linear or exponential temperature decrease over iterations? (y = exponential / n = linear): ")
    coeff = float(input("Sim Annealing: Enter desired temperature decrease coefficient: "))

    # run the annealing
    best_chain, stability_over_time = annealing_bruteforce(protein, repeat_amount, iteration_amount,
                                                           amount_of_pulls_per_iteration, start_temp, coeff,
                                                           exponential)

    # plot stability over time for hillclimb statistics
    data_plot_request = input("Sim Annealing: Do you want to plot the stability over time? (y/n): ")

    if data_plot_request == "y":
        data_plot_annealing(stability_over_time, protein)

    plot_chain_request(best_chain)


# Restart Hill Climb Algorithm
def restart_hill_climb_flow(protein):
    amount_of_reset_checks = int(input("Restart Hillclimb: Enter the amount of chain restart checks: "))
    amt_stab_change_checks = int(
        input("Restart Hillclimb: How many stability change checks before checking if chain should restart?: "))
    amt_pulls_per_stab_change_check = int(
        input("Restart Hillclimb: Enter the amount of nodepulls performed per stability change check: "))

    local_minima_chains, stability_over_time = hill_climber(protein, amount_of_reset_checks, amt_stab_change_checks,
                                                            amt_pulls_per_stab_change_check)

    # plot stability over time for hillclimb statistics
    data_plot_request = input("Restart Hillclimb: Do you want to plot the stability over time? (y/n): ")

    if data_plot_request == "y":
        data_plot_hillclimb(stability_over_time, protein)

    plot_chain_request(local_minima_chains)


def self_learning(proteinSeq):
    agent = None
    reward = 0
    stabilityProgress = []
    proteinSequence = proteinSeq

    # Create environment
    protein = Environment.Protein(proteinSequence)

    # Create agent
    agent = Agents.QLearningAgent()
    if agent.import_values(proteinSequence):
        print("Imported Q values")
    else:
        print("New Sequence")
    # showProtein(protein.state, protein.getStability())
    while not agent.terminate():
        while True:
            action = agent.PerceiveAndAct(protein, reward, protein.actions)
            if action is None: break
            reward = protein.ProcessAction(action)
        # Create new environment and reset agent
        stabilityProgress.append(protein.getStability())
        protein = Environment.Protein(proteinSequence)
        agent.reset()
    while True:
        action = agent.runBestSolution(protein, protein.actions)
        if action is None:
            break
        reward = protein.ProcessAction(action)
        agent.export_values(proteinSequence)
    showProtein(protein.state, protein.getStability(), stabilityProgress)


# Main userflow
def main():
    if len(sys.argv) != 3:
        print("Please use the format: python main.py [protein_string] [optimalization_type]")
        return

    optimalization_type = str(sys.argv[2])
    protein = str(sys.argv[1])

    # testing if only H, P, C combinations exist in protein.
    protein_test = [amino for amino in protein if (amino == "H" or amino == "P" or amino == "C")]

    if len(protein_test) != len(protein):
        print("Please only use H's, C's and P's for your protein.")
        return

    if optimalization_type == "R":
        # make random chain
        grid_object = Grid(protein)
        grid_object = random_chain(grid_object)

        # for plotting compatibility
        chains = {}
        chains[grid_object.stability] = [grid_object.grid, grid_object.grid_chain]

        plot_chain_request(chains)

    if optimalization_type == "SA":
        annealing_flow(protein)

    if optimalization_type == "RHC":
        restart_hill_climb_flow(protein)

    if optimalization_type == "SL":
        self_learning(protein)

    return


if __name__ == "__main__":
    main()
