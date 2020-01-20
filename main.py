from classes.ChainLattice import ChainLattice
from algorithms.chaingeneration.chaingenerate import generate_chain
from algorithms.chaingeneration.multiplechains import multiple_chains
from algorithms.optimizingalgorithms.chainpulling import chain_pulling
from algorithms.optimizingalgorithms.evochainpull import find_best_pulled_chain, find_best_c
from algorithms.upperbound import calc_upperbound
from visualisation.data import get_chain_data, get_plot_data, get_chain_from_file, write_chain_to_csv
from visualisation.plot3D import plot_chain3D, plot_multiple_chains
import os

# Get this shitty greedy_tries out. It slows down the program by a significant deal.
greedy_tries = 20

def clear_terminal():
    os.system('cls' if os.name == 'nt' else 'clear')

def plotting_and_data_handler(Chain):
    plot = input("\nDo you want to plot the chain? (y/n): ")
    if (plot == "y"):
        plot_chain3D(Chain)

    save_data = input("Do you want to save the amino moves of the chain into a .csv file? (Warning: only saves the non-pulled version of the chain) (y/n): ")
    if (save_data == "y"):
        write_chain_to_csv(get_chain_data(Chain))

def get_user_input_for_generating_chain():
    protein = input("Enter desired protein-chain (C's, H's and P's): ")

    for amino in protein:
        if amino != "H" and amino != "C" and amino != "P":
            print("Only C, H and P aminos allowed.")
            exit(0)
    
    # calculate upperbound
    if protein:
        min_stability = calc_upperbound(protein)
        print(f"The naive minimal stability of this protein is {min_stability}\n")
        
    str_greedy = input("Random chain generation or greedy-move algorithm chain generation? (y = greedy / n = random): ")

    str_hillclimb = input("Chain-pulling for better stability (Hill Climb Algorithm) after generation of the chain? (y/n): ")

    if (str_greedy != "y" and str_greedy != "n") or (str_hillclimb != "y" and str_hillclimb != "n"):
        print("Only answer with y or n please.")
        exit(0)

    return protein, str_greedy, str_hillclimb

def main():
    clear_terminal()
    ################### Loading existing chain data or generating own?
    from_csv = input("Plot an existing protein chain from a .csv file, or generate a new one? (y = load / n = generate): ")

    if (from_csv == "y"):
        file = input("Give the filename from the data folder (without .csv extension): ")
        Chain = get_chain_from_file(file)
        
        if Chain:
            plot_chain3D(Chain)

        prompt_rerun()

    ################### Get chain data input for generation
    protein, str_greedy, str_hillclimb = get_user_input_for_generating_chain()

    ################### Chain generation

    # if no chain pulling
    if str_hillclimb == "n":
        # Random generation
        if str_greedy == "n":
            Chain = generate_chain(protein, False, greedy_tries)

        # Optimized (greedymoves) generation
        if str_greedy == "y":
            Chain = generate_chain(protein, True, greedy_tries)

    # Using Mehmet's hillclimb ChainPulling algorithm after generating chains
    if str_hillclimb == "y":
        chain_generations = int(input("How many chains should be generated to pull?: "))
        pull_times_per_chain = int(input("How many random pulls should be executed within each chain?: "))

        if str_greedy == "n":
            pulled_random_chains_list = chain_pulling(protein, False, greedy_tries, chain_generations, pull_times_per_chain)
            # just getting the best Chain
            Chain = pulled_random_chains_list[-1]
        else:
            pulled_greedy_chains_list = chain_pulling(protein, True, greedy_tries, chain_generations, pull_times_per_chain)
            Chain = pulled_greedy_chains_list[-1]

    # Handling plotting and data of chain
    if Chain:
        plotting_and_data_handler(Chain)

    prompt_rerun()

def prompt_rerun():
    rerun = 0
    while rerun != "y" and rerun != "n":
        rerun = input("Do you want to try another protein chain? (y/n):")

    if rerun == "y":
        main()

    exit(0)

if __name__ == "__main__":
    #main()
    find_best_pulled_chain("PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP", 100, True)

    # statistics plotting for quantifying quality of algorithm
    #protein = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
    #iterations = 500
    #use_optimize_algorithm = True
    #chain_nr, chain_stability = multiple_chains(protein, iterations, use_optimize_algorithm, greedy_tries)
    #plot_multiple_chains(chain_nr, chain_stability)
