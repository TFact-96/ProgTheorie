from classes.AminoLattice import AminoLattice
from algorithms.chaingeneration.chaingenerate import generate_chain
from algorithms.chaingeneration.bruteforce import bruteforce_chains
from algorithms.chaingeneration.multiplechains import multiple_chains
from visualisation.data import get_chain_data, get_plot_data, get_chain_from_file, write_chain_to_csv
from visualisation.plot3D import plot_chain2D, plot_chain3D, plot_multiple_chains
import os

# Get this shitty optimalization_tries out. It slows down the program by a significant deal.
optimalization_tries = 20

def clear_terminal():
    os.system('cls' if os.name == 'nt' else 'clear')

def plotting_and_data_handler(lattice):
    plot = input("\nDo you want to plot the chain? (y/n): ")
    if (plot == "y"):
        if ThreeD:
            plot_chain3D(lattice)
        else:
            plot_chain2D(lattice)

    save_data = input("Do you want to save the atom moves of the chain into a .csv file? (y/n): ")
    if (save_data == "y"):
        write_chain_to_csv(get_chain_data(lattice))

def get_user_input_for_generating_chain():
    amino = input("Enter desired amino-chain (C's, H's and P's): ")

    for atom in amino:
        if atom != "H" and atom != "C" and atom != "P":
            print("Only C, H and P atoms allowed.")
            exit(0)

    str_optimize = input("Random generation or greedy-move algorithm optimization? (y = greedy / n = random): ")
    str_brute = input("Brute force generate chains and find best one? (y/n): ")

    if (str_optimize != "y" and str_optimize != "n") or (str_brute != "y" and str_brute != "n"):
        print("Only answer with y or n please.")
        exit(0)

    if str_brute == "y":
        try:
            iterations = int(input("How many chain generations for brute forcing?: "))
        except ValueError:
            print("Only positive integer for brute force iterations please.")
            exit(0)
        
        if iterations < 0:
            print("Only positive integer for brute force iterations please.")
            exit(0)
            
        return iterations, amino, str_optimize, str_brute

    return 0, amino, str_optimize, str_brute

def main():
    clear_terminal()

    ################### Loading existing chain data or generating own?
    from_csv = input("Plot an existing amino chain from a .csv file, or generate a new one? (y = load / n = generate): ")

    if (from_csv == "y"):
        file = input("Give the filename from the data folder (without .csv extension): ")
        lattice = get_chain_from_file(file)
        if lattice:
            if ThreeD:
                plot_chain3D(lattice)
            else:
                plot_chain2D(lattice)

        prompt_rerun()

    ################### Get chain data input for generation
    iterations, amino, str_optimize, str_brute = get_user_input_for_generating_chain()

    ################### Chain generation

    # Random generation
    if str_optimize == "n":
        # bruteforce
        if str_brute == "y":
            lattice = bruteforce_chains(amino, iterations, False, 0, ThreeD)
        # non bruteforce
        else:
            lattice = generate_chain(amino, False, 0, ThreeD)

    # Optimized (greedymoves) generation
    if str_optimize == "y":
        # bruteforce
        if str_brute == "y":
            lattice = bruteforce_chains(amino, iterations, True, optimalization_tries, ThreeD)
        # not bruteforce
        else:
            lattice = generate_chain(amino, True, optimalization_tries, ThreeD)

    # Handling plotting and data of chain
    if lattice:
        plotting_and_data_handler(lattice)
    
    prompt_rerun()

def prompt_rerun():
    rerun = 0
    while rerun != "y" and rerun != "n":
        rerun = input("Do you want to try another amino chain? (y/n):")
    
    if rerun == "y":
        main()
    
    exit(0)
    
if __name__ == "__main__":
    
    # statistics plotting for quantifying quality of algorithm
    #amino = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
    #iterations = 500
    #use_optimize_algorithm = True
    #chain_nr, chain_stability = multiple_chains(amino, iterations, use_optimize_algorithm, optimalization_tries)
    #plot_multiple_chains(chain_nr, chain_stability)
    
    # Mehmet's hill climb
    ThreeD = False
    greedy = True
    greedy_chain_generations = 100
    pull_times_per_chain = 100
    amino = "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP"
    best_stability = 0
    best_chain_list = []
    
    for i in range(greedy_chain_generations):
        # generate_chain(amino, Greedymove=True, optimalization_tries bij het greedy genereren, 3D=True/2D=False)
        k = generate_chain(amino, greedy, optimalization_tries, ThreeD)
        
        if k:
            k.set_stability_and_bonds()
            print(f"Generated greedy chain nr {i}: Stability {k.stability}")

            k.random_pull(pull_times_per_chain)
            k.set_stability_and_bonds()
            print(f"After pulling this chain: Stability {k.stability}\n")
            
            if k.stability < best_stability:
                best_k = k
                best_chain_list.append(k.chain)
                best_stability = k.stability
                print(f"That was a record chain.\n")    
                
    print(f"Plotting best chain")
    best_chain_list[-1].plot_chain()       
    #main()