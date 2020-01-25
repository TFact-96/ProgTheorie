from classes.Grid import Grid
from algorithms.HillClimb import hill_climber
from algorithms.SimAnnealing import simulated_annealing
from visualisation.DataPlots import data_plot_hillclimb, data_plot_annealing
from visualisation.PlotBestChain import plot_best_chain

annealing = False
hill_climb = True
amino = "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP"

def main():
    if annealing:
        best_chains, stability_over_time = simulated_annealing(amino, 1000, 2, False,
                True, 0.001, 0.997)


        # plot stability over time for hillclimb statistics
        data_plot_annealing(stability_over_time, amino)
        
        # 3D plot of the best chain of the hill_climber
        plot_best_chain(best_chains)
            
    if hill_climb:
        # Grid class with node-pulling algorithm built in.
        grid_class = Grid(amino)
        
        # hill_climber finds a local maximum, and resets to a new random chain when it finds it.
        # This makes it a Restart Hillclimbing Algorithm.
        # The more resets, the more it discovers local maximas, the higher the chance
        # to find the global maxima.
        # first argument = amount of times a whole chain should be pulled before the next check
        # second argument = amount of reset checks you want to do after pulling the chain n times
        best_chains, stability_over_time = hill_climber(amino, 10, 20)
        
        # plot stability over time for hillclimb statistics
        data_plot_hillclimb(stability_over_time, amino)
        
        # 3D plot of the best chain of the hill_climber
        plot_best_chain(best_chains)
    return
    
if __name__ == "__main__":
    main()