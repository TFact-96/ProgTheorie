from classes.Grid import Grid
from algorithms.HillClimb import hill_climber
from visualisation.DataPlots import data_plot_hillclimb
from visualisation.PlotBestHillClimb import plot_best_chain_hillclimb

def main():
    # Grid class with node-pulling algorithm built in.
    grid_class = Grid("PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP")
    
    # hill_climber finds a local maximum, and resets to a new random chain when it finds it.
    # This makes it a Restart Hillclimbing Algorithm.
    # The more resets, the more it discovers local maximas, the higher the chance
    # to find the global maxima.
    # first argument = amount of times a whole chain should be pulled before the next check
    # second argument = amount of reset checks you want to do after pulling the chain n times
    grid_class, stability_over_time = hill_climber(50, 10, grid_class)
    
    # plot stability over time for hillclimb statistics
    data_plot_hillclimb(stability_over_time, grid_class.amino)
    
    # 3D plot of the best chain of the hill_climber
    plot_best_chain_hillclimb(grid_class)
    return
    
if __name__ == "__main__":
    main()