from classes.GridClass import Grid
from algorithms.HillClimb import hill_climber, plot_best_c
from visualisation.DataPlots import hill_climb_plot

def main():
    # Grid class with node-pulling algorithm built in.
    grid_class = Grid("PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP")
    
    # hill_climber finds a local maximum, and resets to a new random chain when it finds it.
    # This makes it a Restart Hillclimbing Algorithm.
    # The more resets, the more it discovers local maximas, the higher the chance
    # to find the global maxima.
    # first argument = amount of reset checks you want to do
    # second argument = amount of times a whole chain should be pulled before the next check
    grid_class, stability_over_time = hill_climber(20, 50, grid_class)
    
    # plot stability over time for hillclimb statistics
    hill_climb_plot(stability_over_time)
    
    # 3D plot of the best chain of the hill_climber
    plot_best_c(grid_class)
    return
    
if __name__ == "__main__":
    main()