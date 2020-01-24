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
    # first argument = amount of resets you want to do (doesnt reset if still is a chance for improvement)
    # second argument = amount of times a whole chain should be pulled
    grid_class, stability_over_time = hill_climber(10, 100, grid_class)
    
    # plot stability over time for hillclimb statistics
    hill_climb_plot(stability_over_time)
    
    # 3D plot of the best chain of the hill_climber
    plot_best_c(grid_class)
    return
    
if __name__ == "__main__":
    main()