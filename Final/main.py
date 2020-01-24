from classes.GridClass import Grid
from algorithms.HillClimb import hill_climber, plot_best_c

def main():
    # Grid class with node-pulling algorithm built in.
    grid_class = Grid("PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP")
    
    # hill_climber finds a local maximum, and resets to a new random chain when it finds it.
    # This makes it a Resetting Hillclimbing Algorithm.
    # The more resets, the more it discovers local maximas, the higher the chance
    # to find the global maxima.
    # first argument = amount of resets you want to do.
    grid_class, stability_over_time = hill_climber(10, grid_class)
    
    # 3D plot of the best chain of the hill_climber
    plot_best_c(grid_class)
    return
    
if __name__ == "__main__":
    main()