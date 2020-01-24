from çlasses.GridClass import Grid
from algorithms.HillClimb import hill_climber, find_best_c

def main():
    grid_class = Grid("PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP")
    
    # first number = amount of resets you want to do
    # it makes it a Resetting Hillclimbing Algorithm
    hill_climber(10, grid_class)
    find_best_c(grid_class)
    return
    
if __name__ == "__main__":
    main()