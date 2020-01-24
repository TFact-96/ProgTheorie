from GridClass import Grid
from HillClimb import hill_climber, find_best_c

def main():
    grid_class = Grid("PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP")
    hill_climber(1, grid_class)
    find_best_c(grid_class)
    return
    
if __name__ == "__main__":
    main()