from classes.Grid import Grid

# creates a random chain
def random_chain(amino):
    chains = {}
    
    # set grid class
    grid = Grid(amino)
    
    # generate a chain in this grid
    grid.create_chain()
    
    # update neighbor bonds and set stability
    grid.update_neighbours()
    
    return grid
