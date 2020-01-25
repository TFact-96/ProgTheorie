from visualisation.Plot3D import Plot3D

def plot_best_chain_hillclimb(grid_class):
    best_chain_key = min(grid_class.best_chain.keys())
    best_chain_double = grid_class.best_chain[best_chain_key]

    best_chain = best_chain_double[0]
    best_grid = best_chain_double[1]

    best_stability, best_hh = grid_class.update_neighbours(best_grid, best_chain)
    
    # prepare for plotting
    x = []
    y = []
    z = []
    color = []

    # get coords and node colors of chain
    for key, value in best_chain.items():
        node = best_grid[value[0]].nodes[0]

        x.append(node.x)
        y.append(node.y)
        z.append(node.z)

        # set amino colors
        if node.type == "H":
            color.append('red')
            
        elif node.type == "C":
            color.append('yellow')
        
        else:
            color.append('blue')
            
    Plot3D(x, y, z, best_hh, best_stability, color)
