from visualisation.Plot3D import Plot3D

# 3D plotting the best chain out of a dict
# in the form chain[stability] = grid_object
def plot_best_chain(best_chains):
    # choose best chain out of the dict of chains
    best_stability = min(best_chains.keys())
    best_grid_object = best_chains[best_stability]
    
    # prepare for plotting
    x = []
    y = []
    z = []
    color = []

    # get coords and node colors of chain
    for key, value in best_grid_object.grid_chain.items():
        node = best_grid_object.grid[value[0]].nodes[0]

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
    
    # 3D plot this chain
    Plot3D(x, y, z, best_grid_object.hh_bonds, best_grid_object.ch_bonds, best_grid_object.cc_bonds, best_grid_object.stability, color)