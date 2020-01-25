from visualisation.Plot3D import Plot3D

def plot_best_chain(best_chains):
    best_chain_key = min(best_chains.keys())
    best_chain_class = best_chains[best_chain_key]
    
    # prepare for plotting
    x = []
    y = []
    z = []
    color = []

    # get coords and node colors of chain
    for key, value in best_chain_class.grid_chain.items():
        node = best_chain_class.grid[value[0]].nodes[0]

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
            
    Plot3D(x, y, z, best_chain_class.hh_bonds, best_chain_class.stability, color)