import copy
from visualisation.plot3D import plot3D

def hill_climber(max_iteration, grid_class):
    
    # Current protein chain
    current_hilltop, grid = grid_class.create_chain()        
    current_stability = grid_class.update_neighbours(grid, current_hilltop)[0]
    
    # Save as best chain
    best_c = current_hilltop
    best_stab_c = current_stability
    best_grid = grid
    grid_class.best_chain[best_stab_c] = [best_c, best_grid]
    
    # for statistic plotting
    stability_over_time = []
    
    # try iteration amount of random chains
    for iteration in range(max_iteration):
        best_c_found = False
        print(f"try: {iteration}")

        # amount the >whole< chain should be pulled
        for it in range(100):
            print(f"pulling whole chain iteration: {it}")
            
            # try to pull each node
            for index in range(1, len(current_hilltop) - 1):
                print(f"pullmove on node: {index}")
                temp_grid, temp_chain = grid_class.pull_move(
                    grid[current_hilltop[index][0]].nodes[0], grid, current_hilltop,
                )

                temp_stability = grid_class.update_neighbours(temp_grid, temp_chain)[0]
                
                # if stability is better after pull, save this as best chain
                if temp_stability < best_stab_c and temp_stability < grid_class.pivot_upperbound:
                    print("pullmove made better stab!")
                    best_c = copy.deepcopy(temp_chain)
                    best_grid = copy.deepcopy(temp_grid)
                    best_stab_c = temp_stability
                    best_c_found = True
                    stability_over_time.append(best_stab_c)
        
        # go further with this chain if a better stability was found after 'it'-amount of iterations
        if best_c_found:
            current_hilltop = copy.deepcopy(best_c)
            current_stability = best_stab_c
            grid = copy.deepcopy(best_grid)
            best_c_found = False

        # if no stability change after this amount, must be local maximum. Save in best_chain list
        else:
            grid_class.best_chain[best_stab_c] = [best_c, best_grid]

            current_hilltop, grid = grid_class.create_chain()
            current_stability = grid_class.update_neighbours(grid, current_hilltop)[0]

            best_c = current_hilltop
            best_stab_c = current_stability
            best_grid = grid
    
    return grid_class, stability_over_time


def plot_best_c(grid_class):
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
            
    plot3D(x, y, z, best_hh, best_stability, color)

