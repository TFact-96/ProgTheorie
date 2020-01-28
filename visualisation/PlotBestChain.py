from visualisation.Plot3D import Plot3D
from classes.Grid import Grid

# 3D plotting the best chain out of a dict
# in the form chain[stability] = grid_object
def plot_best_chain(best_chains):
    # choose best chain and grid out of the dict of chains
    best_stability = min(best_chains.keys())
    best_grid, best_chain = best_chains[best_stability]

    # prepare for plotting
    x = []
    y = []
    z = []
    color = []
    protein = []

    # get coords and node colors of chain
    for key, value in best_chain.items():
        node = best_grid[value[0]].nodes[0]

        x.append(node.x)
        y.append(node.y)
        z.append(node.z)

        # set amino colors
        if node.type == "H":
            color.append("red")

        elif node.type == "C":
            color.append("yellow")

        else:
            color.append("blue")

        # for calculating bonds soon
        protein.append(node.type)

    # create best grid object and calculate bonds
    best_grid_object = Grid("".join(protein))
    best_grid_object.grid = best_grid
    best_grid_object.grid_chain = best_chain
    best_grid_object.update_all_bonds()

    # 3D plot this chain
    Plot3D(
        x,
        y,
        z,
        best_grid_object.hh_bonds,
        best_grid_object.ch_bonds,
        best_grid_object.cc_bonds,
        best_grid_object.stability,
        color,
    )
