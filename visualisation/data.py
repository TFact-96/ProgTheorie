###################################### Returns resulting stability of whole chain and a tuple list of fold codes per atom
def get_chain_data(lattice, print_data):
    # cant get data if chain was stuck
    if lattice.chain_stuck:
        return

    # set stability level and bonds of the chain
    lattice.set_stability_and_bonds()

    moves = [[node.type, node.fold_code] for node in lattice.chain]

    if print_data:
        print(f"\nThe stability of this amino-acid is {lattice.stability}.\n")
        print(f"Its moves are:")

        for node in moves:
            print(f"{node[0]} {node[1]}")

    return lattice.stability, moves

###################################### Preparing lists for plotting amino chain
def get_plot_data(lattice):
    # prepare for plotting
    x = []
    y = []
    z = []
    color = []

    # set stability and bonds info
    lattice.set_stability_and_bonds()

    # append all nodes in coord list
    for node in lattice.chain:
        x.append(node.x)
        y.append(node.y)
        z.append(node.z)

        # set node colors
        if node.type == "H":
            node.color = 'red'
        elif node.type == "C":
            node.color = 'orange'

        color.append(node.color)

    return x, y, z, color
