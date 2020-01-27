from classes.Grid import Grid
from classes.Node import Node

# creates a protein chain starting at (0, 0, 0)
# then filling in the grid with Nodes by doing a random move from last added
# Node to the next one until it has length of whole protein string.
def create_chain(grid_object):
    # n'th Node
    index = 0

    # Create grid of protein size
    grid_object.create_grid(int(len(grid_object.protein)))

    # Create first Node and set type as first protein char
    first_node = Node(0, 0, 0)
    first_node.type = grid_object.protein[0]

    # Put in Grid
    grid_object.add_point(first_node, index)

    # repeat until the whole protein has finished
    while index < len(grid_object.protein) - 1:

        # select last added Node
        current_node_key = grid_object.grid_chain[index][0]
        current_node = grid_object.grid[current_node_key].nodes[0]

        # choose random move
        random_move = np.random.randint(len(grid_object.moves))

        # make new coordinates of new Node
        new_x = current_node.x + grid_object.moves[random_move][0]
        new_y = current_node.y + grid_object.moves[random_move][1]
        new_z = current_node.z + grid_object.moves[random_move][2]

        # if the move doesnt overlap own chain -> add to grid
        if not grid_object.overlap(new_x, new_y, new_z):
            index += 1
            new_node = Node(new_x, new_y, new_z)
            new_node.n = index
            new_node.type = grid_object.protein[index]
            grid_object.add_point(new_node, index)

        # start over if chain is stuck
        if grid_object.chain_stuck(new_x, new_y, new_z):
            index = 0
            grid_object.clear_grid()
            grid_object.grid_chain = {}
            grid_object.add_point(first_node, index)

# inputs an protein chain and 
# returns a grid object filled with a randomly generated chain
def random_chain(protein):
    
    # create grid object
    grid_object = Grid(protein)
    
    # generate a chain in this grid
    create_chain(grid_object)
    
    # update bonds and set stability
    grid_object.update_all_bonds()
    
    return grid
