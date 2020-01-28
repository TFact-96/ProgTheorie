import numpy as np
from classes.Grid import Grid
from classes.Node import Node


def make_first_node(grid_object, index):
    """
    Create first Node and set type as first protein character.
    :param grid_object: grid object
    :param index: index
    :return: first node
    """
    first_node = Node(0, 0, 0)
    first_node.type = grid_object.protein[0]

    # Put in Grid
    grid_object.add_point(first_node, index)

    return first_node


def make_new_coords(grid_object, current_node):
    """
    Choose random move for a node.
    :param grid_object: grid object
    :param current_node: current node
    """
    random_move = np.random.randint(len(grid_object.moves))

    return (
        current_node.x + grid_object.moves[random_move][0],
        current_node.y + grid_object.moves[random_move][1],
        current_node.z + grid_object.moves[random_move][2],
    )


def random_chain(grid_object):
    """
    Creates a protein chain starting at (0, 0, 0)
    then filling in the grid with Nodes by doing a random move from last added
    Node to the next one until it has length of whole protein string.
    :param grid_object: grid object
    :return: grid object
    """

    # n'th Node
    index = 0

    # Clear grid if it exists
    if not len(grid_object.grid) != 0:
        grid_object.clear_grid()
        grid_object.grid_chain = {}

    first_node = make_first_node(grid_object, index)

    # repeat until the whole protein has finished
    while index < len(grid_object.protein) - 1:

        # select last added Node
        current_node_key = grid_object.grid_chain[index][0]
        current_node = grid_object.grid[current_node_key].nodes[0]

        new_x, new_y, new_z = make_new_coords(grid_object, current_node)

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

    grid_object.update_all_bonds()

    return grid_object
