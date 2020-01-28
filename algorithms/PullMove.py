import numpy as np
import random
from classes.GridPoint import GridPoint

# check if diagonal coordinates from a point is not filled
# returns all available coordinates
def check_diagonals(grid_object, x, y, z):
    """
    Check if diagonal coordinates from a point is not filled, returns all available coordinates.
    Arguments: grid_object, x, y, z
    """
    # diagonal moves for chain pulling
    available_moves = []

    for move in grid_object.diagonal_moves:
        # cant overflow the grid
        if f"{x + move[0], y + move[1], z + move[2]}" not in grid_object.grid:
            grid_object.grid[f"{x + move[0], y + move[1], z + move[2]}"] = GridPoint(
                False, [x + move[0], y + move[1], z + move[2]]
            )

        # if its not filled
        if not grid_object.grid[f"{x + move[0], y + move[1], z + move[2]}"].filled:
            available_moves.append(move)

    return available_moves


def create_vectors(grid_object, node):
    """
    Returns coords of i, i+1 (or i-1) and the vector between these coords. Arguments: grid_object, node
    """
    node_i_coords = [node.x, node.y, node.z]

    # if node at left end of chain, use i+1 to fold into the middle
    if node.n < (len(grid_object.protein) / 2):
        node_i1 = grid_object.grid_chain[int(node.n) + 1]

    # right end of chain use i-1 to fold into the middle
    else:
        node_i1 = grid_object.grid_chain[int(node.n) - 1]

    # get coords and vector
    node_i1_coords = np.array([node_i1[1][0], node_i1[1][1], node_i1[1][2]])
    vector1 = node_i1_coords - np.array(node_i_coords)

    return node_i_coords, node_i1_coords, vector1


def check_requirements(
    grid_object, available_moves, vector1, node_i_coords, node_i1_coords
):
    """
    Check requirements if L and C are free --> returns coords of L, coords of C, and True if so.
    Ff multiple L and C's are free, returns a random choice of them.
    Arguments: grid_object, available_moves, vector1, node_i_coords, node_i1_coords
    """
    viable_moves = []
    found = False

    for move in available_moves:
        L = node_i_coords + np.array(move)
        C = L - vector1

        if (
            (not grid_object.overlap(L[0], L[1], L[2]))
            and (not grid_object.overlap(C[0], C[1], C[2]))
            and (np.linalg.norm(L - node_i1_coords) == 1.0)
        ):
            viable_moves.append([L, C])
            found = True

    if not found:
        return 0, 0, False

    random_choice = random.choice(viable_moves)
    return random_choice[0], random_choice[1], True


def move_residue_left(index, grid_object):
    """
    Function to move residue in first half of chain. Arguments: index, grid_object
    """

    residue_node_key = grid_object.grid_chain[index][0]
    residue_node = grid_object.grid[residue_node_key].nodes[0]

    residue_node_next_key = grid_object.grid_chain[index + 2][0]
    residue_node_next = grid_object.grid[residue_node_next_key].nodes[0]

    grid_object.transfer_point(
        residue_node, residue_node_next.x, residue_node_next.y, residue_node_next.z,
    )


def move_residue_right(grid_object, node):
    """
    Function to move residue in second half of chain. Arguments: grid_object, node
    """
    index_from_end = len(grid_object.protein)
    # residue of chain follows in footsteps
    while index_from_end > node.n + 2:
        index_from_end -= 1
        residue_node_key = grid_object.grid_chain[index_from_end][0]
        residue_node = grid_object.grid[residue_node_key].nodes[0]

        residue_node_next_key = grid_object.grid_chain[index_from_end - 2][0]
        residue_node_next = grid_object.grid[residue_node_next_key].nodes[0]

        grid_object.transfer_point(
            residue_node, residue_node_next.x, residue_node_next.y, residue_node_next.z,
        )


# pulling a node in the grid_object diagonally. Always towards the middle of the chain.
def pull_move(grid_object, node):
    """
    The main pull move. Arguments: grid_object, node
    """

    node_i_coords, node_i1_coords, vector1 = create_vectors(grid_object, node)

    available_moves = check_diagonals(
        grid_object, node_i_coords[0], node_i_coords[1], node_i_coords[2]
    )

    L, C, check = check_requirements(
        grid_object, available_moves, vector1, node_i_coords, node_i1_coords
    )

    if check:
        # For left side of chain folding towards the middle
        if node.n < (len(grid_object.protein) / 2):
            # residue of chain follows in footsteps
            for index in range(int(node.n - 1)):

                move_residue_left(index, grid_object)

            # Previous node moves to C
            previous_node_key = grid_object.grid_chain[int(node.n) - 1][0]
            previous_node = grid_object.grid[previous_node_key].nodes[0]

            grid_object.transfer_point(previous_node, C[0], C[1], C[2])

        # for right side of chain folding towards the middle
        else:

            move_residue_right(grid_object, node)

            # Previous node moves to C
            previous_node_key = grid_object.grid_chain[int(node.n) + 1][0]
            previous_node = grid_object.grid[previous_node_key].nodes[0]

            grid_object.transfer_point(previous_node, C[0], C[1], C[2])

        # node moves to L
        grid_object.transfer_point(node, L[0], L[1], L[2])
