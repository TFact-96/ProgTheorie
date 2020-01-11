import os
import csv
import time
from classes.AminoLattice import AminoLattice

# clearing the terminal
def clear_terminal():
    os.system('cls' if os.name == 'nt' else 'clear')

###################################### Returns resulting stability of whole chain and a tuple list of (atom, fold_code)
def get_chain_data(lattice):
    # cant get data if chain was stuck
    if lattice.chain_stuck:
        return

    # set final stability level and bonds of the chain in the object
    lattice.set_stability_and_bonds()

    move_data = [(node.type, node.fold_code) for node in lattice.chain]

    return lattice.stability, move_data


##################################### Write move_data to a csv file with datestamp in name
def write_chain_to_csv(move_data):
    timestr = time.strftime("%Y%m%d-%H%M%S")
    filename = "AminoChain" + f"(S={move_data[0]})-" + timestr

    with open(f"data/{filename}.csv", mode='w') as chain:
        atom = csv.writer(chain, delimiter=',')

        atom.writerow(['atom', 'fold_code'])

        for node in move_data[1]:
            atom.writerow(node)

###################################### Generate chain from existing CSV
def get_chain_from_file(file):
    amino = []
    fold_codes = []
    moves = []
    new_x = 0
    new_y = 0
    new_z = 0

    # read datafile and put amino string and fold codes into lists
    with open(f"data/{file}.csv", mode='r') as chain_data:
        chain_reader = csv.DictReader(chain_data)
        line_count = 1
        for atom in chain_reader:
            amino.append(atom["atom"])
            fold_codes.append(atom["fold_code"])
            line_count += 1

    # make lattice object with this amino
    lattice = AminoLattice("".join(amino))

    # get all move coordinates from the fold codes
    for fold in fold_codes:
        if fold != "0":
            moves.append(lattice.moves[lattice.fold_code_to_index[fold]])

    # create atom objects based on moves from the zeroeth and put them in the lattice chain
    for i in range(len(moves)):
        new_x += moves[i][0]
        new_y += moves[i][1]
        new_z += moves[i][2]

        new_coords = [new_x, new_y, new_z]
        last_atom = lattice.chain[-1]
        fold_code = fold_codes[i]

        atom = lattice.create_atom_object(new_coords, last_atom, fold_code)
        lattice.chain.append(atom)

    return lattice

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
            node.color = 'yellow'

        color.append(node.color)

    return x, y, z, color
