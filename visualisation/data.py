import os
import csv
import time
from classes.ChainLattice import ChainLattice

# clearing the terminal
def clear_terminal():
    os.system('cls' if os.name == 'nt' else 'clear')

###################################### Returns resulting stability of whole chain and a tuple list of (amino, fold_code)
def get_chain_data(Chain):
    move_data = [(amino.type, amino.fold_code) for amino_key, amino in Chain.state.items()]
    return Chain.stability, move_data

##################################### Write move_data to a csv file with datestamp in name
def write_chain_to_csv(move_data):
    timestr = time.strftime("%Y%m%d-%H%M%S")
    filename = "AminoChain" + f"(S={move_data[0]})-" + timestr

    with open(f"data/{filename}.csv", mode='w') as chain:
        amino = csv.writer(chain, delimiter=',')

        amino.writerow(['amino', 'fold_code'])

        for amino in move_data[1]:
            amino.writerow(amino)

    print(f"File saved to {filename}.csv in the /data folder.")

###################################### Generate chain from existing CSV
def get_chain_from_file(file, ThreeD):
    protein = []
    fold_codes = []
    moves = []
    new_x = 0
    new_y = 0
    new_z = 0

    # check if file exists
    if not os.path.isfile(f"data/{file}.csv"):
        print("File not found.")
        return

    # read datafile and put amino string and fold codes into lists
    with open(f"data/{file}.csv", mode='r') as chain_data:
        chain_reader = csv.DictReader(chain_data)
        line_count = 1
        for amino in chain_reader:
            protein.append(amino["amino"])
            fold_codes.append(amino["fold_code"])
            line_count += 1

    # make Chain object with this amino
    Chain = ChainLattice("".join(protein), ThreeD)

    # get all move coordinates from the fold codes
    for fold in fold_codes:
        if fold != "0":
            moves.append(Chain.moves[Chain.fold_code_to_index[fold]])

    # create amino objects based on moves from the zeroeth and put them in the Chain chain
    for i in range(len(moves)):
        new_x += moves[i][0]
        new_y += moves[i][1]
        if ThreeD:
            new_z += moves[i][2]

        if ThreeD:
            new_coords = [new_x, new_y, new_z]
        else:
            new_coords = [new_x, new_y]

        last_amino = Chain.state[i]
        fold_code = fold_codes[i]

        amino = Chain.create_amino_object(new_coords, last_amino, fold_code)
        Chain.state[i + 1] = amino

    return Chain

###################################### Preparing lists for plotting amino chain
def get_plot_data(Chain, ThreeD):
    # prepare for plotting
    x = []
    y = []
    z = []
    color = []

    # set stability and bonds info
    Chain.set_stability_and_bonds()

    # append all aminos in coord list
    for amino_key, amino in Chain.state.items():
        x.append(amino.x)
        y.append(amino.y)

        if ThreeD:
            z.append(amino.z)

        # set amino colors
        if amino.type == "H":
            amino.color = 'red'
        elif amino.type == "C":
            amino.color = 'yellow'

        color.append(amino.color)

    if ThreeD:
        return x, y, z, color
    else:
        return x, y, color
