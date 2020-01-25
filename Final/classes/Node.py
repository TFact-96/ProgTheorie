class Node:
    # object for an node in the amino chain
    def __init__(self, x, y, z):

        # coords
        self.x = x
        self.y = y
        self.z = z
        self.n = 0

        # atom type (default P)
        self.type = "P"

        # the fold that this atom makes to the next atom (default 0)
        self.fold_code = 0

        # All neighbours of a single node
        self.neighbours = []

    # print the type if printing the object
    def __repr__(self):
        return f"{self.type}"