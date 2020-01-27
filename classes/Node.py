# Object for an Amino in the Protein Chain (Node in the Grid)
class Node:
    def __init__(self, x, y, z):

        # coords
        self.x = x
        self.y = y
        self.z = z
        self.n = 0

        # amino type (default P)
        self.type = "P"

    # print the type if printing the object
    def __repr__(self):
        return f"{self.type}"