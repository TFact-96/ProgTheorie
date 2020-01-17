class Amino3D:
    # object for an atom in the amino chain
    def __init__(self, x, y, z):
        # nth atom in the chain
        self.n = 0

        # coords
        self.x = x
        self.y = y
        self.z = z

        # amino type (default P)
        self.type = "P"

        # the fold that this atom makes to the next atom (default 0)
        self.fold_code = 0

        # default color (P)
        self.color = 'blue'

        # All neighbours of a single node
        self.neighbours = []

    # print the type if printing the object
    def __str__(self):
        return f"{self.type}"

class Amino2D:
    # object for an atom in the amino chain
    def __init__(self, x, y):
        # nth atom in the chain
        self.n = 0

        # coords
        self.x = x
        self.y = y

        # amino type (default P)
        self.type = "P"

        # the fold that this atom makes to the next atom (default 0)
        self.fold_code = 0

        # default color (P)
        self.color = 'blue'

    # print the type if printing the object
    def __str__(self):
        return f"{self.type}"
