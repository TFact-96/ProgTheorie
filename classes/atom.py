class Atom:
    # object for an atom in the amino chain
    def __init__(self, x, y):
        # nth atom
        self.n = 0

        # coords
        self.x = x
        self.y = y

        # atom type (C, H, P)
        self.type = type

        # the fold that this atom makes to the next atom (default 0)
        self.fold_code = 0

        # default color (P)
        self.color = "blue"

    # print the type if printing the object
    def __str__(self):
        return f"{self.type}"
