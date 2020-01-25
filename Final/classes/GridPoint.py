class GridPoint:
    def __init__(self, filled, coords):
        self.filled = filled
        self.coords = coords
        self.nodes = []

    def add_node(self, node):
        self.nodes.append(node)
        self.filled = True

    def remove_node(self, node):
        self.nodes.remove(node)
        if self.nodes == []:
            self.filled = False

    def __repr__(self):
        return f"{self.nodes}"
