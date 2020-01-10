from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from visualisation.data import get_plot_data

###################################### Plotting the chain
def plot_chain(lattice):
    x, y, z, color = get_plot_data(lattice)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # plot the bond lines
    for hh_bond in lattice.hh_bonds:
        ax.plot3D(hh_bond[0], hh_bond[1], hh_bond[2], "--", markersize=1, color='red', zorder=1)
        #plt.text((hh_bond[0][0] + hh_bond[0][1]) / 2, (hh_bond[1][0] + hh_bond[1][1]) / 2, "-1")

    for ch_bond in lattice.ch_bonds:
        ax.plot3D(ch_bond[0], ch_bond[1], ch_bond[2], "--", markersize=1, color='orange', zorder=1)

    for cc_bond in lattice.cc_bonds:
        ax.plot3D(cc_bond[0], cc_bond[1], cc_bond[2], "--", markersize=1, color='yellow', zorder=1)

    # plot the chain itself
    ax.plot3D(x, y, z, "-", linewidth=3, color='black', zorder=1, label=f"Amino chain (Stability: {lattice.stability})")
    ax.scatter3D(x, y, z, color=color, zorder=2)

    # plot atomtype name at its node coord
    for i in range(len(x)):
        ax.text(x[i] + 0.05, y[i] + 0.05, z[i] + 0.05, lattice.chain[i].type)

    plt.show()

    # rest info and show plot
