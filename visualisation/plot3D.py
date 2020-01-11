from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from visualisation.data import get_plot_data

###################################### Plotting the chain
def plot_chain(lattice):
    x, y, z, color = get_plot_data(lattice)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # plot the bond lines
    for hh_bond in lattice.hh_bonds:
        ax.plot3D(
            hh_bond[0], hh_bond[1], hh_bond[2],
            "--",
            markersize=1,
            color='red',
            zorder=-1
        )

        # Plot how much a bond reduces stability
        #plt.text(
            # (hh_bond[0][0] + hh_bond[0][1]) / 2,
            # (hh_bond[1][0] + hh_bond[1][1]) / 2,
            # (hh_bond[2][0] + hh_bond[2][1]) / 2,
            # "-1"
        # )

    for ch_bond in lattice.ch_bonds:
        ax.plot3D(
            ch_bond[0], ch_bond[1], ch_bond[2],
            "--",
            markersize=1,
            color='orange',
            zorder=-1
        )

    for cc_bond in lattice.cc_bonds:
        ax.plot3D(
            cc_bond[0], cc_bond[1], cc_bond[2],
            "--",
            markersize=1,
            color='yellow',
            zorder=-1
        )

    # plot text at start and ending atom
    ax.text(
         x[0], y[0], z[0] + 0.1,
         "Start"
    )

    ax.text(
         x[-1], y[-1], z[-1] + 0.1,
         "Finish"
    )

    # plot the chain itself
    ax.plot3D(
        x, y, z,
        "-",
        linewidth=3,
        color='black',
        zorder=0,
    )

    ax.scatter3D(
        x, y, z,
        color=color,
        edgecolor='black',
        s=100,
        depthshade=False
    )

    custom_legend = [
                        Line2D([0], [0], color='black', lw=2, label="Amino chain"),
                        Line2D([0], [0], marker='o', color='black', label='H-atom',
                               markerfacecolor='red', markersize=8),
                        Line2D([0], [0], marker='o', color='black', label='P-atom',
                               markerfacecolor='blue', markersize=8),
                        Line2D([0], [0], marker='o', color='black', label='C-atom',
                               markerfacecolor='yellow', markersize=8),
                        Line2D([0], [0], linestyle="--", color='red', lw=1, label="H-H (Stability -1)"),
                        Line2D([0], [0], linestyle="--", color='orange', lw=1, label="H-C (Stability -1)"),
                        Line2D([0], [0], linestyle="--", color='yellow', lw=1, label="C-C (Stability -5)"),

                    ]

    # Hide grid lines
    ax.grid(False)

    # Hide axes ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    # plt axis and grid off
    plt.axis('off')
    plt.grid(b=None)

    ax.legend(handles=custom_legend)
    ax.set_title(f"Amino molecule chain (Stability: {lattice.stability})")
    plt.show()

    # rest info and show plot
