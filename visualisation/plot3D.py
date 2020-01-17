from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from visualisation.data import get_plot_data
import numpy as np

###################################### Plotting the chain
def plot_chain3D(Chain):
    x, y, z, color = get_plot_data(Chain, True)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # plot the bond lines
    for hh_bond in Chain.hh_bonds:
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

    for ch_bond in Chain.ch_bonds:
        ax.plot3D(
            ch_bond[0], ch_bond[1], ch_bond[2],
            "--",
            markersize=1,
            color='orange',
            zorder=-1
        )

    for cc_bond in Chain.cc_bonds:
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
                        Line2D([0], [0], color='black', lw=2, label="Protein chain"),
                        Line2D([0], [0], marker='o', color='black', label='H-amino',
                               markerfacecolor='red', markersize=8),
                        Line2D([0], [0], marker='o', color='black', label='P-amino',
                               markerfacecolor='blue', markersize=8),
                        Line2D([0], [0], marker='o', color='black', label='C-amino',
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
    ax.set_title(f"Protein chain (Stability: {Chain.stability})")
    plt.show()

def plot_chain2D(Chain):
    x, y, color = get_plot_data(Chain, False)


    # plot the bond lines
    for hh_bond in Chain.hh_bonds:
        plt.plot(
            hh_bond[0], hh_bond[1],
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

    for ch_bond in Chain.ch_bonds:
        plt.plot(
            ch_bond[0], ch_bond[1],
            "--",
            markersize=1,
            color='orange',
            zorder=-1
        )

    for cc_bond in Chain.cc_bonds:
        plt.plot(
            cc_bond[0], cc_bond[1],
            "--",
            markersize=1,
            color='yellow',
            zorder=-1
        )

    # plot text at start and ending atom
    plt.text(
         x[0], y[0] + 0.1,
         "Start"
    )

    plt.text(
         x[-1], y[-1] + 0.1,
         "Finish"
    )

    # plot the chain itself
    plt.plot(
        x, y,
        "-",
        linewidth=3,
        color='black',
        zorder=0,
    )

    plt.scatter(
        x, y,
        color=color,
        edgecolor='black',
        s=100,
    )

    custom_legend = [
                        Line2D([0], [0], color='black', lw=2, label="Protein chain"),
                        Line2D([0], [0], marker='o', color='black', label='H-amino',
                               markerfacecolor='red', markersize=8),
                        Line2D([0], [0], marker='o', color='black', label='P-amino',
                               markerfacecolor='blue', markersize=8),
                        Line2D([0], [0], marker='o', color='black', label='C-amino',
                               markerfacecolor='yellow', markersize=8),
                        Line2D([0], [0], linestyle="--", color='red', lw=1, label="H-H (Stability -1)"),
                        Line2D([0], [0], linestyle="--", color='orange', lw=1, label="H-C (Stability -1)"),
                        Line2D([0], [0], linestyle="--", color='yellow', lw=1, label="C-C (Stability -5)"),

                    ]

    # plt axis and grid off
    plt.axis('off')
    plt.grid(b=None)

    plt.legend(handles=custom_legend)
    plt.title(f"Protein chain (Stability: {Chain.stability})")
    plt.show()

def plot_multiple_chains(chain_nr, chain_data):
    x = chain_nr
    y = chain_data

    average = np.average(chain_data)
    standard_dev = float(np.std(chain_data))

    # 95% chance the stability is within this interval
    TwoStdDevAway = average - (2 * standard_dev)
    PositiveTwoStdAway = average + (2 * standard_dev)

    # stability cant be over zero
    if PositiveTwoStdAway > 0:
        PositiveTwoStdAway = 0

    plt.hlines(average, chain_nr[0], chain_nr[-1], zorder=1, colors='red', linewidth=3, label="Average stability: %.2f" % (average))
    plt.hlines(PositiveTwoStdAway, chain_nr[0], chain_nr[-1], zorder=1, colors='orange', linestyles="--", label="95 percent confidence interval (std. dev = %.2f )" % (float(standard_dev)))
    plt.hlines(TwoStdDevAway, chain_nr[0], chain_nr[-1], zorder=1, colors='orange', linestyles="--")
    plt.plot(x, y, zorder=0, label="Stabilities of multiple chains")
    plt.xlabel("Chain number")
    plt.ylabel("Stability")
    plt.ylim(1, -35)
    plt.title("Generations of multiple amino chains")
    plt.legend()
    plt.show()
