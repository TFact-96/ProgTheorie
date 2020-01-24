import matplotlib.pyplot as plt

def hill_climb_plot(stability_over_time):
    x = [i for i in range(len(stability_over_time))]
    
    plt.plot(x, stability_over_time)
    plt.xlabel("Chainlength pulls")
    plt.ylabel("Stability")
    plt.show()
    return