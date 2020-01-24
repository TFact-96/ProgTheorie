import matplotlib.pyplot as plt

def hill_climb_plot(stability_over_time):    
    plt.plot(stability_over_time)
    plt.xlabel("Chainlength pulls")
    plt.ylabel("Stability")
    plt.title("Restart Hill Climb Algorithm")
    plt.show()
    return