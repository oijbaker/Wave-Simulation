import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import os
# print the path
print("Current working directory:", os.getcwd())

#simulation_data = pd.read_csv("simulations/output/waves.csv", header=None)
simulation_data = pd.read_csv("src/output/wave_output.txt", sep=" ",header=None)
size = int(np.sqrt(simulation_data.shape[1]))  # Assuming square matrices
plot_data = np.array([sim[:-1].reshape(size,size) for sim in simulation_data.values])
print(plot_data.shape)

# animate an imshow of each matrix in plot_data

fig, ax = plt.subplots()
im = ax.imshow(plot_data[0], animated=True, cmap='ocean')
plt.colorbar(im, ax=ax, label='Wave Intensity')

# function to update the imshow
def update(frame):
    # print(f"Updating frame {frame}")
    im.set_array(plot_data[frame])
    return [im]

# create the animation
ani = animation.FuncAnimation(fig, update, frames=len(plot_data), interval=1, blit=True)
plt.title("Wave Simulation")
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.show()