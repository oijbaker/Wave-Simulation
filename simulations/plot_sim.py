import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import os
from matplotlib.patches import Circle

# print the path
print("Current working directory:", os.getcwd())

simulator = "cpp" # or "cpp"
object_type = "diffraction" # or "diffraction"
radius = 10
gap = 10

if simulator == "py":
    simulation_data = pd.read_csv("simulations/output/waves.csv", header=None)
    size = int(np.sqrt(simulation_data.shape[1]))  # Assuming square matrices
    plot_data = np.array([sim.reshape(size,size) for sim in simulation_data.values])
elif simulator == "cpp":
    simulation_data = pd.read_csv("src/output/wave_output.txt", sep=" ",header=None)
    size = int(np.sqrt(simulation_data[:-1].shape[1]))
    plot_data = np.array([sim[:-1].reshape(size,size) for sim in simulation_data.values])

print(plot_data.shape)



# animate an imshow of each matrix in plot_data

fig, ax = plt.subplots()
im = ax.imshow(plot_data[0], animated=True, cmap='ocean')
plt.colorbar(im, ax=ax, label='Wave Intensity')

if object_type == "circle":
    # create a circle patch
    obj = Circle((size//2, size//2), radius, color='black', fill=False, lw=2)
    ax.add_patch(obj)

    def update(frame):
        # print(f"Updating frame {frame}")
        im.set_array(plot_data[frame])
        return [im, obj]
    
elif object_type == "diffraction":
    # create a line patch
    obj1 = plt.Line2D((0, size//2-gap//2), (size//2, size//2), color='black', lw=2)
    obj2 = plt.Line2D((size//2+gap//2, size), (size//2, size//2), color='black', lw=2)
    ax.add_line(obj1)
    ax.add_line(obj2)

    def update(frame):
        # print(f"Updating frame {frame}")
        im.set_array(plot_data[frame])
        return [im, obj1, obj2]


# create the animation
ani = animation.FuncAnimation(fig, update, frames=len(plot_data), interval=1, blit=True)
plt.title("Wave Simulation")
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.show()