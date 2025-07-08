import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.ndimage import gaussian_filter

# Load data (columns: x, y, cloud)
data = np.loadtxt("simulation_output.txt")
x = data[:, 0]  # x position
y = data[:, 1]  # y position
cloud = data[:, 2]  # Cloud identifier (optional)

# Constants
num_particles = 100000  # Total particles per frame
num_frames = len(data) // num_particles

# Grid parameters for heatmap
grid_size = 1000
extent = [-3000, 1500, -500, 500]  # Doubled field of view

# Set up plot
fig, ax = plt.subplots(figsize=(6, 6))
img = ax.imshow(np.zeros((grid_size, grid_size)),
                extent=extent,
                origin='lower',
                cmap='coolwarm',
                vmin=0, vmax=2.5)  # Adjust vmax to control brightness

ax.set_xlim(-3000, 1500)  # Adjusted for larger scene
ax.set_ylim(-500, 500)
ax.set_title("Cloud-Cloud Collision Simulation")

def update(frame):
    start = frame * num_particles
    end = start + num_particles
    x_frame = x[start:end]
    y_frame = y[start:end]

    # Create histogram and apply strong blur
    heatmap, _, _ = np.histogram2d(
        x_frame, y_frame,
        bins=grid_size,
        range=[[-3000,1500], [-500, 500]]
    )
    blurred = gaussian_filter(heatmap, sigma=12)  # Increased blur for smoothness

    img.set_data(blurred.T)  # Transpose to match the image layout
    return img,

ani = animation.FuncAnimation(
    fig, update, frames=num_frames,
    blit=True, interval=10
)

# To save as a GIF, uncomment this line:
#ani.save("cloud_collision_blur.gif", writer='pillow', fps=40, dpi=180)

plt.tight_layout()
plt.show()
