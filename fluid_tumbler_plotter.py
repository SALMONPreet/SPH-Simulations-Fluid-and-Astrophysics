import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

# Configuration
num_particles = 600
box_size = 15.0
log_interval = 10  # Data was saved every 10 simulation steps
save_steps = [10,50,200, 1000]
save_frames = [step // log_interval for step in save_steps]

# Load data
with open("simulation_output.txt", "r") as f:
    lines = f.readlines()

# Each frame has num_particles lines
frames = len(lines) // num_particles

# Reshape data into (frames, num_particles, 2)
data = np.array([[float(value) for value in line.split()] for line in lines])
data = data.reshape((frames, num_particles, 2))

# Plot setup
fig, ax = plt.subplots()
scat = ax.scatter([], [], s=500, color='blue', alpha=0.09, edgecolors='none')

ax.set_xlim(-2.5, 12.5)
ax.set_ylim(-2.5, 12.5)
ax.set_title("2D Fluid Simulation")
ax.set_xlabel("X")
ax.set_ylabel("Y")

# Update function for animation
def update(frame):
    positions = data[frame]
    scat.set_offsets(positions)
    ax.set_title(f"2D Fluid Simulation")

    # Save snapshot at specified frames
    #if frame in save_frames:
        #filename = f"fluid_frame_{frame * log_interval}.png"
        #plt.savefig(filename, dpi=300, bbox_inches='tight')
        #print(f"Saved snapshot: {filename}")

    return scat,

# Create animation
ani = animation.FuncAnimation(
    fig, update, frames=frames, interval=30, blit=True
)
#ani.save("cloud_collision_blur.gif", writer='pillow', fps=40, dpi=180)
plt.tight_layout()
plt.show()
