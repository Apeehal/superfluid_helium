import numpy as np
import matplotlib.pyplot as plt

# Define the Node class
class Node:
    def __init__(self, x, y, temperature=1.9, pressure=0.0, velocity=(0.0, 0.0)):
        self.x = x
        self.y = y
        self.temperature = temperature
        self.pressure = pressure
        self.velocity = velocity

# Define grid dimensions and constants
num_nodes_x = 60
num_nodes_y = 40
dx = dy = 1.0  # Spatial step (assuming uniform grid)
alpha = 0.01   # Thermal diffusivity
k = 0.5        # Thermal conductivity (for the Neumann condition)
q_left = -10.0  # Heat flux at the left boundary
q_right = 0 # Heat flux at the right boundary
dt = 0.1       # Time step size
time_steps = 1000  # Number of time steps

# Create a grid of nodes with initial temperature condition (global 1.9K)
nodes = [[Node(x, y) for x in range(num_nodes_x)] for y in range(num_nodes_y)]

# Function to apply FDM and calculate next temperature values with Neumann BCs (specified heat flux)
def fdm_step(nodes, alpha, k, q_left, q_right, dt, dx, dy):
    # Create a copy of current temperatures to store the new values
    new_nodes = [[Node(x, y, node.temperature) for x, node in enumerate(row)] for y, row in enumerate(nodes)]
    
    # Update interior nodes using FDM
    for y in range(1, num_nodes_y - 1):
        for x in range(1, num_nodes_x - 1):
            T_xx = (nodes[y][x+1].temperature - 2 * nodes[y][x].temperature + nodes[y][x-1].temperature) / dx**2
            T_yy = (nodes[y+1][x].temperature - 2 * nodes[y][x].temperature + nodes[y-1][x].temperature) / dy**2
            new_nodes[y][x].temperature = nodes[y][x].temperature + alpha * dt * (T_xx + T_yy)
    
    # Neumann Boundary Conditions with heat flux
    # Left boundary (specified heat flux q_left)
    for y in range(num_nodes_y):
        T1 = nodes[y][1].temperature  # First interior node
        new_nodes[y][0].temperature = T1 - (q_left * dx) / k  # Apply Neumann BC for left boundary

    # Right boundary (specified heat flux q_right)
    for y in range(num_nodes_y):
        T_last = nodes[y][num_nodes_x - 2].temperature  # Last interior node
        new_nodes[y][num_nodes_x - 1].temperature = T_last + (q_right * dx) / k  # Apply Neumann BC for right boundary
    
    # Top and bottom boundaries (we assume no flux for simplicity, can apply similar logic if needed)
    for x in range(num_nodes_x):
        new_nodes[0][x].temperature = new_nodes[1][x].temperature  # Bottom boundary (no flux)
        new_nodes[num_nodes_y - 1][x].temperature = new_nodes[num_nodes_y - 2][x].temperature  # Top boundary (no flux)
    
    # Return the updated grid
    return new_nodes

# Simulation loop for a number of time steps
for step in range(time_steps):
    nodes = fdm_step(nodes, alpha, k, q_left, q_right, dt, dx, dy)

# Extract temperature values for plotting
temperature_data = np.array([[node.temperature for node in row] for row in nodes])

# Plot the temperature
plt.figure(figsize=(8, 6))
plt.imshow(temperature_data, cmap='hot', interpolation='nearest', origin='lower')

# Add color bar
plt.colorbar(label='Temperature (K)')

# Add labels and title
plt.title(f'Temperature Distribution with Neumann Boundary Condition (Heat Flux) after {time_steps} Time Steps')
plt.xlabel('X-axis (nodes)')
plt.ylabel('Y-axis (nodes)')

# Set ticks for better readability
plt.xticks(ticks=np.arange(num_nodes_x), labels=np.arange(num_nodes_x))
plt.yticks(ticks=np.arange(num_nodes_y), labels=np.arange(num_nodes_y))

# Show the plot
plt.show()
