import numpy as np
import matplotlib.pyplot as plt

# Constants
num_nodes_x = 60  # Number of nodes in the x direction
num_nodes_y = 40  # Number of nodes in the y direction
dx = dy = 0.001  # Spatial step in the x direction
q_top = 2  # Heat flux at the top boundary
q_bottom = 2  # Heat flux at the bottom boundary
dt = 0.01  # Time step
num_timesteps = 1 # Total number of time steps

# Node class
class Node:
    def __init__(self, x, y, temperature=1.9, density_superfluid=162.9, density_normal=147.5, entropy=756.95, viscosity=2.5*10**-6):
        self.x = x
        self.y = y
        self.temperature = temperature
        self.density_superfluid = density_superfluid
        self.density_normal = density_normal
        self.entropy = entropy
        self.viscosity = viscosity
        self.k = self.calculate_k()
        self.alpha = self.calculate_alpha()

    def calculate_k(self):
        L = num_nodes_y * dy # Length of the grid in meters
        total_density = self.density_superfluid + self.density_normal
        k = (((L * total_density * self.entropy)**2) / (8 * self.viscosity))
        return k

    def calculate_alpha(self):
        c_p = 5000  # Specific heat constant
        total_density = self.density_superfluid + self.density_normal
        alpha = self.k / (total_density * c_p)
        return alpha




# Main simulation function
def run_simulation():
    # Step 1: Initialize the grid
    nodes = [[Node(x, y) for x in range(num_nodes_x)] for y in range(num_nodes_y)]

    # Step 2: Simulation loop
    for t in range(num_timesteps):
        # Apply Neumann boundary conditions
        for x in range(num_nodes_x):
            # Top boundary (specified heat flux)
            T_interior_top = nodes[1][x].temperature  # First interior node below the boundary
            k_local_top = nodes[0][x].k  # Use local k value for boundary node
            nodes[0][x].temperature = T_interior_top + (q_top * dy) / k_local_top  # Update top boundary temperature

            # Bottom boundary (specified heat flux)
            T_interior_bottom = nodes[num_nodes_y - 2][x].temperature  # First interior node above the boundary
            k_local_bottom = nodes[num_nodes_y - 1][x].k  # Use local k value for boundary node
            nodes[num_nodes_y - 1][x].temperature = T_interior_bottom + (q_bottom * dy) / k_local_bottom  # Update bottom boundary temperature

        # Step 3: Compute temperature gradients
        T_xx = np.zeros((num_nodes_y, num_nodes_x))
        T_yy = np.zeros((num_nodes_y, num_nodes_x))

        for y in range(1, num_nodes_y - 1):
            for x in range(1, num_nodes_x - 1):
                T_xx[y][x] = (nodes[y][x + 1].temperature - 2 * nodes[y][x].temperature + nodes[y][x - 1].temperature) / (dx ** 2)
                T_yy[y][x] = (nodes[y + 1][x].temperature - 2 * nodes[y][x].temperature + nodes[y - 1][x].temperature) / (dy ** 2)

        # Step 4: Update interior node temperatures based on gradients
        for y in range(1, num_nodes_y - 1):
            for x in range(1, num_nodes_x - 1):
                # Update temperature using thermal diffusivity alpha
                nodes[y][x].temperature += nodes[y][x].alpha * (T_xx[y][x] + T_yy[y][x]) * dt

    # Step 5: Visualize the final temperature distribution
    visualize_results(nodes)

# Visualization function
def visualize_results(nodes):
    temperatures = np.array([[node.temperature for node in row] for row in nodes])
    plt.imshow(temperatures, cmap='hot', interpolation='nearest')
    plt.colorbar(label='Temperature (K)')
    plt.title('Temperature Distribution After Simulation')
    plt.xlabel('X Nodes')
    plt.ylabel('Y Nodes')
    plt.show()

# Run the simulation
run_simulation()

"""
# Initialize the grid with uniform temperature
def initialize_grid():
    nodes = []
    for y in range(num_nodes_y):
        row = []
        for x in range(num_nodes_x):
            # All nodes start at 1.9 K
            node = Node(x, y, temperature=1.9)
            row.append(node)
        nodes.append(row)
    return nodes

# Apply Neumann boundary conditions
def apply_neumann_bc(nodes, q_top, q_bottom, dy):
    # Top boundary
    for x in range(num_nodes_x):
        T_interior = nodes[1][x].temperature  # First interior node below the boundary
        k_local = nodes[0][x].k  # Use local k value for boundary node
        nodes[0][x].temperature = T_interior + (q_top * dy) / k_local  # Apply Neumann BC for top boundary

    # Bottom boundary
    for x in range(num_nodes_x):
        T_interior = nodes[num_nodes_y - 2][x].temperature  # First interior node above the boundary
        k_local = nodes[num_nodes_y - 1][x].k  # Use local k value for boundary node
        nodes[num_nodes_y - 1][x].temperature = T_interior + (q_bottom * dy) / k_local  # Apply Neumann BC for bottom boundary

    print("Top boundary temperatures:", [nodes[0][x].temperature for x in range(num_nodes_x)])
    print("Bottom boundary temperatures:", [nodes[num_nodes_y - 1][x].temperature for x in range(num_nodes_x)])

# Compute temperature gradients
def compute_temperature_gradients(nodes):
    T_xx = np.zeros((num_nodes_y, num_nodes_x))
    T_yy = np.zeros((num_nodes_y, num_nodes_x))

    for y in range(1, num_nodes_y - 1):
        for x in range(1, num_nodes_x - 1):
            T_xx[y][x] = (nodes[y][x + 1].temperature - 2 * nodes[y][x].temperature + nodes[y][x - 1].temperature) / (dx ** 2)
            T_yy[y][x] = (nodes[y + 1][x].temperature - 2 * nodes[y][x].temperature + nodes[y - 1][x].temperature) / (dy ** 2)
    return T_xx, T_yy

# Main simulation function
def run_simulation():
    nodes = initialize_grid()  # Initialize the grid

    for t in range(num_timesteps):
        apply_neumann_bc(nodes, q_top, q_bottom, dy)  # Update boundary conditions
        T_xx, T_yy = compute_temperature_gradients(nodes)  # Calculate gradients

        # Update temperatures based on the calculated gradients
        for y in range(1, num_nodes_y - 1):
            for x in range(1, num_nodes_x - 1):
                # Update temperature using thermal diffusivity alpha
                nodes[y][x].temperature += nodes[y][x].alpha * (T_xx[y][x] + T_yy[y][x]) * dt

    # Visualize the final temperature distribution
    visualize_results(nodes)

# Visualization function
def visualize_results(nodes):
    temperatures = np.array([[node.temperature for node in row] for row in nodes])
    plt.imshow(temperatures, cmap='hot', interpolation='nearest')
    plt.colorbar(label='Temperature (K)')
    plt.title('Temperature Distribution After Simulation')
    plt.xlabel('X Nodes')
    plt.ylabel('Y Nodes')
    plt.show()

# Run the simulation
run_simulation()
"""