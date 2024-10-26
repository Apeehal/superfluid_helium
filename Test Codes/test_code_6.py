
import numpy as np
import matplotlib.pyplot as plt

# Constants
num_nodes_x = 60  # Number of nodes in the x direction
num_nodes_y = 40  # Number of nodes in the y direction
dx = 0.001  # Spatial step in the x direction
dy = 0.001  # Spatial step in the y direction
q_top = 20  # Heat flux at the top boundary
q_bottom = 20  # Heat flux at the bottom boundary
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
        L = num_nodes_y * 10**-3  # Length of the grid in meters
        total_density = self.density_superfluid + self.density_normal
        k = (((L * total_density * self.entropy)**2) / (8 * self.viscosity))
        return k

    def calculate_alpha(self):
        c_p = 5000  # Specific heat constant
        total_density = self.density_superfluid + self.density_normal
        alpha = self.k / (total_density * c_p)
        return alpha

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

# Apply Neumann boundary conditions, only updating boundary nodes
def apply_neumann_bc_boundary_only(nodes, q_top, q_bottom, dy):
    # Top boundary
    for x in range(num_nodes_x):
        T_interior = nodes[1][x].temperature  # First interior node below the boundary
        k_local = nodes[0][x].k  # Use local k value for boundary node
        nodes[0][x].temperature = T_interior + (q_top * dy) / k_local  # Apply Neumann BC for top boundary
        print(nodes[0][0].temperature)

    # Bottom boundary
    for x in range(num_nodes_x):
        T_interior = nodes[num_nodes_y - 2][x].temperature  # First interior node above the boundary
        k_local = nodes[num_nodes_y - 1][x].k  # Use local k value for boundary node
        nodes[num_nodes_y - 1][x].temperature = T_interior + (q_bottom * dy) / k_local  # Apply Neumann BC for bottom boundary

# Main simulation function
def run_simulation():
    nodes = initialize_grid()  # Initialize the grid

    # Time-stepping loop, updating only boundary nodes
    for t in range(num_timesteps):
        apply_neumann_bc_boundary_only(nodes, q_top, q_bottom, dy)  # Only update boundary nodes

    # Visualize the final temperature distribution
    visualize_results(nodes)

# Visualization function
def visualize_results(nodes):
    temperatures = np.array([[node.temperature for node in row] for row in nodes])
    plt.imshow(temperatures, cmap='hot', interpolation='nearest')
    plt.colorbar(label='Temperature (K)')
    plt.title('Temperature Distribution After Simulation (Boundary-Only Update)')
    plt.xlabel('X Nodes')
    plt.ylabel('Y Nodes')
    plt.show()

# Run the simulation
run_simulation()



