import numpy as np
import matplotlib.pyplot as plt


# Constants
num_nodes_x = 60  # Number of nodes in the x direction
num_nodes_y = 40  # Number of nodes in the y direction
dx = 0.001  # Spatial step in the x direction
dy = 0.001  # Spatial step in the y direction
q_top = 0.2  # Heat flux at the top boundary
q_bottom = 0.2  # Heat flux at the bottom boundary
dt = 0.000000001  # Time step for FDM update
num_timesteps = 100  # Total number of timesteps

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

    def update_properties(self):
        # Update properties based on the current temperature
        self.density_superfluid = self.calculate_density_superfluid()
        print(self.density_superfluid)
        self.density_normal = self.calculate_density_normal()
        self.entropy = self.calculate_entropy()
        self.viscosity = self.calculate_viscosity()

    def calculate_density_superfluid(self):
        # Example: Simple relationship; modify based on your system
        return (self.temperature)  # Adjust as needed

    def calculate_density_normal(self):
        # Example: Simple relationship; modify based on your system
        return (self.temperature)  # Adjust as needed

    def calculate_entropy(self):
        # Example: Simple relationship; modify based on your system
        return (self.temperature)  # Adjust as needed

    def calculate_viscosity(self):
        # Example: Simple relationship; modify based on your system
        return (self.temperature)  # Adjust as needed

    def calculate_k(self):
        L = num_nodes_y * (dy**2)  # Length of the grid in meters
        total_density = self.density_superfluid + self.density_normal
        k = (((L * total_density * self.entropy)**2) / (8 * self.viscosity))
        return k

    def calculate_alpha(self):
        c_p = 5000  # Specific heat constant
        total_density = self.density_superfluid + self.density_normal
        alpha = self.calculate_k() / (total_density * c_p)
        return alpha

# Initialize the grid with uniform temperature
def initialize_grid():
    nodes = []
    for y in range(num_nodes_y):
        row = []
        for x in range(num_nodes_x):
            node = Node(x, y, temperature=1.9)
            row.append(node)
        nodes.append(row)
    return nodes

# Apply Neumann boundary conditions to update boundary nodes only
def apply_neumann_bc_boundary_only(nodes, q_top, q_bottom, dy):
    # Top boundary
    for x in range(num_nodes_x):
        T_interior = nodes[1][x].temperature  # First interior node below the boundary
        k_local = nodes[0][x].calculate_k()  # Use local k value for boundary node
        nodes[0][x].temperature = T_interior + (q_top * dy) / k_local  # Apply Neumann BC for top boundary

    # Bottom boundary
    for x in range(num_nodes_x):
        T_interior = nodes[num_nodes_y - 2][x].temperature  # First interior node above the boundary
        k_local = nodes[num_nodes_y - 1][x].calculate_k()  # Use local k value for boundary node
        nodes[num_nodes_y - 1][x].temperature = T_interior + (q_bottom * dy) / k_local  # Apply Neumann BC for bottom boundary

# Update interior nodes based on the FDM method
def update_interior_nodes(nodes, dt):
    new_temperatures = np.array([[node.temperature for node in row] for row in nodes])

    for y in range(1, num_nodes_y - 1):
        for x in range(1, num_nodes_x - 1):
            T_xx = (nodes[y][x + 1].temperature - 2 * nodes[y][x].temperature + nodes[y][x - 1].temperature) / (dx ** 2)
            T_yy = (nodes[y + 1][x].temperature - 2 * nodes[y][x].temperature + nodes[y - 1][x].temperature) / (dy ** 2)
            alpha = nodes[y][x].calculate_alpha()
            new_temperatures[y][x] += alpha * (T_xx + T_yy) * dt  # Update using alpha and FDM

    for y in range(1, num_nodes_y - 1):
        for x in range(1, num_nodes_x - 1):
            nodes[y][x].temperature = new_temperatures[y][x]

# Main simulation function
def run_simulation():
    nodes = initialize_grid()  # Initialize the grid

    # Iterate over the specified number of timesteps
    for t in range(num_timesteps):
        apply_neumann_bc_boundary_only(nodes, q_top, q_bottom, dy)  # Update boundary nodes
        update_interior_nodes(nodes, dt)  # Update interior nodes based on FDM

        # Update properties for each node based on the new temperature after updating temperatures
        for row in nodes:
            for node in row:
                node.update_properties()

    # Visualize the final temperature distribution
    visualize_results(nodes)

# Visualization function
def visualize_results(nodes):
    temperatures = np.array([[node.temperature for node in row] for row in nodes])
    plt.imshow(temperatures, cmap='hot', interpolation='nearest')
    plt.colorbar(label='Temperature (K)')
    plt.title(f'Temperature Distribution After {num_timesteps} Timesteps')
    plt.xlabel('X Nodes')
    plt.ylabel('Y Nodes')
    plt.show()

# Run the simulation
run_simulation()
