import numpy as np
import matplotlib.pyplot as plt

# Define the Node class with properties related to superfluid and normal densities, entropy, etc.
class Node:
    def __init__(self, x, y, temperature=1.9, pressure=0.0, velocity=(0.0, 0.0), 
                 density_superfluid=162.9, density_normal=147.5, entropy=756.95, viscosity=2.5*10**-6):
        self.x = x
        self.y = y
        self.temperature = temperature
        self.pressure = pressure
        self.velocity = velocity
        self.density_superfluid = density_superfluid
        self.density_normal = density_normal
        self.entropy = entropy
        self.viscosity = viscosity
        self.k = self.calculate_k()  # Thermal conductivity based on given formula
        self.alpha = self.calculate_alpha()  # Thermal diffusivity based on given formula
    
    # Calculate thermal conductivity k
    def calculate_k(self):
        L = num_nodes_y * dy
        #total_density = self.density_superfluid + self.density_normal
        #print(self.density_superfluid)
        total_density = 162.9 + 147.5
        #k = (((L * total_density * self.entropy)**2) / (8 * self.viscosity))
        k = (((L * total_density * self.entropy)**2) / (8 * self.viscosity))
        return k

    # Calculate thermal diffusivity alpha
    def calculate_alpha(self):
        c_p = 5000  # Specific heat constant
        #total_density = self.density_superfluid + self.density_normal
        total_density = 162.9 + 147.5
        alpha = self.k / (total_density * c_p)
        return alpha

# Define grid dimensions
num_nodes_x = 60
num_nodes_y = 40
dx = dy = 0.01  # Spatial step (assuming uniform grid)
dt = 0.1       # Time step size
time_steps = 200  # Number of time steps
q_top = 2  # Heat flux at the top boundary
q_bottom = 2  # Heat flux at the bottom boundary

# Function to initialize densities, entropy, and viscosity for each node
def initialize_node_properties(x, y):
    # Example: Define varying densities, entropy, and viscosity
    # You can adjust this function to simulate real-world scenarios
    temperature = 1.9
    density_superfluid = 162.9  # Example variation with position
    density_normal = 147.5  # Example variation with position
    entropy = 756.95  # Example entropy variation
    viscosity =  2.5*10**-6  # Example viscosity variation
    return temperature, density_superfluid, density_normal, entropy, viscosity

# Create a grid of nodes with initialized density, entropy, and viscosity
nodes = [[Node(x, y, *initialize_node_properties(x, y)) for x in range(num_nodes_x)] for y in range(num_nodes_y)]

# Function to apply FDM and calculate next temperature values with varying k, alpha, and Neumann BCs
def fdm_step(nodes, dt, dx, dy, q_top, q_bottom):
    # Create a copy of current temperatures to store the new values
    new_nodes = [[Node(x, y, node.temperature, node.pressure, node.velocity, node.density_superfluid, 
                       node.density_normal, node.entropy, node.viscosity) for x, node in enumerate(row)] for y, row in enumerate(nodes)]
    

        # Neumann Boundary Conditions with heat flux at top and bottom
    # Bottom boundary (specified heat flux q_bottom)
    for x in range(num_nodes_x):
        T_interior = nodes[1][x].temperature  # First interior node above the boundary
        k_local = nodes[0][x].k  # Use local k value for boundary node
        new_nodes[0][x].temperature = T_interior + (q_bottom * dy) / k_local  # Apply Neumann BC for bottom boundary

    # Top boundary (specified heat flux q_top)
    for x in range(num_nodes_x):
        T_interior = nodes[num_nodes_y - 2][x].temperature  # Last interior node below the boundary
        k_local = nodes[num_nodes_y - 1][x].k  # Use local k value for boundary node
        new_nodes[num_nodes_y - 1][x].temperature = T_interior + (q_top * dy) / k_local  # Apply Neumann BC for top boundary

    # Left and right boundaries (assuming no flux for simplicity, can be adjusted if needed)
    for y in range(num_nodes_y):
        new_nodes[y][0].temperature = new_nodes[y][1].temperature  # Left boundary (no flux)
        new_nodes[y][num_nodes_x - 1].temperature = new_nodes[y][num_nodes_x - 2].temperature  # Right boundary (no flux)
    
    # Update interior nodes using FDM with varying alpha
    for y in range(1, num_nodes_y - 1):
        for x in range(1, num_nodes_x - 1):
            T_xx = (nodes[y][x+1].temperature - 2 * nodes[y][x].temperature + nodes[y][x-1].temperature) / dx**2
            T_yy = (nodes[y+1][x].temperature - 2 * nodes[y][x].temperature + nodes[y-1][x].temperature) / dy**2
            alpha_local = nodes[y][x].alpha  # Use local alpha value
            new_nodes[y][x].temperature = nodes[y][x].temperature + alpha_local * dt * (T_xx + T_yy)

   


    #for y in range(num_nodes_y):
    #    for x in range(num_nodes_x):

    #        new_nodes[y][x].density_normal = 147.5
    #        new_nodes[y][x].density_superfluid = 162.9
    #        new_nodes[y][x].entropy = 756.95


    # Return the updated grid
    return new_nodes

# Simulation loop for a number of time steps
for step in range(time_steps):
    nodes = fdm_step(nodes, dt, dx, dy, q_top, q_bottom)

# Extract temperature values for plotting
temperature_data = np.array([[node.temperature for node in row] for row in nodes])

# Plot the temperature
plt.figure(figsize=(8, 6))
plt.imshow(temperature_data, cmap='hot', interpolation='nearest', origin='lower')

# Add color bar
plt.colorbar(label='Temperature (K)')

# Add labels and title
plt.title(f'Temperature Distribution with Varying k and alpha (Density/Entropy-Based) after {time_steps} Time Steps')
plt.xlabel('X-axis (nodes)')
plt.ylabel('Y-axis (nodes)')

# Set ticks for better readability
plt.xticks(ticks=np.arange(num_nodes_x), labels=np.arange(num_nodes_x))
plt.yticks(ticks=np.arange(num_nodes_y), labels=np.arange(num_nodes_y))

# Show the plot
plt.show()
