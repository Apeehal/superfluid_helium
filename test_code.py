import numpy as np
import matplotlib.pyplot as plt



# Define grid dimensions and constants
num_nodes_x = 60
num_nodes_y = 40

grid = [(x, y) for y in range(num_nodes_y) for x in range(num_nodes_x)]

dx = dy = 0.01  # Spatial step (assuming uniform grid)
q_top =  0.2  # Heat flux at the left boundary
q_bottom = 0.2 # Heat flux at the right boundary
dt = 1       # Time step size
time_steps = 20 # Number of time steps

viscosity = 2.5*10**-6
lambda_transition = 2.17




# Define the Node class
class Node:
    def __init__(self, x, y, temperature=0.0, temperature_gradient=(0.0,0.0), pressure=0.0, pressure_gradient=(0.0,0.0), velocity_normal=(0.0, 0.0), velocity_superfluid=(0.0,0.0), density_normal=0.0, density_superfluid=0.0, entropy=0.0):
    
        self.x = x
        self.y = y
        self.temperature = temperature
        self.temperature_gradient = temperature_gradient
        self.pressure = pressure
        self.pressure_gradient = pressure_gradient
        self.velocity_normal = velocity_normal
        self.velocity_superfluid = velocity_superfluid
        self.density_normal = density_normal
        self.density_superfluid = density_superfluid
        self.entropy = entropy


    def __repr__(self):
        return f"Node({self.x}, {self.y}, T={self.temperature}, grad_T={self.temperature_gradient}, P={self.pressure}, grad_P={self.pressure_gradient}, v_n={self.velocity_normal}, v_s={self.velocity_superfluid}, rho_n={self.density_normal}, rho_s={self.density_superfluid}, s={self.entropy})"



def initialize_node(x, y):
    # Example initialization: You can modify the logic to reflect realistic conditions
    temperature = 1.9
    density_normal = 147.5
    density_superfluid = (density_normal/((temperature/lambda_transition)**5.6)) - density_normal
    entropy = 756.95  # Entropy can be constant or vary
    return Node(x, y, temperature=temperature,density_normal=density_normal, density_superfluid=density_superfluid, entropy=entropy)



# Create a grid of nodes with initial temperature condition (global 1.9K)
nodes = [[initialize_node(x, y) for x in range(num_nodes_x)] for y in range(num_nodes_y)]

############################
# Assign Initial Conditions
############################
"""
initial_temperature = 1.9
for row in nodes:
    for node in row:
        node.temperature = initial_temperature


initial_pressure = 0.101325
for row in nodes:
    for node in row:
        node.pressure = initial_pressure


initial_entropy = 756.95
for row in nodes:
    for node in row:
        node.entropy = initial_entropy


initial_normal_density = 147.5
for row in nodes:
    for node in row:
        node.density_normal = initial_normal_density


initial_superfluid_density = (initial_normal_density/((initial_temperature/lambda_transition)**5.6)) - initial_normal_density
for row in nodes:
    for node in row:
        node.density_superfluid = initial_superfluid_density



initial_normal_velocity = (0,0)
for row in nodes:
    for node in row:
        node.velocity_normal = initial_normal_density


initial_superfluid_velocity = (0,0)
for row in nodes:
    for node in row:
        node.velocity_superfluid = initial_superfluid_density
"""

############################
# Assign Pressure Boundary conditions
############################
for y in range(num_nodes_y):
    for x in range(num_nodes_x):
           
        if x == 0:  # Left boundary
            nodes[y][x].pressure = 101325





c_p = 5000
def fdm_step(nodes, q_top, q_bottom, dt, dx, dy):
    # Create a copy of current temperatures to store the new values
    new_nodes = [[Node(x, y, node.temperature, node.temperature_gradient, node.pressure, node.pressure_gradient, node.velocity_normal, node.velocity_superfluid, node.density_normal, node.density_superfluid, node.entropy) for x, node in enumerate(row)] for y, row in enumerate(nodes)]

    # Update interior nodes using FDM
    for y in range(1, num_nodes_y - 1):
        for x in range(1, num_nodes_x - 1):
            T_xx = (nodes[y][x+1].temperature - 2 * nodes[y][x].temperature + nodes[y][x-1].temperature) / dx**2
            T_yy = (nodes[y+1][x].temperature - 2 * nodes[y][x].temperature + nodes[y-1][x].temperature) / dy**2
            k = (((num_nodes_y*10**-3)*(nodes[y][x].density_superfluid + nodes[y][x].density_normal)*(nodes[y][x].entropy))**2)/(8*viscosity)
            alpha = k/((nodes[y][x].density_superfluid+nodes[y][x].density_normal)*(c_p))
            new_nodes[y][x].temperature = nodes[y][x].temperature + alpha * dt * (T_xx + T_yy)
    
    # Neumann Boundary Conditions with heat flux at top and bottom
    # Bottom boundary (specified heat flux q_bottom)
    for x in range(num_nodes_x):
        T_interior = nodes[1][x].temperature  # First interior node above the boundary
        k = (((num_nodes_y*10**-3)*(nodes[0][x].density_superfluid + nodes[0][x].density_normal)*(nodes[0][x].entropy))**2)/(8*viscosity)
        new_nodes[0][x].temperature = T_interior + (q_bottom * dy) / k  # Apply Neumann BC for bottom boundary

    # Top boundary (specified heat flux q_top)
    for x in range(num_nodes_x):
        T_interior = nodes[num_nodes_y - 2][x].temperature  # Last interior node below the boundary
        k = (((num_nodes_y*10**-3)*(nodes[num_nodes_y - 1][x].density_superfluid + nodes[num_nodes_y - 1][x].density_normal)*(nodes[num_nodes_y - 1][x].entropy))**2)/(8*viscosity)
        new_nodes[num_nodes_y - 1][x].temperature = T_interior + (q_top * dy) / k  # Apply Neumann BC for top boundary
    
    # Left and right boundaries (assuming no flux for simplicity, can be adjusted if needed)
    for y in range(num_nodes_y):
        new_nodes[y][0].temperature = new_nodes[y][1].temperature  # Left boundary (no flux)
        new_nodes[y][num_nodes_x - 1].temperature = new_nodes[y][num_nodes_x - 2].temperature  # Right boundary (no flux)
    

    

    # Return the updated grid
    return new_nodes

# Simulation loop for a number of time steps
for step in range(time_steps):
    nodes = fdm_step(nodes, q_top, q_bottom, dt, dx, dy)

# Extract temperature values for plotting
temperature_data = np.array([[node.temperature for node in row] for row in nodes])

# Plot the temperature
plt.figure(figsize=(8, 6))
plt.imshow(temperature_data, cmap='hot', interpolation='nearest', origin='lower')

# Add color bar
plt.colorbar(label='Temperature (K)')

# Add labels and title
plt.title(f'Temperature Distribution with Neumann Boundary Condition (Top & Bottom Heat Flux) after {time_steps} Time Steps')
plt.xlabel('X-axis (nodes)')
plt.ylabel('Y-axis (nodes)')

# Set ticks for better readability
plt.xticks(ticks=np.arange(num_nodes_x), labels=np.arange(num_nodes_x))
plt.yticks(ticks=np.arange(num_nodes_y), labels=np.arange(num_nodes_y))

# Show the plot
plt.show()