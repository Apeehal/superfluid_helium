import pandas as pd
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# Define the grid dimensions
num_nodes_x = 60  # Number of nodes in the x-direction
num_nodes_y = 40  # Number of nodes in the y-direction

# Generate the grid of nodes
grid = [(x, y) for y in range(num_nodes_y) for x in range(num_nodes_x)]

lambda_transition = 2.17

class Node:
    def __init__(self, x, y, heat_flux=(0.0,0.0), temperature=0.0, temperature_gradient=(0.0,0.0), pressure=0.0, pressure_gradient=(0.0,0.0), velocity_normal=(0.0, 0.0), velocity_superfluid=(0.0,0.0), density_normal=0.0, density_superfluid=0.0, entropy=0.0):
    
        self.x = x
        self.y = y
        self.heat_flux = heat_flux
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
        return f"Node({self.x}, {self.y}, q={self.heat_flux}, T={self.temperature}, grad_T={self.temperature_gradient}, P={self.pressure}, grad_P={self.pressure_gradient}, v_n={self.velocity_normal}, v_s={self.velocity_superfluid}, rho_n={self.density_normal}, rho_s={self.density_superfluid}, s={self.entropy})"



# Create a grid of nodes
nodes = [[Node(x, y) for x in range(num_nodes_x)] for y in range(num_nodes_y)]


############################
# Assign Initial Conditions
############################

initial_temperature = 100
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



initial_heat_flux = 0
for row in nodes:
    for node in row:
        node.heat_flux = initial_heat_flux


initial_normal_velocity = (0,0)
for row in nodes:
    for node in row:
        node.velocity_normal = initial_normal_density


initial_superfluid_velocity = (0,0)
for row in nodes:
    for node in row:
        node.velocity_superfluid = initial_superfluid_density


############################
# Assign boundary conditions
############################
for y in range(num_nodes_y):
    for x in range(num_nodes_x):
        if y == 0:  # Bottom boundary
            nodes[y][x].heat_flux = 0.2  
            nodes[y][x].temperature = 300

        elif y == num_nodes_y - 1:  # Top boundary
            nodes[y][x].heat_flux = 0.2  
            nodes[y][x].temperature = 300

           

        if x == 0:  # Left boundary
            nodes[y][x].pressure = 101325
        elif x == num_nodes_x - 1:  # Right boundary
            nodes[y][x].pressure = 0.101325






############################
############################
# Start Computation ########
############################
############################


############################
# Solve for Temperatures ###
############################

if x != 0 and x != num_nodes_x-1 and y != range(1, num_nodes_y):
    print(y)



def temperature_averaging_step(nodes):

    new_nodes = [[Node(x, y, node.temperature) for x, node in enumerate(row)] for y, row in enumerate(nodes)]

    for y in range(num_nodes_y):
        for x in range(num_nodes_x):

            if x != 0 and x != num_nodes_x-1 and y != 0 and y != num_nodes_y-1:  
                T_avg = (nodes[y][x+1].temperature + nodes[y][x-1].temperature + nodes[y+1][x].temperature + nodes[y-1][x].temperature)/4
                #k = 
                #Q_component = 
                #T = 
            
                new_nodes[y][x].temperature = T_avg
            
            if x == num_nodes_x - 1 and (y != 0 and y != num_nodes_y - 1):
                T_avg = (nodes[y][x-1].temperature + nodes[y+1][x].temperature + nodes[y-1][x].temperature) / 3
                new_nodes[y][x].temperature = T_avg

            if x == 0 and (y != 0 and y != num_nodes_y - 1):
                T_avg = (nodes[y][x+1].temperature + nodes[y+1][x].temperature + nodes[y-1][x].temperature) / 3
                new_nodes[y][x].temperature = T_avg
    
    return new_nodes



# Apply the temperature averaging step for multiple iterations
iterations = 2 # Number of times we want to average temperatures
for i in range(iterations):
    nodes = temperature_averaging_step(nodes)

# Extract temperature values for plotting
temperature_data = np.array([[node.temperature for node in row] for row in nodes])
# Plot the temperature distribution
plt.figure(figsize=(8, 6))
plt.imshow(temperature_data, cmap='hot', interpolation='nearest', origin='lower')

# Add color bar and labels
plt.colorbar(label='Temperature (K)')
plt.title(f'Temperature Distribution after {iterations} Averaging Steps')
plt.xlabel('X-axis (nodes)')
plt.ylabel('Y-axis (nodes)')

# Show the plot
plt.show()






"""

# Function to apply FDM and calculate next temperature values
def fdm_step(nodes, alpha, dt, dx, dy):
    # Create a copy of current temperatures to store the new values
    new_nodes = [[Node(x, y, node.temperature) for x, node in enumerate(row)] for y, row in enumerate(nodes)]
    
    for y in range(1, num_nodes_y - 1):  # Avoid boundary nodes
        for x in range(1, num_nodes_x - 1):  # Avoid boundary nodes
            T_xx = (nodes[y][x+1].temperature - 2 * nodes[y][x].temperature + nodes[y][x-1].temperature) / dx**2
            T_yy = (nodes[y+1][x].temperature - 2 * nodes[y][x].temperature + nodes[y-1][x].temperature) / dy**2
            new_nodes[y][x].temperature = nodes[y][x].temperature + alpha * dt * (T_xx + T_yy)
    
    # Return the updated grid
    return new_nodes

# Simulation loop for a number of time steps
for step in range(time_steps):
    nodes = fdm_step(nodes, alpha, dt, dx, dy)


# Extract temperature values for plotting
temperature_data = np.array([[node.temperature for node in row] for row in nodes])

# Plot the temperature
plt.figure(figsize=(8, 6))
plt.imshow(temperature_data, cmap='hot', interpolation='nearest', origin='lower')

# Add color bar
plt.colorbar(label='Temperature (K)')

# Add labels and title
plt.title(f'Temperature Distribution after {time_steps} Time Steps')
plt.xlabel('X-axis (nodes)')
plt.ylabel('Y-axis (nodes)')

# Set ticks for better readability
plt.xticks(ticks=np.arange(num_nodes_x), labels=np.arange(num_nodes_x))
plt.yticks(ticks=np.arange(num_nodes_y), labels=np.arange(num_nodes_y))

# Show the plot
plt.show()




















# Extract temperature values into a 2D array
pressure_data = np.array([[node.pressure for node in row] for row in nodes])

# Plot the temperature
plt.figure(figsize=(8, 6))
plt.imshow(pressure_data, cmap='hot', interpolation='nearest', origin='lower')

# Add color bar
plt.colorbar(label='Pressure (Pa)')

# Add labels and title
plt.title('Temperature Distribution with Global Initial Condition')
plt.xlabel('X-axis (nodes)')
plt.ylabel('Y-axis (nodes)')

# Set ticks for better readability
plt.xticks(ticks=np.arange(num_nodes_x), labels=np.arange(num_nodes_x))
plt.yticks(ticks=np.arange(num_nodes_y), labels=np.arange(num_nodes_y))

# Show the plot
plt.show()
"""
