





import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

# Constants
num_nodes_x = 10  # Number of nodes in the x-direction
num_nodes_y = 10  # Number of nodes in the y-direction
dx = 0.01         # Space step in meters
dy = 0.01         # Space step in meters
time_step = 0.01  # Time step in seconds
num_timesteps = 100  # Number of time steps to simulate

# Node class as provided
class Node:
    def __init__(self, x, y, temperature=1.9, density_superfluid=162.9, density_normal=147.5, 
                 csv_file='helium_Property_Table.csv', entropy=756.95, viscosity=2.5e-6, 
                 pressure=0.876, v_s=(0, 0), v_n=(0, 0)):
        self.x = x
        self.y = y
        self.temperature = temperature
        self.density_superfluid = density_superfluid
        self.density_normal = density_normal
        self.entropy = entropy
        self.viscosity = viscosity
        self.pressure = pressure
        self.v_s = v_s
        self.v_n = v_n
        
        # Read the CSV file for properties
        df = pd.read_csv(csv_file, header=None)
        self.temperatures = df[0].values
        self.densities = df[1].values
        self.entropies = df[2].values
        self.density_interpolator = interp1d(self.temperatures, self.densities, kind='linear', fill_value='extrapolate')

    def update_properties(self):
        self.density_superfluid = self.calculate_density_superfluid()
        self.density_normal = self.calculate_density_normal()
        self.entropy = self.calculate_entropy()

    def calculate_density_normal(self):
        return self.density_interpolator(self.temperature)

    def calculate_density_superfluid(self):
        rho_n = self.density_interpolator(self.temperature)
        rho_s = (rho_n / ((self.temperature / 2.17) ** 5.6)) - rho_n
        return rho_s

    def calculate_entropy(self):
        return self.entropy

    def calculate_k(self):
        L = num_nodes_y * (dy ** 2)
        total_density = self.density_superfluid + self.density_normal
        k = (((L * total_density * self.entropy) ** 2) / (8 * self.viscosity))
        return k

    def calculate_alpha(self):
        c_p = 5000  # Specific heat constant
        total_density = self.density_superfluid + self.density_normal
        alpha = self.calculate_k() / (total_density * c_p)
        return alpha

def initialize_grid():
    nodes = []
    for y in range(num_nodes_y):
        row = []
        for x in range(num_nodes_x):
            pressure = 101325 if x == 0 else 10132.5  # Pressure in Pascals
            node = Node(x, y, temperature=1.9, pressure=pressure)
            row.append(node)
        nodes.append(row)
    return nodes

def update_grid(nodes):
    for row in nodes:
        for node in row:
            node.update_properties()

def run_simulation():
    nodes = initialize_grid()
    for t in range(num_timesteps):
        update_grid(nodes)
        # Add logic for updating the grid based on FDM equations, if necessary
        # This can involve calculating the next state based on current states
        print(f'Timestep {t+1}/{num_timesteps} completed.')

# Run the simulation
run_simulation()








def update_temperature(nodes):
    new_temperatures = np.zeros((num_nodes_y, num_nodes_x))

    for y in range(1, num_nodes_y - 1):  # Avoid boundaries for now
        for x in range(1, num_nodes_x - 1):
            # Get the current node
            current_node = nodes[y][x]

            # Apply the finite difference scheme for temperature
            T_center = current_node.temperature
            T_left = nodes[y][x - 1].temperature
            T_right = nodes[y][x + 1].temperature
            T_up = nodes[y - 1][x].temperature
            T_down = nodes[y + 1][x].temperature

            # Calculate Laplacian
            laplacian_T = (T_left + T_right + T_up + T_down - 4 * T_center) / (dx ** 2)

            # Update temperature using explicit finite difference
            alpha = current_node.calculate_alpha()  # Assuming alpha is defined in the Node class
            new_temperature = T_center + alpha * laplacian_T * time_step
            new_temperatures[y, x] = new_temperature

    # Update the temperatures in the original nodes
    for y in range(1, num_nodes_y - 1):
        for x in range(1, num_nodes_x - 1):
            nodes[y][x].temperature = new_temperatures[y, x]

def run_simulation():
    nodes = initialize_grid()
    for t in range(num_timesteps):
        update_grid(nodes)  # Update physical properties based on current temperature
        update_temperature(nodes)  # Update temperature using finite difference logic
        # Add logic for updating the grid based on FDM equations, if necessary
        print(f'Timestep {t+1}/{num_timesteps} completed.')

# Run the simulation
run_simulation()
