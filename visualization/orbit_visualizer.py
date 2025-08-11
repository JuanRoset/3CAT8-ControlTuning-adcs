import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.basemap import Basemap
from matplotlib import rc

plt.rc('font', **{'family':'serif', 'serif':['Times'], 'size': 12})
rc('text', usetex=True)

# Function to generate points on a sphere in 3D
def generate_sphere_points(radius, num_points):
    phi = np.linspace(0, np.pi, num_points)
    theta = np.linspace(0, 2*np.pi, num_points)
    phi, theta = np.meshgrid(phi, theta)
    x = radius * np.sin(phi) * np.cos(theta)
    y = radius * np.sin(phi) * np.sin(theta)
    z = radius * np.cos(phi)
    return x, y, z

# Read the data from CSV
data = pd.read_csv('../data/orbit_file.csv')

# Extracting position data
time = data['Time']
x = data['PositionX']
y = data['PositionY']
z = data['PositionZ']

longitude = data['Longitude']
latitude = data['Latitude']

q_orbit_0 = data['QuaternionW']
q_orbit_1 = data['QuaternionX']
q_orbit_2 = data['QuaternionY']
q_orbit_3 = data['QuaternionZ']

# Identify the index where the longitude crosses +180 to -180
cross_idxs = [0]
for i in range(len(longitude) - 1):
    if abs(longitude[i + 1] - longitude[i]) > 180:
        cross_idxs.append(i + 1)
cross_idxs.append(len(longitude) - 1)

# Create 3D plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw={'width_ratios': [0.5, 1]})

ax1.axis('off')
ax2.axis('off')

ax1 = fig.add_subplot(121, projection='3d')

# Plot orbit
ax1.plot(x, y, z, label='Trajectory')

# Set aspect ratio to 'equal'
ax1.set_aspect('equal')

# Generate points for Earth's surface
radius = 6371 * 10**3  # Earth's radius in meters
num_points = 25
x_earth, y_earth, z_earth = generate_sphere_points(radius, num_points)

# Plot Earth as wireframe sphere
ax1.plot_wireframe(x_earth, y_earth, z_earth, color='g', alpha=0.3, label='Earth')

# Set labels and title
ax1.set_xlabel('$X_{ECI} (m)$')
ax1.set_ylabel('$Y_{ECI} (m)$')
ax1.set_zlabel('$Z_{ECI} (m)$')
ax1.set_title('Orbit trajectory')
ax1.legend(frameon=False)

# Add Basemap on ax2
ax2 = fig.add_subplot(122)
m = Basemap(ax=ax2)
m.drawcoastlines()
m.drawcountries()
m.drawparallels(np.arange(-90., 91., 30.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 60.), labels=[0, 0, 0, 1])

# Plot ground track by segments
for i in range(len(cross_idxs) - 1):
    # Split the data into two segments
    longitude_segment = longitude[cross_idxs[i]:cross_idxs[i + 1]]
    latitude_segment = latitude[cross_idxs[i]:cross_idxs[i + 1]]

    # print x and y (longitude, latitude)
    x, y = m(longitude_segment, latitude_segment)
    m.plot(x, y, marker=None, color='r', linestyle='-')

# Set title
ax2.set_title('Ground track')

# Show plot
plt.show()