import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from scipy.integrate import cumulative_trapezoid

plt.style.use("../styles/pretty_plots.mplstyle")

# Read the data from CSV
data = pd.read_csv('../data/simulation_output.csv')

# Specify orbital period
period = 1.41 * 3600

plot_type = ["detumbling", "nadir"]
plot_type = plot_type[1]

# Extracting angular rate data
time = [data['Time'], data['Time'] / 3600, data['Time'] / period]
timelabel = ["Time (s)", "Time (h)", "Time (orbits)"]

time_type = ["seconds", "hours", "orbits"]
time_type_idx = 1

time = time[time_type_idx]
timelabel = timelabel[time_type_idx]

# Extract magnetorquer current and power data
i_mtq_x = data['Current_mtq_X']
i_mtq_y = data['Current_mtq_Y']
i_mtq_z = data['Current_mtq_Z']
i_total_mtq = abs(i_mtq_x) + abs(i_mtq_y) + abs(i_mtq_z)

p_mtq_x = data['Power_mtq_X']
p_mtq_y = data['Power_mtq_Y']
p_mtq_z = data['Power_mtq_Z']
p_total_mtq = abs(p_mtq_x) + abs(p_mtq_y) + abs(p_mtq_z)

# Extract reaction wheels current and power data
rate_RW = data['AngRate_RW']
i_RW = data['Current_RW']
p_RW = data['Power_RW']

# Compute the cumulative integrals
e_total_mtq = cumulative_trapezoid(p_total_mtq, time, initial=0)
e_total_RW = cumulative_trapezoid(p_RW, time, initial=0)

# Sum the cumulative integrals
e_total = e_total_mtq + e_total_RW

# Compute total power and energy being drawn
i_total = i_total_mtq + i_RW
p_total = p_total_mtq + p_RW

# Extract the disturbance torques
gg_torque_x = data['Gg_Torque_X']
gg_torque_y = data['Gg_Torque_Y']
gg_torque_z = data['Gg_Torque_Z']

air_torque_x = data['Air_Torque_X']
air_torque_y = data['Air_Torque_Y']
air_torque_z = data['Air_Torque_Z']

srp_torque_x = data['SRP_Torque_X']
srp_torque_y = data['SRP_Torque_Y']
srp_torque_z = data['SRP_Torque_Z']

# Extract the control torques
mtq_torque_x = data['Torque_mtq_X']
mtq_torque_y = data['Torque_mtq_Y']
mtq_torque_z = data['Torque_mtq_Z']

rw_torque_z = data['Torque_RW']

# Extracting angular rate data
omega_x = data['AngularRateX']
omega_y = data['AngularRateY']
omega_z = data['AngularRateZ']

# Extracting error quaternion data
error_0 = data['Quaternion_error_W']
error_1 = data['Quaternion_error_X']
error_2 = data['Quaternion_error_Y']
error_3 = data['Quaternion_error_Z']

# Process data
omega_x = omega_x * 180 / np.pi
omega_y = omega_y * 180 / np.pi
omega_z = omega_z * 180 / np.pi

error_0 = abs(np.arccos(error_0) * 2) * 180 / np.pi

# Create plot with specified colors, grid, and font
fig0, (ax00, ax01, ax02) = plt.subplots(3, 1, figsize=(16, 14))

# Plotting on ax00
ax00.plot(time, omega_x, label='X body axis')
ax00.plot(time, omega_y, label='Y body axis')
ax00.plot(time, omega_z, label='Z body axis')
ax00.set_ylabel('Angular rate (º/s)')
ax00.legend(loc = 'lower right')
ax00.grid(True)
ax00.set_xlim(min(time), max(time))
ax00.tick_params(axis='x', labelbottom=False)
#ax00.tick_params(axis='both', which='major')

# Plotting on ax01
ax01.plot(time, error_1, label='X body axis')
ax01.plot(time, error_2, label='Y body axis')
ax01.plot(time, error_3, label='Z body axis')
ax01.set_ylabel('Error quaternion')
ax01.legend(loc = 'lower right')
ax01.grid(True)
ax01.set_xlim(min(time), max(time))
ax01.tick_params(axis='x', labelbottom=False)
#ax01.tick_params(axis='both', which='major')

# Plotting on ax02
ax02.plot(time, error_0, label='Total')  # Plotting only one line for clarity
ax02.set_xlabel(timelabel)
ax02.set_ylabel('Attitude error (º)')
ax02.legend(loc = 'lower right')
ax02.grid(True)
if plot_type == "detumbling":
    ax02.set_ylim(0, 180)
ax02.set_xlim(min(time), max(time))
ax02.tick_params(axis='both', which='major')

# Adjusting vertical space between subplots and reducing external margins
fig0.subplots_adjust(hspace=0.1, left=0.1, right=0.9, top=0.95, bottom=0.1)




# Plot the actuator's currents, power and energy
fig1, (ax10, ax11, ax12) = plt.subplots(3, 1, figsize=(16, 14))

# Plotting on ax10
ax10.plot(time, i_total, label='Total', color='black')
ax10.plot(time, i_mtq_x, label='MTQ-X')
ax10.plot(time, i_mtq_y, label='MTQ-Y')
ax10.plot(time, i_mtq_z, label='MTQ-Z')
ax10.plot(time, i_RW, label='RW-Z')
ax10.set_ylabel('Current (A)')
ax10.legend(loc = 'upper right')
ax10.grid(True)
ax10.set_xlim(min(time), max(time))
ax10.tick_params(axis='x', labelbottom=False)

# Plotting on ax11
ax11.plot(time, p_total, label='Total', color='black')
ax11.plot(time, p_mtq_x, label='MTQ-X')
ax11.plot(time, p_mtq_y, label='MTQ-Y')
ax11.plot(time, p_mtq_z, label='MTQ-Z')
ax11.plot(time, p_RW, label='RW-Z')
ax11.set_ylabel('Power (W)')
ax11.legend(loc = 'upper right')
ax11.grid(True)
ax11.set_xlim(min(time), max(time))
ax11.tick_params(axis='x', labelbottom=False)

ax12.plot(time, e_total, label='Total', color='black')
ax12.set_xlabel(timelabel)
ax12.set_ylabel('Energy consumed (J)')
ax12.legend(loc = 'right')
ax12.grid(True)
ax12.set_xlim(min(time), max(time))
ax12.tick_params(axis='both', which='major')

# Adjusting vertical space between subplots and reducing external margins
fig1.subplots_adjust(hspace=0.1, left=0.1, right=0.9, top=0.95, bottom=0.1)




# Plot the disturbance torques and currents
fig2, (ax20, ax21, ax22) = plt.subplots(3, 1, figsize=(16, 14))

ax20.plot(time, gg_torque_x, label='Gravitational torque X')
ax20.plot(time, gg_torque_y, label='Gravitational torque Y')
ax20.plot(time, gg_torque_z, label='Gravitational torque Z')
ax20.set_ylabel('Torque (Nm)')
ax20.legend(loc = 'upper right')
ax20.grid(True)
ax20.set_xlim(min(time), max(time))
ax20.tick_params(axis='x', labelbottom=False)
 
ax21.plot(time, air_torque_x, label='Aerodynamic torque X')
ax21.plot(time, air_torque_y, label='Aerodynamic torque Y')
ax21.plot(time, air_torque_z, label='Aerodynamic torque Z')
ax21.set_ylabel('Torque (Nm)')
ax21.legend(loc = 'upper right')
ax21.grid(True)
ax21.set_xlim(min(time), max(time))
ax21.tick_params(axis='x', labelbottom=False)

ax22.plot(time, srp_torque_x, label='SRP torque X')
ax22.plot(time, srp_torque_y, label='SRP torque Y')
ax22.plot(time, srp_torque_z, label='SRP torque Z')
ax22.set_xlabel(timelabel)
ax22.set_ylabel('Torque (Nm)')
ax22.legend(loc = 'upper right')
ax22.grid(True)
ax22.set_xlim(min(time), max(time))
ax22.tick_params(axis='both', which='major')

# Adjusting vertical space between subplots and reducing external margins
fig2.subplots_adjust(hspace=0.1, left=0.1, right=0.9, top=0.95, bottom=0.1)




# Plot the actuator's currents, power and energy
fig3, (ax30, ax31) = plt.subplots(2, 1, figsize=(16, 8))

# Plotting on ax10
ax30.plot(time, mtq_torque_x, label='MTQ-X')
ax30.plot(time, mtq_torque_y, label='MTQ-Y')
ax30.plot(time, mtq_torque_z, label='MTQ-Z')
ax30.plot(time, rw_torque_z, label='RW-Z')
ax30.set_ylabel('Torque (Nm)')
ax30.legend(loc = 'upper right')
ax30.grid(True)
ax30.set_xlim(min(time), max(time))
if plot_type == "detumbling":
    ax30.tick_params(axis='both', which='major')
    ax30.set_xlabel(timelabel)
else:
    ax30.tick_params(axis='x', labelbottom=False)

# Plotting on ax11
ax31.plot(time, rate_RW)  # Plotting only one line for clarity
ax31.set_xlabel(timelabel)
ax31.set_ylabel('RW angular rate (rad/s)')
ax31.grid(True)
ax31.set_xlim(min(time), max(time))
ax31.tick_params(axis='both', which='major')

# Adjusting vertical space between subplots and reducing external margins
if plot_type == "detumbling":
    fig3.subplots_adjust(hspace=0.4, left=0.1, right=0.9, top=0.95, bottom=0.1)
else:
    fig3.subplots_adjust(hspace=0.1, left=0.1, right=0.9, top=0.95, bottom=0.1)

 
plt.show()