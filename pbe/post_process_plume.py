import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
plt.rcParams['font.size'] = 14



# Read csv and create df
script_folder = Path(__file__).resolve().parent
csv_path = script_folder / "plume_variables.csv"

plume_df = pd.read_csv(csv_path)

print(plume_df.columns)
# Plot Temperature vs time
fig, ax = plt.subplots(
    figsize=(8,5))

ax.plot(
    plume_df['time'],
    plume_df[' Temperature (K)'],
    c='black')

ax.set_xlabel("Time (s)")
ax.set_ylabel("Temperature (K)")
ax.grid()

plt.show()
fig.savefig('Temp_plot.png', format='png')


# Plot smw, sw, si vs time
fig2, ax2 = plt.subplots(
    figsize=(8,5))

ax2.plot(
    plume_df['time'],
    plume_df['smw'],
    c='black')

ax2.plot(
    plume_df['time'],
    plume_df['sw'],
    c='blue')

ax2.plot(
    plume_df['time'],
    plume_df['si'],
    c='cyan')

ax2.set_xlabel("Time (s)")
ax2.set_ylabel("Super saturation ratio")
ax2.grid()

plt.show()
fig2.savefig('Saturation_plot.png', format='png')


def pw_sat(T):
    return (np.exp(54.842763 - 6763.22/T - 4.210 * np.log(T) + 0.000367 * T + np.tanh(0.0415 * (T - 218.8)) * 
    (53.878 - 1331.22/T - 9.44523 * np.log(T) + 0.014025 * T)))

def pi_sat(T):
    return (np.exp(9.550426 - 5723.265/T + 3.53068 * np.log(T) -  0.00728332 * T ))

pwsat = pw_sat(plume_df[' Temperature (K)'])
pisat = pi_sat(plume_df[' Temperature (K)'])
# Plot Pvap vs Temp
fig3, ax3 = plt.subplots(
    figsize=(8,5))

ax3.plot(
    plume_df[' Temperature (K)'],
    pwsat,
    c='blue')

ax3.plot(
    plume_df[' Temperature (K)'],
    pisat,
    c='cyan')

ax3.plot(
    plume_df[' Temperature (K)'],
    plume_df['Vapour Pressure (Pa)'],
    c='black'
)

ax3.set_xlim(220,260)
ax3.set_xlabel('Temperature (K)')
ax3.set_ylabel('Pressure (Pa)')
ax3.set_ylim(0,100)
plt.show()
fig3.savefig('Vapour_Pressure_plot.png', format='png')
