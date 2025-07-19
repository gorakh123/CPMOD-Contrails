import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path

# Read csv and create df
script_folder = Path(__file__).resolve().parent
csv_path = script_folder / "plume_variables.csv"

plume_df = pd.read_csv(
    csv_path, 
    header=None, 
    sep=r'\s+')

plume_df.columns = [
    'T', 
    'smw', 
    'sw', 
    'si', 
    'Pw', 
    'time']


# Plot Temperature vs time
fig, ax = plt.subplots(
    figsize=(8,5))

ax.plot(
    plume_df['time'],
    plume_df['T'],
    c='black')

ax.set_xlabel("Time (s)")
ax.set_ylabel("Temperature (K)")
ax.grid()

plt.show()


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
