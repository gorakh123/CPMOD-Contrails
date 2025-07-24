import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
plt.rcParams['font.size'] = 14

# Read csv and create df
script_folder = Path(__file__).resolve().parent
csv_path = script_folder / "moments.csv"

moments_df = pd.read_csv(csv_path)
print(moments_df.columns)

for col in ['n_ice(#/m^3)', 'n_water(#/m^3)', 'n_soot(#/m^3)']:
    moments_df[col] = (
        moments_df[col]
          .astype(str)                     # ensure it’s a string
          .str.strip()                     # strip any whitespace
          .replace('', pd.NA)              # empty → NA
          .pipe(pd.to_numeric, errors='coerce')  # parse floats (and sci‑notation)
          .fillna(0)                       # decide how to handle missing → here 0
          .astype(int)                     # finally cast to int
    )

fig, ax = plt.subplots(figsize=(8,5))

ax.plot(
    moments_df[' time(s)'],
    moments_df['n_soot(#/m^3)'],
    label='soot',
    c='black')

ax.plot(
    moments_df[' time(s)'],
    moments_df['n_water(#/m^3)'],
    label='water',
    c='blue')

ax.plot(
    moments_df[' time(s)'],
    moments_df['n_ice(#/m^3)'],
    label='ice',
    c='cyan')

ax.set_xlabel("Time (s)")
ax.set_ylabel("number concentration (#/m^3)")
ax.set_ylim(0,1e12)
ax.grid()

plt.legend()
plt.show()

fig.savefig("number_conc_plot.png", format="png")


moments_df['r_i(nm)'] = moments_df['r_i(nm)'].str.strip()
moments_df['r_i(nm)'] = pd.to_numeric(moments_df['r_i(nm)'], errors='coerce')
moments_df['r_i(nm)'] = moments_df['r_i(nm)'].fillna(0)

moments_df['r_act(nm)'] = moments_df['r_act(nm)'].str.strip()
moments_df['r_act(nm)'] = pd.to_numeric(moments_df['r_act(nm)'], errors='coerce')
moments_df['r_act(nm)'] = moments_df['r_act(nm)'].fillna(0)

fig2, ax2 = plt.subplots(figsize=(8,5))
ax2.plot(
    moments_df[' time(s)'],
    moments_df['r_i(nm)'],
    c='cyan',
    label='ice'
    )

ax2.plot(
    moments_df[' time(s)'],
    moments_df['r_act(nm)'],
    c='blue',
    label='ice and water',
    linestyle='--')

ax2.set_xlabel("Time (s)")
ax2.set_ylabel("Radius (nm)")
plt.legend()
ax2.grid()

plt.show()
fig2.savefig("moments_plot.png", format='png')



