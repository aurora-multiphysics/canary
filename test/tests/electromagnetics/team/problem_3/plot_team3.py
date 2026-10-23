import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read CSV
df = pd.read_csv("clf_Aformsolve_out_line_sample_0001.csv")

# Compute magnitude of complex component 2
b_mag_2 = np.sqrt(
    df["b_field_real_2"]**2 +
    df["b_field_imag_2"]**2
)

# Plot magnitude along x_1
plt.figure(figsize=(8, 5))
plt.plot(df["x_1"], b_mag_2, "o-", linewidth=1.5)

plt.xlabel("x_1")
plt.ylabel(r"$|b_2|$")
plt.title("Magnitude of Complex Field Component 2 Along x_1")
plt.grid(True)
plt.tight_layout()

plt.show()
