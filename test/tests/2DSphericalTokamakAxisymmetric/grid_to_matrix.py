import pandas as pd

src = "plasma_jtheta.csv"
dst = "plasma_jtheta_bilinear.csv"

df = pd.read_csv(src)

pivot = (
    df.pivot(
        index="Z_m",
        columns="R_m",
        values="Jtheta_A_per_m2",
    )
    .sort_index(axis=0)
    .sort_index(axis=1)
)

R = pivot.columns.to_numpy()
Z = pivot.index.to_numpy()

with open(dst, "w") as f:
    # First row = R coordinates
    f.write(",".join(f"{r:.17g}" for r in R) + "\n")

    # Remaining rows:
    # Z_j, J(R0,Zj), J(R1,Zj), ...
    for z, row in zip(Z, pivot.to_numpy()):
        f.write(
            f"{z:.17g},"
            + ",".join(f"{j:.17g}" for j in row)
            + "\n"
        )

print(f"Created {dst}")
print(f"Grid size = {len(R)} x {len(Z)}")
print(f"R = {R.min()} -> {R.max()}")
print(f"Z = {Z.min()} -> {Z.max()}")