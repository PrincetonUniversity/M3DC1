import numpy as np
import matplotlib.pyplot as plt

base = 'bz_50'

r = np.loadtxt(f'{base}_r.csv', delimiter=',')
z = np.loadtxt(f'{base}_z.csv', delimiter=',')

with open(f'{base}_bz.csv') as f:
    lines = f.read().strip().splitlines()
bz = np.array([list(map(float, line.split(','))) for line in lines])

print(f"r:  {r.shape}  [{r[0]:.4f}, {r[-1]:.4f}]")
print(f"z:  {z.shape}  [{z[0]:.4f}, {z[-1]:.4f}]")
print(f"bz: {bz.shape}")

plt.figure(figsize=(7, 5))
plt.pcolormesh(r, z, bz, shading='auto')
plt.colorbar(label='B_Z (External)')
plt.xlabel('R')
plt.ylabel('Z')
plt.title('B_Z (External)')
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.show()
