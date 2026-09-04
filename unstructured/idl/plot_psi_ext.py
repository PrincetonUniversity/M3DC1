import numpy as np
import matplotlib.pyplot as plt
from scipy.io import readsav

infile = 'external_field.sav'

sav = readsav(infile)
r = sav['rgrid']
z = sav['zgrid']
psi_ext = sav['psi_ext']
phi0 = sav['phi0']

print(f"r:       {r.shape}  [{r.min():.4f}, {r.max():.4f}]")
print(f"z:       {z.shape}  [{z.min():.4f}, {z.max():.4f}]")
print(f"psi_ext: {psi_ext.shape}")

plt.figure(figsize=(7, 5))
plt.pcolormesh(r, z, psi_ext, shading='auto')
plt.colorbar(label='psi_ext')
plt.xlabel('R')
plt.ylabel('Z')
plt.title(f'psi_ext (phi = {float(phi0):.4g})')
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.show()
