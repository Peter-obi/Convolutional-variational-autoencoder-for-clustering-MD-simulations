import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances

# Set up the universe object with the topology from the PDB files
universe_traj1 = mda.Universe("final_pam.pdb", "pam.nc")
universe_traj2 = mda.Universe("final_nam.pdb", "nam.nc")

# Define the selection for Cα atoms
ca_selection = universe_traj1.select_atoms("name CA")

# Define the cutoff distance for contacts
cutoff = 8.0  # in Angstrom

# Extract contact maps for each trajectory
contact_maps_traj1 = []
for ts in universe_traj1.trajectory[::5]:
    dists = distances.distance_array(ca_selection.positions, ca_selection.positions)
    contacts = dists < cutoff
    contact_maps_traj1.append(contacts.astype(np.float32))

contact_maps_traj2 = []
for ts in universe_traj2.trajectory[::5]:
    dists = distances.distance_array(ca_selection.positions, ca_selection.positions)
    contacts = dists < cutoff
    contact_maps_traj2.append(contacts.astype(np.float32))

# Save the contact maps to files
np.save("contact_maps_traj1.npy", contact_maps_traj1)
np.save("contact_maps_traj2.npy", contact_maps_traj2)

