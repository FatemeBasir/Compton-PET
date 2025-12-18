import os
import sys
import glob
import uproot
import numpy as np
import matplotlib.pyplot as plt

from lor_mlem import calculate_lor_tof


def read_root_files():
    desired_pdg_codes = [22, -22]

    # Use glob to match file patterns
    file_paths = glob.glob("../build/bin/*.root")

    n_files = len(file_paths)

    if len(sys.argv) > 1:
        print(sys.argv[1])
        n_files = int(sys.argv[1])
        

    # List to hold valid files
    valid_files_with_trees = []
    valid_files_with_trees2 = []

    for path in file_paths[0:n_files]:
        print(f'Reading file {path}')
        try:
            # Try to open the file and access the "Hits" tree
            with uproot.open(path) as file:
                if "Hits" in file:
                    valid_files_with_trees.append(f"{path}:Hits")
                    valid_files_with_trees2.append(f"{path}:Header")
        except (uproot.exceptions.KeyInFileError, ValueError, OSError) as e:
            # Handle cases where the file is damaged or does not contain the "Hits" tree
            print(f"Skipping file {path} due to error: {e}")

    # Concatenate valid files
    if valid_files_with_trees:
        arrays = uproot.concatenate(valid_files_with_trees, library="np")
        arrays2 = uproot.concatenate(valid_files_with_trees2, library="np")
        
        # Extracting the relevant branches from the concatenated arrays
        fX = arrays["fX"]
        fY = arrays["fY"]
        fZ = arrays["fZ"]
        fEvent = arrays["fEvent"]
        fTime = arrays["fTime"]
        fPDG = arrays["fPDG"]
        fEnergy = arrays["fEnergy"]

        # Create a mask to filter the desired PDG codes
        mask = np.isin(fPDG, desired_pdg_codes)

        # Apply the mask to filter the data
        fX = fX[mask]
        fY = fY[mask]
        fZ = fZ[mask]
        fEvent = fEvent[mask]
        fTime = fTime[mask]
        fEnergy = fEnergy[mask]

        # Extracting the relevant branches from Header tree
        mcX = arrays2["mcX"]
        mcY = arrays2["mcY"]
        mcZ = arrays2["mcZ"]
        mcEvent = arrays2["mcEvent"]

        return fX, fY, fZ, fEvent, fTime, fEnergy, mcX, mcY, mcZ, mcEvent
    else:
        print("No valid ROOT files found.")
        return None, None, None, None, None


def plot_reco(points):
    fig = plt.figure(figsize=(7, 6))
    plt.hist2d(points[:, 0], points[:, 1], bins=100, cmap='gray')
    
    plt.ylim(-30, 30)
    plt.xlim(-30, 30)
    plt.xlabel('X [mm]')
    plt.ylabel('X [mm]')
    plt.colorbar()
    plt.tight_layout()
    #plt.savefig('plot.pdf')
    plt.savefig('reco.png')
    #plt.show()

    fig = plt.figure(figsize=(7, 6))
    plt.hist2d(points[:, 0], points[:, 1], bins=100, cmap='rainbow')
    
    plt.ylim(-30, 30)
    plt.xlim(-30, 30)
    plt.xlabel('X [mm]')
    plt.ylabel('X [mm]')
    plt.colorbar()
    plt.tight_layout()
    #plt.savefig('plot.pdf')
    plt.savefig('reco_col.png')
    #plt.show()


def plot_header(points):
    fig = plt.figure(figsize=(6, 6))
    plt.plot(points[0], points[1], 'o')
    plt.ylim(-30, 30)
    plt.xlim(-30, 30)
    plt.xlabel('X [mm]')
    plt.ylabel('X [mm]')
    plt.tight_layout()
    #plt.savefig('plot.pdf')
    plt.savefig('header.png')
    #plt.show()


# Read the data from the ROOT file
fX, fY, fZ, fEvent, fTime, fEnergy, mcX, mcY, mcZ, mcEvent = read_root_files()

# Calculate the LOR for reconstruction
interaction_points = calculate_lor_tof(fX, fY, fZ, fEvent, fTime, fEnergy, energy_window=(0.4, 0.6), time_window = 1.0, angle_range=(85, 95), azimuthal_angle_cut=5.0, c=299792458.0)

header = [mcX, mcY]

print(header)

# Plot the reconstructed data
plot_reco(interaction_points)
plot_header(header)

print(interaction_points)
