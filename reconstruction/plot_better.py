import os
import glob
import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from lor_mlem import calculate_lor_tof


def read_root_files_in_chunks(chunk_size):
    desired_pdg_codes = [22, -22]

    # Use glob to match file patterns
    file_paths = glob.glob("../build/*.root")

    # List to hold valid files
    valid_files_with_trees = []
    valid_files_with_trees2 = []

    for path in file_paths:
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
        for arrays in uproot.iterate(valid_files_with_trees, library="np", step_size=chunk_size):
            # Extracting the relevant branches from the concatenated arrays
            fX = arrays["fX"]
            fY = arrays["fY"]
            fZ = arrays["fZ"]
            fEvent = arrays["fEvent"]
            fTime = arrays["fTime"]
            fPDG = arrays["fPDG"]
            fEnergy = arrays["fEnergy"]
            
            arrays2 = uproot.concatenate(valid_files_with_trees2, library="np", step_size=chunk_size)

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

            yield fX, fY, fZ, fEvent, fTime, fEnergy, mcX, mcY, mcZ, mcEvent
        else:
            print("No valid ROOT files found.")
            yield None, None, None, None, None



def plot_reco(points):
    fig = plt.figure(figsize=(6, 6))
    plt.hist2d(points[:, 0], points[:, 1], bins=250, cmap="rainbow")
    plt.ylim(-50, 50)
    plt.xlim(-50, 50)
    plt.tight_layout()
    #plt.savefig('plot.pdf')
    plt.savefig('reco.png')
    #plt.show()


def plot_header(points):
    fig = plt.figure(figsize=(6, 6))
    plt.plot(points[0], points[1], 'o')
    plt.ylim(-50, 50)
    plt.xlim(-50, 50)
    plt.tight_layout()
    #plt.savefig('plot.pdf')
    plt.savefig('header.png')
    #plt.show()


# Read the data from the ROOT file
#fX, fY, fZ, fEvent, fTime, fEnergy, mcX, mcY, mcZ, mcEvent = read_root_files()

chunk_size = 100000  # Number of events to process in each chunk

# Read and process the data in chunks
all_interaction_points = []

for iteration, (fX, fY, fZ, fEvent, fTime, fEnergy, mcX, mcY, mcZ, mcEvent) in enumerate(read_root_files_in_chunks(chunk_size)):
    if fX.size > 0 and fY.size > 0 and fZ.size > 0 and fEvent.size > 0 and fTime.size > 0:
        # Calculate the LOR for reconstruction
        interaction_points = calculate_lor_tof(fX, fY, fZ, fEvent, fTime, fEnergy, energy_window=(0.49, 0.52), time_window = 1.0, distance_cut=100.0, c=299792458.0)

        # Ensure interaction_points is not empty
        if interaction_points.size > 0:

            # Append to the list of all interaction points
            all_interaction_points.append(interaction_points)

            if len(all_interaction_points) > 0:
                # Convert to numpy array for plotting
                all_interaction_points_array = np.vstack(all_interaction_points)

                header = [mcX, mcY]

                print(header)

                # Plot the reconstructed data
                plot_reco(all_interaction_points_array)
                plot_header(header)

                print(interaction_points)
