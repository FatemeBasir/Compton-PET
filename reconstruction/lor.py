import uproot
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm
import glob
import os
from tqdm import tqdm

all_interaction_points = []
all_lor_lengths = []
all_delta_ts = []
all_angles = []
all_azimuth_angles = []
all_energy = []

def plot_data_cuts(filename, data, cuts, name="", axis="", limits=(0.0, 1.0), bins=50):
    filename += "_" + name.replace(" ", "_").lower()

    print("Plotting data...", filename)

    fig = plt.figure()
    plt.hist(data, bins=50, alpha=0.75, edgecolor="black", range=limits)
    plt.xlabel(axis)
    plt.ylabel("Frequency")
    plt.title(name)
    plt.grid(True)
    plt.axvline(cuts[0], color="r", linestyle="dashed", linewidth=1)
    plt.axvline(cuts[1], color="r", linestyle="dashed", linewidth=1)
    plt.savefig(filename + ".pdf")
    plt.savefig(filename + ".png")
    plt.close()


def read_root_file(file_name):
    print(f"Reading file: {file_name}")
    # Open the ROOT file
    file = uproot.open(file_name)

    # Access the tree named "Hits" for reconstructed data
    hits_tree = file["Hits"]

    # Access the tree named "Header" for MC truth data
    header_tree = file["Header"]

    # Extracting the relevant branches from Hits tree
    fX = hits_tree["fX"].array(library="np")
    fY = hits_tree["fY"].array(library="np")
    fZ = hits_tree["fZ"].array(library="np")
    fEvent = hits_tree["fEvent"].array(library="np")
    fTime = hits_tree["fTime"].array(library="np")
    fEnergy = hits_tree["fEnergy"].array(library="np")

    # Extracting the relevant branches from Header tree
    mcX = header_tree["mcX"].array(library="np")
    mcY = header_tree["mcY"].array(library="np")
    mcZ = header_tree["mcZ"].array(library="np")
    mcEvent = header_tree["mcEvent"].array(library="np")

    print("Events", fEvent, fX, fY, fZ)
    print("MC Truth", mcEvent, mcX, mcY, mcZ)

    all_energy.append(fEnergy)

    return fX, fY, fZ, fEvent, fTime, fEnergy, mcX, mcY, mcZ, mcEvent


def calculate_lor_tof(
    file_prefix,
    fX,
    fY,
    fZ,
    fEvent,
    fTime,
    fEnergy,
    time_window=0.5,
    energy_window=(0.45, 0.55),
    distance_cut=1200.0,
    angle_range=(85, 95),
    azimuthal_angle_cut=10.0,
    c=299792458.0,
    no_select=False,
):
    interaction_points = []
    lor_lengths = []
    delta_ts = []
    angles = []
    azimuth_angles = []

    unique_events = np.unique(fEvent)

    for event in unique_events:
        indices = np.where(fEvent == event)[0]
        if len(indices) != 2:
            continue

        pos1 = np.array([fX[indices[0]], fY[indices[0]], fZ[indices[0]]])
        pos2 = np.array([fX[indices[1]], fY[indices[1]], fZ[indices[1]]])
        time1 = fTime[indices[0]]
        time2 = fTime[indices[1]]
        energy1 = fEnergy[indices[0]]
        energy2 = fEnergy[indices[1]]

        midpoint = (pos1 + pos2) / 2
        delta_t = time1 - time2

        delta_ts.append(delta_t)
        all_delta_ts.append(delta_t)

        # Convert TOF to distance in mm (c is in m/s, so convert to mm/ns)
        c_mm_ns = c * 1e-6
        distance_from_midpoint = (delta_t * c_mm_ns) / 2

        lor_length = np.linalg.norm(pos2 - pos1)
        lor_lengths.append(lor_length)
        all_lor_lengths.append(lor_length)

        # Calculate the angle between the LOR and the detector plane
        detector_vector = pos2 - pos1
        z_axis = np.array([0, 0, 1])
        angle = np.degrees(
            np.arccos(
                np.dot(detector_vector, z_axis)
                / (np.linalg.norm(detector_vector) * np.linalg.norm(z_axis))
            )
        )

        angles.append(angle)
        all_angles.append(angle)

        # Calculate azimuthal angle difference
        phi1 = np.arctan2(pos1[1], pos1[0])
        phi2 = np.arctan2(pos2[1], pos2[0])
        azimuthal_angle_diff = np.degrees(np.abs(phi2 - phi1))

        if azimuthal_angle_diff > 180:
            azimuthal_angle_diff = 360 - azimuthal_angle_diff

        azimuth_angles.append(azimuthal_angle_diff)
        all_azimuth_angles.append(azimuthal_angle_diff)

        if not no_select:
            # Apply energy windowing
            if not (
                energy_window[0] <= energy1 <= energy_window[1]
                and energy_window[0] <= energy2 <= energy_window[1]
            ):
                continue

            # Apply coincidence timing window
            delta_t = abs(time2 - time1)
            if delta_t > time_window:
                continue

            # Apply distance cut
            if lor_length == 0 or lor_length > distance_cut:
                continue

            # Apply polar angle cut
            if not (angle_range[0] <= angle <= angle_range[1]):
                continue

            # Apply azimuth angle cut
            if abs(180 - azimuthal_angle_diff > azimuthal_angle_cut):
                continue

        interaction_point = (
            midpoint + distance_from_midpoint * (pos2 - pos1) / lor_length
        )
        interaction_points.append(interaction_point)

    # Plotting data with cuts
    plot_data_cuts(
        file_prefix,
        fEnergy,
        energy_window,
        "Energy Distribution",
        "Energy [MeV]",
        (0.0, 1.0),
    )
    plot_data_cuts(
        file_prefix,
        lor_lengths,
        [0.0, distance_cut],
        "LOR Distances",
        "Distance [mm]",
        (1.0, 2000.0),
    )
    plot_data_cuts(
        file_prefix,
        delta_ts,
        [-time_window, time_window],
        "Time Window",
        "Time [ns]",
        (-3.0, 3.0),
    )
    plot_data_cuts(
        file_prefix, angles, angle_range, "Polar Angles", "Angle [deg]", (0.0, 180.0)
    )
    plot_data_cuts(
        file_prefix,
        azimuth_angles,
        [175.0, 185.0],
        "Azimuth Angles",
        "Angle [deg]",
        (90.0, 270.0),
        100.0,
    )

    return np.array(interaction_points)


def plot_interaction_points(interaction_points, mc_points, file_prefix):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    ax.scatter(
        interaction_points[:, 0],
        interaction_points[:, 1],
        interaction_points[:, 2],
        c="r",
        marker="o",
        label="Reconstructed Points",
    )
    ax.scatter(
        mc_points[:, 0],
        mc_points[:, 1],
        mc_points[:, 2],
        c="b",
        marker="^",
        label="MC Truth Points",
    )

    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    ax.set_zlabel("Z (mm)")
    ax.set_title("Interaction Points Reconstruction")
    ax.legend()

    plt.savefig(f"{file_prefix}_3d.png")
    plt.savefig(f"{file_prefix}_3d.pdf")
    plt.close()


def calculate_signed_residuals(interaction_points, mc_points):
    # Calculate the signed residuals (distances between reconstructed points and MC truth points)
    signed_residuals = np.linalg.norm(interaction_points - mc_points, axis=1) * np.sign(
        interaction_points[:, 0] - mc_points[:, 0]
    )
    return signed_residuals


def plot_residuals_histogram(residuals, file_prefix):
    plt.hist(
        residuals,
        bins=50,
        alpha=0.75,
        edgecolor="black",
        range=(-10.0, 10.0),
        density=True
    )
    plt.xlabel("Signed Residual (mm)")
    plt.ylabel("Frequency")
    plt.title("Signed Residuals Histogram")
    plt.grid(True)

    # Fit the central part of the residuals to a Gaussian distribution
    filtered_residuals = residuals[(residuals > -10) & (residuals < 10)]
    mu, std = norm.fit(filtered_residuals)
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)

    plt.plot(x, p, "k", linewidth=2)
    title = f"Signed Residuals Histogram\nFit results: mu = {mu:.2f},  std = {std:.2f}"
    plt.title(title)

    plt.savefig(f"{file_prefix}_residuals.png")
    plt.savefig(f"{file_prefix}_residuals.pdf")
    plt.close()

    return mu, std


def process_file(file_name):

    file_prefix = os.path.splitext(os.path.basename(file_name))[0]

    # Read the data from the ROOT file
    fX, fY, fZ, fEvent, fTime, fEnergy, mcX, mcY, mcZ, mcEvent = read_root_file(
        file_name
    )

    # Calculate the interaction points using LOR and TOF
    interaction_points = calculate_lor_tof(
        file_prefix, fX, fY, fZ, fEvent, fTime, fEnergy, no_select=False
    )

    # Prepare the MC truth points with matching event numbers
    mc_points = []
    for event in np.unique(fEvent):
        mc_index = np.where(mcEvent == event)[0]
        if len(mc_index) > 0:
            mc_point = np.array([mcX[mc_index[0]], mcY[mc_index[0]], mcZ[mc_index[0]]])
            mc_points.append(mc_point)

    # Ensure mc_points and interaction_points are of the same length
    if len(interaction_points) > len(mc_points):
        interaction_points = interaction_points[: len(mc_points)]
    elif len(mc_points) > len(interaction_points):
        mc_points = mc_points[: len(interaction_points)]

    interaction_points = np.array(interaction_points)
    mc_points = np.array(mc_points)

    # Calculate the signed residuals
    signed_residuals = calculate_signed_residuals(interaction_points, mc_points)

    # Plot the 3D interaction points
    plot_interaction_points(interaction_points, mc_points, file_prefix)

    # Plot the residuals histogram and fit with Gaussian
    mu, std = plot_residuals_histogram(signed_residuals, file_prefix)

    return signed_residuals, std


def main():
    # Specify the pattern for the output files
    file_pattern = "../build/bin/output*_t*.root"

    # Get a list of all the output files
    files = glob.glob(file_pattern)[0:100]

    all_residuals = []
    individual_resolutions = []

    for file_name in tqdm(files, position=0, leave=True):
        print(f"Processing file: {file_name}")
        residuals, resolution = process_file(file_name)
        all_residuals.extend(residuals)
        individual_resolutions.append(resolution)
        print(
            f"Spatial Resolution (Gaussian fit std) for {file_name}: {resolution:.2f} mm"
        )

    time_window=0.5
    energy_window=(0.45, 0.55)
    distance_cut=1200.0
    angle_range=(85, 95)
    azimuthal_angle_cut=10.0

    full_energy = np.concatenate(all_energy)

    # Plotting data with cuts
    plot_data_cuts(
        'combined',
        full_energy,
        energy_window,
        "Combined Energy Distribution",
        "Energy [MeV]",
        (0.1, 1.0),
        bins=100
    )
    plot_data_cuts(
        'combined',
        all_lor_lengths,
        [0.0, distance_cut],
        "Combined LOR Distances",
        "Distance [mm]",
        (1.0, 2000.0),
    )
    plot_data_cuts(
        'combined',
        all_delta_ts,
        [-time_window, time_window],
        "Combined Time Window",
        "Time [ns]",
        (-1.0, 1.0),
    )
    plot_data_cuts(
        'combined', all_angles, angle_range, "Combined Polar Angles", "Angle [deg]", (0.0, 180.0)
    )
    plot_data_cuts(
        'combined',
        all_azimuth_angles,
        [175.0, 185.0],
        "Combined Azimuth Angles",
        "Angle [deg]",
        (90.0, 200.0),
        100.0,
    )


    # Convert to numpy array
    all_residuals_np = np.array(all_residuals)

    # Plot the combined residuals histogram and fit with Gaussian
    mu_combined, std_combined = plot_residuals_histogram(
        all_residuals_np, file_prefix="combined"
    )

    # Print out the combined spatial resolution
    print(f"Combined Spatial Resolution (Gaussian fit std): {std_combined:.2f} mm")
    print(f"Mean value from combined Gaussian fit: {mu_combined:.2f} mm")

if __name__ == "__main__":
    main()
