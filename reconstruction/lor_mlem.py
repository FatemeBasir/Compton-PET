import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import glob
import os
from scipy.spatial import cKDTree

def read_root_file(file_name):
    print(f"Reading file: {file_name}")
    # Open the ROOT file
    file = uproot.open(file_name)
    
    # Access the tree named "Hits" for reconstructed data
    hits_tree = file["Hits"]
    
    # Access the tree named "Header" for MC truth data
    header_tree = file["Header"]

    # Extracting the relevant branches from Hits tree
    fX = hits_tree["fX"].array()
    fY = hits_tree["fY"].array()
    fZ = hits_tree["fZ"].array()
    fEvent = hits_tree["fEvent"].array()
    fTime = hits_tree["fTime"].array()
    fEnergy = hits_tree["fEnergy"].array()

    # Extracting the relevant branches from Header tree
    mcX = header_tree["mcX"].array()
    mcY = header_tree["mcY"].array()
    mcZ = header_tree["mcZ"].array()
    mcEvent = header_tree["mcEvent"].array()

    return fX, fY, fZ, fEvent, fTime, fEnergy, mcX, mcY, mcZ, mcEvent


def select_two_main_photons(indices, fEnergy):
    if len(indices) < 2:
        return indices
    
    energy_diffs = np.abs(fEnergy[indices] - 0.511)
    sorted_indices = indices[np.argsort(energy_diffs)]
    return sorted_indices[:2]


def calculate_lor_tof(fX, fY, fZ, fEvent, fTime, fEnergy, energy_window=(0.4, 0.6), time_window=1.0, distance_cut=1000.0, angle_range=(85, 95), azimuthal_angle_cut=5.0, c=299792458.0):
    interaction_points = []
    unique_events = np.unique(fEvent)

    counter = 0

    for event in unique_events:
        try:
            if(counter%(int(len(unique_events)/100)) == 0):
                print(str(int(counter/len(unique_events)*100))+'%')
        except Exception as e:
            print(e)

        counter += 1

        indices = np.where(fEvent == event)[0]
        if len(indices) < 2:
            continue

        indices = select_two_main_photons(indices, fEnergy)

        pos1 = np.array([fX[indices[0]], fY[indices[0]], fZ[indices[0]]])
        pos2 = np.array([fX[indices[1]], fY[indices[1]], fZ[indices[1]]])
        time1 = fTime[indices[0]]
        time2 = fTime[indices[1]]
        energy1 = fEnergy[indices[0]]
        energy2 = fEnergy[indices[1]]

        # Apply energy windowing
        if not (energy_window[0] <= energy1 <= energy_window[1] and energy_window[0] <= energy2 <= energy_window[1]):
            continue

        # Apply coincidence timing window
        delta_t = abs(time2 - time1)
        if delta_t > time_window:
            continue

        midpoint = (pos1 + pos2) / 2
        delta_t = time1 - time2

        # Convert TOF to distance in mm (c is in m/s, so convert to mm/ns)
        c_mm_ns = c * 1e-6
        distance_from_midpoint = (delta_t * c_mm_ns) / 2

        lor_length = np.linalg.norm(pos2 - pos1)

        if lor_length == 0 or lor_length > distance_cut:
            continue

        # Calculate the angle between the LOR and the detector plane
        detector_vector = pos2 - pos1
        z_axis = np.array([0, 0, 1])
        angle = np.degrees(np.arccos(np.dot(detector_vector, z_axis) / (np.linalg.norm(detector_vector) * np.linalg.norm(z_axis))))

        if not (angle_range[0] <= angle <= angle_range[1]):
            continue

        # Calculate azimuthal angle difference
        phi1 = np.arctan2(pos1[1], pos1[0])
        phi2 = np.arctan2(pos2[1], pos2[0])
        azimuthal_angle_diff = np.degrees(np.abs(phi2 - phi1))

        if azimuthal_angle_diff > 180:
            azimuthal_angle_diff = 360 - azimuthal_angle_diff
        
        if abs(180 - azimuthal_angle_diff > azimuthal_angle_cut):
            continue

        #print(azimuthal_angle_diff)

        interaction_point = midpoint + distance_from_midpoint * (pos2 - pos1) / lor_length
        #interaction_point = midpoint
        interaction_points.append(interaction_point)

        #print(fTime[indices[0]], fTime[indices[1]], pos1, pos2, interaction_point)

    return np.array(interaction_points)

def plot_interaction_points(interaction_points, mc_points, file_prefix):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    ax.scatter(interaction_points[:, 0], interaction_points[:, 1], interaction_points[:, 2], c='r', marker='o', label='Reconstructed Points')
    ax.scatter(mc_points[:, 0], mc_points[:, 1], mc_points[:, 2], c='b', marker='^', label='MC Truth Points')
    
    ax.set_xlabel('X (mm)')
    ax.set_ylabel('Y (mm)')
    ax.set_zlabel('Z (mm)')
    ax.set_title('Interaction Points Reconstruction')
    ax.legend()
    
    plt.savefig(f"{file_prefix}_3d.png")
    plt.savefig(f"{file_prefix}_3d.pdf")
    plt.close()

def calculate_signed_residuals(interaction_points, mc_points):
    # Use nearest neighbor approach to find the closest MC point for each interaction point
    tree = cKDTree(mc_points)
    distances, indices = tree.query(interaction_points)
    matched_mc_points = mc_points[indices]
    signed_residuals = np.linalg.norm(interaction_points - matched_mc_points, axis=1) * np.sign(interaction_points[:, 0] - matched_mc_points[:, 0])
    return signed_residuals

def plot_residuals_histogram(residuals, file_prefix):
    plt.hist(residuals, bins=50, alpha=0.75, edgecolor='black', density=True)
    plt.xlabel('Signed Residual (mm)')
    plt.ylabel('Frequency')
    plt.title('Signed Residuals Histogram')
    plt.grid(True)

    # Fit the central part of the residuals to a Gaussian distribution
    filtered_residuals = residuals[(residuals > -100) & (residuals < 100)]
    mu, std = norm.fit(filtered_residuals)
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)
    
    plt.plot(x, p, 'k', linewidth=2)
    title = f"Signed Residuals Histogram\nFit results: mu = {mu:.2f},  std = {std:.2f}"
    plt.title(title)
    
    plt.savefig(f"{file_prefix}_residuals.png")
    plt.savefig(f"{file_prefix}_residuals.pdf")
    plt.close()
    
    return mu, std

def initialize_image(interaction_points, image_size, voxel_size):
    """Initialize the image using the LOR-derived interaction points."""
    image = np.zeros(image_size)
    for point in interaction_points:
        voxel = tuple((point / voxel_size).astype(int))
        if all(0 <= v < s for v, s in zip(voxel, image_size)):
            image[voxel] += 1
    return image

def forward_project(image, system_matrix):
    """Forward project the image to the detector space."""
    return system_matrix @ image.flatten()

def back_project(ratio, system_matrix, image_shape):
    """Back project the ratio to the image space."""
    back_proj = system_matrix.T @ ratio
    return back_proj.reshape(image_shape)

def normalize_image(image):
    """Normalize the image to ensure non-negative values."""
    return np.maximum(image, 0)

def mlem_reconstruction(projections, system_matrix, initial_image, num_iterations=10):
    """Perform MLEM reconstruction."""
    image_shape = initial_image.shape
    image = initial_image.copy()

    for iteration in range(num_iterations):
        print(f"Iteration {iteration + 1}/{num_iterations}")
        estimated_projections = forward_project(image, system_matrix)
        ratio = projections / (estimated_projections + 1e-10)  # Avoid division by zero
        back_proj = back_project(ratio, system_matrix, image_shape)
        image *= back_proj
        image = normalize_image(image)

    return image

def plot_image(image, file_prefix, title="Reconstructed Image"):
    """Plot the reconstructed image."""
    plt.imshow(image.sum(axis=-1), cmap='gray')
    plt.title(title)
    plt.colorbar()
    plt.savefig(f"{file_prefix}_mlem_reconstructed.png")
    plt.savefig(f"{file_prefix}_mlem_reconstructed.pdf")
    plt.close()

def extract_interaction_points_from_image(image, voxel_size):
    """Extract interaction points from the reconstructed image."""
    points = []
    threshold = image.mean() + 2 * image.std()  # Simple thresholding
    indices = np.argwhere(image > threshold)
    for index in indices:
        point = index * voxel_size
        points.append(point)
    return np.array(points)

def process_file(file_name):
    file_prefix = os.path.splitext(os.path.basename(file_name))[0]

    # Read the data from the ROOT file
    fX, fY, fZ, fEvent, fTime, fEnergy, mcX, mcY, mcZ, mcEvent = read_root_file(file_name)
    
    # Calculate the interaction points using LOR and TOF
    interaction_points = calculate_lor_tof(fX, fY, fZ, fEvent, fTime, fEnergy)
    
    # Prepare the MC truth points with matching event numbers
    mc_points = []
    for event in np.unique(fEvent):
        mc_index = np.where(mcEvent == event)[0]
        if len(mc_index) > 0:
            mc_point = np.array([mcX[mc_index[0]], mcY[mc_index[0]], mcZ[mc_index[0]]])
            mc_points.append(mc_point)
    
    # Ensure mc_points and interaction_points are of the same length
    if len(interaction_points) > len(mc_points):
        interaction_points = interaction_points[:len(mc_points)]
    elif len(mc_points) > len(interaction_points):
        mc_points = mc_points[:len(interaction_points)]
    
    interaction_points = np.array(interaction_points)
    mc_points = np.array(mc_points)

    # Define the voxel size and image size
    voxel_size = 2 # Adjust voxel size as needed
    image_size = (64, 64, 64)  # Image dimensions

    # Create a simple system matrix (for illustration purposes)
    # The number of rows in the system matrix should match the number of projections
    # The number of columns should match the number of voxels in the image
    num_bins = 100
    num_voxels = np.prod(image_size)
    system_matrix = np.random.rand(num_bins, num_voxels)

    # Create simulated projections (for illustration purposes)
    # Replace this with your actual projections
    projections = np.random.rand(num_bins)

    print(interaction_points, image_size, voxel_size)

    # Initialize the image using LOR-derived interaction points
    initial_image = initialize_image(interaction_points, image_size, voxel_size)

    # Perform MLEM reconstruction
    num_iterations = 20
    reconstructed_image = mlem_reconstruction(projections, system_matrix, initial_image, num_iterations)

    # Extract interaction points from the MLEM reconstructed image
    mlem_interaction_points = extract_interaction_points_from_image(reconstructed_image, voxel_size)

    # Calculate the signed residuals based on the MLEM reconstructed interaction points
    signed_residuals = calculate_signed_residuals(mlem_interaction_points, mc_points)
    
    # Plot the 3D interaction points
    plot_interaction_points(interaction_points, mc_points, file_prefix)
    
    # Plot the residuals histogram and fit with Gaussian
    mu, std = plot_residuals_histogram(signed_residuals, file_prefix)

    # Plot the reconstructed image
    plot_image(reconstructed_image, file_prefix, title="MLEM Reconstructed Image")

    return signed_residuals, std

def main():
    # Specify the pattern for the output files
    file_pattern = "../build/output*_t*.root"
    
    # Get a list of all the output files
    files = glob.glob(file_pattern)
    
    all_residuals = []
    individual_resolutions = []
    
    for file_name in files:
        print(f"Processing file: {file_name}")
        residuals, resolution = process_file(file_name)
        all_residuals.extend(residuals)
        individual_resolutions.append(resolution)
        print(f"Spatial Resolution (Gaussian fit std) for {file_name}: {resolution:.2f} mm")
    
    # Convert to numpy array
    all_residuals = np.array(all_residuals)
    
    # Plot the combined residuals histogram and fit with Gaussian
    mu_combined, std_combined = plot_residuals_histogram(all_residuals, "combined")
    
    # Print out the combined spatial resolution
    print(f"Combined Spatial Resolution (Gaussian fit std): {std_combined:.2f} mm")
    print(f"Mean value from combined Gaussian fit: {mu_combined:.2f} mm")

if __name__ == "__main__":
    main()

