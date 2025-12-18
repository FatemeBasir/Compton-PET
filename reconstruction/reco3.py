import ROOT
#from ROOT import *
import numpy as np
from scipy.ndimage import rotate
import matplotlib.pyplot as plt

# Open the ROOT file
file = ROOT.TFile("../MF3-build/output0_t0.root")
tree = file.Get("Hits")

# Create arrays to store the data
fX = []
fY = []
fZ = []

# Loop over the entries in the tree and extract the data
for event in tree:
    fX.append(event.fX)
    fY.append(event.fY)
    fZ.append(event.fZ)

fX = np.array(fX)
fY = np.array(fY)
fZ = np.array(fZ)


def create_sinogram(fX, fY, num_bins, num_angles):
    max_radius = np.sqrt(np.max(fX**2 + fY**2))
    sinogram = np.zeros((num_angles, num_bins))
    
    for angle_idx in range(num_angles):
        theta = np.pi * angle_idx / num_angles
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        
        # Calculate the projections
        projection = fX * cos_theta + fY * sin_theta
        
        # Create histogram for this angle
        hist, bin_edges = np.histogram(projection, bins=num_bins, range=(-max_radius, max_radius))
        sinogram[angle_idx, :] = hist
    
    return sinogram


num_bins = 180
num_angles = 180
sinogram = create_sinogram(fX, fY, num_bins, num_angles)


def mlem_reconstruction(sinogram, num_iterations, num_bins, num_angles):
    reconstructed_image = np.ones((num_bins, num_bins))
    projection_bins = np.linspace(-num_bins // 2, num_bins // 2, num_bins)
    
    for iteration in range(num_iterations):
        back_projection = np.zeros_like(reconstructed_image)
        
        for angle_idx in range(num_angles):
            theta = 180 * angle_idx / num_angles
            projection = sinogram[angle_idx, :]
            
            # Forward projection
            forward_proj = np.sum(reconstructed_image, axis=1)
            
            # Avoid division by zero
            forward_proj[forward_proj == 0] = 1e-6
            
            # Compute correction ratio
            correction_ratio = projection / forward_proj
            
            # Expand correction_ratio to 2D
            correction_ratio_2d = np.tile(correction_ratio, (num_bins, 1))
            
            # Backward projection
            correction = rotate(correction_ratio_2d, theta, reshape=False)
            back_projection += correction
        
        reconstructed_image *= back_projection / num_angles
        
    return reconstructed_image

num_iterations = 20
reconstructed_image = mlem_reconstruction(sinogram, num_iterations, num_bins, num_angles)



# Display sinogram
# plt.figure(figsize=(10, 8))
plt.imshow(sinogram, cmap='gray', aspect='auto')
plt.title('Sinogram')
plt.xlabel('Bin')
plt.ylabel('Angle')
plt.colorbar()
plt.tight_layout()
plt.savefig('sinogram.pdf')
plt.savefig('sigogram.png')
plt.show()

# Display reconstructed image
# plt.figure(figsize=(10, 8))
plt.imshow(reconstructed_image, cmap='gray')
plt.title('Reconstructed Image')
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar()
plt.tight_layout()
plt.savefig('image.png')
plt.savefig('image.pdf')
plt.show()
