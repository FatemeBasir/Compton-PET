import os
import glob
import ROOT
import numpy as np
import matplotlib.pyplot as plt
from skimage.transform import radon, iradon, iradon_sart

# Step 1: Load the ROOT file and extract the data
#file = ROOT.TFile("../build/output0_t1.root")
#tree = file.Get("Hits")  # Replace "tree_name" with the actual name of the tree

print("Loading ROOT files...")

n_files = 10

os.chdir("../build/bin/")

tree_name = "Hits"
chain = ROOT.TChain(tree_name)

files = glob.glob("*.root")

for file in files[0:n_files]:
	chain.Add(file)

fX = []
fY = []
fZ = []

for entry in chain:
    fX.append(entry.fX)
    fY.append(entry.fY)
    fZ.append(entry.fZ)

fX = np.array(fX)
fY = np.array(fY)
fZ = np.array(fZ)

# Step 2: Preprocess the data (if needed)

# Step 3: Initialize the image grid
print("Initializing image grid...")

image_size = (256, 256)  # Adjust as per your data
image = np.ones(image_size)

# Create a circular mask
radius = image_size[0] // 2
Y, X = np.ogrid[:image_size[0], :image_size[1]]
dist_from_center = np.sqrt((X - radius) ** 2 + (Y - radius) ** 2)
mask = dist_from_center <= radius

# Apply the mask to the image
image[~mask] = 0

# Step 4: Compute the initial sinogram based on detector positions
print("Creating initial sinogram based on detector positions...")
num_angles = 180
theta = np.linspace(0., 180., num_angles, endpoint=False)
initial_sinogram = radon(image, theta=theta, circle=True)
sinogram = np.zeros_like(initial_sinogram)

# Normalize detector coordinates to the sinogram space
def normalize_to_sinogram_space(value, min_val, max_val, size):
    return int((value - min_val) / (max_val - min_val) * (size - 1))

# Populate the sinogram based on the detector positions
for (x, y, z) in zip(fX, fY, fZ):
    bin_x = normalize_to_sinogram_space(x, np.min(fX), np.max(fX), sinogram.shape[0])
    bin_y = normalize_to_sinogram_space(y, np.min(fY), np.max(fY), sinogram.shape[1])
    
    if 0 <= bin_x < sinogram.shape[0] and 0 <= bin_y < sinogram.shape[1]:
        sinogram[bin_x, bin_y] += 1

# Step 5: Implement the MLEM algorithm
print("Using MLEM algorithm...")
num_iterations = 50  # Number of MLEM iterations
epsilon = 1e-6  # Small value to prevent division by zero

for iteration in range(num_iterations):
    projection = radon(image, theta=theta, circle=True)
    ratio = sinogram / (projection + epsilon)
    correction = iradon_sart(ratio, theta=theta)
    image *= correction
    image /= np.sum(image) / np.sum(sinogram)
    print(f"Iteration {iteration+1}/{num_iterations} complete")

# Step 6: Visualize the reconstructed image
print("Plotting reconstructed image...")
plt.imshow(image, cmap='gray')
plt.title("Reconstructed Image using MLEM")
plt.xlabel("X position (pixels)")
plt.ylabel("Y position (pixels)")
plt.colorbar()
plt.savefig("image.png")
plt.savefig("image.pdf")
plt.show()
