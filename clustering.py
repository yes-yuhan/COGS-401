#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 10:34:48 2024

@author: liuyuh
"""

import numpy as np 
from sklearn.cluster import KMeans  
import nibabel as nib 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os
import subprocess


"""
Set the working directory
"""
os.chdir("/data/chamal/projects/yuhan/projects/2023-07-thalamus-cortex-covariance")


"""
Set the variables
"""
sex = "female"
atlas = "AHRA_annotation_full"
# Available age groups include: 10_to_13, 12_to_15, 14_to_17, 16_to_19
age_group = "10_to_13"
# Available measures include: cortex_area_left, cortex_gyrification, cortex_jd,
# cortex_mc_mid_left, cortex_thickness_left, cortex_volume_left, thalamus_jd
measure = "cortex_area_left"

"""
Define the helper functions
"""
def map_tha_vector_to_mask(tha_mask, correlation_matrix): 
    mask = tha_mask 
    tha_vector_to_mask = {} 
    row = 0 
    # For each voxel in the mask
    for i in range(mask.shape[0]): 
        for j in range(mask.shape[1]): 
            for k in range(mask.shape[2]): 
                # If this voxel is part of the mask
                if mask[i,j,k] > 0.5: 
                     tha_vector_to_mask[row]=(i,j,k) 
                     row += 1 
    return tha_vector_to_mask 

# This function maps every column of the row vector to the cortical mask
# It takes a cortical mask file and a correlation matrix as parameters
# and returns a dictionary of cortical coordinates corresponding to each row in the correlation matrix
def map_cor_vector_to_mask(cor_mask_data, correlation_matrix): 
    mask = cor_mask 
    cor_vector_to_mask = {} 
    col = 0 
    for i in range(mask.shape[0]): 
            for j in range(mask.shape[1]): 
                    for k in range(mask.shape[2]): 
                            if mask[i,j,k] > 0.5: 
                                     cor_vector_to_mask[col]=(i,j,k) 
                                     col += 1 
    return cor_vector_to_mask
# This function creates a new cortical file with the correlation values corresponding to
# the designated row number in the correlation matrix
def create_new_cortical_mask(mask_path, correlation_row, output_path):
    # Load the cortical mask
    mask = nib.load(mask_path)
    mask_data = mask.get_fdata()

    # Initialize a counter for the columns in the correlation matrix
    col = 0

    # For each voxel in the mask
    for i in range(mask_data.shape[0]):
        for j in range(mask_data.shape[1]):
            for k in range(mask_data.shape[2]):
                # If this voxel is part of the mask
                if mask_data[i, j, k] > 0.5:
                    # Set the value of the voxel to the correlation value
                    mask_data[i, j, k] = correlation_row[col]

                    # Increment the column counter
                    col += 1

    # Save the new cortical mask
    new_mask_img = nib.Nifti1Image(mask_data, affine=mask.affine)
    nib.save(new_mask_img, output_path)


"""
Set the paths
"""
path_cormat = "derivatives/matrices/"+atlas+"/population_correlations/"+sex+"/"+age_group+"/thalamus_jd_X_"+measure+".npy"
path_tha_mask = "resources/yohan/derivatives/PNC/masks/AHRA_annotation_full/thalamus_mask_highres.mnc"
path_cor_mask = "resources/yohan/derivatives/PNC/masks/AHRA_annotation_full/cortex_mask_lowres.mnc"
path_outdir = "analyses/clustering/k_means/"+sex+"/"+age_group+"/thalamus_jd_X_"+measure+"/"
os.makedirs(path_outdir, exist_ok=True)


"""
Load the correlation matrice and the mask files
"""
cormat = np.load(path_cormat)
# Get the data of the thalamus mask
tha_mask = nib.load(path_tha_mask).get_fdata()
cor_mask = nib.load(path_cor_mask).get_fdata()


"""
Perform K-means clustering
"""
kmeans = KMeans(n_clusters=6, random_state=0, n_init=10, max_iter=300, algorithm="elkan").fit(cormat)

clusters = np.zeros(tha_mask.shape, dtype=int)


"""
Map the vector indix back to the thalamic mask
"""
tha_vector_to_mask = map_tha_vector_to_mask(tha_mask, cormat)


"""
Assign voxels to clusters
"""
for i in range(6):
    voxel_indices = np.where(kmeans.labels_ == i)[0]
    mask_indices = [tha_vector_to_mask[v] for v in voxel_indices]
    for index in mask_indices:
        clusters[index] = i + 1
         

"""
Save clusters to NIfTI files
"""
for i in range(1,7):
    print(f"Processing cluster {i} ...")
    cluster_mask = (clusters == i).astype(np.float32)
    nii_filename = os.path.join(path_outdir, f'cluster_{i}.nii')
    
    img = nib.Nifti1Image(cluster_mask, np.eye(4))
    
    nib.save(img, nii_filename)
    print(f"NIfTI file created: {nii_filename}")


"""
Visualisation as a sanity check
"""
# Initialize the figure for 3D plotting
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Define colors
colors = ['r', 'g', 'b', 'y', 'c', 'm']

# Iterate over each cluster label
for i in range(6):
    # Get the voxel indices for the current cluster
    voxel_indices = np.where(clusters == i+1)
    x, y, z = voxel_indices
    ax.scatter(x,y,z, color=colors[i], label=f'Cluster {i+1}', s=70)

ax.set_xlabel('X axis')               
ax.set_ylabel('Y axis')    
ax.set_zlabel('Z axis')

# Show the legend

plt.show()

# Save the plot to a file
plot_path = os.path.join(path_outdir, f'thalamic_clusters')
plt.savefig(plot_path)
"""

    
cluster_mask = (clusters == i)
x, y, z = np.where(cluster_mask)
ax.scatter(x,y,z, c=colors[i-1], label=f'Cluster {i}', s=1)


# Map the voxel indices to mask indices
mask_indices = [tha_vector_to_mask[v] for v in voxel_indices if v in tha_vector_to_mask]

# Check if there are any indices to plot
if not mask_indices:
    print(f"No indices for cluster {i+1}.")
    continue

# Extract the coordinates for plotting
x_coords, y_coords, z_coords = zip(*mask_indices)

# Plot the cluster points
ax.scatter(x_coords, y_coords, z_coords, c=colors[i-1], label=f'Cluster {i+1}', s=1)  # Added 's=1' to reduce point size

"""

"""
(Aborted) Save clusters to MINC files 

for i in range(1,7):
    cluster_mask = (clusters == i).astype(np.float32)
    npy_temp_filename = os.path.join(path_outdir, f'cluster_{i}.npy')
    raw_temp_filename = os.path.join(path_outdir, f'cluster_{i}.raw')
    mnc_filename = os.path.join(path_outdir, f'cluster_{i}.mnc')
    np.save(npy_temp_filename, cluster_mask)
    cluster_mask_loaded = np.load(npy_temp_filename).astype(np.float32)
    cluster_mask_loaded.tofile(raw_temp_filename)
    if not os.path.exists(raw_temp_filename):
        print(f"Failed to create raw binary file: {raw_temp_filename}")
        continue
    try:
        subprocess.run(['mincconvert', raw_temp_filename, mnc_filename], check=True)
        print(f"MINC file created: {mnc_filename}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while converting {raw_temp_filename}: {e}")
    os.remove(npy_temp_filename)
    os.remove(raw_temp_filename)
""" 
