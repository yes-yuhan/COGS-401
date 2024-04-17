'''
These helper functions are created to aid the conversion between the numpy array for the correlation matrix, the mask index, and the world coordinates.
'''

import nibabel as nib 
import numpy as np
import pandas as pd


### This function merges a thalamic mask file with a cortical mask file
def merge_masks(tha_mask_path, cor_mask_path, output_path): 

    # Load the target thalamic and cortical mask files
    tha_mask = nib.load(tha_mask_path).get_fdata() 
    cor_mask = nib.load(cor_mask_path).get_fdata() 

    # Create a new nii file with the masks merged
    merged_mask = tha_mask + cor_mask 
    merged_mask_img = nib.Nifti1Image(merged_mask, affine=nib.load(tha_mask_path).affine)

    # Save the merged mask
    nib.save(merged_mask_img, output_path)


### This function maps a set of MNI coordinates (x, y, z) that points to a voxel in the thalamic mask to a row of the correlation matrix
def create_mask_to_voxel_mapping(mask_path, correlation_matrix):
    # Load the mask
    mask = nib.load(mask_path).get_fdata()

    # Initialize the mapping as an empty dictionary
    mask_to_voxel = {}

    # Initialize a counter for the rows in the correlation matrix
    row = 0

    # For each voxel in the mask
    for i in range(mask.shape[0]):
        for j in range(mask.shape[1]):
            for k in range(mask.shape[2]):
                # If this voxel is part of the mask
                if mask[i, j, k] > 0.5:
                    # Add the mapping from the 3D index to the row
                    mask_to_voxel[(i, j, k)] = row

                    # Increment the row counter
                    row += 1

    return mask_to_voxel


### This function converts the MNI coordinates to voxel indices
def mni_to_voxel(mni_coord, affine):
    # Convert the MNI coordinates to voxel indices
    voxel_coord = np.round(np.linalg.inv(affine).dot(mni_coord + [1]))[:3].astype(int)
    # Returns the voxel indices
    return voxel_coord


### This function finds the row number of the correlation matrix for a given set of MNI coordinates
### If the set of MNI coordinates does not point to a voxel on the thalamic mask, this function will return None
def find_row_for_mni(mask_path, correlation_matrix, mni_coord):
    # Load the mask
    mask = nib.load(mask_path)
    mask_data = mask.get_fdata()

    # Convert the MNI coordinate to voxel indices
    voxel_coord = mni_to_voxel(mni_coord, mask.affine)

    # Initialize a counter for the rows in the correlation matrix
    row = 0

    # For each voxel in the mask
    for i in range(mask_data.shape[0]):
        for j in range(mask_data.shape[1]):
            for k in range(mask_data.shape[2]):
                # If this voxel is part of the mask
                if mask_data[i, j, k] > 0.5:
                    # If this voxel is the one we're looking for
                    if np.array_equal([i, j, k], voxel_coord):
                        return row

                    # Increment the row counter
                    row += 1

    # If the voxel was not found in the mask
    return None


### This function creates an excel file with the MNI coordinates for all the voxels in the mask
def create_mni_coordinate_excel(mask_path, output_path):
    # Load the mask
    mask = nib.load(mask_path)
    mask_data = mask.get_fdata()

    # Initialize a list to store the MNI coordinates
    mni_coords = []

    # For each voxel in the mask
    for i in range(mask_data.shape[0]):
        for j in range(mask_data.shape[1]):
            for k in range(mask_data.shape[2]):
                # If this voxel is part of the mask
                if mask_data[i, j, k] > 0.5:
                    # Convert the voxel indices to MNI coordinates
                    mni_coord = mask.affine.dot(np.array([i, j, k, 1]))[:3]

                    # Add the MNI coordinate to the list
                    mni_coords.append(mni_coord)

    # Create a DataFrame from the list of MNI coordinates
    df = pd.DataFrame(mni_coords, columns=['X', 'Y', 'Z'])

    # Write the DataFrame to an Excel file
    df.to_excel(output_path, index=False)


### This function creates a new thalamic file with a single voxel highlighted in an extreme value
def create_new_thalamic_mask(mask_path, voxel_coord, output_path):
    # Load the thalamic mask
    mask = nib.load(mask_path)
    mask_data = mask.get_fdata()

    # Change the value of the voxel of interest to an extreme value
    mask_data[voxel_coord[0], voxel_coord[1], voxel_coord[2]] = np.max(mask_data) + 1

    # Save the new thalamic mask
    new_mask_img = nib.Nifti1Image(mask_data, affine=mask.affine)
    nib.save(new_mask_img, output_path)


### This function creates a new cortical file with the correlation values corresponding to
### the designated row number in the correlation matrix
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

### This function maps every column of the row vector to the cortical mask
### It takes a cortical mask file and a correlation matrix as parameters
### and returns a dictionary of cortical coordinates corresponding to each row in the correlation matrix
def map_cor_vector_to_mask(cor_mask, correlation_matrix): 
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

### This function maps rows of the correlation matrix to the indices in a thalamic mask
### It takes a thalamic mask file and a correlation matrix as parameters
### and returns a dictionary of thalamic coordinates corresponding to each row in the correlation matrix
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
