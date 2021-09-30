

# Purpose: to generate a volumetric display of the voxel shape in 2D-MRI and 3D-MRI
# This program uses a formula for the voxel sensitivity function that will be discussed in class
# First objective is that you type the codem aking sure the indentations of each line are correct
# This program requires a specific folder structure where the result tiff images will be stored, specifically
# Voxel Sensitivity with the two subfolders; 1_Results 2D and 1_Results 3D


# ===> Required Python libraries ========================================

import sys
import os
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
from skimage.io import imsave

# ==> Global variables ==========================================

Nx = 256  # ==> matrix size in x-direction
Ny = 256  # ==> matrix size in y-direction
Nz = 256  # ==> matrix size in z-direction (slice direction)

delta_xy = 20  # ==> voxel size in transverse plane
SLTH = 20  # ==> voxel size in slice encoding direction

temp = np.zeros([Nx, Ny, Nz], 'int16')  # initialize the arrays
image = np.zeros_like(temp)

scan_mode = "3D"  # Main control variable

# ==> Queries user for working folders ================================

SubjectFolder = filedialog.askdirectory(title='Select the Subject Folder', initialdir=os.getcwd())

# Checks that the Subject folders exist
if SubjectFolder:
    print("'" + SubjectFolder + "set as the source folder")
else:
    sys.exit("An invalid path was selected. Exiting...")
    
Folder1 = "1_Results_2D"
Folder2 = "1_Results_3D"

# Checks that the sub folders exist
if not os.path.exists(os.path.join(SubjectFolder, Folder1)):
    sys.exit("The selected Subject folder does not contain a " + Folder1 + " sub folder")
if not os.path.exists(os.path.join(SubjectFolder, Folder2)):
    sys.exit("The selected Subject folder does not contain a " + Folder2 + " sub folder")
os.chdir(SubjectFolder)
# ==> Functions =========================================================


def pv_scale_2d(A, Min, Max, r, c):
    A1 = np.empty_like(A)
    Amax = A.max()
    Amin = A.min()
    for j in range(0, r-1):
        for k in range(0, c-1):
            if np.abs(Amax-Amin) == 0:
                A1[j, k] = 0
            else:
                A1[j, k] = (Max-Min)*((A[j, k] - Amin)/(Amax-Amin)) + Min
    return A1


def q(x, n):
    if x == 0:
        f = 1
    else:
        f = (1/n)*np.sin(np.pi*x/delta_xy)/(np.sin(np.pi*x/(n*delta_xy)))
    return f


def slice_profile(x, sig):
    return np.exp(-(x/sig)**2)


def vsf_2D(x, y, z):
    return q(x, Nx)*q(y, Ny)*slice_profile(z, SLTH/2)


def vsf_3D(x, y, z):
    return q(x, Nx)*q(y, Ny)*q(z, Nz)


def VSF(nx, ny, nz):
    array = np.empty_like(temp)
    for j in range(0, nx):
            for k in range(0, ny):
                for sl in range(0, nz):
                    if scan_mode == "2D":
                        array[j, k, sl] = (2**8)*vsf_2D(j-nx/2, k-ny/2, sl-nz/2)
                    if scan_mode == "3D":
                        array[j, k, sl] = (2**8)*vsf_3D(j-nx/2, k-ny/2, sl-nz/2)
            print("outer loop index ==> "+str(j))
    return np.abs(array)

# ===> Main Program =============================================


image = VSF(Nx, Ny, Nz)

for sl1 in range(0, Nz-1):
    if scan_mode == "2D":
        imsave("1_Results_2D/2DVSF"+str(sl1) + ".tif", np.uint16(image[:, :, sl1]))
        print("Saving VSF =>", sl1)
    if scan_mode == "3D":
        imsave("1_Results_3D/3DVSF"+str(sl1) + ".tif", np.uint16(image[:, :, sl1]))
        print("Saving VSF =>", sl1)
