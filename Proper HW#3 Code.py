

import sys
import os
from tkinter import filedialog
import pydicom as dicom
import numpy as np
from skimage.transform import radon, iradon
from collections import defaultdict
import matplotlib.pyplot as plt

# Control variables ===================================

test_mode = 0         # 0 = test mode yes & 1 = all slices
zf = 1                      # zoom factor
dp = 0.25                  # display power factor for Fourier domain
sl_disp = 40            # slice number )test mode only) 

# ==================================================


class DICOMImageInfo:
    def __init__(self, filePath, imageNumber, echoTime, repetitionTime):
        self.filePath = filePath
        self.imageNumber = imageNumber
        self.echoTime = echoTime
        self.repetitionTime = repetitionTime

    def getPath(self):
        return self.filePath

    def getImageNumber(self):
        return self.imageNumber

    def getEchoTime(self):
        return self.echoTime

    def getRepetitionTime(self):
        return self.repetitionTime

#   Open DICOM files in folder


def getDICOMFromFolder(folder):
    echoTimeDict = defaultdict(list)
    for filename in os.listdir(folder):
        dicomFile = os.path.join(folder, filename)
        try:
            ds = dicom.read_file(dicomFile)                         # Tries to read DICOM files and tags
            instanceNumber = int(ds.InstanceNumber)
            echoTime = float(ds.EchoTime)
            repetitionTime = int(ds.RepetitionTime)
            dicomData = DICOMImageInfo(dicomFile, instanceNumber, echoTime, repetitionTime)
            echoTimeDict[dicomData.getEchoTime()].append(dicomData)
        except:
            continue                                                # Skipping invalid files
            print("WARNING: Mismatched Image Skipped. Check Folders")
            # Displays warning if skipped - will update to be more precise
    return echoTimeDict

# ======================================================
# ================= 1. Ask for Subject =================
# ======================================================
#   Queries user for Subject folder


SubjectFolder = filedialog.askdirectory(title='Select the Subject Folder', initialdir=os.getcwd())
os.chdir(SubjectFolder)
#   Checks that the Subject folder exists
if SubjectFolder:
    print("'" + SubjectFolder + "'set as the source folder")
else:
    sys.exit("An invalid path was selected. Exiting...")

Folder1 = "1_DAdcm"

#   Checks that the sub folders exist
if not os.path.exists(os.path.join(SubjectFolder, Folder1)):
    sys.exit("The selected Subject folder does not contain a " + Folder1 + " sub folder")

# ==Check dicom folder 1a for reading====

echoTimeDict = getDICOMFromFolder(Folder1)

# gets dict of list of DICOM objects
'''if len(echoTimeDict) != 1:          # error check number of unique Echo Times
    print("The folder must have exactly one TE in the DICOM files")
    sys.exit("This folder has " + str(len(echoTimeDict)))'''
firstKey = list(echoTimeDict.keys())[0]
lenFirstList1 = len(echoTimeDict[firstKey])
for key in list(echoTimeDict.keys()):             # Sorts all lists
    echoTimeDict[key].sort(key=lambda tempData1: tempData1.getImageNumber(), reverse=False)
Ns = int(len(echoTimeDict[firstKey]))

if test_mode < 1:
    sl_one = sl_disp
    sl_two = sl_one+1

    print("test mode=yes")
else:
    sl_one = 0
    sl_two = Ns

    print("full run")
print(Ns, "slices")


# Main loop of the program=====================================================


for sl in range(sl_one, sl_two):
    img1 = dicom.read_file(echoTimeDict[firstKey][sl].getPath())
    data1 = img1.pixel_array

    if 'PixelData' in img1:
        img1Row = int(img1.Rows)
        img1Col = int(img1.Columns)
    else:
        print("No dimensional data in " + img1)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), facecolor='white')

    image = data1
    theta = np.linspace(0., 90., max(image.shape), endpoint=False)
    sinogram = radon(image, theta=theta, circle=True)


    image1 = ax1.imshow(image, cmap='gray')
    image2 = ax2.imshow(np.real(sinogram), cmap='gray')

    ax1.set_title('Image in anatomic space', color='blue')
    ax2.set_title('Radon transform or sinogram', color='red')

    plt.show()

    reconstruction_fbp = np.abs(iradon(sinogram, theta=tom_angle, circle=True))
    error = reconstruction_fbp - image
    print(f"FBP rms reconstruction error: {np.sqrt(np.mean(error**2)):.3g}")

    imkwargs = dict(vmin=-0.2, vmax=0.2)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharex=True, sharey=True)
    ax1.set_title("Reconstruction\nFiltered back projection")
    ax1.imshow(reconstruction_fbp, cmap=plt.cm.Greys_r)
    ax2.set_title("Reconstruction error\nFiltered back projection")
    ax2.imshow(reconstruction_fbp - image, cmap=plt.cm.Greys_r, **imkwargs)
    plt.show()
