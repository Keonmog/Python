# Imported libs

import sys
import os
from tkinter import filedialog
import pydicom as dicom
import numpy as np
from skimage.transform import rescale, radon, iradon
from collections import defaultdict
import matplotlib.pyplot as plt

# Control variables ===================================

test_mode = 0         # 0 = test mode yes & 1 = all slices
zf = 1                      # zoom factor
dp = 0.1                  # display power factor for Fourier domain
sl_disp = 30            # slice number )test mode only) 

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
            ds = dicom.read_file(dicomFile)         # Tries to read DICOM files and tags
            instanceNumber = int(ds.InstanceNumber)
            echoTime = float(ds.EchoTime)
            repetitionTime = int(ds.RepetitionTime)
            dicomData = DICOMImageInfo(dicomFile, instanceNumber, echoTime, repetitionTime)
            echoTimeDict[dicomData.getEchoTime()].append(dicomData)
        except:
            continue                 # Skipping invalid files
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
for key in list(echoTimeDict.keys()):                                             # Sorts all lists
    echoTimeDict[key].sort(key=lambda tempData1: tempData1.getImageNumber(), reverse=False)
Ns = int(len(echoTimeDict[firstKey]))

if test_mode == 0:
    sl_one = sl_disp
    sl_two = sl_one+1

    print("test mode=yes")
else:
    sl_one = 0
    sl_two = Ns

    print("full run")
print(Ns, "slices")

# ==============================================================================


def resize_2d(A, zf, r, c):
    Azoom = np.uint16(rescale(A, zf, order=3, mode='edge', preserve_range=True))
    return Azoom


def flex_mask(R, r, c, pw):
    cmask = np.empty((r, c))
    for j in range(0, r):
        for k in range(0, c):
            if ((j-0.5*r)**pw+(k-0.5*c)**pw) <= R**pw:
                cmask[j, k] = 0
            else:
                cmask[j, k] = 1
    return np.uint16(cmask)


def line_mask(r, c):
    cmask = np.empty((r, c))
    for j in range(0, r):
        for k in range(0, c):
            if j == np.int(img1Row/2):
                cmask[j, k] = 1
            else:
                cmask[j, k] = 1
    return np.uint16(cmask)

# Main loop of the program=====================================================


for sl in range(sl_one, sl_two):
    img1 = dicom.read_file(echoTimeDict[firstKey][sl].getPath())
    data1 = resize_2d(img1.pixel_array, zf, img1.Rows, img1.Columns)

    if 'PixelData' in img1:
        img1Row = int(img1.Rows)
        img1Col = int(img1.Columns)
    else:
        print("No dimensional data in " + img1)
    
    image = data1
    IMAGE = np.fft.fftshift(np.fft.fft2(image))
    # theta = np.linspace(0., 180., max(image.shape), endpoint=False)
    # sinogram = radon(image, theta=theta, circle=True)
    # IMAGE = (np.fft.fft2(image))

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(12, 6), facecolor='white')


    image1 = ax1.imshow(image, cmap='gray')
    image2 = ax2.imshow(np.abs(IMAGE**dp), cmap='gray')
    image3 = ax3.imshow(np.real(IMAGE**dp), cmap='gray')
    image4 = ax4.imshow(np.imag(IMAGE**dp), cmap='gray')

    ax1.set_title('Image in anatomic space', color='blue')
    ax2.set_title('k-space: modulus', color='red')
    ax3.set_title('k-space: real', color='red')
    ax4.set_title('k-space: imaginary', color='red')

    # plt.colorbar(image1, ax=ax1)
    # plt.colorbar(image2, ax=ax2)
    fig.tight_layout()
    plt.show()




