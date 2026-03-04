import glob
import os
import re
import tifffile
import numpy as np
from skimage.io import imread
from scipy.fft import fft, ifft, fftfreq

### Parameters
D = 30.87 # measured source-object distance (cm)
d = 14.9 # measured object-detector distance (cm)
pixel_size = 12.7 / 343 # cm   (we measured on our images 343px -> 12.7cm, so 1px = 12.7/343 cm)
datapath = "C:/Users\\prete\\Documents\\M2\\Practical courses\\CT\\data\\Projections\\360" # path to the folder that contains the projections

### Functions
def get_sorted_files(directory, ext="*.png"):
    """Load and sort files in natural order (img1, img2, ..., img100) to suit our projection data."""
    files = glob.glob(os.path.join(directory, ext))
    return sorted(files, key=lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split('([0-9]+)', s)])


def weighted_img(img, d, D, px_size):
    """Compute the weighted image according to the geometry given by d, D, and px_size."""
    # Create centered 2D grid
    h, w = img.shape
    u_px = np.arange(w) - w / 2
    v_px = np.arange(h) - h / 2
    U_px, V_px = np.meshgrid(u_px, v_px)

    # Pixel to cm conversion
    U_cm = U_px * px_size
    V_cm = V_px * px_size
    
    # u -> s coordinates conversion
    factor = 1 + (d / D)
    S_k = U_cm / factor
    T_k = V_cm / factor
    
    # Calculate weighting for each pixel in the image
    distance = (S_k / D)**2 + (T_k / D)**2
    weight = 1.0 / np.sqrt(1 + distance)

    return img * weight


def filtering(sinogram, pixel_size):
    """Apply Ram-Lak filter to the sinogram in frequency domain."""
    n_angles, n_detectors = sinogram.shape
    freq = fftfreq(n_detectors, d=pixel_size)
    ram_lak_filter = np.abs(2 * freq)
    sinogram_fft = fft(sinogram, axis=1)
    filtered_fft = sinogram_fft * ram_lak_filter
    sino_filt = np.real(ifft(filtered_fft, axis=1)).astype(np.float32)

    return sino_filt[:, :n_detectors].astype(np.float32)


def backprojection(filtered_sinogram, theta, D, h, px_size):
    """Perform back-projection to reconstruct the image from filtered sinogram."""
    n_angles, n_detectors = filtered_sinogram.shape
    image = np.zeros((h, h), dtype=np.float32)
    
    # Image grid in cm that covers the same area as the detector
    grid_cm = (np.arange(h) - h/2) * px_size
    Y, Z = np.meshgrid(grid_cm, grid_cm)
    
    theta_rad = np.deg2rad(theta)
    det_center_idx = n_detectors / 2
    
    for i, angle in enumerate(theta_rad):
        cos_a, sin_a = np.cos(angle), np.sin(angle)
        
        # Coordinates rotation
        t = Y * cos_a - Z * sin_a # transverse axis (parallel to the detector)
        p = Y * sin_a + Z * cos_a # depth axis (towards source)
        
        # Compute magnification factor U
        dist_src = D - p # distance from source to point (t,p) along optical axis
        dist_src[dist_src < 0.01] = 0.01
        U = dist_src / D
        
        # Projection onto detector
        s_proj_cm = t / U
        s_proj_idx = (s_proj_cm / px_size) + det_center_idx
        
        # Interpolation of the sinogram values
        profile = filtered_sinogram[i, :]
        values = np.interp(s_proj_idx, np.arange(n_detectors), profile, left=0, right=0)
        
        # Weigh on back-projection
        bp_weigh = 1 / (U**2)
        image += values * bp_weigh

    return image * (np.pi / (2 * n_angles))


### Main function
def reconstruction():
    # Load the projections and the first image for geometry
    imgs = get_sorted_files(datapath)
    n_angles = len(imgs)
    img0 = imread(imgs[0], as_gray=True).astype(np.float32)
    h, w = img0.shape

    # Calculate the mean intensity on left border (20 pixels) inside the first image
    # We use it to normalize the sinograms later (otherwise we don't see anything)
    I0 = np.mean(img0[:, :20], axis=1)
    I0[I0 <= 0] = np.max(I0)

    # Compute weighted images and store them in a table
    stack = np.zeros((n_angles, h, w), dtype=np.float32)
    for idx, img in enumerate(imgs):
        C_beta_n = imread(img, as_gray=True).astype(np.float32)
        stack[idx] = weighted_img(C_beta_n, d, D, pixel_size)

    # Initialize output volume for reconstruction
    volume = np.zeros((w, h, h), dtype=np.float32)
    theta = np.linspace(0., 360., n_angles, endpoint=False)
    
    for x in range(w):
        if x % 50 == 0:
            print(f"Slice {x}/{w}")

        # The sinogram is the stack of projections at position x
        sinogram = stack[:, :, x]

        # Conversion into density using Beer-Lambert law
        norm = sinogram / I0[np.newaxis, :]
        norm[norm <= 0] = 1e-6
        density = -np.log(norm)
        density[density < 0] = 0
        
        # Filtering and back-projection
        filtered = filtering(density, pixel_size)
        reconstructed = backprojection(filtered, theta, D, h, pixel_size)
        
        volume[x, :, :] = reconstructed

    tifffile.imwrite("reconstruction.tif", volume.astype(np.float32))
    print("Reconstruction saved as reconstruction.tif")
    return


reconstruction()