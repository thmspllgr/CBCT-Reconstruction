This project was started as part of my second year of Master. It aims at reconstructing 3D volumes from 2D projections of a rotating sample obtained via x-rays imaging. It is an extension of the classical CT Filtered Back Projection (FBP) algorithm: instead of considering x-rays parallel to the detector (with a x-ray source at infinity), we consider a ponctual source that emits x-rays in all directions with a detector of finite size, thus the x-rays arrive following a cone shape. This is the principle of Cone Beam Computerized Tomography (CBCT).

## Structure

```
cbct-reconstruction
├── data
│   ├── 360                  # Folder containing 360 projections of a sample object (cylindric object with some internal and external structure)
│   ├── 15                   # Folder containing 15 projections of the same object
│   └── 4                    # Folder containing 4 projections of the same object
├── reconstruction.py        # Entry point
└── README.md                # Documentation
```

## Usage

1. Inside `reconstruction.py`, change the `datapath` variable to the path to the data you want to use.

2. Then run:
    ```
    python reconstruction.py
    ```

## Dependencies

- `numpy`
- `matplotlib`
- `tifffile`
- `scikit-image`
- `scipy`

## Notes

1. If you want to get a sample sinogram, just add inside the reconstruction loop:
        ```
        if x == 125:    # or any other slice index
            tifffile.imwrite("sinogram.tif", norm.astype(np.float32))
        ```

2. If you want to get a slice from the 3D reconstruction, you can also add inside the loop:
        ```
        if x == 125:
            tifffile.imwrite("slice.tif", reconstructed.astype(np.float32))
        ```

3. If you want to visualize the whole reconstructed 3D volume, you can use ```napari``` which supports tiff stacks. In your terminal just run:
        ```
        pip install napari[all]
        napari reconstruction.tif
        ```

4. If you want to compare our implementation with the `iradon` function from `skimage`, you can replace the filtering and backprojection calls with:
        ```
        reconstructed = iradon(density.T, theta=theta, circle=True)
        ```
And add `from skimage.transform import iradon` at the top of `reconstruction.py`. By default `iradon` uses the Ram-Lak filter so we can use it as is.