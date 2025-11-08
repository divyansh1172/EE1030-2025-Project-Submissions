# SVD Image Compression

A C implementation of image compression using Singular Value Decomposition (SVD). This project decomposes grayscale images using SVD and allows reconstruction with varying levels of compression.

## Features

- Supports both JPEG and PNG input images
- Converts color images to grayscale automatically
- Computes SVD decomposition using Jacobi eigenvalue algorithm
- Generates compressed images with different k-rank approximations
- Produces error analysis and visualization plots
- Calculates Frobenius norm errors for quality assessment

## Prerequisites

- C compiler (gcc/clang)
- LibJPEG library
- LibPNG library 
- Python 3.x with NumPy and Matplotlib (for plotting)

## Project Structure
```
├── codes/
│   ├── c_libs/
│   │   └── jacobi.h           # Header file with function declarations
│   ├── c_main/
│   │   ├── jacobi.c           # Implementation of SVD and image processing
│   │   └── svd.c              # Main program
│   └── python/
│       └── plot.py            # Plotting script for error analysis
├── figs/                      # Output directory for plots
└── tables/                    # LaTeX tables with error metrics
```

## Installation

### Installing Dependencies

**macOS (Homebrew):**
```bash
brew install jpeg libpng
```

**Ubuntu/Debian:**
```bash
sudo apt-get install libjpeg-dev libpng-dev
```

**Python packages:**
```bash
pip install numpy matplotlib
```

### Building the Project
```bash
# Navigate to the c_main directory
cd codes/c_main

# Compile the program
gcc svd.c jacobi.c -o svd -I/opt/homebrew/include -L/opt/homebrew/lib -ljpeg -lpng -lm -O3
```

**Note:** Adjust the include and library paths based on your system configuration. For Linux systems, you may use:
```bash
gcc svd.c jacobi.c -o svd -ljpeg -lpng -lm -O3
```

## Usage

### Running SVD Compression
```bash
./svd input_image.jpg
# or
./svd input_image.png
```

### Generating Plots
```bash
python codes/python/plot.py input_image
```

## Output Files

The program generates several output files:

- `*_k{N}.jpg/png` - Compressed images with k singular values
- `*_singular_vals.txt` - List of computed singular values
- `*_frob_norm.txt` - Frobenius norm error data
- `*_spectrum_plot.png` - Visualization of singular value spectrum
- `*_absolute_error_plot.png` - Absolute error vs k graph
- `*_relative_error_plot.png` - Relative error vs k graph

## Algorithm Overview

The compression process follows these steps:

1. Read input image and convert to grayscale
2. Compute SVD decomposition: **A = UΣV^T**
   - Use Jacobi eigenvalue algorithm for **A^T A**
   - Sort singular values in descending order
   - Compute corresponding **U** matrix
3. Generate k-rank approximations for various k values
4. Calculate error metrics for each approximation
5. Save compressed images and analysis data

## Error Metrics

The program calculates two types of errors:

- **Absolute Frobenius norm error:** ||A - A_k||_F
- **Relative Frobenius norm error:** ||A - A_k||_F / ||A||_F

where A is the original image matrix and A_k is the rank-k approximation.

## Example
```bash
# Compress an image
./svd photo.jpg

# Generate analysis plots
python codes/python/plot.py photo
```

This will produce compressed versions of the image with different compression ratios and corresponding error analysis plots.


## Author

Divyansh Maheshwari 

EE25BTECH11037

## Acknowledgments

- SVD implementation based on the Jacobi eigenvalue algorithm
- Image processing using LibJPEG and LibPNG libraries