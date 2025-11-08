# SVD Image Compression — EE1030 Matrix Theory (Dr. GVV Sharma)

Project for EE1030 (Matrix Theory) by Dr. GVV Sharma: a compact C implementation that demonstrates image compression via Singular Value Decomposition (SVD). The program reads an image, computes the SVD (Jacobi-based on AᵀA), and writes rank-k reconstructions and error metrics.

## Quick summary
- Language: C (pure C implementation)
- Goal: Demonstrate low-rank image approximation using SVD; compute/visualize error metrics
- Course: EE1030 Matrix Theory (assignment)

## Requirements

### macOS
- C compiler: gcc/clang (Xcode command line tools)
- math library: -lm
- Optional image libs: libjpeg, libpng
- Package manager: Homebrew (optional)
Install example:
```bash
brew install jpeg libpng
```

### Linux (Debian/Ubuntu)
- C compiler and build tools:
```bash
sudo apt update
sudo apt install build-essential pkg-config
```
- Image libraries:
```bash
sudo apt install libjpeg-dev libpng-dev
```
- Optional (if using cmake or other tooling):
```bash
sudo apt install cmake
```

### Linux (Fedora/RHEL)
```bash
sudo dnf install @development-tools pkgconfig
sudo dnf install libjpeg-turbo-devel libpng-devel
```

### Linux (Arch)
```bash
sudo pacman -S base-devel pkgconf libjpeg-turbo libpng
```

Notes:
- If you prefer single-file headers, you can use stb_image.h and stb_image_write.h (no external image libs needed).
- On Linux, package names may vary slightly by distribution (libjpeg-dev vs libjpeg-turbo-devel).

## Build
If Makefile exists:
```bash
make
```
Manual compile example (linking libjpeg/libpng if used):
```bash
gcc -o svd src/main.c src/jacobi_svd.c src/image_io.c -lm -ljpeg -lpng
```
If using stb headers, omit -ljpeg -lpng.

## Usage
Basic:
```bash
./svd <input_image> 
# Example: keep top 50 singular values
./svd examples/einstein.jpg 
```

Options
- input_image: JPEG/PNG path
- k: number of singular values to keep (1 .. min(width,height))
- -o out_dir: output directory (default: results/)

## Outputs
- Reconstructed images: <name>_k<N>.(png|jpg)
- Singular values: <name>_singular_vals.txt
- Error metrics: <name>_errors.txt
- LaTeX table example: tables/einstein_norm.tex


## How it works (concise)
1. Read image and convert to grayscale matrix A (m×n)
2. Form ATA = AᵀA (n×n)
3. Compute eigenpairs of ATA via Jacobi method → V, λ
4. Singular values σ = sqrt(λ); sort descending
5. Compute U columns from A·V / σ
6. Reconstruct A_k = U_k Σ_k V_kᵀ for requested k
7. Save reconstructed image(s) and compute metrics:
   - Absolute Frobenius error ||A − A_k||_F
   - Relative Frobenius error ||A − A_k||_F / ||A||_F
   - RMSE and PSNR

## Tips and limitations
- Work in grayscale to reduce memory/time.
- Jacobi SVD is educational and slower than LAPACK for large images.
- For large images, downsample before SVD or use smaller k.
- Ensure correct memory management; free allocations.

## Academic notes
- This repository contains work for EE1030 assignment by Dr. GVV Sharma.
- Include your name and roll number when submitting per course instructions.
- Follow course honor code when reusing or sharing code.

## Author
Divyansh Maheshwari 

EE25BTECH11037