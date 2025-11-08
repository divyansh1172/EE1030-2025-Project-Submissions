# Image Compression using SVD
EE1030 Matrix Theory (Dr. GVV Sharma)  
By: Divyansh Maheshwari (EE25BTECH11037)

## Overview
This project implements image compression using Singular Value Decomposition (SVD) in C. It converts images to grayscale, performs SVD using the Jacobi method, and creates compressed versions with different ranks.

## Features
- Pure C implementation
- Grayscale image processing
- SVD-based compression
- Error metrics calculation (Frobenius norm)


## Technical Details
The compression works by:
1. Converting image to matrix A
2. Computing SVD: A = UΣVᵀ
3. Creating rank-k approximations
4. Measuring compression quality

## Performance Tips
- Use smaller images for faster processing
- Lower rank k gives higher compression but lower quality
- Program works best with grayscale images

## References
- Gilbert Strang's MIT Courseware on SVD