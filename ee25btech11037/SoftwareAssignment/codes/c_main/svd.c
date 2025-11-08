#include <stdlib.h>
#include <stdio.h>
// #include "jacobi.h"
#include "/Users/divy/Desktop/Coding/software_proj/codes/c_libs/jacobi.h"

int main(int arg, char *argv[])
{
    if (arg != 2)
    {
        fprintf(stderr, "Usage: %s <input.jpg>\n", argv[0]);
        return 1;
    }

    const char *input_file = argv[1];

    int m, n;
    unsigned char *image_data = NULL;
    if (is_png(input_file))
    {
        image_data = read_png(input_file, &m, &n);
    }
    else
    {
        image_data = read_jpeg(input_file, &m, &n);
    }

    if (image_data == NULL)
    {
        fprintf(stderr, "Error reading input file\n");
        return 1;
    }

    //  allocating memory for all the matrices
    double **a = allocate_matrix(m, n);   // a matrix
    double **at = allocate_matrix(n, m);  // transpose of a
    double **ata = allocate_matrix(n, n); // matrix prod of ata
    double **V = allocate_matrix(n, n);   // eigenvector matrix
    double **U = allocate_matrix(m, n);   // left multiplication matrix
    double *singular_vals = malloc(n * sizeof(double));
    double **A_reconstructed = allocate_matrix(m, n);

    // copying image data into a
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            a[i][j] = (double)image_data[i * n + j] / 255.0;
        }
    }
    free(image_data);

    // running all the functions and doing the svd decomposition (ONCE)
    transpose_mat(m, n, a, at);              // taking transpose of a
    mult_at_a(m, n, ata, a, at);             // multiplying at and a
    jacobi_eigen(n, ata, V, singular_vals);  // computing singular values and V
    sort_svd(n, singular_vals, V);           // sorting descending
    compute_U(m, n, a, U, singular_vals, V); // computing U from the formula SVD decomposition

    save_singular_vals(input_file, singular_vals, n);
    // Define the k values you want to save
    int k_values[] = {1, 5, 10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800,850, 900, 950, 1000, n}; // Add or modify as needed
    int num_k_values = sizeof(k_values) / sizeof(k_values[0]);

    char name[200];
    read_name(input_file, name);

    char error_filename[200];
    // char rel_error_filename[200];
    snprintf(error_filename, sizeof(error_filename), "%s_frob_norm.txt", name);
    // snprintf(rel_error_filename, sizeof(rel_error_filename), "%s_rel_frob_norm.txt", name);

    // double A_frob = frob(a, m, n);
    double A_frob = frob_norm(singular_vals, 0, n);
    FILE *err_fp = fopen(error_filename, "w");
    // FILE *rel_err_fp = fopen(rel_error_filename, "w");
    // Generate and save images for each k value
    for (int i = 0; i < num_k_values; i++)
    {
        int k = k_values[i];

        // Skip if k exceeds the number of singular values
        if (k > n)
            continue;

        // Reconstruct the matrix with k components
        reconstruct_matrix(m, n, k, U, singular_vals, V, A_reconstructed);

        // Create output filename
        char output_file[256];
        if (is_png(input_file))
        {
            snprintf(output_file, sizeof(output_file), "%s_k%d.png", name, k);
            // Save the reconstructed image
            write_png(output_file, m, n, A_reconstructed);
        }
        else
        {
            snprintf(output_file, sizeof(output_file), "%s_k%d.jpg", name, k);
            // Save the reconstructed image
            write_jpeg(output_file, m, n, A_reconstructed, 90);
        }

        double err = frob_norm(singular_vals, k, n);
        fprintf(err_fp, "%-5d %.6lf %.6lf\n", k, err, err / A_frob);
        // fprintf(rel_err_fp, "%-5d %.4lf\n", k, err / A_frob);

        printf("Saved: %s (k=%d) ; Frobenius error is %.4lf\n", output_file, k, err);
    }

    // Clean up
    free_matrix(a);
    free_matrix(at);
    free_matrix(ata);
    free_matrix(V);
    free_matrix(U);
    free(singular_vals);
    free_matrix(A_reconstructed);
    // fclose(rel_err_fp);
    fclose(err_fp);
    return 0;
}