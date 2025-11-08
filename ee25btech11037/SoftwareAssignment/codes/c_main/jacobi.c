// THIS IS FOR FULL C IMPLEMENTATION
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <jpeglib.h>
#include <png.h>
// #include "jacobi.h"
#include "/Users/divy/Desktop/Coding/software_proj/codes/c_libs/jacobi.h"
// this function is correct and working properly
// allocates memory for a matrix in a heap
double **allocate_matrix(int m, int n)
{
    double **mat = malloc(m * sizeof(double *));
    mat[0] = malloc(m * n * sizeof(double));
    for (int i = 0; i < m; i++)
    {
        mat[i] = mat[0] + i * n;
    }
    return mat;
}

// this function is correct and working properly
// free a allocated matrix
void free_matrix(double **mat)
{
    free(mat[0]);
    free(mat);
}

//  this function is correct and working properly
// Transpose a matrix a
void transpose_mat(int m, int n, double **a, double **at)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            at[j][i] = a[i][j];
        }
    }
    return;
}

// this function is correct and working properly
// Multiply two mats at and a with a[m][n]
void mult_at_a(int m, int n, double **res, double **a, double **at)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            res[i][j] = 0;
            for (int k = 0; k < m; k++)
            {
                res[i][j] += at[i][k] * a[k][j];
            }
        }
    }
}

// implementing jacobi for singualr values and V matrix
//  working fine
void jacobi_eigen(int n, double **a, double **V, double *singular_vals)
{
    // initializing V as an identity matrix
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            V[i][j] = (i == j) ? 1.00 : 0.00;
        }
    }

    // changing the matrix a till all the non-diagonal elems get close to or zero
    int comp = 1;
    while (comp)
    {
        comp = 0;

        // getting only non-diagonal elems
        for (int i = 0; i < n - 1; i++) // iterating on a row
        {
            for (int j = i + 1; j < n; j++) // iterating over each column
            {
                double Aij = a[i][j];
                if (fabs(Aij) > 1e-9)
                {
                    double Aii = a[i][i], Ajj = a[j][j];
                    double theta = 0.5 * atan2(2 * Aij, Ajj - Aii);
                    double c = cos(theta), s = sin(theta);
                    // applying to rows
                    for (int k = 0; k < n; k++)
                    {
                        double Aki = c * a[i][k] - s * a[j][k];
                        double Akj = s * a[i][k] + c * a[j][k];
                        a[i][k] = Aki;
                        a[j][k] = Akj;
                    }

                    // apllying to columns
                    for (int k = 0; k < n; k++)
                    {
                        double Aki = a[k][i];
                        double Akj = a[k][j];
                        a[k][i] = c * Aki - s * Akj;
                        a[k][j] = s * Aki + c * Akj;
                    }

                    // updating V
                    for (int k = 0; k < n; k++)
                    {
                        double Vki = V[k][i];
                        double Vkj = V[k][j];
                        V[k][i] = c * Vki - s * Vkj;
                        V[k][j] = s * Vki + c * Vkj;
                    }
                    comp = 1;
                }
            }
        }
    }
    // storing the singular vals
    for (int i = 0; i < n; i++)
    {
        singular_vals[i] = sqrt(fabs(a[i][i]));
    }
    return;
}

// sorting SVD result to descending order
void sort_svd(int n, double *S, double **V)
{
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            if (S[j] > S[i])
            {
                // Swap singular values
                double temp = S[i];
                S[i] = S[j];
                S[j] = temp;

                // Swap corresponding columns in V
                for (int k = 0; k < n; k++)
                {
                    temp = V[k][i];
                    V[k][i] = V[k][j];
                    V[k][j] = temp;
                }
            }
        }
    }
}

// computing U=1/(simga[i])*A*v_i
void compute_U(int m, int n, double **a, double **U, double *S, double **V)
{
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            U[i][j] = 0.0;
    for (int i = 0; i < n; i++)
    {
        if (S[i] < 1e-9)
            continue;

        for (int row = 0; row < m; row++)
        {
            U[row][i] = 0.0;
            for (int col = 0; col < n; col++)
                U[row][i] += a[row][col] * V[col][i];
            U[row][i] /= S[i];
        }
    }
}

// reconstructing A_k matrix for different values of k
void reconstruct_matrix(int m, int n, int k, double **U, double *S,
                        double **V, double **A_reconstructed)
{
    // Initialize result to zero
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A_reconstructed[i][j] = 0.0;
        }
    }

    // Use only the top k singular values
    int rank = (k < n) ? k : n;

    // A = sum from i=0 to k-1 of: sigma_i * u_i * v_i^T
    for (int i = 0; i < rank; i++)
    {
        for (int row = 0; row < m; row++)
        {
            for (int col = 0; col < n; col++)
            {
                A_reconstructed[row][col] += S[i] * U[row][i] * V[col][i];
            }
        }
    }
}

// reading the file
// working perfectly
unsigned char *read_jpeg(const char *filename, int *m, int *n)
{
    // setup
    struct jpeg_decompress_struct cinfo; // control structure for image handling which holds all the decompression state
    struct jpeg_error_mgr jerr;          // for error handling
    FILE *fp;                            // standard file pointer
    JSAMPROW row_pointer[1];             // jsamprow is an unsigned char* // array that points to a row that reads one line at a time

    unsigned char *buffer;             // input image array to read the data using this file
    cinfo.err = jpeg_std_error(&jerr); // connect the jerr error handler to cinfo

    fp = fopen(filename, "rb"); // opening the file
    if (fp == NULL)             // if the file does not open then returning NULL
    {
        return NULL;
    }

    jpeg_create_decompress(&cinfo); // intializing the decopression object
    jpeg_stdio_src(&cinfo, fp);     // telling libjpeg to get data from fp

    (void)jpeg_read_header(&cinfo, TRUE); // read the file header to get image data

    cinfo.out_color_space = JCS_GRAYSCALE; // forcing the image to be greyscale
    cinfo.output_components = 1;           // ensuring only one component

    (void)jpeg_start_decompress(&cinfo); // Start the decompression. After this, output dimensions are final.

    *m = cinfo.output_height; // getting the final dimensions of the inputted image
    *n = cinfo.output_width;

    int r = *n;
    buffer = (unsigned char *)malloc((*m) * r * sizeof(unsigned char));
    if (buffer == NULL) // if buffer doesnt get allocated data then returning NULL
    {
        jpeg_destroy_decompress(&cinfo);
        fclose(fp);
        return NULL;
    }

    while (cinfo.output_scanline < *m) // reading each line in the image using the rowpointer
    {
        row_pointer[0] = &buffer[cinfo.output_scanline * r]; // setting the pointer to first element of the current readig a line
        (void)jpeg_read_scanlines(&cinfo, row_pointer, 1);   // telling lib to read a line and write it directly to address in row_pointer
    }

    (void)jpeg_finish_decompress(&cinfo); // finish decompression
    jpeg_destroy_decompress(&cinfo);      // destroy the decompressor
    fclose(fp);                           // closing the file
    return buffer;
}

void write_jpeg(const char *filename, int m, int n, double **A_k, int quality)
{
    struct jpeg_compress_struct cinfo; // control structure for image handling which holds all the decompression state
    struct jpeg_error_mgr jerr;        // for error handling
    FILE *of;                          // pointer to a file
    JSAMPROW row_pointer[1];           // row pointer
    unsigned char *buffer;             // preparing to write file

    buffer = (unsigned char *)malloc(m * n * sizeof(unsigned char)); // allocating memory for image
    if (buffer == NULL)
    {
        return;
    }

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            double val = A_k[i][j] * 255.0;

            // Clamp to [0, 255] range
            if (val > 255.0)
                val = 255.0;
            if (val < 0.0)
                val = 0.0;

            buffer[i * n + j] = (unsigned char)(val + 0.5); // +0.5 for rounding
        }
    }

    cinfo.err = jpeg_std_error(&jerr); // error handler
    jpeg_create_compress(&cinfo);      // initializing the compression object

    of = fopen(filename, "wb");
    if (of == NULL)
    {
        free(buffer);
        jpeg_destroy_compress(&cinfo);
        return;
    }

    jpeg_stdio_dest(&cinfo, of); // Tell libjpeg to write to our opened file.

    // setting image parameters
    cinfo.image_height = m;
    cinfo.image_width = n;
    cinfo.input_components = 1;
    cinfo.in_color_space = JCS_GRAYSCALE;

    // Set default compression parameters
    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, quality, TRUE); // setting output quality

    jpeg_start_compress(&cinfo, TRUE);               // start the compressor
    while (cinfo.next_scanline < cinfo.image_height) // Loop while there are still scanlines to write
    {
        row_pointer[0] = &buffer[cinfo.next_scanline * n];  // address of row where data is to be written
        (void)jpeg_write_scanlines(&cinfo, row_pointer, 1); // writing one line to the output file
    }

    jpeg_finish_compress(&cinfo); // finish the decompressor

    fclose(of);                    // closing file
    jpeg_destroy_compress(&cinfo); // detroying the compressor
    free(buffer);                  // freeing buffer
}

unsigned char *read_png(const char *filename, int *m, int *n)
{
    FILE *fp = fopen(filename, "rb");
    if (!fp)
        return NULL;

    png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png)
    {
        fclose(fp);
        return NULL;
    }

    png_infop info = png_create_info_struct(png);
    if (!info)
    {
        png_destroy_read_struct(&png, NULL, NULL);
        fclose(fp);
        return NULL;
    }

    if (setjmp(png_jmpbuf(png)))
    {
        png_destroy_read_struct(&png, &info, NULL);
        fclose(fp);
        return NULL;
    }

    png_init_io(png, fp);
    png_read_info(png, info);

    *m = png_get_image_height(png, info);
    *n = png_get_image_width(png, info);
    png_byte color_type = png_get_color_type(png, info);
    png_byte bit_depth = png_get_bit_depth(png, info);

    // Convert to grayscale
    if (color_type == PNG_COLOR_TYPE_PALETTE)
        png_set_palette_to_rgb(png);

    if (color_type == PNG_COLOR_TYPE_RGB || color_type == PNG_COLOR_TYPE_RGB_ALPHA)
        png_set_rgb_to_gray(png, 1, 0.299, 0.587);

    if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
        png_set_expand_gray_1_2_4_to_8(png);

    if (png_get_valid(png, info, PNG_INFO_tRNS))
        png_set_tRNS_to_alpha(png);

    if (bit_depth == 16)
        png_set_strip_16(png);

    if (color_type == PNG_COLOR_TYPE_GRAY_ALPHA || color_type == PNG_COLOR_TYPE_RGB_ALPHA)
        png_set_strip_alpha(png);

    png_read_update_info(png, info);

    unsigned char *buffer = (unsigned char *)malloc((*m) * (*n) * sizeof(unsigned char));
    if (!buffer)
    {
        png_destroy_read_struct(&png, &info, NULL);
        fclose(fp);
        return NULL;
    }

    png_bytep *row_pointers = (png_bytep *)malloc(sizeof(png_bytep) * (*m));
    for (int y = 0; y < *m; y++)
        row_pointers[y] = buffer + y * (*n);

    png_read_image(png, row_pointers);

    free(row_pointers);
    png_destroy_read_struct(&png, &info, NULL);
    fclose(fp);

    return buffer;
}

void write_png(const char *filename, int m, int n, double **A_k)
{
    FILE *fp = fopen(filename, "wb");
    if (!fp)
        return;

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png)
    {
        fclose(fp);
        return;
    }

    png_infop info = png_create_info_struct(png);
    if (!info)
    {
        png_destroy_write_struct(&png, NULL);
        fclose(fp);
        return;
    }

    if (setjmp(png_jmpbuf(png)))
    {
        png_destroy_write_struct(&png, &info);
        fclose(fp);
        return;
    }

    png_init_io(png, fp);

    png_set_IHDR(png, info, n, m, 8, PNG_COLOR_TYPE_GRAY,
                 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT);

    png_write_info(png, info);

    unsigned char *buffer = (unsigned char *)malloc(m * n * sizeof(unsigned char));
    if (!buffer)
    {
        png_destroy_write_struct(&png, &info);
        fclose(fp);
        return;
    }

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            double val = A_k[i][j] * 255.0;
            if (val > 255.0)
                val = 255.0;
            if (val < 0.0)
                val = 0.0;
            buffer[i * n + j] = (unsigned char)(val + 0.5);
        }
    }

    png_bytep *row_pointers = (png_bytep *)malloc(sizeof(png_bytep) * m);
    for (int y = 0; y < m; y++)
        row_pointers[y] = buffer + y * n;

    png_write_image(png, row_pointers);
    png_write_end(png, NULL);

    free(row_pointers);
    free(buffer);
    png_destroy_write_struct(&png, &info);
    fclose(fp);
}

// Froebenius norm of A-A_k matrix
double frob_norm(double *s, int k, int n)
{
    double error = 0;
    for (int i = k; i < n; i++)
    {
        error += s[i] * s[i];
    }
    return sqrt(error) * 255.00;
}

// checking for png or jpg file
int is_png(const char *filename)
{
    int len = strlen(filename);
    if (len < 4)
        return 0;
    if ((strcmp(filename + len - 4, ".png") == 0) || (strcmp(filename + len - 4, ".PNG") == 0))
    {
        return 1;
    }
    return 0;
}

void read_name(const char *input, char *name)
{
    int len = strlen(input), i;
    for (i = 0; i < len; i++)
    {
        if (input[i] != '.')
        {
            name[i] = input[i];
        }
        else
        {
            name[i] = '\0';
            return;
        }
    }
    name[i] = '\0';
    return;
}

void save_singular_vals(const char *filename, double *s, int n)
{
    char base_name[200];
    read_name(filename, base_name);

    char txt_name[250];
    snprintf(txt_name, sizeof(txt_name), "%s_singular_vals.txt", base_name);

    FILE *fp = fopen(txt_name, "w");
    for (int i = 0; i < n; i++)
    {
        fprintf(fp, "%-5d %.4lf\n", i + 1, s[i] * 255.0);
    }
    fclose(fp);
}
