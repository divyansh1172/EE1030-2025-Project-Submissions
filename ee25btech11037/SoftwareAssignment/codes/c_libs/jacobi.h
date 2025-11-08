#ifndef JACOBI_H
#define JACOBI_H

// Added stdio.h for FILE* type used in is_png
#include <stdio.h>

double **allocate_matrix(int m, int n);
void free_matrix(double **mat);
void transpose_mat(int m, int n, double **a, double **at);
void mult_at_a(int m, int n, double **res, double **a, double **at);
void jacobi_eigen(int n, double **a, double **V, double *singular_vals);
void sort_svd(int n, double *S, double **V);
void compute_U(int m, int n, double **a, double **U, double *S, double **V);
void reconstruct_matrix(int m, int n, int k, double **U, double *S, double **V, double **A_reconstructed);

unsigned char *read_jpeg(const char *filename, int *m, int *n);
void write_jpeg(const char *filename, int m, int n, double **A_k, int quality);
unsigned char *read_png(const char *filename, int *m, int *n);
void write_png(const char *filename, int m, int n, double **A_k);
double frob_norm(double *s, int k, int n);
int is_png(const char *filename);
void read_name(const char *input, char *name);
void save_singular_vals(const char *filename, double *s, int n);
#endif