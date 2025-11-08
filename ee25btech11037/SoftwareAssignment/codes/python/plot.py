import sys
import numpy as np
import matplotlib.pyplot as plt

# --- 1. Get the base filename from the command line ---
if len(sys.argv) != 2:
    print("Usage: python plot_svd.py <base_filename>")
    print("Example: python plot_svd.py my_image")
    sys.exit(1)  # Exit with an error code

base_filename = sys.argv[1]

# --- 2. Define all filenames ---
singular_vals_file = f"{base_filename}_singular_vals.txt"
frob_norm_file = f"{base_filename}_frob_norm.txt"

spectrum_plot_file = f"{base_filename}_spectrum_plot.png"
abs_error_plot_file = f"{base_filename}_absolute_error_plot.png"
rel_error_plot_file = f"{base_filename}_relative_error_plot.png"


# --- 3. Load data with error checking ---
try:
    # Load singular values (Index, SingularValue)
    sv_data = np.loadtxt(
        singular_vals_file, 
        delimiter=None  # Use whitespace
    )

    # Load frobenius norm error (k, Absolute_Error, Relative_Error)
    error_data = np.loadtxt(
        frob_norm_file,
        delimiter=None # Use whitespace
        # MODIFICATION: Removed skiprows=1
    )

except FileNotFoundError as e:
    print(f"Error: Could not find a required data file.")
    print(f"Details: {e}")
    print(f"Please make sure '{singular_vals_file}' and '{frob_norm_file}' exist.")
    print("Hint: Run the C program first (the one that creates .txt files).")
    sys.exit(1)
except Exception as e:
    print(f"An error occurred while reading the files: {e}")
    # This can happen if frob_norm_file is empty or has wrong column count
    print("Hint: Check that frob_norm.txt has 3 columns and no header.")
    sys.exit(1)

# --- 4. Create Plot 1: Singular Value Spectrum ---
try:
    print(f"Generating {spectrum_plot_file}...")
    plt.figure(figsize=(10, 6))

    # sv_data[:, 0] is the index column
    # sv_data[:, 1] is the SingularValue column
    plt.plot(sv_data[:, 0], sv_data[:, 1], 'b-')

    plt.title(f'Singular Value Spectrum for "{base_filename}"', fontsize=16)
    plt.xlabel('k (Singular Value Index)', fontsize=12)
    plt.ylabel('Singular Value (Linear Scale)', fontsize=12)
    plt.grid(True, which="both", linestyle='--', linewidth=0.5)

    plt.savefig(spectrum_plot_file)
    plt.close()  # Close the plot to free memory

except Exception as e:
    print(f"Error creating spectrum plot: {e}")

# --- 5. Create Plot 2: Absolute Frobenius Norm Error ---
try:
    print(f"Generating {abs_error_plot_file}...")
    plt.figure(figsize=(10, 6))

    # error_data[:, 0] is the k column
    # error_data[:, 1] is the Absolute_Error column
    plt.plot(error_data[:, 0], error_data[:, 1], 'r-o', markersize=4)

    plt.title(f'Absolute Frobenius Norm Error vs. k for "{base_filename}"', fontsize=16)
    plt.xlabel('k (Number of Singular Values Used)', fontsize=12)
    plt.ylabel('Absolute Frobenius Norm Error', fontsize=12)
    
    plt.grid(True, linestyle='--', linewidth=0.5)

    plt.savefig(abs_error_plot_file)
    plt.close()

except Exception as e:
    print(f"Error creating absolute error plot: {e}")

# --- 6. Create Plot 3: Relative Frobenius Norm Error ---
try:
    print(f"Generating {rel_error_plot_file}...")
    plt.figure(figsize=(10, 6))

    # error_data[:, 0] is the k column
    # error_data[:, 2] is the Relative_Error column
    plt.plot(error_data[:, 0], error_data[:, 2], 'g-o', markersize=4)

    # Format Y-axis as percentage
    from matplotlib.ticker import PercentFormatter
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1.0))

    plt.title(f'Relative Frobenius Norm Error vs. k for "{base_filename}"', fontsize=16)
    plt.xlabel('k (Number of Singular Values Used)', fontsize=12)
    plt.ylabel('Relative Frobenius Norm Error (Loss)', fontsize=12)
    
    plt.grid(True, linestyle='--', linewidth=0.5)

    plt.savefig(rel_error_plot_file)
    plt.close()

except Exception as e:
    print(f"Error creating relative error plot: {e}")

print("\nDone. All plotting complete.")