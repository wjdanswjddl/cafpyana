# Covariance Analysis Script - Python Conversion

This is a Python conversion of the original C++ covariance analysis script for systematic uncertainty calculations in particle physics analysis.

## Key Changes from C++ to Python

### 1. **Modern Python Libraries**
- **numpy**: For efficient array operations and matrix calculations
- **matplotlib**: For plotting covariance matrices and histograms
- **uproot**: For reading ROOT files (replaces ROOT C++ libraries)
- **pathlib**: For cross-platform path handling

### 2. **Simplified Structure**
- Consolidated constants into a `Constants` class
- Streamlined function signatures with type hints
- Removed verbose ROOT-specific code while maintaining functionality
- Used Python's built-in features for cleaner code

Current `main` branch is based on python v3.9.15.
It is for running the repository without any issue at gpvm servers with the `spack`.
There is no need to open an SL7 image.
If python version is updated in gpvm servers, compatibility issue should be revisited.
