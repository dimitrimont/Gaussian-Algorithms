# Linear Algebra Methods


## Overview
This Python program provides custom implementations of key linear algebra methods used for solving systems of equations and matrix operations:

Gauss-Jordan Elimination
Matrix Inversion using Gauss-Jordan method
Gaussian Elimination with back-substitution
Matrix Determinant calculation

## Key Features

Interactive Interface: Command-line menu for selecting different calculation methods
Partial Pivoting: Implementation includes pivot selection for numerical stability
Row Operation Tracking: Counts row swaps for accurate determinant calculation
Custom Implementations: All algorithms are implemented from first principles without relying on library functions

Requirements

Python 3.x
NumPy

Usage

Run the program:
Copypython Gaussian_Algorithms_MontgomeryDimitri.py

Follow the interactive menu prompts to select a calculation method:

Option 1: Gauss-Jordan Elimination
Option 2: Gauss-Jordan Inversion
Option 3: Gaussian Elimination
Option 4: Matrix Determinant
Option 5: Exit



Implementation Details
Gauss-Jordan Elimination

Transforms a matrix to reduced row echelon form
Uses partial pivoting to select the largest magnitude pivot
Returns the solution to a system of linear equations

Matrix Inversion

Appends an identity matrix to the original matrix
Applies Gauss-Jordan elimination to transform the left side to identity
The right side becomes the inverse matrix

Gaussian Elimination

Performs forward elimination to transform the matrix to upper triangular form
Uses back-substitution to solve for variables
Includes partial pivoting for stability

Determinant Calculation

Based on the product of diagonal elements after Gaussian elimination
Accounts for sign changes due to row swaps

Author
Dimitri Montgomery
