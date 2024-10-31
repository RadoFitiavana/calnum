
# Calnum: C++ Numerical Calculus and Signal Processing Library

**Calnum** is a collection of C++ implementations for a range of numerical calculus and signal processing techniques, focused on applications in numerical linear algebra, partial differential equations, and discrete transforms. This library includes both fundamental algorithms and innovative approaches, such as optimized storage and advanced matrix manipulation, designed for efficiency in large-scale computations.

## Features

### 1. Numerical Linear Algebra
- Efficient algorithms for matrix operations, vector manipulation, and linear system solutions.
- Optimized storage techniques for handling large, sparse matrices using a modified Cuthill-McKee ordering.
- Specialized algorithms for solving large and sparse systems, reducing memory requirements and computational cost.

### 2. Solving Partial Differential Equations (PDEs)
- Implementation of second-order PDE solvers in 2D using numerical linear algebra.
- Utilizes the discrete sine transform (DST) to solve boundary value problems efficiently.

### 3. Discrete Sine Transform (DST) and Signal Processing
- Functions for computing the discrete sine transform, useful in both signal processing and PDE solvers.
- Algorithms designed to enhance performance in signal filtering and transformation tasks.

## Key Algorithms and Methods

- **Modified Cuthill-McKee (MCM) Ordering**: An optimization to standard MCM for sparse matrix storage, significantly improving memory usage in large-scale linear systems.
- **Sparse Matrix Solvers**: Efficient routines for solving large, sparse systems arising in numerical simulations.
- **Numerical Methods for PDEs**: Leveraging discrete transforms and numerical linear algebra to provide stable and accurate solutions for 2D second-order PDEs.

