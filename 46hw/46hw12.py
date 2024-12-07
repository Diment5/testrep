import numpy as np

def gaussian_elimination_no_pivoting(A, b):
    """
    Gaussian elimination without pivoting, using 4-digit floating-point arithmetic.
    """
    n = len(b)
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)

    # Forward elimination
    for i in range(n):
        for j in range(i + 1, n):
            if A[i, i] == 0:
                raise ValueError("Zero pivot encountered!")
            m = round(A[j, i] / A[i, i], 4)  # Compute multiplier with rounding
            A[j, i:] = np.round(A[j, i:] - m * A[i, i:], 4)  # Eliminate with rounding
            b[j] = round(b[j] - m * b[i], 4)

    # Back substitution
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = round((b[i] - np.dot(A[i, i + 1:], x[i + 1:])) / A[i, i], 4)

    return x

def gaussian_elimination_with_pivoting(A, b):
    """
    Gaussian elimination with partial pivoting, using 4-digit floating-point arithmetic.
    """
    n = len(b)
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)

    # Forward elimination with partial pivoting
    for i in range(n):
        # Pivoting: find the row with the largest absolute value in column i
        max_row = np.argmax(np.abs(A[i:, i])) + i
        if A[max_row, i] == 0:
            raise ValueError("Zero pivot encountered!")

        # Swap rows in A and b
        if max_row != i:
            A[[i, max_row]] = A[[max_row, i]]
            b[[i, max_row]] = b[[max_row, i]]

        # Elimination
        for j in range(i + 1, n):
            m = round(A[j, i] / A[i, i], 4)  # Compute multiplier with rounding
            A[j, i:] = np.round(A[j, i:] - m * A[i, i:], 4)  # Eliminate with rounding
            b[j] = round(b[j] - m * b[i], 4)

    # Back substitution
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = round((b[i] - np.dot(A[i, i + 1:], x[i + 1:])) / A[i, i], 4)

    return x

# Define the system
A = [
    [6, 2, 2],
    [2, 2/3, 1/3],
    [1, 2, -1]
]
b = [-2, 1, 0]

# Solve using Gaussian elimination without pivoting
x_no_pivot = gaussian_elimination_no_pivoting(A, b)
print("Solution without pivoting:", x_no_pivot)

# Solve using Gaussian elimination with partial pivoting
x_with_pivot = gaussian_elimination_with_pivoting(A, b)
print("Solution with partial pivoting:", x_with_pivot)



def power_method(A, tol=1e-8, max_iter=1000):
    """
    Power Method to find the dominant eigenvalue and eigenvector.
    """
    n = A.shape[0]
    x = np.random.rand(n)  # Start with a random vector
    x = x / np.linalg.norm(x)
    
    for i in range(max_iter):
        x_new = np.dot(A, x)
        x_new = x_new / np.linalg.norm(x_new)
        
        # Check convergence
        if np.linalg.norm(x_new - x) < tol:
            break
        x = x_new
    
    eigenvalue = np.dot(x.T, np.dot(A, x))
    return eigenvalue, x, i + 1

def hilbert_matrix(n):
    """
    Generate the Hilbert matrix of size n x n.
    """
    return np.array([[1 / (i + j - 1) for j in range(1, n + 1)] for i in range(1, n + 1)])

# Part (a): Apply the power method to Hilbert matrices of size 4:4:20
results_a = []
for n in range(4, 21, 4):
    H = hilbert_matrix(n)
    eigenvalue, eigenvector, iterations = power_method(H)
    results_a.append((n, eigenvalue, iterations))


for n in range(4, 21, 4):
    H = hilbert_matrix(n)
    eigenvalue, eigenvector, iterations = power_method(H)
    print(f"n = {n}: Dominant Eigenvalue = {eigenvalue:.6f}, Iterations = {iterations}")
def modified_power_method_smallest(A, tol=1e-8, max_iter=1000):
    """
    Power Method modified to find the smallest eigenvalue using inverse iteration.
    """
    n = A.shape[0]
    x = np.random.rand(n)
    x = x / np.linalg.norm(x)

    A_inv = np.linalg.inv(A)  # Compute the inverse of A

    for i in range(max_iter):
        x_new = np.dot(A_inv, x)
        x_new = x_new / np.linalg.norm(x_new)

        # Check convergence
        if np.linalg.norm(x_new - x) < tol:
            break
        x = x_new

    smallest_eigenvalue = 1 / np.dot(x.T, np.dot(A, x))
    return smallest_eigenvalue, x, i + 1


# Part (b): Find the smallest eigenvalue for n = 16
n = 16
H = hilbert_matrix(n)
smallest_eigenvalue, eigenvector, iterations = modified_power_method_smallest(H)

print(f"Part (b):\nSmallest Eigenvalue = {smallest_eigenvalue:.6f}, Iterations = {iterations}")

# Part (c): Error estimation consistency
# Perturbation E
E = np.random.rand(n, n) * 1e-3  # Small perturbation matrix
A_perturbed = H + E
smallest_eigenvalue_perturbed, _, _ = modified_power_method_smallest(A_perturbed)

# Theoretical bound using Bauer-Fiske theorem
P = np.linalg.eig(H)[1]  # Eigenvectors of H
P_norm = np.linalg.norm(P, ord=2)
P_inv_norm = np.linalg.norm(np.linalg.inv(P), ord=2)
E_norm = np.linalg.norm(E, ord=2)

error_bound = P_norm * P_inv_norm * E_norm
actual_error = abs(smallest_eigenvalue_perturbed - smallest_eigenvalue)

print("\nPart (c):")
print(f"Actual Error = {actual_error:.6e}")
print(f"Theoretical Error Bound = {error_bound:.6e}")
print(f"Consistency: {'Yes' if actual_error <= error_bound else 'No'}")
