import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import scipy.linalg as scila
import time

def driver():
    ''' Create matrix for testing different ways of solving a square linear system '''

    sizes = [100, 500, 1000, 2000, 4000, 5000]
    right_hand_sides_threshold = 0
    
    for N in sizes:
        print(f"\nN = {N}")

        # Right hand side
        b = np.random.rand(N, 1)
        A = np.random.rand(N, N)

        # Measure time for standard solver
        start_time = time.time()
        x = scila.solve(A, b)
        end_time = time.time()
        standard_solve_time = end_time - start_time
        
        r = la.norm(np.matmul(A, x) - b)
        print(f"Residual (standard): {r}")
        print(f"Time (standard): {standard_solve_time:.4f} s")

        # Measure time for LU factorization and LU solver
        start_time = time.time()
        P, L, U = scila.lu(A)
        lu_factorization_time = time.time() - start_time
        
        start_time = time.time()
        y = scila.solve_triangular(L, np.matmul(P, b), lower=True)
        x_lu = scila.solve_triangular(U, y)
        lu_solve_time = time.time() - start_time
        
        r_lu = la.norm(np.matmul(A, x_lu) - b)
        print(f"Residual (LU): {r_lu}")
        print(f"LU factor time: {lu_factorization_time:.4f} s")
        print(f"LU solve time: {lu_solve_time:.4f} s")

        total_lu_time = lu_factorization_time + lu_solve_time
        if total_lu_time < standard_solve_time:
            right_hand_sides_threshold += 1

    print(f"\nLU faster in {right_hand_sides_threshold} cases")
    interpret_results(right_hand_sides_threshold, len(sizes))


def interpret_results(threshold_count, total_cases):
    ''' Interpret the results of LU solver performance compared to the standard solver '''
    print("\nResults:")
    if threshold_count > total_cases / 2:
        print("LU was faster in most cases, better for multiple systems with the same matrix.")
    else:
        print("Standard solver faster for most cases, better for single right-hand side problems.")


def create_rect(N, M):
    ''' This subroutine creates an ill-conditioned rectangular matrix '''
    a = np.linspace(1, 10, M)
    d = 10 ** (-a)

    D2 = np.zeros((N, M))
    for j in range(0, M):
        D2[j, j] = d[j]

    # Create matrices needed to manufacture the low rank matrix
    A = np.random.rand(N, N)
    Q1, R = la.qr(A)
    A = np.random.rand(M, M)
    Q2, R = la.qr(A)
    
    B = np.matmul(Q1, D2)
    B = np.matmul(B, Q2)
    return B

if __name__ == "__main__":
    driver()
