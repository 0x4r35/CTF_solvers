#!/usr/bin/env python3
"""
Solve script for challenge (n1)^3 over GF(257).

Given:
- output.txt contains two printed objects from Sage:
  1) The vector of homogeneous cubic polynomials f(x) = ((x U)^{∘3}) T evaluated at symbolic x = (x0,...,x_{n-1}).
  2) The ciphertext vector y = f(flag_digits) for the unknown plaintext digits over base 257.

Approach (Jennrich's algorithm in GF(257)):
- Parse the polynomial tuple into per-coordinate lists of monomials with coefficients.
- Build slice matrices A = M^(k1)(I,I,r) and B = M^(k2)(I,I,r) using a random r in GF(257)^n, where
  M^(k) is the symmetric 3-tensor with entries M_{abc} such that f_k(x) = sum_{a,b,c} M_{abc} x_a x_b x_c.
- Diagonalize C = A * B^{-1} over GF(257) by scanning eigenvalues in the field and computing eigenvectors via nullspaces.
- Columns of U_hat are the eigenvectors (up to scale/permutation).
- Generate n random inputs x_i, evaluate Y_i = f(x_i), compute S_hat = ((x_i U_hat)^{∘3}), then solve S_hat T_hat = Y to obtain T_hat.
- Decrypt y: s = y T_hat^{-1}; take component-wise cube root via exponent 171 (3*171 ≡ 1 mod 256); then x = s^{1/3} U_hat^{-1}.
- Reconstruct bytes from base-257 digits and print the flag.
"""

import re
import sys
import random
from typing import List, Tuple

P = 257

def mod(x: int) -> int:
    return x % P

INV3 = pow(3, -1, P)
INV6 = pow(6, -1, P)


class PolyData:
    """
    Holds parsed data for one polynomial f_k, with monomials separated by pattern and
    precomputed M-coefficients per term (M_{abc}).
    """

    __slots__ = (
        'n',
        'terms_iii',   # List[Tuple[int, int, int]]: (i, coeff_poly, m)
        'terms_iij',   # List[Tuple[int, int, int, int]]: (i, j, coeff_poly, m)
        'terms_ijk',   # List[Tuple[int, int, int, int, int]]: (i, j, k, coeff_poly, m)
    )

    def __init__(self, n: int):
        self.n = n
        self.terms_iii: List[Tuple[int, int, int]] = []
        self.terms_iij: List[Tuple[int, int, int, int]] = []
        self.terms_ijk: List[Tuple[int, int, int, int, int]] = []

    def add_term(self, coeff_poly: int, vars_exp: List[Tuple[int, int]]):
        """
        Add a monomial term with coefficient from polynomial representation.
        vars_exp: list of (var_index, exponent). Sum of exponents must be 3.
        """
        coeff_poly = mod(coeff_poly)
        # Expand into multi-set of indices with multiplicity equal to exponent
        idxs: List[int] = []
        for v, e in vars_exp:
            if e <= 0:
                continue
            idxs.extend([v] * e)
        if len(idxs) != 3:
            raise ValueError(f"Monomial degree !=3: {vars_exp}")
        i, j, k = idxs[0], idxs[1], idxs[2]
        # Classify by multiplicity pattern
        if i == j == k:
            # x_i^3 term: M[i,i,i] = coeff
            m = coeff_poly
            self.terms_iii.append((i, coeff_poly, m))
        else:
            # Sort to detect pattern, but preserve actual structure for classification
            # For iij: exactly two equal
            if i == j != k:
                irep, other = i, k
            elif i == k != j:
                irep, other = i, j
            elif j == k != i:
                irep, other = j, i
            else:
                irep = other = None  # all distinct
            if irep is not None:
                # x_i^2 * x_j
                m = mod(coeff_poly * INV3)
                self.terms_iij.append((irep, other, coeff_poly, m))
            else:
                # all distinct
                a, b, c = sorted((i, j, k))
                m = mod(coeff_poly * INV6)
                self.terms_ijk.append((a, b, c, coeff_poly, m))

    def slice_matrix(self, r: List[int]) -> List[List[int]]:
        """
        Compute A = M^(k)(I, I, r) as an n x n symmetric matrix over GF(257).
        r is length-n vector over GF(257).
        """
        n = self.n
        A = [[0] * n for _ in range(n)]
        # x_i^3: A[i,i] += m * r[i]
        for i, _c, m in self.terms_iii:
            A[i][i] = mod(A[i][i] + m * r[i])
        # x_i^2 * x_j: contributions
        for i, j, _c, m in self.terms_iij:
            # (a,b,c) permutations: (i,i,j), (i,j,i), (j,i,i)
            A[i][i] = mod(A[i][i] + m * r[j])
            tmp = mod(m * r[i])
            A[i][j] = mod(A[i][j] + tmp)
            A[j][i] = mod(A[j][i] + tmp)
        # x_i * x_j * x_k (i<j<k)
        for i, j, k, _c, m in self.terms_ijk:
            # Contribute to pairs (i,j),(i,k),(j,k)
            tmp = mod(m * r[k])
            A[i][j] = mod(A[i][j] + tmp)
            A[j][i] = mod(A[j][i] + tmp)
            tmp = mod(m * r[j])
            A[i][k] = mod(A[i][k] + tmp)
            A[k][i] = mod(A[k][i] + tmp)
            tmp = mod(m * r[i])
            A[j][k] = mod(A[j][k] + tmp)
            A[k][j] = mod(A[k][j] + tmp)
        return A

    def eval_at(self, x: List[int]) -> int:
        """
        Evaluate this polynomial at vector x over GF(257).
        """
        n = self.n
        if len(x) != n:
            raise ValueError("x length mismatch")
        res = 0
        # x_i^3
        for i, c, _m in self.terms_iii:
            res = mod(res + c * pow(x[i], 3, P))
        # x_i^2 * x_j
        for i, j, c, _m in self.terms_iij:
            res = mod(res + c * mod(pow(x[i], 2, P) * x[j]))
        # x_i * x_j * x_k
        for i, j, k, c, _m in self.terms_ijk:
            res = mod(res + c * mod(x[i] * x[j] * x[k]))
        return res


def parse_output_file(path: str) -> Tuple[List[PolyData], List[int]]:
    """
    Parse output.txt to obtain list of PolyData for each coordinate and the ciphertext vector y.
    """
    txt = open(path, 'r', encoding='utf-8', errors='ignore').read()

    # Extract two top-level parenthesis blocks
    tuples = []
    i = 0
    while i < len(txt):
        if txt[i] == '(': 
            depth = 1
            j = i + 1
            while j < len(txt) and depth > 0:
                if txt[j] == '(':
                    depth += 1
                elif txt[j] == ')':
                    depth -= 1
                j += 1
            if depth != 0:
                raise ValueError("Unbalanced parentheses in output.txt")
            tuples.append(txt[i:j])
            if len(tuples) == 2:
                break
            i = j
        else:
            i += 1
    if len(tuples) < 2:
        raise ValueError("Could not find two tuples in output.txt")

    poly_tuple_str = tuples[0]
    y_tuple_str = tuples[1]

    # Split tuple items at top-level commas (depth 1)
    def split_tuple_items(s: str) -> List[str]:
        assert s[0] == '(' and s[-1] == ')'
        content = s[1:-1]
        items = []
        cur = []
        depth = 0
        for ch in content:
            if ch == '(':
                depth += 1
            elif ch == ')':
                depth -= 1
            if ch == ',' and depth == 0:
                items.append(''.join(cur))
                cur = []
            else:
                cur.append(ch)
        if cur:
            items.append(''.join(cur))
        return [it.strip() for it in items if it.strip()]

    poly_items = split_tuple_items(poly_tuple_str)
    n = len(poly_items)
    # Parse y vector entries
    y_items = split_tuple_items(y_tuple_str)
    y = [mod(int(it)) for it in y_items]
    if len(y) != n:
        # Fallback: sometimes Sage prints vector of GF elements differently; try splitting by commas directly
        raise ValueError(f"Ciphertext length {len(y)} doesn't match poly dimension {n}")

    # Parse each polynomial string into PolyData
    polys: List[PolyData] = []
    for k, pstr in enumerate(poly_items):
        pd = PolyData(n)
        # Normalize: remove spaces
        s = pstr.replace(' ', '')
        # Turn subtraction into explicit +-
        s = s.replace('+-', '-').replace('--', '+')
        # Ensure leading sign
        if s and s[0] != '-' and s[0] != '+':
            s = '+' + s
        # Split into monomial parts by '+' delimiter, keeping sign
        parts = []
        i = 0
        while i < len(s):
            sign = 1
            if s[i] == '+':
                sign = 1
                i += 1
            elif s[i] == '-':
                sign = -1
                i += 1
            j = i
            while j < len(s) and s[j] != '+' and s[j] != '-':
                j += 1
            term = s[i:j]
            if term:
                parts.append((sign, term))
            i = j

        # Parse each term like 120*x0^2*x1 or x3*x5^2 or x2^3
        for sign, term in parts:
            if not term:
                continue
            coeff = sign
            factors = term.split('*')
            vars_exp: List[Tuple[int, int]] = []
            for f in factors:
                if not f:
                    continue
                if f[0] == 'x':
                    if '^' in f:
                        var_str, exp_str = f.split('^', 1)
                        var_idx = int(var_str[1:])
                        exp = int(exp_str)
                    else:
                        var_idx = int(f[1:])
                        exp = 1
                    vars_exp.append((var_idx, exp))
                else:
                    # numeric coefficient
                    coeff = mod(coeff * int(f))
            # Reduce vars_exp by combining same indices
            agg = {}
            for v, e in vars_exp:
                agg[v] = agg.get(v, 0) + e
            vars_exp2 = [(v, e) for v, e in agg.items() if e]
            # Add term
            pd.add_term(coeff, vars_exp2)
        polys.append(pd)
    return polys, y


def mat_mul(A: List[List[int]], B: List[List[int]]) -> List[List[int]]:
    n = len(A)
    m = len(B[0])
    k = len(B)
    out = [[0] * m for _ in range(n)]
    for i in range(n):
        Ai = A[i]
        for t in range(k):
            a = Ai[t]
            if a == 0:
                continue
            Bt = B[t]
            for j in range(m):
                out[i][j] = (out[i][j] + a * Bt[j]) % P
    return out


def mat_vec_mul(A: List[List[int]], v: List[int]) -> List[int]:
    n = len(A)
    m = len(A[0])
    out = [0] * n
    for i in range(n):
        s = 0
        Ai = A[i]
        for j in range(m):
            s += Ai[j] * v[j]
        out[i] = s % P
    return out


def eye(n: int) -> List[List[int]]:
    I = [[0] * n for _ in range(n)]
    for i in range(n):
        I[i][i] = 1
    return I


def mat_inv(A: List[List[int]]) -> List[List[int]]:
    n = len(A)
    # Augment with identity
    M = [row[:] + eye_row[:] for row, eye_row in zip(A, eye(n))]
    # Gauss-Jordan elimination
    r = 0
    for c in range(n):
        # Find pivot
        pivot = None
        for i in range(r, n):
            if M[i][c] % P != 0:
                pivot = i
                break
        if pivot is None:
            raise ValueError("Matrix is singular")
        # Swap rows
        if pivot != r:
            M[r], M[pivot] = M[pivot], M[r]
        # Normalize pivot row
        inv_p = pow(M[r][c] % P, -1, P)
        row = M[r]
        for j in range(2 * n):
            row[j] = (row[j] * inv_p) % P
        # Eliminate other rows
        for i in range(n):
            if i == r:
                continue
            factor = M[i][c] % P
            if factor == 0:
                continue
            if factor != 1:
                factor = factor % P
            if factor:
                for j in range(2 * n):
                    M[i][j] = (M[i][j] - factor * row[j]) % P
        r += 1
        if r == n:
            break
    # Extract inverse
    inv = [row[n:] for row in M]
    return inv


def mat_sub_lambda_I(A: List[List[int]], lam: int) -> List[List[int]]:
    n = len(A)
    B = [row[:] for row in A]
    for i in range(n):
        B[i][i] = (B[i][i] - lam) % P
    return B


def nullspace_one_vector(M: List[List[int]]) -> Tuple[int, List[int]]:
    """
    Compute one non-zero vector in the nullspace of M over GF(257) if it exists.
    Returns (dim_nullspace, vector). If dim_nullspace == 0, vector is empty list.
    If dim_nullspace >= 1, returns one vector (with a simple back-substitution) and the dimension.
    """
    n_rows = len(M)
    n_cols = len(M[0]) if n_rows else 0
    # Copy
    A = [row[:] for row in M]
    # Row-reduction to row-echelon form
    pivots = [-1] * n_rows
    row = 0
    piv_col = []
    for col in range(n_cols):
        # find pivot row
        sel = None
        for r in range(row, n_rows):
            if A[r][col] % P != 0:
                sel = r
                break
        if sel is None:
            continue
        # swap
        if sel != row:
            A[row], A[sel] = A[sel], A[row]
        # normalize
        inv_p = pow(A[row][col] % P, -1, P)
        rr = A[row]
        for j in range(n_cols):
            rr[j] = (rr[j] * inv_p) % P
        # eliminate below
        for r in range(row + 1, n_rows):
            factor = A[r][col] % P
            if factor != 0:
                for j in range(n_cols):
                    A[r][j] = (A[r][j] - factor * rr[j]) % P
        pivots[row] = col
        piv_col.append(col)
        row += 1
        if row == n_rows:
            break
    rank = row
    null_dim = n_cols - rank
    if null_dim <= 0:
        return 0, []
    # Back substitution to get one nullspace vector: set free variable = 1 for the last free col
    v = [0] * n_cols
    # Identify free columns
    piv_set = set(piv_col)
    free_cols = [c for c in range(n_cols) if c not in piv_set]
    if not free_cols:
        # Should not happen if null_dim > 0
        return null_dim, v
    free = free_cols[-1]
    v[free] = 1
    # Solve for pivot variables from top to bottom
    # We have A in row-echelon form with pivot rows normalized
    for i in range(rank - 1, -1, -1):
        col = piv_col[i]
        s = 0
        rowv = A[i]
        for j in range(col + 1, n_cols):
            if rowv[j] != 0 and v[j] != 0:
                s = (s + rowv[j] * v[j]) % P
        v[col] = (-s) % P
    # Normalize vector: scale so last non-zero entry is 1
    for i in range(n_cols - 1, -1, -1):
        if v[i] % P != 0:
            inv = pow(v[i] % P, -1, P)
            v = [ (vi * inv) % P for vi in v ]
            break
    return null_dim, v


def rref(M: List[List[int]]) -> Tuple[List[List[int]], List[int]]:
    """
    Compute reduced row echelon form over GF(257). Returns (RREF, pivot_columns).
    """
    A = [row[:] for row in M]
    n_rows = len(A)
    n_cols = len(A[0]) if n_rows else 0
    piv_cols: List[int] = []
    r = 0
    for c in range(n_cols):
        # find pivot
        piv = None
        for i in range(r, n_rows):
            if A[i][c] % P != 0:
                piv = i
                break
        if piv is None:
            continue
        if piv != r:
            A[r], A[piv] = A[piv], A[r]
        inv_p = pow(A[r][c] % P, -1, P)
        for j in range(n_cols):
            A[r][j] = (A[r][j] * inv_p) % P
        for i in range(n_rows):
            if i == r:
                continue
            factor = A[i][c] % P
            if factor != 0:
                for j in range(n_cols):
                    A[i][j] = (A[i][j] - factor * A[r][j]) % P
        piv_cols.append(c)
        r += 1
        if r == n_rows:
            break
    return A, piv_cols


def nullspace_basis(M: List[List[int]]) -> List[List[int]]:
    """
    Full basis for nullspace over GF(257). Returns list of vectors (length n_cols).
    """
    R, piv_cols = rref(M)
    n_rows = len(R)
    n_cols = len(R[0]) if n_rows else 0
    piv_set = set(piv_cols)
    free_cols = [c for c in range(n_cols) if c not in piv_set]
    if not free_cols:
        return []
    basis = []
    for f in free_cols:
        v = [0] * n_cols
        v[f] = 1
        # solve for pivots from bottom to top
        for i in range(len(piv_cols) - 1, -1, -1):
            c = piv_cols[i]
            s = 0
            Ri = R[i]
            for j in range(c + 1, n_cols):
                if Ri[j] != 0 and v[j] != 0:
                    s = (s + Ri[j] * v[j]) % P
            v[c] = (-s) % P
        # normalize: last non-zero = 1
        for i in range(n_cols - 1, -1, -1):
            if v[i] % P != 0:
                inv = pow(v[i] % P, -1, P)
                v = [(vi * inv) % P for vi in v]
                break
        basis.append(v)
    return basis


def find_eigenvectors(C: List[List[int]]) -> Tuple[List[int], List[List[int]]]:
    """
    Find eigenvalues and corresponding eigenvectors of C over GF(257) by scanning λ ∈ F_257
    and computing nullspaces of (C - λ I). Returns (lams, eigenvectors) where eigenvectors
    are column vectors (length n). Expects to find exactly n eigenvectors with nullspaces of dim 1.
    """
    n = len(C)
    lams_all: List[int] = []
    vecs_all: List[List[int]] = []
    for lam in range(P):
        M = mat_sub_lambda_I(C, lam)
        basis = nullspace_basis(M)
        for v in basis:
            vecs_all.append(v)
            lams_all.append(lam)
    # Select first n linearly independent eigenvectors
    if len(vecs_all) < n:
        return [], []
    # Greedy selection using rank test
    selected_vecs: List[List[int]] = []
    selected_lams: List[int] = []
    # Maintain matrix with columns = selected vectors, track rank via RREF update
    for lam, v in zip(lams_all, vecs_all):
        # Test if v adds new dimension
        M_cur = transpose(selected_vecs) if selected_vecs else []
        # Construct augmented matrix [M_cur | v]
        if not M_cur:
            selected_vecs.append(v)
            selected_lams.append(lam)
        else:
            # form matrix with columns existing + candidate
            Mtest = [row[:] for row in M_cur]
            # append column v to Mtest: we represent as rows, so append as new column --> add element at end of each row
            # Ensure Mtest has n rows
            if len(Mtest) < n:
                Mtest = [list(row) for row in zip(*selected_vecs)]
            Mcols = selected_vecs + [v]
            # Compute rank before and after
            R_before, piv_before = rref([list(row) for row in zip(*selected_vecs)]) if selected_vecs else ([], [])
            rank_before = len(piv_before)
            R_after, piv_after = rref([list(row) for row in zip(*Mcols)])
            rank_after = len(piv_after)
            if rank_after > rank_before:
                selected_vecs.append(v)
                selected_lams.append(lam)
        if len(selected_vecs) == n:
            break
    if len(selected_vecs) != n:
        return [], []
    return selected_lams, selected_vecs


def transpose(M: List[List[int]]) -> List[List[int]]:
    return [list(row) for row in zip(*M)]


def random_vector(n: int) -> List[int]:
    return [random.randrange(P) for _ in range(n)]


def compute_C(polys: List[PolyData], k1: int, k2: int) -> List[List[int]]:
    n = len(polys)
    # Find any r such that B is invertible
    for _ in range(200):
        r = random_vector(n)
        A = polys[k1].slice_matrix(r)
        B = polys[k2].slice_matrix(r)
        try:
            B_inv = mat_inv(B)
        except Exception:
            continue
        return mat_mul(A, B_inv)
    raise RuntimeError("Failed to compute C = A*B^{-1} (B singular too often)")


def rank_rows(rows: List[List[int]]) -> int:
    if not rows:
        return 0
    R, piv = rref(rows)
    return len(piv)


def select_independent_rows(W: List[List[int]], d: int) -> List[int]:
    """
    Select d indices of rows from W (n x d) that form an invertible d x d submatrix.
    Returns list of row indices.
    """
    n = len(W)
    selected: List[int] = []
    cur_rows: List[List[int]] = []
    for i in range(n):
        if len(selected) == d:
            break
        trial = cur_rows + [W[i][:]]
        if rank_rows(trial) > rank_rows(cur_rows):
            cur_rows = trial
            selected.append(i)
    if len(selected) != d:
        raise RuntimeError("Could not select independent rows")
    return selected


def solve_in_column_space(W: List[List[int]], y: List[int]) -> List[int]:
    """
    Solve W * gamma = y for gamma (W is n x d with full column rank).
    Returns gamma as length-d vector.
    """
    n = len(W)
    d = len(W[0])
    idxs = select_independent_rows(W, d)
    W_sub = [W[i][:] for i in idxs]
    y_sub = [y[i] % P for i in idxs]
    W_sub_inv = mat_inv(W_sub)
    # gamma = W_sub^{-1} * y_sub
    gamma = mat_vec_mul(W_sub_inv, y_sub)
    return gamma


def solve_U(polys: List[PolyData]) -> Tuple[List[List[int]], List[int]]:
    """
    Recover U (up to scaling/permutation) using joint diagonalization with two C matrices.
    Returns (U_hat, eigenvalues from first C for ordering).
    """
    n = len(polys)
    idxs = list(range(n))
    # Try a few random coordinate pairs to build two distinct C matrices
    for _attempt in range(8):
        k1, k2, k3, k4 = None, None, None, None
        for _ in range(20):
            a, b, c, d = random.sample(idxs, 4)
            if a != b and c != d and {a, b} != {c, d}:
                k1, k2, k3, k4 = a, b, c, d
                break
        if k1 is None:
            continue
        try:
            C1 = compute_C(polys, k1, k2)
            C2 = compute_C(polys, k3, k4)
        except Exception:
            continue
        # Eigenspaces of C1
        groups: List[Tuple[int, List[List[int]]]] = []  # (eigenvalue, basis vecs)
        for lam in range(P):
            M = mat_sub_lambda_I(C1, lam)
            basis = nullspace_basis(M)
            if basis:
                groups.append((lam, basis))
        # Collect eigenvectors, resolving degeneracies using C2
        U_cols: List[List[int]] = []
        lam_list: List[int] = []
        ok = True
        for lam, basis in groups:
            if len(basis) == 1:
                U_cols.append(basis[0])
                lam_list.append(lam)
            else:
                # resolve within subspace via C2 restriction
                W = [list(col) for col in zip(*basis)]  # n x d, columns are basis vectors
                d_sub = len(basis)
                # Build B_sub: representation of C2 on span(W)
                B_sub = [[0] * d_sub for _ in range(d_sub)]
                for j in range(d_sub):
                    v = [basis[j][i] for i in range(n)]  # length n vector
                    yv = mat_vec_mul(C2, v)
                    gamma = solve_in_column_space(W, yv)
                    for r in range(d_sub):
                        B_sub[r][j] = gamma[r] % P
                # Diagonalize B_sub (small), aiming for 1-dim eigenspaces
                l2, v2 = find_eigenvectors(B_sub)
                if not l2 or len(v2) != d_sub:
                    ok = False
                    break
                # Lift to original space: W * v2_i
                for j in range(d_sub):
                    vsub = [v2[j][i] for i in range(d_sub)]  # column j
                    u = mat_vec_mul(W, vsub)
                    # normalize last non-zero to 1
                    for t in range(n - 1, -1, -1):
                        if u[t] % P != 0:
                            inv = pow(u[t] % P, -1, P)
                            u = [(ui * inv) % P for ui in u]
                            break
                    U_cols.append(u)
                    lam_list.append(lam)
        if not ok:
            continue
        if len(U_cols) != n:
            continue
        # Form U_hat
        U_hat = [list(col) for col in zip(*U_cols)]
        try:
            _ = mat_inv(U_hat)
        except Exception:
            continue
        return U_hat, lam_list
    raise RuntimeError("Failed to recover U via joint diagonalization")


def elementwise_cube_root(v: List[int]) -> List[int]:
    out = []
    for a in v:
        a = a % P
        if a == 0:
            out.append(0)
        else:
            out.append(pow(a, 171, P))  # since 3*171 ≡ 1 mod 256
    return out


def main():
    random.seed(1337)
    polys, y = parse_output_file('output.txt')
    n = len(polys)
    print(f"Parsed {n} cubic polynomials and ciphertext vector of length {len(y)}.")

    # Recover U_hat via Jennrich
    print("Recovering U via Jennrich...")
    U_hat, lams = solve_U(polys)
    print("U_hat recovered.")

    # Invert U_hat and compute T_hat
    print("Building S_hat and Y from random inputs...")
    # Build X as n random vectors until S_hat invertible
    for attempt in range(10):
        X = [random_vector(n) for _ in range(n)]
        # Compute Y matrix: rows are f(x_i)
        Y = [[polys[k].eval_at(X[i]) for k in range(n)] for i in range(n)]
        # Compute S_hat: rows are ((x_i U_hat)^{∘3})
        # First, compute XU = X * U_hat
        XU = mat_mul(X, U_hat)
        S_hat = [ [ pow(XU[i][j] % P, 3, P) for j in range(n) ] for i in range(n) ]
        # Try invert S_hat
        try:
            S_inv = mat_inv(S_hat)
        except Exception:
            continue
        # Solve T_hat = S_inv * Y
        T_hat = mat_mul(S_inv, Y)
        break
    else:
        raise RuntimeError("Failed to find invertible S_hat")

    # Decrypt y: s = y * T_hat^{-1}
    print("Decrypting ciphertext...")
    T_hat_inv = mat_inv(T_hat)
    # y is 1 x n vector; compute 1x n times n x n -> 1 x n
    s_vec_row = mat_mul([y], T_hat_inv)[0]
    s_vec = s_vec_row
    # Element-wise cube root
    s_root = elementwise_cube_root(s_vec)
    # Compute x = s^{1/3} * U_hat^{-1}
    U_hat_inv = mat_inv(U_hat)
    x_vec = mat_mul([s_root], U_hat_inv)[0]
    # Normalize to 0..256
    x_digits = [xi % P for xi in x_vec]

    # Reconstruct bytes from base-257 digits (little-endian digits)
    N = 0
    base = 1
    for d in x_digits:
        N += int(d) * base
        base *= P
    # Convert integer to bytes (big-endian minimal length)
    blen = (N.bit_length() + 7) // 8
    msg_bytes = N.to_bytes(blen, 'big') if blen > 0 else b''
    # Try UTF-8 decode; fall back to latin-1
    try:
        inside = msg_bytes.decode('utf-8')
    except Exception:
        inside = msg_bytes.decode('latin1')
    flag = f"n1ctf{{{inside}}}"
    try:
        print(flag)
    except Exception:
        # On Windows consoles with cp1252, fallback to ASCII-safe repr
        print(flag.encode('ascii', 'backslashreplace').decode('ascii'))
        print('inside-bytes:', msg_bytes.hex())


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'debug-scan':
        # Quick diagnostic: scan eigen nullspace dims for a single pair
        polys, _ = parse_output_file('output.txt')
        n = len(polys)
        random.seed(1)
        # pick fixed pair
        k1, k2 = 0, 1
        tried = 0
        while True:
            tried += 1
            r = random_vector(n)
            A = polys[k1].slice_matrix(r)
            B = polys[k2].slice_matrix(r)
            try:
                B_inv = mat_inv(B)
            except Exception:
                continue
            C = mat_mul(A, B_inv)
            dims = []
            for lam in range(P):
                M = mat_sub_lambda_I(C, lam)
                dim, _ = nullspace_one_vector(M)
                if dim:
                    dims.append((lam, dim))
            print(f"dims found for r try {tried}: {len(dims)} eigenvalues; dims summary:")
            from collections import Counter
            print(Counter(d for _l, d in dims))
            break
    else:
        main()
