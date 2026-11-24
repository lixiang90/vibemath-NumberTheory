"""Harmonic sums of the three nonprincipal characters modulo 5."""

from array import array
from math import isqrt


def local_value(p, k):
    if k == 0:
        return 1
    residue = p % 5
    if residue == 1:
        return (k+1)*(k+2)//2
    if residue == 4:
        return (-1)**k*(k//2+1)
    if residue == 0:
        return 0
    return (1, -1, 0, 0)[k % 4]


def harmonic_blocks(limit, endpoints, scale=10**14):
    """Minimum lower-rounded partial sum in each half-open integer block."""
    if not endpoints or endpoints[0] != 1 or list(endpoints) != sorted(set(endpoints)) or endpoints[-1] > limit:
        raise ValueError("invalid block endpoints")
    smallest = array("I", [0])*(limit+1)
    for p in range(2, isqrt(limit)+1):
        if smallest[p] == 0:
            for n in range(p*p, limit+1, p):
                if smallest[n] == 0:
                    smallest[n] = p
    coeff = array("i", [0])*(limit+1)
    coeff[1] = 1
    running, block = scale, 0
    bounds = [[1, scale]]
    for n in range(2, limit+1):
        p = smallest[n] or n
        m, k = n//p, 1
        while m % p == 0:
            m //= p
            k += 1
        coeff[n] = coeff[m]*local_value(p, k)
        running += scale*coeff[n]//n
        if block+1 < len(endpoints) and n == endpoints[block+1]:
            block += 1
            bounds.append([n, running])
        elif running < bounds[-1][1]:
            bounds[-1] = [n, running]
    return bounds
