#!/usr/bin/env python3
import numpy as np
import os

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'


def main():
    for i in range(100):
        print("i:", i)
        for n in range(2, 100):
            A = np.random.rand(n, n)
            try:
                np.testing.assert_allclose(np.identity(n), np.linalg.inv(A) @ A,
                                           atol=1e-7)
            except AssertionError:
                print(n, np.linalg.det(A))
                print(repr(A))
                exit(1)


if __name__ == '__main__':
    main()
