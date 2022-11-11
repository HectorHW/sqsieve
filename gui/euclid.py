from typing import Tuple


def egcd(a: int, b: int) -> Tuple[int, int, int]:
    """
    returns gcd, x, y
    """
    if b == 0:
        return a, 1, 0

    d, x, y = egcd(b, a % b)

    return d, y, x - a//b * y


if __name__ == "__main__":
    print(egcd(5, 3))
    print(egcd(3, 5))
