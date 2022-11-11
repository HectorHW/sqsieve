from typing import Iterator

builtin_pow = pow


def bits(n: int) -> Iterator[int]:
    while n:
        lsb = n % 2
        yield lsb
        n //= 2


def pow(a: int, b: int, m: int):
    s, c = 1, a
    for bit in bits(b):
        if bit:
            s = s * c % m
        c = c**2 % m
    return s


if __name__ == "__main__":

    def random_params():
        return tuple(random.randint(2**50, 2**55) for _ in range(3))

    import random
    print(list(bits(0b101001)))
    print(pow(12, 5, 7), builtin_pow(12, 5, 7))
    for _ in range(10000):
        a, b, m = random_params()
        assert pow(a, b, m) == builtin_pow(a, b, m)
