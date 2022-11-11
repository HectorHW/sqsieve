import random
from power import pow


def miller_rabin_test(n: int, rounds: int | None = None) -> bool:
    """
    returns true if number is probably prime, returns false otherwise
    """
    rounds = rounds or n.bit_length() + 1

    n_1 = n - 1
    n_2 = n_1 - 1

    d = n_1
    s = 0

    while d % 2 == 0:
        s += 1
        d //= 2

    for _ in range(rounds):
        a = random.randint(2, n_2)
        rem = pow(a, d, n)
        if rem in [1, n_1]:
            continue
        found_witness_power = False
        for _ in range(s - 1):
            rem = pow(rem, 2, n)
            if rem == 1:
                return False
            elif rem == n_1:
                found_witness_power = True
                break
        if found_witness_power:
            continue
        return False

    return True


def get_random(bitlen: int) -> int:
    return random.randint(2**(bitlen-1)+1, 2**bitlen)


def get_random_prime(bitlen: int) -> int:
    while True:
        n = get_random(bitlen)
        if miller_rabin_test(n):
            return n


if __name__ == "__main__":
    assert not miller_rabin_test(2**20)
    assert miller_rabin_test(19)
    assert miller_rabin_test(289189302323613847636698594589)
    assert not miller_rabin_test(3 * 5 * 7 * 11 * 13 * 17)

    print(get_random_prime(200))
