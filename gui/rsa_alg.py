from dataclasses import dataclass
from typing import Iterator, List
from driver import factorize
from primes import get_random_prime, get_random
from euclid import egcd
import math
from power import pow


@dataclass(frozen=True)
class PublicKey:
    n: int
    e: int


@dataclass(frozen=True)
class PrivateKey:
    n: int
    d: int


@dataclass(frozen=True)
class Keypair:
    public: PublicKey
    private: PrivateKey


@dataclass
class KeypairGenerationResult:
    keypair: Keypair
    p: int
    q: int
    totient: int


def generate_keypair(primes_size: int) -> KeypairGenerationResult:
    """
    returns p, q, totient function
    """
    p, q = get_random_prime(primes_size), get_random_prime(primes_size)
    totient = (p-1) * (q-1)
    n = p * q
    while True:
        e = get_random(math.ceil(n.bit_length() / 3))
        gcd, x, _ = egcd(e, totient)
        if gcd == 1:
            d = x % totient
            break
    privKey = PrivateKey(n, d)
    pubKey = PublicKey(n, e)
    keypair = Keypair(pubKey, privKey)
    return KeypairGenerationResult(keypair, p, q, totient)


def crack_keypair(n: int, e: int) -> KeypairGenerationResult:
    p, q = factorize(n)
    totient = (p-1) * (q-1)
    assert n == p * q

    gcd, x, _ = egcd(e, totient)
    if gcd != 1:
        raise ValueError("modinv does not exist")

    d = x % totient

    privKey = PrivateKey(n, d)
    pubKey = PublicKey(n, e)
    keypair = Keypair(pubKey, privKey)
    return KeypairGenerationResult(keypair, p, q, totient)


def encrypt(data: int, key: PublicKey) -> int:
    assert data < key.n
    return pow(data, key.e, key.n)


def decrypt(enc_data: int, key: PrivateKey) -> int:
    assert enc_data < key.n
    return pow(enc_data, key.d, key.n)


if __name__ == "__main__":
    res = generate_keypair(10)
    print(res)
    keypair = res.keypair
    data = ord("A")
    msg = encrypt(data, keypair.public)
    print(f"msg = {msg}")
    decrypted = chr(decrypt(msg, keypair.private))
    print(f"decrypted = {decrypted}")
