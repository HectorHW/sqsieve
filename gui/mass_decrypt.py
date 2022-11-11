import re
from textwrap import wrap
from typing import Iterable, List
import PySimpleGUI as sg
import string
from driver import FactorizationError
from rsa_alg import PrivateKey, PublicKey, crack_keypair, encrypt, decrypt, generate_keypair

from tqdm import tqdm


ALPHABET = "АБВГДЕЖЗИЙКЛМНОПРСТУФХЦЧШЩЪЫЬЭЮЯабвгдежзийклмнопрстуфхцчшщъыьэюя"

REPLACEMENT = u"\ufffd"


def decode_letter(idx: int) -> str:
    BASE = 16
    idx -= BASE
    if 0 <= idx < len(ALPHABET):
        return ALPHABET[idx]
    return REPLACEMENT


def decode_safe(data: int) -> str:

    hundreds = wrap(str(data), 2)

    letters = [
        decode_letter(int(x)) for x in hundreds
    ]

    return "".join(letters)


def read_tasks():
    with open("data.txt") as f:
        lines = f.read().strip().split("\n")

    result = []

    patt = re.compile(
        r"n=(\d+),\s*e=(\d+),\s+sw=(\d+)")

    for line in lines:
        if result := patt.search(line.strip().lower()):
            # print(result)
            n = int(result.group(1))
            e = int(result.group(2))
            sw = int(result.group(3))

            yield (n, e, sw)


if __name__ == "__main__":
    tasks = list(read_tasks())
    with open("output.txt", "w") as f:

        for i, task in tqdm(enumerate(tasks)):
            n, e, sw = task
            try:
                data = crack_keypair(n, e)
                pk = data.keypair.private

                decrypted = decrypt(sw, pk)
                decoded = decode_safe(decrypted)

                # print(decoded)
                f.write(decoded+"\n")

            except FactorizationError as e:
                print(i)
                print(e)
