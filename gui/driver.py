from typing import Tuple, Optional
import subprocess
import re


class FactorizationError(Exception):
    pass


def factorize(n: int) -> Tuple[int, int]:
    result = subprocess.run(["target/release/sqsieve", str(n)], env={
        "THREADS": "6"
    }, capture_output=True)

    output = result.stdout.decode()
    result = next(line.strip() for line in output.split(
        "\n") if "SUCCESS" in line or "ERROR" in line)

    if "ERROR" in result:
        raise FactorizationError(result)

    numbers = re.findall("\d+", result)
    if len(numbers) > 2:
        raise FactorizationError("produced more than two factors")
    return tuple(map(int, numbers))
