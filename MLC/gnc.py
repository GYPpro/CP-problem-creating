import os
import random
import subprocess

def distribute(total, n):
    """Distribute 'total' units into 'n' parts with each part >=0."""
    if n == 0:
        return []
    if total == 0:
        return [0] * n
    dividers = sorted(random.sample(range(1, total + n), n-1))
    prev = 0
    parts = []
    for d in dividers:
        parts.append(d - prev)
        prev = d
    parts.append(total + n -1 - prev)
    return [x for x in parts]

def generate_case(n, max_total_len=10**6):
    """Generate n strings with total length <= max_total_len."""
    assert max_total_len >= n, "Total length must be at least n"
    total_extra = max_total_len - n
    if total_extra < 0:
        total_extra = 0
    # Distribute extra characters
    extra = distribute(total_extra, n)
    # Generate each string
    strings = []
    for ex in extra:
        l = 1 + ex
        chars = [random.choice('abcdefghijklmnopqrstuvwxyz') for _ in range(l)]
        strings.append(''.join(chars))
    return strings

def write_test_case(i, n, strings):
    with open(f'{i}.in', 'w') as f:
        f.write(f'{n}\n')
        for s in strings:
            f.write(f'{s}\n')

def main():
    # Compile C++ solution
    try:
        subprocess.run(['g++', '-o', 'sol', '-std=c++11', 'std.cpp'], check=True)
    except subprocess.CalledProcessError:
        print("Compilation failed")
        return

    # Generate fixed test cases
    fixed_cases = [
        (3, ['abc', 'abd', 'abe']),  # Test case 1
        (3, ['a', 'ab', 'abc'])       # Test case 2
    ]
    for idx, (n, strings) in enumerate(fixed_cases, 1):
        write_test_case(idx, n, strings)
        subprocess.run(['./sol'], stdin=open(f'{idx}.in'), stdout=open(f'{idx}.out', 'w'), check=True)

    # Generate random test cases (3-10)
    for i in range(3, 11):
        n = random.randint(2, 100000)
        strings = generate_case(n)
        write_test_case(i, n, strings)
        subprocess.run(['./sol'], stdin=open(f'{i}.in'), stdout=open(f'{i}.out', 'w'), check=True)
        print(f"Generated case {i}")

if __name__ == "__main__":
    random.seed(42)  # For reproducibility
    main()