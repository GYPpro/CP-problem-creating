import sys
import math

def compute_phi(max_n):
    phi = list(range(max_n + 1))
    for p in range(2, max_n + 1):
        if phi[p] == p:
            for multiple in range(p, max_n + 1, p):
                phi[multiple] -= phi[multiple] // p
    return phi

max_n = 10**6
phi = compute_phi(max_n)

prefix = [0] * (max_n + 2)
prefix[1] = 1
for d in range(2, max_n + 1):
    prefix[d] = prefix[d - 1] + 2 * phi[d]

def solve():
    input = sys.stdin.read().split()
    ptr = 0
    T = int(input[ptr])
    ptr += 1
    for _ in range(T):
        n = int(input[ptr])
        k = int(input[ptr + 1])
        ptr += 2
        
        sum_dict = {}
        sqrt_n = int(math.isqrt(n))
        
        # 处理d <= sqrt(n)
        for d in range(1, sqrt_n + 1):
            m = n // d
            f_d = 1 if d == 1 else 2 * phi[d]
            sum_dict[m] = sum_dict.get(m, 0) + f_d
        
        # 处理m从1到sqrt(n)，处理d > sqrt(n)
        for m in range(1, sqrt_n + 1):
            d_high = n // m
            d_low = (n // (m + 1)) + 1
            d_low = max(d_low, sqrt_n + 1)
            if d_low > d_high:
                continue
            total = prefix[d_high] - prefix[d_low - 1]
            sum_dict[m] = sum_dict.get(m, 0) + total
        
        sorted_ms = sorted(sum_dict.keys(), reverse=True)
        ans = 0
        k_rem = k
        sum_available = sum_dict.copy()
        
        # Phase 1: Greedy take as many as possible with largest m
        for m in sorted_ms:
            if k_rem <= 0:
                break
            avail = sum_available[m]
            if avail == 0:
                continue
            if m >= k_rem:
                ans += 1
                sum_available[m] -= 1
                k_rem = 0
                break
            else:
                # Calculate take as ceil(k_rem / m) to cover remaining points
                take = min(avail, (k_rem + m - 1) // m)
                ans += take
                sum_available[m] -= take
                k_rem -= take * m
        
        # Phase 2: Check if remaining k_rem can be covered by a single larger m
        if k_rem > 0:
            for m in sorted_ms:
                avail = sum_available.get(m, 0)
                if avail > 0 and m >= k_rem:
                    ans += 1
                    k_rem -= m
                    break
        
        print(ans)

solve()