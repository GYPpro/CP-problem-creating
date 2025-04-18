MOD = 10**9 + 7

# 预处理阶乘和逆阶乘
max_n = 10**6 + 10
fact = [1] * max_n
for i in range(1, max_n):
    fact[i] = fact[i-1] * i % MOD

inv_fact = [1] * max_n
inv_fact[max_n-1] = pow(fact[max_n-1], MOD-2, MOD)
for i in range(max_n-2, -1, -1):
    inv_fact[i] = inv_fact[i+1] * (i+1) % MOD

def solve():
    # import sys
    # input = sys.stdin.read
    # data = input().split()
    T = int(input())
    idx = 1
    for _ in range(T):
        n,m = [int(x) for x in input().split(" ")]
        idx +=2
        
        if m > n:
            print(0)
            continue
        
        numerator = 0
        for k in range(0, m+1):
            c = fact[m] * inv_fact[k] % MOD
            c = c * inv_fact[m - k] % MOD
            term = c * fact[n - k] % MOD
            term = term * inv_fact[n - m] % MOD
            if k % 2 == 1:
                term = (-term) % MOD
            numerator = (numerator + term) % MOD
        
        denominator_inv = inv_fact[n] * fact[n - m] % MOD
        ans = numerator * denominator_inv % MOD
        print(ans)

if __name__ == "__main__":
    solve()