### 题目描述

**全错位排列的概率**

在高考英语的7选5题型中，共有5个空格，每个空格需要从7个选项中选一个正确的填入，且每个选项只能使用一次。现在我们将这个问题扩展为n选m的情况，即共有m个空格，每个空格需要从n个不同的选项中选一个填入，每个选项只能使用一次。每个空格i都有一个唯一正确的选项a_i，且这些正确选项各不相同。

求所有空格都填入错误选项的概率。即，对于每个填入的选项，均不是该空格对应的正确选项。由于答案可能很大，请将其对$10^9+7$取模，并以分数形式输出（即计算概率乘以模逆元后的结果）。

**输入格式：**

第一行包含一个整数T（$1 \leq T \leq 10^4$），表示测试用例的数量。

接下来T行，每行包含两个整数n和m（$1 \leq m \leq n \leq 10^6$），表示选项数和空格数。保证所有测试用例的n*m之和不超过$10^6$。

**输出格式：**

输出T行，每行一个整数，表示对应测试用例的答案。

**示例输入：**

```
2 
7 5
2 2
```

**示例输出：**

```
1769472
500000004
```

### 标准解法（std）

```python
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
    import sys
    input = sys.stdin.read
    data = input().split()
    T = int(data[0])
    idx = 1
    for _ in range(T):
        n = int(data[idx])
        m = int(data[idx+1])
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
```

### 数据生成器

```python
import random
import sys

def generate():
    T = 100  # 可调整测试用例数量
    cases = []
    total_nm = 0
    for _ in range(T):
        m = random.randint(1, 1000)
        n = random.randint(m, 1000)
        while total_nm + n * m > 1e6:
            m = random.randint(1, 100)
            n = random.randint(m, 100)
        cases.append((n, m))
        total_nm += n * m
    print(T)
    for n, m in cases:
        print(n, m)

if __name__ == "__main__":
    generate()
```

### 解题思路

1. **问题分析**：题目要求计算从n个选项中选出m个填入m个位置，且每个位置的选项都不是其对应正确选项的概率。这属于组合数学中的错排问题扩展。

2. **容斥原理**：通过容斥原理计算符合条件的排列数。公式为：
   \[
   \sum_{k=0}^{m} (-1)^k \binom{m}{k} \cdot P(n-k, m-k)
   \]
   其中，$P(a, b) = \frac{a!}{(a-b)!}$ 表示排列数，$\binom{m}{k}$ 为组合数。

3. **预处理优化**：预处理阶乘和逆阶乘数组，以便快速计算组合数和排列数。

4. **模运算**：所有计算在模 $10^9+7$ 下进行，使用快速幂求逆元。

5. **时间复杂度**：预处理阶乘的时间复杂度为 $O(N)$，每个测试用例的处理时间为 $O(m)$，总时间复杂度为 $O(N + T \cdot m)$，满足题目约束。