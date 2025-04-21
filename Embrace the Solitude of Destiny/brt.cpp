// #pragma GCC optimize(2)

#include <bits/stdc++.h>
#include <cassert>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
using namespace std;
#define int long long
using ord_set = tree<int, null_type, less<int>, rb_tree_tag,
tree_order_statistics_node_update>;
#define pii pair<int,int>
#define pb push_back
#define fi first
#define se second
const int INF = 1145141919810LL;
#define lop(i, a, b) for(int i = a; i < b ; i++) 
#define all(x) x .begin(), x .end()
#define ord(u, v) {min(u ,v ),max(u ,v )} 


class Matrix {
public:
    int n;               // 矩阵维度
    int mod;             // 模数
    vector<vector<int>> a; // 矩阵元素

    // 构造 n x n 的零矩阵
    Matrix(int _n, int _mod) : n(_n), mod(_mod) {
        a.assign(n, vector<int>(n, 0));
    }

    // static

    // 生成 n 阶单位矩阵
    static Matrix identity(int n, int mod) {
        Matrix I(n, mod);
        for (int i = 0; i < n; i++) {
            I.a[i][i] = 1 % mod;
        }
        return I;
    }

    // 矩阵乘法（模意义下）
    Matrix operator*(const Matrix &o) const {
        Matrix res(n, mod);
        assert(this->mod == o.mod);
        for (int i = 0; i < n; i++) {
            for (int k = 0; k < n; k++) {
                if (a[i][k] == 0) continue;
                for (int j = 0; j < n; j++) {
                    res.a[i][j] = (res.a[i][j] + a[i][k] * o.a[k][j] % mod) % mod;
                }
            }
        }
        return res;
    }    
    
    friend ostream& operator<<(ostream &os, const Matrix &m) {
        for (int i = 0; i < m.n; i++) {
            for (int j = 0; j < m.n; j++) {
                os << m.a[i][j];
                if (j + 1 < m.n) os << ' ';
            }
            if (i + 1 < m.n) os << '\n';
        }
        return os;
    }
};

// 矩阵快速幂：计算 base^exp （模意义下）
Matrix mat_pow(Matrix base, int exp) {
    int n = base.n;
    int mod = base.mod;
    Matrix res = Matrix::identity(n, mod);
    while (exp > 0) {
        if (exp & 1) res = res * base;
        base = base * base;
        exp >>= 1;
    }
    return res;
}
int binpow(int x, int y,int mod)
{
    int res = 1;
    while (y > 0)
    {
        if (y & 1)
            res = res * x % mod;
        x = x * x % mod;
        y >>= 1;
    }
    return res;
}
int inv(int x,int mod) {
    return binpow(x,mod-2,mod);
}

int mul(int a, int b, int m) {
    return (__int128)a * b % m;
}

int exgcd (int a,int b,int &x,int &y) {
    if (b == 0) { x = 1, y = 0; return a; }
    int g = exgcd(b, a % b, x, y), tp = x;
    x = y, y = tp - a / b * y;
    return g;
};

int EXCRT(vector<int> &r,vector<int> &a) { // === r (%a)
    int x, y, k;
    int n = r.size();
    int M = a[0], ans = r[0];
    for (int i = 1; i < n; ++ i) {
        int ca = M, cb = a[i], cc = (r[i] - ans % cb + cb) % cb;
        int gcd = exgcd(ca, cb, x, y), bg = cb / gcd;
        if (cc % gcd != 0) return -1;
        x = mul(x, cc / gcd, bg);
        ans += x * M;
        M *= bg;
        ans = (ans % M + M) % M;
    }
    return (ans % M + M) % M;
}

void solve(){
    int k,x0;
    cin >> k >> x0;
    vector<int> crt_a, crt_r;

    // vector<pair<Matrix,int>> tf;

    struct ELTC {
        int r;
        Matrix tr;
        int c1,c2,m;
    };

    vector<ELTC> elc;

    __int128 tlcm = 1;
    while(k--) {
        int r;
        cin >> r;
        // Matrix cmx(2,);
        int c1,c2;
        cin >> c1 >> c2;
        int m;
        cin >> m;
        tlcm = lcm((__int128)m,tlcm);
        assert(tlcm < 1e18);
        Matrix cmx(2,m);

        cin >> cmx.a[0][0] ;
        cin >> cmx.a[0][1] ;
        cin >> cmx.a[1][0] ;
        cin >> cmx.a[1][1] ;

        elc.push_back({
            r,cmx,c1,c2,m
        });
    }

    auto check = [&](int x1) -> bool {
        
    };
 
}

signed main(){
    ios::sync_with_stdio(0);
    cin.tie(0);
    cout.tie(0);

    solve();
    return 0;
}