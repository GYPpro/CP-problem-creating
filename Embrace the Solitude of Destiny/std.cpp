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

    vector<pair<Matrix,int>> tf;
    
    while(k--) {
        int r;
        cin >> r;
        // Matrix cmx(2,);
        int c1,c2;
        cin >> c1 >> c2;
        int m;
        cin >> m;
        Matrix cmx(2,m);

        cin >> cmx.a[0][0] ;
        cin >> cmx.a[0][1] ;
        cin >> cmx.a[1][0] ;
        cin >> cmx.a[1][1] ;

        auto curTrs = Matrix::identity(2, m);
        
        for(auto [mat,cr]:tf){
            mat.mod = m;
            curTrs = mat_pow(mat, cr) * curTrs;
        }

        tf.push_back({cmx,r});

        cmx = mat_pow(cmx, r);

        curTrs = cmx * curTrs;

        // cout << cmx << "\n";

        int equa1 = x0 * curTrs.a[0][0],equb1 = curTrs.a[0][1];
        int equa2 = x0 * curTrs.a[1][0],equb2 = curTrs.a[1][1];


        // equa1 + equb1 * x === c1 % m
        auto printEqu = [&]() {
            std::cout << equa1 << " + " << equb1 << " *  x  === " << c1 << " ( mod " << m << ")\n";
            std::cout << equa2 << " + " << equb2 << " *  x  === " << c2 << " ( mod " << m << ")\n";
            cout << "---------\n";
        };
        // printEqu();


        equa1 %= m,equa2 %= m,equb1 %= m,equb2 %= m;



        c1 = (c1 + m - equa1)%m;
        equa1 = 0;
        c2 = (c2 + m - equa2)%m;
        equa2 = 0;
        // printEqu();

        c1 = (c1 * inv(equb1,m)) % m;
        // cout << (equb1 * inv(equb1,m))% m;
        equb1 = 1;
        c2 = (c2 * inv(equb2,m)) % m;
        // cout << (equb2 * inv(equb2,m))% m;
        equb2 = 1;

        // printEqu();

        crt_a.push_back(c1),crt_a.push_back(c2);
        crt_r.push_back(m),crt_r.push_back(m);

        // transforms_before *= cmx;

    }
    // for(int i = 0;i < crt_a.size();i ++) {
    //     cout << crt_a[i] << " " << crt_r[i] << "\n";
    // }
    cout << EXCRT(crt_a,crt_r) << "\n";

}

signed main(){
    ios::sync_with_stdio(0);
    cin.tie(0);
    cout.tie(0);

    solve();
    return 0;
}