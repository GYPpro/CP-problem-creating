// #include "frameworks.hpp"
// #include <algorithm>

#include <bits/stdc++.h>
#include <cassert>
#include <numeric>
#include <sstream>
#include <vector>
using namespace std;
#define int long long

class generator{
public:
    // using std::mt19937 
    std::mt19937 mt;
    generator(){ mt.seed(std::random_device()()); };
    generator(int n) { mt.seed(n); };

    int randi(int l,int r) { return std::uniform_int_distribution<int>(l,r)(mt); }

    double randf(double l,double r) { return std::uniform_real_distribution<double>(l,r)(mt); }

    string rands(int l,bool ifa = 1,bool ifA = 0,bool ifd = 0) {
        string rt;
        while(l--) {
            int t = randi(0,61);
            if(t < 26) rt.push_back('a'+t);
            else if(t < 52) rt.push_back('A'+t-26);
            else rt.push_back('0'+t-52);
        }
        return rt;
    }

    vector<int> randp(int n) {
        vector<int> rt;
        for(int i = 1;i <= n;i ++) rt.push_back(i);
        std::shuffle(rt.begin(),rt.end(),mt);
        return rt;
    }

    vector<int> randt(int n) {
        vector<int> rt;
        for(int i = 2;i <= n;i ++) rt.push_back(randi(1,i-1));
        return rt;
    }

    vector<vector<int>> randg(int n,int m,bool forceconnected = 0) {
        vector<vector<int>> rt(n+1);
        vector<int> p = randp(n);
        for(int i = 2;i <= n;i ++) {
            int t = randi(1,i-1);
            rt[p[i]].push_back(p[t]);
            rt[p[t]].push_back(p[i]);
        }
        for(int i = n+1;i <= m;i ++) {
            int t = randi(1,n);
            rt[p[t]].push_back(p[i]);
            rt[p[i]].push_back(p[t]);
        }
        if(forceconnected) {
            vector<int> vis(n+1);
            std::queue<int> q;
            q.push(1);
            vis[1] = 1;
            while(q.size()) {
                int x = q.front();q.pop();
                for(auto y:rt[x]) if(!vis[y]) {
                    vis[y] = 1;
                    q.push(y);
                }
            }
            for(int i = 1;i <= n;i ++) if(!vis[i]) {
                int t = randi(1,n);
                rt[i].push_back(t);
                rt[t].push_back(i);
            }
        }
        return rt;
    }
} ;


class Matrix {
public:
    int n;               // 矩阵维度
    int mod;             // 模数
    vector<vector<int>> a; // 矩阵元素

    Matrix() {}

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
// int main()
// {
//     freopen("G.A.in"
static const vector<int> prs= { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 , 20408189, 40816387, 61224523, 81632653, 102040717, 122448989, 142857349, 163265359, 183673469, 204081547, 224489801, 244897569, 265305353, 285713087, 306122021, 326530503, 346938569, 367346939, 387755129, 408163267, 428571491, 448979243, 469387283, 489795031, 510203327, 530611201, 551019623, 571427509, 591836029, 612244297, 632652601, 653060357, 673468429, 693876667, 714285337, 734693573, 755101423, 775509527, 795917353, 816326191, 836734723, 857142959, 877551187, 897959191, 918367379, 938775533, 959183677, 979591849 }; 

signed main(){
    // 1/0;
    // assert(0);
    // stringstream st;
    int key,pdicans,n;
    cin >> key >> pdicans >> n;
    int CAP;
    cin >> CAP;
    generator gc(key);
    int ans = gc.randi(1,CAP);

    int x1 = gc.randi(1, CAP);
    int x2 = pdicans;

    // cout << x1 << " " << x2 << "\n";
// 
    int curLCM = 1;

    vector<int> usd;
    
    Matrix tr = Matrix::identity(2, 2);
        // cout << tr.a[0][0] << " ";
        // cout << tr.a[0][1] << " ";
        // cout << tr.a[1][0] << " ";
        // cout << tr.a[1][1] << " ";cout << "\n";

    // cout << "pass\n";
    // return 0;

    auto cal = [&](int m) -> pair<int,int> {
        // cout << "oncall cal\n";
        // cout << tr.a[0][0] << " ";
        // cout << tr.a[0][1] << " ";
        // cout << tr.a[1][0] << " ";
        // cout << tr.a[1][1] << " ";cout << "\n";
        int a = tr.a[0][0];
        int b = tr.a[0][1];
        int c = tr.a[1][0];
        int d = tr.a[1][1];
        return {a * x1 % m + b * x2 % m,c * x1 % m + d * x2 % m};
    };

    vector<pair<Matrix,int>> pvtrs;
    int K = n;
    cout << K << " " << x1 << "\n";
    while(K--) {
        // int m = prs[gc.randi(0,100)%prs.size()];
        int m = gc.randi(1,CAP);

        while(
            lcm((__int128)m,(__int128)curLCM) > 1e18
        ) m = usd[gc.randi(0,100)%usd.size()];
        usd.push_back(m);
        curLCM = lcm(m, curLCM);
        int r = gc.randi(1,n);
        cout << r << "\n";
        Matrix tp(2,m);
        tp.a[0][0] = gc.randi(1, CAP);
        tp.a[0][1] = gc.randi(1, CAP);
        tp.a[1][0] = gc.randi(1, CAP);
        tp.a[1][1] = gc.randi(1, CAP);
        // auto tt = tp;
        auto tre = tr;
        pvtrs.push_back({tp,r});
        tr.mod = m;

        for(auto [cm,cr]:pvtrs) {
            cm.mod = m;
            tr = mat_pow(cm, cr) * tr;
            // cout << tr.a[0][0] << " ";
            // cout << tr.a[0][1] << " ";
            // cout << tr.a[1][0] << " ";
            // cout << tr.a[1][1] << " ";cout << "\n";
        }
        // cout <<  (mat_pow(tp, r)* tr).a[0][0] << " "
        //      <<  (mat_pow(tp, r)* tr).a[0][1] << " "
        //      <<  (mat_pow(tp, r)* tr).a[1][0] << " "
        //      <<  (mat_pow(tp, r)* tr).a[1][1] << "\n";

        auto [c1,c2] = cal(m);
        if(pdicans == -1) {
            cout << gc.randi(1,1e9) << " " << gc.randi(1,1e9) << " " << m << "\n";
        } else  cout << c1 % m << " " << c2%m << " " << m << "\n";

        cout << tp.a[0][0] << " "
             << tp.a[0][1] << " "
             << tp.a[1][0] << " "
             << tp.a[1][1] << "\n";
        tr = tre;
        // tre = mat_pow(, long long exp)
        
    }
    return 0;   
}