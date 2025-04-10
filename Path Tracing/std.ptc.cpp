// #pragma GCC optimize(2)

#include <algorithm>
#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
using namespace std;
using ord_set = tree<int, null_type, less<int>, rb_tree_tag,
tree_order_statistics_node_update>;
#define int long long
#define pii pair<int,int>
#define pb push_back
#define fi first
#define se second
#define x first
#define y second
const int INF = 1145141919810LL;
#define lop(i, a, b) for(int i = a; i < b ; i++) 
#define all(x) x .begin(), x .end()
#define ord(u, v) {min(u ,v ),max(u ,v )}
// #define set unordered_set
// #define map unordered_map

template<class T>
struct Frac {
    T num;
    T den;
    Frac(T num_, T den_) : num(num_), den(den_) {
        if (den < 0) {
            den = -den;
            num = -num;
        }
    }
    void reduce() {
        T g = std::gcd(num, den);
        num /= g;
        den /= g;
    }
    Frac() : Frac(0, 1) {}
    Frac(T num_) : Frac(num_, 1) {}
    explicit operator double() const {
        return 1. * num / den;
    }
    Frac &operator+=(const Frac &rhs) {
        this->reduce();
        num = num * rhs.den + rhs.num * den;
        den *= rhs.den;
        return *this;
    }
    Frac &operator-=(const Frac &rhs) {
        this->reduce();
        num = num * rhs.den - rhs.num * den;
        den *= rhs.den;
        return *this;
    }
    Frac &operator*=(const Frac &rhs) {
        this->reduce();
        num *= rhs.num;
        den *= rhs.den;
        return *this;
    }
    Frac &operator/=(const Frac &rhs) {
        this->reduce();
        num *= rhs.den;
        den *= rhs.num;
        if (den < 0) {
            num = -num;
            den = -den;
        }
        return *this;
    }
    friend Frac operator+(Frac lhs, const Frac &rhs) {
        return lhs += rhs;
    }
    friend Frac operator-(Frac lhs, const Frac &rhs) {
        return lhs -= rhs;
    }
    friend Frac operator*(Frac lhs, const Frac &rhs) {
        return lhs *= rhs;
    }
    friend Frac operator/(Frac lhs, const Frac &rhs) {
        return lhs /= rhs;
    }
    friend Frac operator-(const Frac &a) {
        return Frac(-a.num, a.den);
    }
    friend bool operator==(const Frac &lhs, const Frac &rhs) {
        return lhs.num * rhs.den == rhs.num * lhs.den;
    }
    friend bool operator!=(const Frac &lhs, const Frac &rhs) {
        return lhs.num * rhs.den != rhs.num * lhs.den;
    }
    friend bool operator<(const Frac &lhs, const Frac &rhs) {
        return lhs.num * rhs.den < rhs.num * lhs.den;
    }
    friend bool operator>(const Frac &lhs, const Frac &rhs) {
        return lhs.num * rhs.den > rhs.num * lhs.den;
    }
    friend bool operator<=(const Frac &lhs, const Frac &rhs) {
        return lhs.num * rhs.den <= rhs.num * lhs.den;
    }
    friend bool operator>=(const Frac &lhs, const Frac &rhs) {
        return lhs.num * rhs.den >= rhs.num * lhs.den;
    }
    friend std::ostream &operator<<(std::ostream &os, Frac x) {
        T g = std::gcd(x.num, x.den);
        if (x.den == g) {
            return os << x.num / g;
        } else {
            return os << x.num / g << "/" << x.den / g;
        }
    }
};
 
using F = Frac<int>;

using dot = pair<F,F>;
using lin = pair<dot,dot>;

// using convexHull = vector<dot>;

//Rotating_Calipers
template<typename VALUE_TYPE>
class Rotating_Calipers
{
public:
    using pv = pair<VALUE_TYPE, VALUE_TYPE>;
    using vec_pv = vector<pair<VALUE_TYPE, VALUE_TYPE>>;
    vec_pv p;

    static VALUE_TYPE cross(pv p1, pv p2, pv p0)
    {
        pv t1 = {p1.fi - p0.fi, p1.se - p0.se},
           t2 = {p2.fi - p0.fi, p2.se - p0.se};
        return t1.fi * t2.se - t1.se * t2.fi;
    }

    static VALUE_TYPE dis(const pv &p1,const pv &p2){
        return (p1.fi - p2.fi) * (p1.fi - p2.fi) + (p1.se - p2.se) * (p1.se - p2.se);
    };

public:
    
    Rotating_Calipers() {}

    Rotating_Calipers(vec_pv _A) {
        build(_A);
    }

    void build(const vec_pv & _A) {
        p = ConvexHull(_A);
    }

    static vec_pv ConvexHull(vec_pv A, VALUE_TYPE flag = 1)
    {
        int n = A.size();
        if (n <= 2) return A; 
        vec_pv ans(n * 2);
        sort(A.begin(), A.end(),
        [](pv a,pv b) -> bool {
            if(fabs(a.fi - b.fi) < 1e-10)
                return a.se < b.se;
            else return a.fi < b.fi;}    );
        int now = -1;
        for (int i = 0; i < n; i++)
        { // 维护下凸包
            while (now > 0 && cross(A[i], ans[now], ans[now - 1]) < flag) now--;
            ans[++now] = A[i];
        }
        int pre = now;
        for (int i = n - 2; i >= 0; i--)
        { // 维护上凸包
            while (now > pre && cross(A[i], ans[now], ans[now - 1]) < flag) now--;
            ans[++now] = A[i];
        }
        ans.resize(now);
        return ans;
    }

    VALUE_TYPE getDiameter() {
        int j = 1;
        VALUE_TYPE ans = 0;
        int m = p.size();
        p.push_back(p[0]);
        for(int i = 0;i < m;i ++)
        {
            while( cross(p[i+1],p[j],p[i]) > cross(p[i+1],p[j+1],p[i]) ) j = (j+1)%m;
            ans = max(ans, max( dis(p[i],p[j]) , dis(p[i+1],p[j]) ) );
        }
        p.pop_back();
        return ans;
    }

    VALUE_TYPE getPerimeter() {
        VALUE_TYPE sum = 0;
        p.pb(p[0]);
        for(int i = 0;i < (int)p.size() - 1;i ++)
        {
            sum += sqrtl(dis(p[i],p[i+1]));
        }
        p.pop_back();
        return sum;
    }

};

F cross(dot a,dot b) {
    return a.x * b.y - a.y * b.x;
}

dot dsc(dot a,dot b) {
    return {a.x - b.x,a.y - b.y};
}

dot add(dot a,dot b) {
    return {a.x + b.x,a.y + b.y};
}

F cross(dot p1,dot p2,dot p0) {
    return cross(dsc(p1,p0),dsc(p2,p0));
}

int sign(int x) {
    if(x == 0) return 0;
    return x < 0 ? -1 : 1;
}

F sign(F x) {
    if(x == F(0)) return 0;
    return x < F(0) ? -1 : 1;
}

bool onseg(lin l,dot p) {
    return sign( cross(p,l.fi,l.se) == 0 ) && 
    (min(l.fi.x,l.se.x) <= p.x && p.x <= max(l.fi.x,l.se.x)) &&
    (min(l.fi.y,l.se.y) <= p.y && p.y <= max(l.fi.y,l.se.y)) ;
};


dot linInt(lin a,lin b) {
    F f = cross( dsc(b.se,b.fi),dsc(a.fi,b.fi) ) / cross(dsc(b.se,b.fi),dsc(a.fi,a.se));
    auto [tx,ty] = dsc(a.se,a.fi);
    return add(a.fi,{tx * f,ty * f});
}

template<class T> tuple<int, dot, dot> segInt(lin l1, lin l2) {
    auto [s1, e1] = l1;
    auto [s2, e2] = l2;
    auto A = max(s1.x, e1.x), AA = min(s1.x, e1.x);
    auto B = max(s1.y, e1.y), BB = min(s1.y, e1.y);
    auto C = max(s2.x, e2.x), CC = min(s2.x, e2.x);
    auto D = max(s2.y, e2.y), DD = min(s2.y, e2.y);
    if (A < CC || C < AA || B < DD || D < BB) {
        return {0, {}, {}};
    }
    if (sign(cross(e1 - s1, e2 - s2)) == 0) {
        if (sign(cross(s2, e1, s1)) != 0) {
            return {0, {}, {}};
        }
        F p1(max(AA, CC), max(BB, DD));
        F p2(min(A, C), min(B, D));
        if (!pointOnSegment(p1, l1)) {
            swap(p1.y, p2.y);
        }
        if (p1 == p2) {
            return {3, p1, p2};
        } else {
            return {2, p1, p2};
        }
    }
    auto cp1 = cross(s2 - s1, e2 - s1);
    auto cp2 = cross(s2 - e1, e2 - e1);
    auto cp3 = cross(s1 - s2, e1 - s2);
    auto cp4 = cross(s1 - e2, e1 - e2);
    if (sign(cp1 * cp2) == 1 || sign(cp3 * cp4) == 1) {
        return {0, {}, {}};
    }
    dot p = linInt(l1, l2);
    if (sign(cp1) != 0 && sign(cp2) != 0 && sign(cp3) != 0 && sign(cp4) != 0) {
        return {1, p, p};
    } else {
        return {3, p, p};
    }
}

bool fastSegIntTest(lin a,lin b) {
    auto [s1,e1] = a;
    auto [s2,e2] = b;
    auto A = max(s1.x,e1.x),AA = min(s1.x,e1.x);
    auto B = max(s1.y,e1.y),BB = min(s1.y,e1.y);
    auto C = max(s2.x,e2.x),CC = min(s2.x,e2.x);
    auto D = max(s2.y,e2.y),DD = min(s2.y,e2.y);

    bool flag_cross = (sign(cross(s1,s2,e1)) * sign(cross(s1,e1,e2))) == 1 &&
                      (sign(cross(s2,s1,e2)) * sign(cross(s2,e2,e1))) == 1;
    bool flag_onseg = onseg(a,s2) || onseg(a,e2) ||
                      onseg(a,s2) || onseg(a,e2);

    return A >= CC && B >= DD && C >= AA && D >= BB && (flag_cross || flag_onseg); 
}

void input(lin &a) {
    int fx,fy;
    cin >> fx >> fy;
    a.fi = dot(F(fx),F(fy));
    cin >> fx >> fy;
    a.se = dot(F(fx),F(fy));
    // cin >> a.fi.x >> a.fi.y >> a.se.x >> a.se.y;
}

dot mirr(lin m,dot a) {
    F A = m.fi.y - m.se.y;
    F B = m.se.x - m.fi.x;
    F C = m.fi.x * m.se.y - m.se.x * m.fi.y;
    F t = (A * a.x + B * a.y + C)/(A * A + B * B);
    t = -t;
    F xm = a.x + A * t;
    F ym = a.y + B * t;
    return { xm * F(2) - a.x , ym * F(2) - a.y};
}

bool onLeft(dot p,lin vec) {
    return cross(vec.se,p,vec.fi) > 0;
}

bool inside(lin vec1,lin vec2,dot p) {
    if(onLeft(vec2.se,vec1)) swap(vec1,vec2);
    // v1在v2的左边
    return onLeft(p,vec2) && ((!onLeft(p,vec1)) && !onseg(vec1,p));
}

lin ext(lin lc) {
    auto vec = dsc(lc.fi,lc.se);
    vec = {vec.fi * (1e9),vec.se * (1e9)};
    return {lc.se,dsc(lc.se,vec)};
}

bool ifAvil(lin vec1,lin vec2,lin seg) {
    if(inside(vec1, vec2, seg.fi) || inside(vec1, vec2, seg.se)) return 1;
    if(fastSegIntTest(ext(vec1), seg) || fastSegIntTest(ext(vec2), seg)) return 1;
    return 0;
    // if(fastSegIntTest( a, lin b))
}

bool ifFstl(lin vec1,lin vec2,lin seg) {
    if(inside(vec1, vec2, seg.fi) && inside(vec1, vec2, seg.se)) return 1;

}

dot convexInt(vector<dot> ch,lin a) {
    auto tmp = a;
    tmp = ext(tmp);
    swap(tmp.fi,tmp.se);
    tmp = ext(tmp);

    vector<dot> ints;
    for(int i = 0;i < ch.size();i ++) {
        auto [status,p1,p2] = segInt(tmp,lin(ch[i],ch[(i+1)%ch.size()] ) );
        if(status == 1) {
            return p1;
        } else if(status == 2) {
            ints.pb(p1);
            ints.pb(p2);
        } else if(status == 3) {
            return p1;
        }
    }
    // for(auto x:ch) {
    //     auto [status,p,_] = segInt(tmp,x);

    // }
}

void solve()
{
    array<lin,3> as;
    for(auto &t:as) input(t);
    vector<int> p(3);
    iota(all(p),1);
    do{
        for(auto tx:p) cout << tx << " ";cout << "\n====\n";

        lin m0 = as[p[0]-1];
        lin m1 = as[p[1]-1];
        lin m2 = as[p[2]-1];

        lin m2t = {mirr(m1,m2.fi),mirr(m1,m2.se)}; 

        // 规定直线m2t参数方程形式为 D = fi + k V; V = se-fi;

        // 先判断m0和m1的位置关系
        // 目标：m2t所在直线。可行域：排除：线段m2，包含：线段m1

        // auto m2_T = {mirr(m1,m2.fi),mirr(m1,m2.se)};         // 
        // auto m1_TT = {mirr(m0,m1.fi),mirr(m0,m1.se)};
        // auto m2_TT = {mirr(m0,m2.fi),mirr(m0,m2.se)};
        // auto 

        // cout << m0.fi.x << " " << m0.fi.y << " " << m0.se.x << " " << m0.se.y << "\n";
        // cout << m1.fi.x << " " << m1.fi.y << " " << m1.se.x << " " << m1.se.y << "\n";
        // cout << m2.fi.x << " " << m2.fi.y << " " << m2.se.x << " " << m2.se.y << "\n";

        // lin t1 = {m0.fi,m1.fi};
        // lin t2 = {m0.fi,m1.se};
        // lin t3 = {m0.se,m1.fi};
        // lin t4 = {m0.se,m1.se};

        // if( ifAvil(t1, t2, m2) ||
        //     ifAvil(t1, t3, m2) ||
        //     ifAvil(t1, t4, m2) ||
        //     ifAvil(t2, t3, m2) ||
        //     ifAvil(t2, t4, m2) ||
        //     ifAvil(t3, t4, m2)
        // ) {
        //     cout << "Yes\n";
        //     return;
        // }
        // if( 
        //     // ifAvil(t1, t2, m2) ||
        //     // ifAvil(t1, t3, m2) ||
        //     ifAvil(t1, t4, m2) ||
        //     ifAvil(t2, t3, m2)
        //     // ifAvil(t2, t4, m2) ||
        //     // ifAvil(t3, t4, m2)
        // ) {
        //     cout << ifAvil(t1, t4, m2) << " " << ifAvil(t2, t3, m2) << "\n";
        //     cout << "Yes\n";
        //     return;
        // }
    } while(next_permutation(all(p)));
    cout << "No\n";
}

void test_mirr(){ //passed
    lin a;
    input(a);
    dot p;
    int fx,fy;
    cin >> fx >> fy;
    p = {F(fx),F(fy)};
    // cin >> p.x >> p.y;
    auto [tx,ty] = mirr(a,p);
    cout << tx << " " << ty << "\n";
}

void test_ext(){
    lin a;
    input(a);
    auto [f,s] = ext(a);
    cout << f.x << " " << f.y << "\n";
    cout << s.x << " " << s.y << "\n";
    // lin a;
    // input(a);
    // auto [tx,ty] = ext(a).fi;
    // cout << tx << " " << ty << "\n";
}

signed main()
{
#ifdef FC
    freopen("G.A.in","r",stdin);
    freopen("G.A.ptc","w",stdout);
#endif
#ifndef FC
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.tie(0);
#endif
    int T = 1;
    cin >> T;
    while(T--) solve();
    // while(T--) test_ext();
    return 0;
}
