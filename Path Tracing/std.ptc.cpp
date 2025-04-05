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
    Frac() : Frac(0, 1) {}
    Frac(T num_) : Frac(num_, 1) {}
    explicit operator double() const {
        return 1. * num / den;
    }
    Frac &operator+=(const Frac &rhs) {
        num = num * rhs.den + rhs.num * den;
        den *= rhs.den;
        return *this;
    }
    Frac &operator-=(const Frac &rhs) {
        num = num * rhs.den - rhs.num * den;
        den *= rhs.den;
        return *this;
    }
    Frac &operator*=(const Frac &rhs) {
        num *= rhs.num;
        den *= rhs.den;
        return *this;
    }
    Frac &operator/=(const Frac &rhs) {
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

bool segIntTest(lin a,lin b) {
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

dot linInt(lin a,lin b) {
    F f = cross( dsc(b.se,b.fi),dsc(a.fi,b.fi) ) / cross(dsc(b.se,b.fi),dsc(a.fi,a.se));
    auto [tx,ty] = dsc(a.se,a.fi);
    return add(a.fi,{tx * f,ty * f});
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
    F C = m.fi.x * m.se.y - m.se.x + m.fi.y;
    F t = (A * a.x + B * a.y + C)/(A * A + B * B);
    t = -t;
    F xm = a.x + A * t;
    F ym = a.y + B * t;
    return {xm*F(2)-a.x,ym * F(2)-a.y};
}

bool onLeft(dot p,lin vec) {
    return cross(vec.se,p,vec.fi) > 0;
}

bool inside(lin vec1,lin vec2,dot p) {
    return onLeft(p, vec1) ^ onLeft(p, vec2);
}

lin ext(lin lc) {
    auto vec = dsc(lc.fi,lc.se);
    vec = {vec.fi * (1e9),vec.se * (1e9)};
    return {dsc(lc.fi,vec),add(lc.se,vec)};
}

bool ifAvil(lin vec1,lin vec2,lin seg) {
    if(inside(vec1, vec2, seg.fi) || inside(vec1,  vec2, seg.se)) return 1;
    if(segIntTest(ext(vec1), seg) || segIntTest(ext(vec2), seg)) return 1;
    return 0;
    // if(segIntTest( a, lin b))
}

void solve()
{
    array<lin,3> as;
    for(auto &t:as) input(t);
    vector<int> p(3);
    iota(all(p),1);
    do{
        lin m0 = as[p[0]-1];
        lin m1 = as[p[1]-1];
        lin m2 = as[p[2]-1];

        m2 = {mirr(m1,m2.fi),mirr(m1,m2.se)};
        m1 = {mirr(m0,m1.fi),mirr(m0,m1.se)};
        m2 = {mirr(m0,m2.fi),mirr(m0,m2.se)};

        lin t1 = {m0.fi,m1.fi};
        lin t2 = {m0.fi,m1.se};
        lin t3 = {m0.se,m1.fi};
        lin t4 = {m0.se,m1.se};

        if(
            ifAvil(t1, t2, m2) ||
            ifAvil(t1, t3, m2) ||
            ifAvil(t1, t4, m2) ||
            ifAvil(t2, t3, m2) ||
            ifAvil(t2, t4, m2) ||
            ifAvil(t3, t4, m2)
        ) {
            cout << "Yes\n";
        } else cout << "No\n";

    } while(next_permutation(all(p)));
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
    return 0;
}
