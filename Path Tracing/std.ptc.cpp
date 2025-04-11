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
        assert(abs(num) <= 1e18);
        assert(abs(den) <= 1e18);
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
 
using F = Frac<__int128>;

F max(F a,F b) {
    return a < b ? b : a;
}

F min(F a,F b) {
    return a < b ? a : b;
}

F abs(F a) {
    return a < F(0) ? -a : a;
}

namespace std {
template<typename T>
class Interval {
public:
    // 构造函数，默认为闭区间
    Interval(T a, T b, bool a_closed = true, bool b_closed = true)
        : l(a), r(b), lc(a_closed), rc(b_closed) {}

    // 判断区间是否非空（当 l < r 或 l==r 且两端都闭合时有效）
    bool valid() const {
        return (l < r) || (l == r && lc && rc);
    }

    // 判断两个区间是否有交集（考虑“接触”情况：若接触端点至少一侧闭合，则认为交集存在）
    bool overlaps(const Interval &o) const {
        if(r < o.l || (r == o.l && (!rc || !o.lc))) return false;
        if(o.r < l || (o.r == l && (!o.rc || !lc))) return false;
        return true;
    }

    // 交集运算：返回两个区间的交集（为空时返回一个无效区间）
    Interval intersect(const Interval &o) const {
        T L = max(l, o.l);
        T R = min(r, o.r);
        bool Lc = (l == o.l ? (lc && o.lc) : (l < o.l ? o.lc : lc));
        bool Rc = (r == o.r ? (rc && o.rc) : (r < o.r ? rc : o.rc));
        Interval ret(L, R, Lc, Rc);
        // 交集无效，则返回一个无效区间（这里用闭合全 false 表示无效）
        return ret.valid() ? ret : Interval(L, L, false, false);
    }

    // 并集运算：若两个区间重叠或相邻（接触时至少一侧闭合），返回合并后的单一区间，
    // 否则返回两个区间
    vector<Interval> unite(const Interval &o) const {
        // 判断“相邻”，如 a.r == b.l 且至少一侧闭合
        bool adjacent = (r == o.l && (rc || o.lc)) || (o.r == l && (o.rc || lc));
        if(!overlaps(o) && !adjacent) return {*this, o};
        T L = (l < o.l ? l : o.l);
        T R = (r > o.r ? r : o.r);
        bool Lc = (l < o.l ? lc : (l > o.l ? o.lc : (lc || o.lc)));
        bool Rc = (r > o.r ? rc : (r < o.r ? o.rc : (rc || o.rc)));
        return { Interval(L, R, Lc, Rc) };
    }

    // 差集运算：计算当前区间减去 o 后剩余的部分，最多返回两个区间
    vector<Interval> subtract(const Interval &o) const {
        // 若没有重叠，则整个区间保留
        if(!overlaps(o)) return {*this};
        Interval inter = intersect(o);
        vector<Interval> res;
        // 左侧剩余区间存在条件：当前区间左端小于交集左端，
        // 或相等但当前左端闭合而交集左端不闭合
        if(l < inter.l || (l == inter.l && lc && !inter.lc))
            res.push_back(Interval(l, inter.l, lc, false));
        // 右侧剩余区间
        if(r > inter.r || (r == inter.r && rc && !inter.rc))
            res.push_back(Interval(inter.r, r, false, rc));
        return res;
    }

    // 方便输出，格式如 [l,r) 表示区间的开闭情况
    void print() const {
        cout << (lc ? "[" : "(") << l << "," << r << (rc ? "]" : ")");
    }

private:
    T l, r;       // 区间左右端点
    bool lc, rc;  // 分别表示左、右端点是否闭合（true：闭区间，false：开区间）
};

template <typename TYPE_NAME,std::size_t N>
class Vector {
public:
    std::array<TYPE_NAME, N> data;

    TYPE_NAME & first = data[0];
    TYPE_NAME & second = data[1];

    // 默认构造函数：初始化为 0
    Vector() {
        data.fill(0.0);
    }

    Vector(TYPE_NAME first,TYPE_NAME second) {
        data[0] = first;
        data[1] = second;
    }

    // 使用 initializer_list 构造
    Vector(std::initializer_list<TYPE_NAME> list) {
        assert(list.size() == N && "Initializer list size must equal vector dimension.");
        std::copy(list.begin(), list.end(), data.begin());
    }

    //拷贝构造
    Vector(const Vector &other) {
        std::copy(other.data.begin(), other.data.end(), data.begin());
    }

    // 移动构造函数
    Vector(Vector &&other) noexcept {
        std::move(other.data.begin(), other.data.end(), data.begin());
    }
    // 移动赋值运算符
    Vector& operator=(Vector &&other) noexcept {
        if (this != &other) {
            std::move(other.data.begin(), other.data.end(), data.begin());
        }
        return *this;
    }

    // 下标运算符
    TYPE_NAME& operator[](std::size_t index) {
        assert(index < N);
        return data[index];
    }
    const TYPE_NAME& operator[](std::size_t index) const {
        assert(index < N);
        return data[index];
    }



    // 向量加法
    Vector operator+(const Vector &other) const {
        Vector result;
        for (std::size_t i = 0; i < N; ++i)
            result.data[i] = data[i] + other.data[i];
        return result;
    }

    // 向量减法
    Vector operator-(const Vector &other) const {
        Vector result;
        for (std::size_t i = 0; i < N; ++i)
            result.data[i] = data[i] - other.data[i];
        return result;
    }

    // 数乘：向量 * 标量（右侧为标量）
    Vector operator*(TYPE_NAME scalar) const {
        Vector result;
        for (std::size_t i = 0; i < N; ++i)
            result.data[i] = data[i] * scalar;
        return result;
    }

    // 相等
    bool operator==(const Vector &other) const {
        for (std::size_t i = 0; i < N; ++i) {
            if (data[i] != other.data[i]) return false;
        }
        return true;
    }

    // 数乘：标量 * 向量（友元函数）
    friend Vector operator*(TYPE_NAME scalar, const Vector &vec) {
        return vec * scalar;
    }

    // 点乘（内积）
    TYPE_NAME dot(const Vector &other) const {
        TYPE_NAME result = 0.0;
        for (std::size_t i = 0; i < N; ++i)
            result += data[i] * other.data[i];
        return result;
    }

    // 求向量的欧几里得范数
    TYPE_NAME norm() const {
        return std::sqrt(this->dot(*this));
    }

    // 输出运算符重载
    friend std::ostream& operator<<(std::ostream &os, const Vector &vec) {
        os << "[";
        for (std::size_t i = 0; i < N; ++i) {
            os << vec.data[i];
            if (i != N - 1)
                os << ", ";
        }
        os << "]";
        return os;
    }
};

template <std::size_t I, typename TYPE_NAME,std::size_t N>
constexpr TYPE_NAME& get(Vector<TYPE_NAME,N>& vec) noexcept {
    static_assert(I < N, "Index out of bounds in Vector::get");
    return vec.data[I];
}

template <std::size_t I,typename TYPE_NAME, std::size_t N>
constexpr const TYPE_NAME& get(const Vector<TYPE_NAME,N>& vec) noexcept {
    static_assert(I < N, "Index out of bounds in Vector::get (const)");
    return vec.data[I];
}

template <std::size_t I,typename TYPE_NAME, std::size_t N>
constexpr TYPE_NAME&& get(Vector<TYPE_NAME,N>&& vec) noexcept {
    static_assert(I < N, "Index out of bounds in Vector::get (rvalue)");
    return std::move(vec.data[I]);
}

    template <typename TYPE_NAME,std::size_t N>
    struct tuple_size<Vector<TYPE_NAME,N>> : std::integral_constant<std::size_t, N> { };

    template <typename TYPE_NAME,std::size_t I, std::size_t N>
    struct tuple_element<I, Vector<TYPE_NAME,N>> {
        using type = TYPE_NAME;
    };
}

// using dot = pair<F,F>;
// using lin = pair<dot,dot>;

using dot = Vector<F, 2>;
using lin = pair<dot, dot>;

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

// dot dsc(dot a,dot b) {
//     return {a.x - b.x,a.y - b.y};
// }

// dot add(dot a,dot b) {
//     return {a.x + b.x,a.y + b.y};
// }

F cross(dot p1,dot p2,dot p0) {
    // return cross(dsc(p1,p0),dsc(p2,p0));
    cross(p1-p0,p2-p0);
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
    // F f = cross( dsc(b.se,b.fi),dsc(a.fi,b.fi) ) / cross(dsc(b.se,b.fi),dsc(a.fi,a.se));


    F f = cross( b.se-b.fi,a.fi-b.fi ) / cross(b.se-b.fi,a.fi-a.se);
    // auto [tx,ty] = dsc(a.se,a.fi);
    // auto tmp = a.se;
    auto [tx,ty] = a.se-a.fi;

    return a.fi + dot({tx * f,ty * f});
}

tuple<int, dot, dot> segInt(lin l1, lin l2) {
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
        dot p1(max(AA, CC), max(BB, DD));
        dot p2(min(A, C), min(B, D));
        if (!onseg(l1,p1)) {
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
    // auto vec = dsc(lc.fi,lc.se);
    auto vec =lc.fi-lc.se;
    // int esp = min()
    F esp(1e18);
    if(vec.fi != F(0)) esp = min(F(esp),F(2e10)/abs(vec.fi));
    if(vec.se != F(0)) esp = min(F(esp),F(2e10)/abs(vec.se));
    assert(esp != 1e18);
    vec = dot({vec.fi * esp,vec.se * esp});
    return {lc.se,lc.se-vec};
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
        if(status != )
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
