#include "frameworks.hpp"
#include <algorithm>

// #include <bits/stdc++.h>
// using namespace std;

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

// int main()
// {
//     freopen("G.A.in","w",stdout);
    
// }


using namespace std;
// using namespace geo;

// T 暴击 (x) 爆伤 (T - x)  c 2c  (0.05 + xc) * (0.5+(T-X)*2c) + (1-(0.05 + xc) )
// 

// import 

lin reExt(lin x) {
    // dl( "for ext:" << x);
    x = purExt(x);
    swap(x.fi,x.se);
    x = purExt(x);
    // ext(x);
    // dl(x);
    return x;
}
bool onLine(dot a, dot b, dot c) {
    return sign(cross(b, a, c)) == 0;
}

F pdot(dot a,dot b) {
    return a.x * b.x + a.y * b.y;
}

int contains(dot p, vector<dot> A) {
    int n = A.size();
    bool in = false;
    for (int i = 0; i < n; i++) {
        dot a = A[i] - p, b = A[(i + 1) % n] - p;
        if (a.y > b.y) {
            swap(a, b);
        }
        if (a.y <= 0 && 0 < b.y && cross(a, b) < 0) {
            in = !in;
        }
        if (cross(a, b) == 0 && pdot(a, b) <= 0) {
            return 1;
        }
    }
    return in ? 2 : 0;
 }

bool fastSegIntTest(lin a,lin b) {
    // dl("on FSI:" << a << " " << b << "\n");
    auto [s1,e1] = a;
    auto [s2,e2] = b;
    auto A = max(s1.x,e1.x),AA = min(s1.x,e1.x);
    auto B = max(s1.y,e1.y),BB = min(s1.y,e1.y);
    auto C = max(s2.x,e2.x),CC = min(s2.x,e2.x);
    auto D = max(s2.y,e2.y),DD = min(s2.y,e2.y);

    bool flag_cross = (sign(cross(s1,s2,e1)) * sign(cross(s1,e1,e2))) == 1 &&
                      (sign(cross(s2,s1,e2)) * sign(cross(s2,e2,e1))) == 1;
    bool flag_onseg = onseg(a,s2) || onseg(a,e2) || 
                      onseg(b,s1) || onseg(b,e1);

    bool flag_final = A >= CC && B >= DD && C >= AA && D >= BB && (flag_cross || flag_onseg); 

    // dl("FSIT END flag_onseg:" << flag_onseg);
    return flag_final;
}


bool intCheck(const calipers & ax, lin tx) {
    auto & tp =  ax.p;
    int n = tp.size();
    for(int i = 0;i < n;i ++) {
        if(fastSegIntTest({tp[i],tp[(i+1)%n]},tx)) return 1;
    }
    return 0;
};

signed main(){
    int cps;
    
    cin >> cps;
    string stus;
    cin >> stus;
    int seed;
    cin >> seed;
    dl(seed);
    generator gc(seed);
    auto randDot = [&]() -> dot {
        return {gc.randi(0,cps),gc.randi(0,cps)};
    };
    auto randLin = [&]() -> lin {
        auto tx = randDot();
        auto ty = randDot();
        while(ty == tx) ty = randDot();
        return {tx,ty};
    };

    auto len = [](lin x) -> int {
        auto [tx,ty] = x.se - x.fi;
        return (double)(tx * tx + ty * ty);
    };

    auto randLinFlag = [&](int minLen) -> lin{
        auto tl = randLin();
        while(len(tl) < minLen) tl = randLin();
        return tl;
    };

    auto randNDot = [&](calipers ax) -> dot{
        auto tl = randDot();
        while(contains(tl, ax.p) != 0)tl = randDot();
        return tl;
    };

    auto randNLin = [&](calipers ax) -> lin {


        auto tx = randNDot(ax);
        auto ty = randNDot(ax);
        while(ty == tx || intCheck(ax,{tx,ty})) ty = randNDot(ax);
        return {tx,ty};
    };

    int falsdkjfa = gc.randi(1,100);
    if(falsdkjfa >= 95) stus = "NO";
    else stus = "YES";

    if(stus == "YES") {
        auto lt = randLinFlag(cps/2);
        auto m1 = randLinFlag(0);
        while(fastSegIntTest(m1,lt) == 0) m1 = randLinFlag(0);
        auto m2 = randLinFlag(0);
        while(
            fastSegIntTest(m2, lt) == 0 ||
            (
                fastSegIntTest(m1, reExt(m2)) != 0 ||
                fastSegIntTest(m2, reExt(m1)) != 0
            ) || 
            fastSegIntTest(m1, m2)
        ) m2 = randLinFlag(0);

        dl(fastSegIntTest(m1, reExt(m2)));
        dl(fastSegIntTest(m2, reExt(m1)));
        dl(" ");
        dl((fastSegIntTest(m2, lt) == 0 ||
            (
                fastSegIntTest(m1, reExt(m2)) != 0 ||
                fastSegIntTest(m2, reExt(m1)) != 0
            ) || 
            fastSegIntTest(m1, m2)));
        
        dl(onLine(m1.fi,m1.se,reExt(m1).fi));
        dl(onLine(m1.fi,m1.se,reExt(m1).se));
        dl(onLine(m2.fi,m2.se,reExt(m2).fi));
        dl(onLine(m2.fi,m2.se,reExt(m2).se));

        dl("THIS");

        dl(onLine(m2.fi, m2.se , m1.fi));
        // dl(onseg(m1,purExt(m2)));
        dl(onseg(purExt(m2),m1.se));
        dl(fastSegIntTest(purExt(m2),m1));
        dl(onLine(m2.fi, m2.se , m1.se));
        dl(onLine(m1.fi, m1.se , m2.fi));
        dl(onLine(m1.fi, m1.se , m2.se));

        dl(" ");
        auto [_,p1,__] = segInt(lt, m1);
        auto [___,p2,____] = segInt(lt, m2);
        lt = ext(lin({p1,p2}));

        calipers cal({
            m1.fi,m1.se,
            m2.fi,m2.se}
        );
        auto m3 = randNLin(cal);
        while(fastSegIntTest(m3, lt) != 0) m3 = randNLin(cal);
        auto m4 = randNLin(cal);
        while(fastSegIntTest(m4, lt) == 0 || fastSegIntTest(m3, m4)) m4 = randNLin(cal);

        // cout << m1.fi.x << " " << m1.fi.y 
        auto print = [&](lin cx) -> void {
            cout << cx.fi.x << " " << cx.fi.y << " " << cx.se.x << " " << cx.se.y << "\n";
        };
        // cout << 1 << "\n";
        print(m1);
        print(m2);
        print(m3);
        print(m4);
    } else if(stus == "NO") {
    // stus == "COMMON";
        
        auto lt = randLinFlag(cps/2);
        auto m1 = randLinFlag(0);
        while(fastSegIntTest(m1,lt) == 0) m1 = randLinFlag(0);
        auto m2 = randLinFlag(0);
        while(
            fastSegIntTest(m2, lt) == 0 ||
            (
                fastSegIntTest(m1, reExt(m2)) != 0 ||
                fastSegIntTest(m2, reExt(m1)) != 0
            ) || 
            fastSegIntTest(m1, m2)
        ) m2 = randLinFlag(0);

        dl(fastSegIntTest(m1, reExt(m2)));
        dl(fastSegIntTest(m2, reExt(m1)));
        dl(" ");
        dl((fastSegIntTest(m2, lt) == 0 ||
            (
                fastSegIntTest(m1, reExt(m2)) != 0 ||
                fastSegIntTest(m2, reExt(m1)) != 0
            ) || 
            fastSegIntTest(m1, m2)));
        
        dl(onLine(m1.fi,m1.se,reExt(m1).fi));
        dl(onLine(m1.fi,m1.se,reExt(m1).se));
        dl(onLine(m2.fi,m2.se,reExt(m2).fi));
        dl(onLine(m2.fi,m2.se,reExt(m2).se));

        dl("THIS");

        dl(onLine(m2.fi, m2.se , m1.fi));
        // dl(onseg(m1,purExt(m2)));
        dl(onseg(purExt(m2),m1.se));
        dl(fastSegIntTest(purExt(m2),m1));
        dl(onLine(m2.fi, m2.se , m1.se));
        dl(onLine(m1.fi, m1.se , m2.fi));
        dl(onLine(m1.fi, m1.se , m2.se));

        dl(" ");
        auto [_,p1,__] = segInt(lt, m1);
        auto [___,p2,____] = segInt(lt, m2);
        lt = ext(lin({p1,p2}));

        calipers cal({
            m1.fi,m1.se,
            m2.fi,m2.se}
        );
        auto m3 = randNLin(cal);
        while(
            fastSegIntTest(m3, lt) != 0 || (
                fastSegIntTest(m3, reExt(m1)) != 0 ||
                fastSegIntTest(m1, reExt(m3)) != 0
            )
        ) m3 = randNLin(cal);
        dl(m1);
        dl(m3);
        // return 0;
        auto checkIfBlock = [&](dot p) -> bool {
            dl("check:" << p);
            F eps = F(1) / 1e3;
            F t = F(0);
            auto vec = m1.se-m1.fi;
            while(t <= F(1)) {
                auto fr = m1.fi + t * vec;
                dl("fr:" << fr << " >< if blocker by m3: "<< (fastSegIntTest({fr,p}, m3) == 1));
                if(
                    fastSegIntTest({fr,p}, m3) == 0
                ) return 0;
                t += eps;
            }
            return 1;
        };

        // auto m4 = randNLin(cal);
        // while(fastSegIntTest(m4, lt) == 0 || fastSegIntTest(m3, m4)) m4 = randNLin(cal);
        dot m4d1 = {0,0},m4d2 = {0,0};
        lin m4 = {m4d1,m4d2};
        do{
            // dl("m4d1 fst failed\n");
            dl((fastSegIntTest(lin(m4), m3) != 0) << " " << intCheck(cal, m4));
            dl(intCheck(cal, m4));
            m4d1 = randNDot(cal);
            while(!checkIfBlock(m4d1))  m4d1 = randNDot(cal);
            dl("SELECTED m4d1 : " << m4d1);
            m4d2 = randNDot(cal);
            while(!checkIfBlock(m4d2) || m4d2 == m4d1)  m4d2 = randNDot(cal);
            dl("SELECTED m4d2 : " << m4d2);
            m4 = {m4d1,m4d2};
        } while (
            // fastSegIntTest(lin(m4), lt) == 0 || 
            fastSegIntTest(lin(m4), m3) != 0 ||
            intCheck(cal, m4)
        );

        // cout << m1.fi.x << " " << m1.fi.y 
        auto print = [&](lin cx) -> void {
            cout << cx.fi.x << " " << cx.fi.y << " " << cx.se.x << " " << cx.se.y << "\n";
        };
        // cout << 1 << "\n";
        print(m1);
        print(m2);
        print(m3);
        print(m4);
    } else {

    }

}