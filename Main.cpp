#include <bits/stdc++.h>

using namespace std;
//structure of a point
struct point
{
    double x;
    double y;
};
//Interval
struct interval
{

    double bottom;
    double top;
    double bound()
    {
        return top - bottom;
    }
    bool operator==(const interval &i)
    {
        return top == i.top && bottom == i.bottom;
    }
    bool operator<(const interval &i)
    {
        return (bottom < i.bottom || (bottom == i.bottom && top < i.top));
    }
};
//Line segment
struct lineSegment
{
    double x1;
    interval span;
};
//Rectangle
struct rectangle
{
    interval x_interval;
    interval y_interval;
};
//Edgetype
enum class edgeType
{
    Left,
    Right,
    Top,
    Bottom
};
//Edge
struct edge
{
    lineSegment line;
    edgeType side;
    bool operator<(const edge &i)
    {
        return (line.x1 < i.line.x1) || (line.x1 == i.line.x1 && side < i.side);
    }
};
//lru
enum class lru
{
    Left,
    Right,
    Undef
};
//ctree
struct cTree
{
    double x;
    lru side;
    cTree *lson;
    cTree *rson;
};
//stripe
struct stripe
{
    interval x_interval;
    interval y_interval;
    struct cTree *tree;
    double x_measure;
};
//Stripe Return
struct stripeRet
{
    vector<interval> l;
    vector<interval> r;
    vector<double> p;
    vector<stripe> s;
};
//partition of y intervals
vector<interval> partition(vector<double> &p)
{
    vector<interval> result;
    sort(p.begin(), p.end());
    for (auto i = p.begin() + 1; i != p.end(); i++)
    {
        result.push_back({*(i - 1), *i});
    }
    return result;
}

void blacken(vector<stripe> &S, vector<interval> &J)
{
    for (auto &s : S)
    {
        for (auto i : J)
        {
            if (s.y_interval.bottom >= i.bottom && s.y_interval.top <= i.top)
            {
                s.x_measure = s.x_interval.top - s.x_interval.bottom;
                s.tree = nullptr;
            }
        }
    }
}

vector<stripe> copy(vector<stripe> &s, vector<double> &p, interval x_int)
{
    vector<stripe> s_new;
    for (auto y_int : partition(p))
    {
        stripe Stripe = {x_int, y_int, nullptr, 0};
        s_new.push_back(Stripe);
    }
    for (auto &Stripe_new : s_new)
    {
        for (auto &Stripe : s)
        {
            if ((Stripe.y_interval.top >= Stripe_new.y_interval.top) && (Stripe.y_interval.bottom <= Stripe_new.y_interval.bottom))
            {
                Stripe_new.x_measure = Stripe.x_measure;
                Stripe_new.tree = Stripe.tree;
            }
        }
    }
    return s_new;
}

vector<stripe> concat(vector<stripe> &S1, vector<stripe> &S2, vector<double> &P, interval x_int)
{
    vector<stripe> S;
    for (auto y_int : partition(P))
    {
        S.push_back({x_int, y_int, nullptr, 0});
    }
    for (auto &s : S)
    {
        bool flag = false;
        for (auto s1 : S1)
        {
            if (s1.y_interval == s.y_interval)
            {
                for (auto s2 : S2)
                {
                    if (s2.y_interval == s.y_interval)
                    {
                        s.x_measure = s1.x_measure + s2.x_measure;
                        if (s1.tree != nullptr && s2.tree != nullptr)
                            s.tree = new cTree{s1.x_interval.top, lru::Undef, s1.tree, s2.tree};
                        else if (s1.tree != nullptr && s2.tree == nullptr)
                            s.tree = s1.tree;
                        else if (s1.tree == nullptr && s2.tree != nullptr)
                            s.tree = s2.tree;
                        else if (s1.tree == nullptr && s2.tree == nullptr)
                            s.tree = nullptr;
                        flag = true;
                        break;
                    }
                }
            }
            if (flag)
                break;
        }
    }
    return S;
}
/*Computes the set of stripes
Input: Set of vertical edges,x_interval 
Output:<L ,R ,P S>*/
stripeRet stripes(vector<edge> &V, interval x_int)
{
    vector<interval> l;
    vector<interval> r;
    vector<double> p;
    vector<stripe> s;
    if (V.size() == 1)
    {
        auto v = V[0];
        if (v.side == edgeType::Left)
        {
            l.push_back(v.line.span);
        }
        else
        {
            r.push_back(v.line.span);
        }
        p = {-INFINITY, v.line.span.bottom, v.line.span.top, INFINITY};
        for (auto y_int : partition(p))
        {
            stripe Stripe = {x_int, y_int, nullptr, 0};
            s.push_back(Stripe);
        }
        for (auto &Stripe : s)
        {
            if (Stripe.y_interval == v.line.span)
            {
                if (v.side == edgeType::Left)
                {
                    Stripe.x_measure = x_int.top - v.line.x1;
                    Stripe.tree = new cTree{v.line.x1, lru::Left, nullptr, nullptr};
                }
                else
                {
                    Stripe.x_measure = v.line.x1 - x_int.bottom;
                    Stripe.tree = new cTree{v.line.x1, lru::Right, nullptr, nullptr};
                }
            }
        }
    }
    else
    {
        //Divide
        nth_element(V.begin(), V.begin() + V.size() / 2, V.end());
        vector<edge> V1, V2;
        for (auto itr = V.begin(); itr != V.begin() + V.size() / 2; itr++)
        {
            V1.push_back(*itr);
        }

        for (auto itr = V.begin() + V.size() / 2; itr != V.end(); itr++)
        {
            V2.push_back(*itr);
        }
        auto xm = (V2[0].line.x1 + max_element(V1.begin(), V1.end())->line.x1) / 2;
        //Conquer
        auto P1 = stripes(V1, {x_int.bottom, xm});
        auto P2 = stripes(V2, {xm, x_int.top});
        //Merge
        vector<interval> LR;
        set_intersection(P1.l.begin(), P1.l.end(), P2.r.begin(), P2.r.end(), back_inserter(LR));
        vector<interval> LTemp;
        set_difference(P1.l.begin(), P1.l.end(), LR.begin(), LR.end(), back_inserter(LTemp));
        set_union(LTemp.begin(), LTemp.end(), P2.l.begin(), P2.l.end(), back_inserter(l));
        vector<interval> RTemp;
        set_difference(P2.r.begin(), P2.r.end(), LR.begin(), LR.end(), back_inserter(RTemp));
        set_union(P1.r.begin(), P1.r.end(), RTemp.begin(), RTemp.end(), back_inserter(r));
        set_union(P1.p.begin(), P1.p.end(), P2.p.begin(), P2.p.end(), back_inserter(p));
        vector<stripe> Sleft = copy(P1.s, p, {x_int.bottom, xm});
        vector<stripe> Sright = copy(P2.s, p, {xm, x_int.top});

        blacken(Sleft, RTemp);
        blacken(Sright, LTemp);
        s = concat(Sleft, Sright, p, x_int);
    }

    stripeRet RETURNS = {l, r, p, s};
    return RETURNS;
}
/*Returns the measure of the set of rectangles
Input :Set of stripes 
Output :measure*/
double measure(vector<stripe> &S)
{
    double sum = 0;
    for (auto s : S)
    {
        if (s.y_interval.top - s.y_interval.bottom == INFINITY)
        {
            sum += 0;
        }
        else
        {
            sum += s.x_measure * (s.y_interval.top - s.y_interval.bottom);
        }
    }
    return sum;
}
/*DFS of the contracte segment Tree
Input: Ctree and xmeasure
*/
std::vector<double> leaf;
void dfs(cTree *tree, double x1)
{
    if (not tree)
        return;
    if (tree->side != lru::Undef)
    {
        if (leaf.empty())
            if (tree->side == lru::Left)
                leaf.push_back(x1);
        leaf.push_back(tree->x);
    }
    cTree *ldaughter = tree->lson, *rdaughter = tree->rson;
    dfs(ldaughter, x1);
    dfs(rdaughter, x1);
}

std::vector<std::tuple<double, double, double, bool>> cont;
/*
Input:x_interval,stripe,yccord,flag 
*/
void query(interval xi, stripe s, double ycoord, bool flag)
{
    cTree *tree = s.tree;
    leaf.clear();
    dfs(tree, xi.bottom);

    int n = int(leaf.size());
    if (n == 0)
        cont.push_back({ycoord, xi.bottom, xi.top, flag});
    if (n % 2)
        leaf.push_back(xi.top), ++n;

    for (int i = 0; i < n; i += 2)
    {
        double &cur = leaf.at(i);
        double &nxt = leaf.at(i + 1);
        if (cur < nxt)
        {
            if (nxt <= xi.bottom or cur >= xi.top)
                continue;
            cont.push_back({ycoord, max(cur, xi.bottom), min(nxt, xi.top), flag});
        }
    }
}
/*
Calculates the contour of the set of rectangles
Input set of stripes and horizontal edges
*/
double contour(vector<stripe> &strips, vector<edge> &hrx)
{
    double out = 0;
    auto it = hrx.begin();
    auto jt = strips.begin();
    while (it != hrx.end() && jt != strips.end())
    {
        if (it->side == edgeType::Bottom)
        {
            if (jt->y_interval.top < it->line.x1)
                ++jt;
            else if (jt->y_interval.top == it->line.x1)
                query(it->line.span, *jt, it->line.x1, 0), ++it;
            else
                ++it;
        }
        else
        {
            if (jt->y_interval.bottom < it->line.x1)
                ++jt;
            else if (jt->y_interval.bottom == it->line.x1)
                query(it->line.span, *jt, it->line.x1, 1), ++it;
            else
                ++it;
        }
    }

    ofstream file;
    file.open("contour_edges.txt");
    int n = int(cont.size());
    sort(cont.begin(), cont.end());
    vector<tuple<double, double>> ep;

    auto &[y1, xl1, xr1, f1] = cont.front();
    for (auto &[y2, xl2, xr2, f2] : cont)
    {
        if (y2 == y1 and f2 == f1 and xr1 >= xl2)
        {
            xr1 = max(xr1, xr2);
            continue;
        }
        file << xl1 << ' ' << y1 << ' ' << ' ' << xr1 - xl1 << ' ' << 0 << '\n';
        ep.push_back({xl1, y1});
        ep.push_back({xr1, y1});
        out += (xr1 - xl1);
        y1 = y2, xl1 = xl2, xr1 = xr2, f1 = f2;
    }

    file << xl1 << ' ' << y1 << ' ' << ' ' << xr1 - xl1 << ' ' << 0 << '\n';
    ep.push_back({xl1, y1});
    ep.push_back({xr1, y1});
    out += (xr1 - xl1);

    std::sort(ep.begin(), ep.end());
    n = int(ep.size());
    for (int i = 0; i + 1 < n; ++i)
    {
        auto [cx, cy] = ep.at(i);
        auto [nx, ny] = ep.at(i + 1);
        if (cx == nx)
        {
            if (cy == ny)
                continue;
            else
            {
                out += (ny - cy);
                file << cx << ' ' << cy << ' ' << 0 << ' ' << ny - cy << '\n';
                ++i;
            }
        }
    }
    file.close();
    return out;
}

/*
Input: set of rectangles
Output: Returns a pair of measure and Contour of the set of rectangles
*/
pair<double, double> RectangleDac(vector<rectangle> &rect)
{
    vector<edge> vrx;
    for (auto r : rect)
    {
        edge e = {{r.x_interval.bottom, r.y_interval}, edgeType::Left};
        vrx.push_back(e);
        e = {{r.x_interval.top, r.y_interval}, edgeType::Right};
        vrx.push_back(e);
    }
    auto res = stripes(vrx, {-INFINITY, INFINITY});
    vector<edge> hrx;
    for (auto r : rect)
    {
        edge e = {{r.y_interval.bottom, r.x_interval}, edgeType::Bottom};
        hrx.push_back(e);
        e = {{r.y_interval.top, r.x_interval}, edgeType::Top};
        hrx.push_back(e);
    }
    std::sort(hrx.begin(), hrx.end());
    pair<double, double> out = {measure(res.s), contour(res.s, hrx)};
    return out;
}
//Driver Function
int main()
{
    // cout << "N of Rectangles:" << endl;
    int n;
    vector<rectangle> Rect;
    // cout << "Enter Points in format x1 x2 y1 y2 (Bottom left and Top Right Points):" << endl;
    freopen("input.txt", "r", stdin);
    cin >> n;
    for (int i = 0; i < n; i++)
    {
        double x1, y1, x2, y2;
        cin >> x1 >> x2 >> y1 >> y2;
        rectangle rect = {{x1, x2}, {y1, y2}};
        Rect.push_back(rect);
    }
    // visualise(x1,x2,y1,y2);
    pair<double, double> S = RectangleDac(Rect);
    cout << S.first << " " << S.second << endl;
}