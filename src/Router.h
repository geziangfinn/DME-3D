#pragma once
#include "global.h"
using namespace std;
#define c_w 0.0002
#define c_v 15.48
#define r_w 0.0001
#define r_v 0.035
#define c_w_standard c_w/1000000000000000 // 单位F/nm
#define c_v_standard c_v/1000000000000000 // 单位F
#define l_w 0.0000000000000018 // 单位H/nm
#define l_v 0.000000000013827 // 单位H
#define c_constraint 300 // 单位 fF
#define LINEAR_DELAY 0
#define ELMORE_DELAY 1
const double eps = 1e-3;// 1e-4 seems to large
const double skewModifyStep = 1;

class PointPair {
    // V_{i,k,n}
public:
    double x1,x2,y1,y2;  // from:(x1,y1), to: (x2,y2)
    PointPair(){}
    PointPair(double _x1,double _x2,double _y1,double _y2) : x1(_x1), y1(_y1), x2(_x2),y2(_y2) {}
    bool operator==(PointPair const& var) const { return (x1 == var.x1 && y1 == var.y1 && x2 == var.x2 && y2 == var.y2); }
};

struct PointPairHasher {
    std::size_t operator()(const PointPair& var) const {
        using boost::hash_combine;
        using boost::hash_value;
        // Start with a hash value of 0.
        std::size_t value = 0;
        // Modify 'seed' by XORing and bit-shifting in
        // one member of 'node' after the other:
        hash_combine(value, hash_value(var.x1));
        hash_combine(value, hash_value(var.y1));
        hash_combine(value, hash_value(var.x2));
        hash_combine(value, hash_value(var.y2));
        // Return the result.
        return value;
    }
};


class GridPoint {
public:
    double x, y;
    int layer;
    GridPoint() {}
    GridPoint(double _x, double _y) : x(_x), y(_y){
        layer=-1;// not necessary to assign a layer
    }
    GridPoint(double _x, double _y, int _z) : x(_x), y(_y), layer(_z) {}
    bool operator==(const GridPoint& rhs) const { return (x == rhs.x && y == rhs.y && layer==rhs.layer); }
    friend inline std::ostream& operator<<(std::ostream& os, const GridPoint& gp) {
        os << "(" << gp.x << "," << gp.y << ")";
        return os;
    }
};

class TAP : public GridPoint {
public:
    int id;
    double capacitance;
    TAP() {}
    TAP(int _id, double _x, double _y) : id(_id) {
        x = _x;
        y = _y;
        layer=-1;
    }
    TAP(int _id, double _x, double _y, int _layer, double _capacitance) : id(_id) {
        x = _x;
        y = _y;
        layer=_layer;
        capacitance=_capacitance;
    }

    string str_xy() {
        stringstream s;
        s << "[" << x << "," << y << "]";
        return s.str();
    }
    friend inline std::ostream& operator<<(std::ostream& os, const TAP& tap) {
        os << "TAP " << tap.id << ", x: " << tap.x << ", y: " << tap.y;
        return os;
    }
};
class Segment {
public:
    GridPoint p1, p2;  // p1 has lower y
    // GridPoint center;
    shared_ptr<Segment> ch[2];
    shared_ptr<Segment> par;
    int id = 0;  // unique id for each segment
    Segment() {}
    Segment(GridPoint u, GridPoint v) : p1(u), p2(v) {
        if (p1.y > p2.y) {
            swap(p1, p2);
        }
    }

    Segment(TAP leaf) : p1(leaf), p2(leaf) {
        if (p1.y > p2.y) {
            swap(p1, p2);
        }
    }

    bool isLeaf() { return p1 == p2; }
    bool operator==(const Segment& rhs) const { return (p1 == rhs.p1 && p2 == rhs.p2); }

    friend inline std::ostream& operator<<(std::ostream& os, const Segment& seg) {
        os << "Seg: (" << seg.p1 << "," << seg.p2 << ")";
        return os;
    }

    double slope() {
        if (isLeaf()) {
            return 0;
        }
        if (p2.x == p1.x) {
            return 0;
        }
        return 1.0 * (p1.y - p2.y) / (p1.x - p2.x);
    }

    Segment intersect(Segment& rhs) {
        double cur_slope = slope();
        double rhs_slope = rhs.slope();
        // cout<<cur_slope<<" "<<rhs_slope;
        // check if 4 points same line
        // if (abs(cur_slope - rhs_slope) < eps) {
        //     if (abs((rhs.p1.y - p1.y) * (p2.x - p1.x) - (p2.y - p1.y) * (rhs.p1.x - p1.x)) < eps) {
        if (rhs.isLeaf()) {  // if current segment is intersecting a single grid point
            Segment ret = rhs;
            if (abs((rhs.p1.y - p1.y) * (p2.x - p1.x) - (p2.y - p1.y) * (rhs.p1.x - p1.x)) < eps) {// check if 4 points same line
                // if ((rhs.p1.y - p1.y) * (p2.x - p1.x) == (p2.y - p1.y) * (rhs.p1.x - p1.x)) {
                if (p1.y - eps <= rhs.p1.y && rhs.p1.y <= p2.y + eps) {  // valid intersection
                    Segment ret = rhs;
                    ret.id = -2;  // return single point intersection
                    return ret;
                }
            }
            ret.id = -1;
            return ret;
        }
        if (abs(cur_slope - rhs_slope) < eps) {// equal slope
            // cout<<" "<<(rhs.p1.y - p1.y) * (p2.x - p1.x)<<" "<<(p2.y - p1.y) * (rhs.p1.x - p1.x)<<endl;;
            if (abs((rhs.p1.y - p1.y) * (p2.x - p1.x) - (p2.y - p1.y) * (rhs.p1.x - p1.x)) < eps) {// check if 4 points same line
                assert(rhs.p1.y <= rhs.p2.y && p1.y <= p2.y);
                GridPoint upper, lower;
                if (rhs.p2.y < p2.y) {
                    upper = rhs.p2;
                } else {
                    upper = p2;
                }
                if (rhs.p1.y > p1.y) {
                    lower = rhs.p1;
                } else {
                    lower = p1;
                }
                if (upper.y < lower.y) {
                    Segment ret;
                    ret.id = -1;
                    return ret;

                    // cout << "No overlap between two segs on the line" << endl;
                    // exit(1);
                }
                return Segment(lower, upper);
            } else {
                Segment ret;
                ret.id = -1;
                return ret;
            }
        } else {
            // might be 1 point or 0
            double A1 = p2.y - p1.y;
            // double B1 = p2.x - p1.x;
            double B1 = p1.x - p2.x;
            double C1 = A1 * p1.x + B1 * p1.y;
            double A2 = rhs.p2.y - rhs.p1.y;
            // double B2 = rhs.p2.x - rhs.p1.x;
            double B2 = rhs.p1.x - rhs.p2.x;
            double C2 = A2 * rhs.p1.x + B2 * rhs.p1.y;
            double det = A1 * B2 - A2 * B1;
            double x = (B2 * C1 - B1 * C2) / det;
            double y = (A1 * C2 - A2 * C1) / det;

            Segment ret;
            // if (p1.y-eps <= y && y <= p2.y +eps) {  // valid intersection
            // if(abs(x-p1.x) <= 0.5 && abs(y-p1.y) <= 0.5){
            //     ret.p1 = p1;
            //     ret.p2 = p1;
            //     ret.id = -3;
            // }else if(abs(x-p2.x) <= 0.5 && abs(y-p2.y) <= 0.5){

            // }
            // else if(abs(x-rhs.p1.x) <= 0.5 && abs(y-rhs.p1.y) <= 0.5){

            // }else if(abs(x-rhs.p2.x) <= 0.5 && abs(y-rhs.p2.y) <= 0.5){

            // }
            if (p1.y - eps <= y && y <= p2.y + eps && rhs.p1.y - eps <= y &&
                y <= rhs.p2.y + eps) {  // valid intersection

                ret.p1 = GridPoint(x, y);
                ret.p2 = GridPoint(x, y);
                ret.id = -2;  // return single point intersection
            } else {
                ret.id = -1;
            }
            return ret;
        }
        // if return with id=-1, means no intersection
    }
};
class TRR {
public:
    Segment core;
    double radius;
    TRR() {}
    TRR(Segment seg, double radi) : core(seg), radius(radi) {}
    friend inline std::ostream& operator<<(std::ostream& os, const TRR& trr) {
        os << trr.core << "; radius:" << trr.radius;
        return os;
    }
    Segment intersect(Segment& seg) {
        vector<GridPoint> trr_boundary_grid;
        vector<Segment> trr_Sides;
        //cout<<"seg slope: "<<seg.slope()<<" p1: "<<seg.p1<<" p2: "<<seg.p2<<endl;
        trr_boundary_grid.emplace_back(core.p1.x, core.p1.y - radius);
        trr_boundary_grid.emplace_back(core.p1.x + radius, core.p1.y);
        trr_boundary_grid.emplace_back(core.p1.x, core.p1.y + radius);
        trr_boundary_grid.emplace_back(core.p1.x - radius, core.p1.y);  // clock-wise
        for (int i = 0; i < 3; i++) {
            trr_Sides.emplace_back(trr_boundary_grid[i], trr_boundary_grid[i + 1]);
        }
        trr_Sides.emplace_back(trr_boundary_grid[3], trr_boundary_grid[0]);
        // for (auto& seg1 : trr_Sides) {
        //     cout << seg1 << endl;
        // }
        for (auto& side : trr_Sides) {
            Segment intersection = side.intersect(seg);
            if (intersection.id != -1) {
                return intersection;
            }
        }
        Segment ret;
        ret.id = -1;
        return ret;
    }

    void draw_TRR(ofstream& stream);
    void draw_core(ofstream& stream);

    bool insideTRR(GridPoint point);

};

class TreeNode {
public:
    int id;
    int layer;
    int tree_layer;
    double load_capacitance;// load from taps
    double delay;
    TRR trr;
    shared_ptr<TreeNode> lc;
    shared_ptr<TreeNode> rc;
    shared_ptr<TreeNode> par;
    TreeNode(int _id) {
        id = _id;
        lc = NULL;
        rc = NULL;
        par = NULL;
        layer=-1;
        delay = 0;
    }
    void set_lc(shared_ptr<TreeNode> child) { lc = child; }
    void set_rc(shared_ptr<TreeNode> child) { rc = child; }
    void set_par(shared_ptr<TreeNode> p) { par = p; }
    void set_layer(int layer){layer=layer;}
    bool isSink(){return lc==NULL&&rc==NULL;}
};
class TreeTopology {
public:
    shared_ptr<TreeNode> root;
    shared_ptr<TreeNode> tmp_root;  // used for refinement evaluation
    int leafNumber;
    int size;
    //alglib::integer_2d_array& HC_result;
    unordered_map<int, shared_ptr<TreeNode>> id_treeNode;

    //TreeTopology(alglib::integer_2d_array& HC_res) : HC_result(HC_res) {}
    TreeTopology(){}
    void inittree(int leafNum, int sz, vector<int> preOrder, vector<int> inOrder);
    // construct binary tree structure based on current Hierarchical clustering result
    void constructTree_old(bool modifyCurrentTree = false);
    void constructTree(vector<int> preOrder, vector<int> inOrder);
    void layerassignment(vector<pair<int,int>> IdAndLayer);
    void treeLayerCal();
    int getSumOfDiameter();
    // randomly switch leaf nodes to reduce sum of diameter
    void refineStructure(int iter = 10000);
    shared_ptr<TreeNode> buildTree(vector<int> pre,vector<int> in, int preStart, int preEnd, int inStart, int inEnd);
};

class GrSteiner : public GridPoint {
public:
    shared_ptr<GrSteiner> lc;
    shared_ptr<GrSteiner> rc;
    shared_ptr<GrSteiner> par;

    GrSteiner(GridPoint p){
        x = p.x;
        y= p.y;
        lc = NULL;
        rc = NULL;
        par = NULL;
    }

    void set_lc(shared_ptr<GrSteiner> child) { lc = child; }
    void set_rc(shared_ptr<GrSteiner> child) { rc = child; }
    void set_par(shared_ptr<GrSteiner> p) { par = p; }
};

class Router {
public:
    int MAX_RUNTIME = 3600;  // test 1
    int NUM_TAPS;        // test 1
    int delay_model=1;//0 for linear and 1 for elmore 
    vector<TAP> taps;
    vector<Segment> vertexMS;  // set of segments
    vector<TRR> vertexTRR;
    vector<double> vertexDistE;
    shared_ptr<TreeTopology> topo;

    string _RUNDIR = "../run_tmp/";
    string _DRAWDIR= "../draw/";

    // Structures to store the final routing result
    GridPoint clockSource;
    vector<GridPoint> pl;
    vector<shared_ptr<GrSteiner>> sol;
        // vector<vector<GridPoint>> sol;

    void init();
    void readInput();
    // To generate topology(try L1 based metric and L2 based ones)
    void NS();          // nearest neighbor topology
    void HC();          // Agglomerative  Hierarchical Clustering
    void CL();          // clustering based algorithm
    void Refinement();  // clustering with refinement 4<= m <=6 in the original paper
    void RGM();         // Recursive Geometric Matching

    // Generate embedding
    void DME();  // Deferred-Merge Embedding
    void route();
    void buildSolution();
    void reportTotalWL();
    void writeSolution();
    void buildTopology();
    void setdelay_model(int);
    void draw_bottom_up();
    void draw_top_down();
    void draw_TRR_pair(TRR trr1,TRR trr2);
    Segment TRRintersect(TRR& trr1,TRR& trr2);
};

double calc_x_RC(shared_ptr<TreeNode> nodeLeft, shared_ptr<TreeNode> nodeRight, shared_ptr<TreeNode> nodeMerge, double L);

double calc_L2_RC(shared_ptr<TreeNode> nodeLeft, shared_ptr<TreeNode> nodeRight, shared_ptr<TreeNode> nodeMerge, int tag);

void update_merge_Capacitance(shared_ptr<TreeNode> nodeMerge, shared_ptr<TreeNode> nodeLeft, shared_ptr<TreeNode> nodeRight, double ea, double eb);

void update_merge_Delay(shared_ptr<TreeNode> nodeMerge, shared_ptr<TreeNode> nodeLeft, shared_ptr<TreeNode> nodeRight, double ea, double eb);

double calc_standard_Capacitance(double capacitance_in_fF);

float calc_delay_RLC(shared_ptr<TreeNode> nodeMerge, shared_ptr<TreeNode> nodeChild, float WL);

void RLC_calculation(shared_ptr<TreeNode> nodeMerge, shared_ptr<TreeNode> nodeLeft, shared_ptr<TreeNode> nodeRight, double &ea, double &eb);

void modify_coord_by_L1(double &ea, double &eb, shared_ptr<TreeNode> nodeLeft, shared_ptr<TreeNode> nodeRight, float x);

void modify_coord_by_L2(double &ea, double &eb, shared_ptr<TreeNode> nodeLeft, shared_ptr<TreeNode> nodeRight, float x);

double L1Dist(GridPoint p1, GridPoint p2);

double min_manhattan_dist(shared_ptr<TreeNode> nodeLeft, shared_ptr<TreeNode> nodeRight);


