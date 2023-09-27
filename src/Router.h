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
const double eps = 1e-4;

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
    TAP(int _id, double _x, double _y, int _layer) : id(_id) {
        x = _x;
        y = _y;
        layer=_layer;
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
};

class TreeNode {
public:
    int id;
    int layer;
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
    int NUM_TAPS = 4;        // test 1
    int delay_model=1;//0 for linear and 1 for elmore 
    vector<TAP> taps;
    vector<Segment> vertexMS;  // set of segments
    vector<TRR> vertexTRR;
    vector<double> vertexDistE;
    shared_ptr<TreeTopology> topo;

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
};

double calc_x_RC(shared_ptr<TreeNode> nodeLeft, shared_ptr<TreeNode> nodeRight, shared_ptr<TreeNode> nodeMerge, double L){
    double beta = r_v * (abs(nodeRight->layer - nodeMerge->layer)*nodeRight->load_capacitance -
                        abs(nodeMerge->layer - nodeLeft->layer)*nodeLeft->load_capacitance +
                        0.5*c_v*(pow((nodeRight->layer - nodeMerge->layer),2) -
                        pow((nodeMerge->layer - nodeLeft->layer),2)));
    double up = (nodeRight->delay - nodeLeft->delay) + r_w * L * (nodeRight->load_capacitance + 0.5 * c_w * L) + beta + r_v * c_w * abs(nodeRight->layer - nodeMerge->layer) * L;
    double down = r_w * (nodeLeft->load_capacitance + nodeRight->load_capacitance  + c_w * L) + r_v * c_w * abs(nodeRight->layer - nodeLeft->layer);
    double x = up/down;
    return x;
}

double calc_L2_RC(shared_ptr<TreeNode> nodeLeft, shared_ptr<TreeNode> nodeRight, shared_ptr<TreeNode> nodeMerge, int tag){
    double alpha, beta, up;
    //tag = 0: |eb| = L'
    //tag = 1: |ea| = L'
    if(tag == 0){
        alpha = r_w * nodeRight->load_capacitance + r_v * c_w * abs(nodeRight->layer - nodeMerge->layer);
        beta = r_v * (abs(nodeRight->layer - nodeMerge->layer)*nodeRight->load_capacitance -
                      abs(nodeMerge->layer - nodeLeft->layer)*nodeLeft->load_capacitance +
               0.5*c_v*((nodeRight->layer - nodeMerge->layer)^2 -
                                                              (nodeMerge->layer - nodeLeft->layer)^2));
        up = sqrt(2 * r_w * c_w * (nodeLeft->delay - nodeRight->delay) + alpha * alpha - 2 * r_w * c_w * beta) - alpha;
    }

    else{
        alpha = r_w * nodeLeft->load_capacitance + r_v * c_w * abs(nodeLeft->layer - nodeMerge->layer);
        beta = r_v * (abs(nodeMerge->layer - nodeLeft->layer)*nodeLeft->load_capacitance -
                     abs(nodeRight->layer - nodeMerge->layer)*nodeRight->load_capacitance +
                     0.5*c_v*(pow((nodeRight->layer - nodeMerge->layer),2) -
                     pow((nodeMerge->layer - nodeLeft->layer),2)));
        up = sqrt(2 * r_w * c_w * (nodeRight->delay - nodeLeft->delay) + alpha * alpha - 2 * r_w * c_w * beta) - alpha;
    }     
    return up/(r_w * c_w);
}

void update_merge_Capacitance(shared_ptr<TreeNode> nodeMerge, shared_ptr<TreeNode> nodeLeft, shared_ptr<TreeNode> nodeRight, float ea, float eb){
    float delta_C = nodeLeft->load_capacitance + nodeRight->load_capacitance + c_w*(ea + eb) + c_v*(abs(nodeRight->layer - nodeLeft->layer));
    //考虑了buffer insertion的电容update
    if(delta_C > c_constraint){
        nodeMerge->load_capacitance = 300;
        //nodeMerge->needBuffer = 1;
        return;
    }
    nodeMerge->load_capacitance = delta_C;
}

void update_merge_Delay(shared_ptr<TreeNode> nodeMerge, shared_ptr<TreeNode> nodeLeft, shared_ptr<TreeNode> nodeRight, double ea, double eb){
    float delta_delay = nodeLeft->delay + 0.695 * (
            0.5 * r_w * c_w * ea * ea +
            r_w * (nodeLeft->load_capacitance + c_v * abs(nodeMerge->layer - nodeLeft->layer)) * ea +
            r_v * (nodeLeft->load_capacitance * abs(nodeMerge->layer - nodeLeft->layer) + 0.5 * c_v * (nodeMerge->layer - nodeLeft->layer) * (nodeMerge->layer - nodeLeft->layer)) -
            (r_w * c_v - r_v * c_w) * abs(nodeMerge->layer - nodeLeft->layer) * ea);
    float delta_delay_right = nodeRight->delay + 0.695 * (
            0.5 * r_w * c_w * eb * eb +
            r_w * (nodeRight->load_capacitance + c_v * abs(nodeMerge->layer - nodeRight->layer)) * eb +
            r_v * (nodeRight->load_capacitance * abs(nodeMerge->layer - nodeRight->layer) + 0.5 * c_v * (nodeMerge->layer - nodeRight->layer) * (nodeMerge->layer - nodeRight->layer)) -
            (r_w * c_v - r_v * c_w) * abs(nodeMerge->layer - nodeRight->layer) * eb);

    /*float delta_delay_right = nodeRight->delay +
            r_w * (ea + eb) * (nodeRight->C + 0.5 * c_w * (ea + eb)) +
            0.5 * r_w * c_w * ea * ea -
            r_w * (nodeRight->C + c_v * abs(nodeMerge->layer - nodeRight->layer) + c_w * (ea + eb)) * ea +
            r_v * (nodeRight->C * abs(nodeMerge->layer - nodeRight->layer) + 0.5 * c_v * abs(nodeMerge->layer - nodeRight->layer) * abs(nodeMerge->layer - nodeRight->layer)) +
            r_w * c_w * abs(nodeMerge->layer - nodeRight->layer) * (ea + eb) -
            (r_w * c_v - r_v * c_w) * abs(nodeMerge->layer - nodeRight->layer) * ea;*/
    printf("/****Test Merge Delay****/\n"
           "Delay before merge: Left: %f, Right: %f, after: Left: %f,  Right: %f\n", nodeLeft->delay,nodeRight->delay, delta_delay, delta_delay_right);
    //按最大delay进行录入
    if(delta_delay >= delta_delay_right)
        nodeMerge->delay = delta_delay;
    else
        nodeMerge->delay = delta_delay_right;
}

double calc_standard_Capacitance(double capacitance_in_fF){
    return (capacitance_in_fF/1000000000000000);
}

float calc_delay_RLC(shared_ptr<TreeNode> nodeMerge, shared_ptr<TreeNode> nodeChild, float WL){// ta calculation
    float t_pdi, theta, omega, numerator, denominator, elmore;
    float wireLength = WL;
    if(wireLength == 0)
        return 0;
    //如果两节点中间没有TSV
    assert(nodeChild->load_capacitance<= c_constraint);
    numerator = wireLength * r_w * (0.5 * wireLength * c_w_standard + calc_standard_Capacitance(nodeChild->load_capacitance));
    //这里分母还没算完，后续要开根
    denominator = wireLength * l_w * (0.5 * wireLength * c_w_standard + calc_standard_Capacitance(nodeChild->load_capacitance));
    if(nodeMerge->layer != nodeChild->layer){
        numerator += r_v * (0.5 * c_v_standard + calc_standard_Capacitance(nodeChild->load_capacitance + wireLength * c_w));
        //denominator += l_v * (0.5 * c_v_standard + calc_standard_Capacitance(nodeChild->C + wireLength * c_w));
    }
    //给分母开根
    denominator = sqrt(denominator);
    //154953990144.000000
    theta = 0.5 * (numerator/denominator);
    omega = 1/denominator;
    elmore = 0.695 * numerator;
    //将单位换算回ps
    t_pdi = roundf(1000000000000000 * (  (1.047 * exp((-1)*theta/0.85))/omega + elmore  ));
    //printf("before : %f\n", 1000000000000000 * (  (1.047 * exp((-1)*theta/0.85))/omega));
    return t_pdi;
}

void RLC_calculation(shared_ptr<TreeNode> nodeMerge, shared_ptr<TreeNode> nodeLeft, shared_ptr<TreeNode> nodeRight, double ea_pointer, double eb_pointer){
    float t_a, t_b;
    t_a = calc_delay_RLC(nodeMerge, nodeLeft, ea_pointer);
    t_b = calc_delay_RLC(nodeMerge, nodeRight, eb_pointer);
    //需要考虑extra wirelength存在的情况。但是尽量不要引入额外的线长。
    while(fabs(t_a + nodeLeft->delay - t_b - nodeRight->delay) >= 1){//? >=1?
        if(t_a + nodeLeft->delay > t_b + nodeRight->delay)
            modify_coord_by_L1(ea_pointer, eb_pointer, nodeLeft, nodeRight, skewModifyStep);
        else
            modify_coord_by_L2(ea_pointer, eb_pointer, nodeLeft, nodeRight, skewModifyStep);

        t_a = calc_delay_RLC(nodeMerge, nodeLeft, ea_pointer);
        t_b = calc_delay_RLC(nodeMerge, nodeRight, eb_pointer);
    }
    //printf("原本的 delay: left: %f, right: %f\n", nodeLeft->delay, nodeRight->delay);
    printf("left: %f, right: %f\n t_a: %f, t_b: %f\n total left: %f, total right: %f\n", nodeLeft->delay, nodeRight->delay, t_a, t_b, t_a+nodeLeft->delay, t_b + nodeRight->delay);
    t_a = calc_delay_RC(nodeMerge, nodeLeft, *ea_pointer);
    t_b = calc_delay_RC(nodeMerge, nodeRight, *eb_pointer);
    printf("RC delay: %f, %f\n", t_a, t_b);

    update_merge_Capacitance(nodeMerge, nodeLeft, nodeRight, *ea_pointer, *eb_pointer);
    //update_merge_Delay(nodeMerge, nodeLeft, nodeRight, ea, eb);
    if(t_a + nodeLeft->delay > t_b + nodeRight->delay)
        nodeMerge->delay = t_a + nodeLeft->delay;
    else
        nodeMerge->delay = t_b + nodeRight->delay;
}

