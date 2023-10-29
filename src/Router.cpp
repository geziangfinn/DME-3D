#include "Router.h"
#include "Setting.h"
#include "Parser.h"
#include "Drawer.h"
using std::cout;
using std::endl;
using std::setprecision;

#define COMPLETE_LINKAGE 0
#define SINGLE_LINKAGE 1
#define L1 1
#define L2 2
#define eps 1e-6
const string padding(30, '=');

void TreeTopology::inittree(int leafNum, int sz, vector<int> preOrder, vector<int> inOrder) {
    leafNumber = leafNum;
    size = sz;
    cout << "Initialize topo: " << leafNum << " leaves and " << sz << " nodes in total" << endl;
    //constructTree(true);
    constructTree(preOrder,inOrder);
}
void TreeTopology::constructTree_old(bool modifyCurrentTree) {
    vector<vector<int>> DAG(size);
    // cout << "rows " << HC_result.rows() << endl;
    //int n_merges = HC_result.rows();
    // int cur_internal_node_idx = leafNumber;
    // for (int i = 0; i < n_merges; i++) {
    //     int p1 = HC_result[i][0];
    //     int p2 = HC_result[i][1];
    //     DAG[cur_internal_node_idx].push_back(p1);
    //     DAG[cur_internal_node_idx].push_back(p2);
    //     cur_internal_node_idx++;
    //     cout << "Merge " << p1 << " and " << p2 << endl;
    // }
    // Construct tree
    tmp_root = make_shared<TreeNode>(size - 1);  // last merge point as root
    std::function<void(shared_ptr<TreeNode>)> buildTree = [&](shared_ptr<TreeNode> curNode) {
        int curId = curNode->id;
        id_treeNode[curId] = curNode;
        if (DAG[curId].size() != 0) {
            auto lc = make_shared<TreeNode>(DAG[curId][0]);
            auto rc = make_shared<TreeNode>(DAG[curId][1]);
            curNode->set_lc(lc);
            curNode->set_rc(rc);
            lc->set_par(curNode);
            rc->set_par(curNode);
            buildTree(curNode->lc);
            buildTree(curNode->rc);
            return;
        }
    };
    buildTree(tmp_root);

    std::function<void(shared_ptr<TreeNode>)> postOrderTraversal = [&](shared_ptr<TreeNode> curNode) {
        int curId = curNode->id;
        if (curNode->lc != NULL && curNode->rc != NULL) {
            postOrderTraversal(curNode->lc);
            postOrderTraversal(curNode->rc);
            cout << "Vis: " << curId << endl;
            return;
        } else {
            cout << "Vis: " << curId << endl;
        }
    };
    // postOrderTraversal(tmp_root);

    if (modifyCurrentTree == true) {
        root = tmp_root;
    }
}

void TreeTopology::constructTree(vector<int> preOrder, vector<int> inOrder)
{
    root=buildTree(preOrder,inOrder,0,preOrder.size()-1,0,inOrder.size()-1);

    std::function<void(shared_ptr<TreeNode>)> postOrderTraversal = [&](shared_ptr<TreeNode> curNode) {
         if(curNode!=nullptr){
            postOrderTraversal(curNode->lc);
            postOrderTraversal(curNode->rc);
            int curId = curNode->id;
            cout<<"Vis: "<<curId<<endl;
            //cout<<"preing "<<curId<<" at layer "<<curNode->layer<<endl;
        }
    };
    postOrderTraversal(root);
}

void TreeTopology::layerassignment(vector<pair<int, int>> IdAndLayer){

    std::function<void(shared_ptr<TreeNode>, int&)> preOrderTraversal = [&](shared_ptr<TreeNode> curNode, int& index) {
        if(curNode!=nullptr){
            int curId = curNode->id;
            assert(IdAndLayer[index].first==curId);
            curNode->layer=IdAndLayer[index].second;
            index++;
            cout<<"preing "<<curId<<" at layer "<<curNode->layer<< " at tree-layer "<<curNode->tree_layer<<endl;
            preOrderTraversal(curNode->lc,index);
            preOrderTraversal(curNode->rc,index);
        }
    };
    int index=0;
    preOrderTraversal(root,index);
}

void TreeTopology::treeLayerCal()
{
    std::function<void(shared_ptr<TreeNode>)> levelTraversal=[&](shared_ptr<TreeNode> curNode){
        queue<shared_ptr<TreeNode>> nodeQueue;
        int count=0;
        assert(curNode!=NULL);
        nodeQueue.emplace(curNode);
        while(!nodeQueue.empty())
        {
            int size=nodeQueue.size();
            while(size>0){   //如果当前队列的结点大小减去1仍然大于0的话
                shared_ptr<TreeNode> cur = nodeQueue.front();
                nodeQueue.pop();     //获取队列的第一个结点，出队列
                cur->tree_layer=count;           //访问当前结点，可以任意改成其他打印操作等等
                if(cur->lc!=NULL){
                    nodeQueue.push(cur->lc);  //左孩子不为空，左孩子入队列
                }
                if(cur->rc!=NULL){
                    nodeQueue.push(cur->rc);  //右孩子不为空，右孩子入队列
                }
                size--;
            }
            count++;                   //如果size不满足条件，退出内层循环，层数加一;
        }
    };
    levelTraversal(root);
}

shared_ptr<TreeNode> TreeTopology::buildTree(vector<int> pre, vector<int> in, int preStart, int preEnd, int inStart, int inEnd)
{
        if(preStart>preEnd)
        {
            return NULL;
        }
        auto curRoot=make_shared<TreeNode>(pre[preStart]);
        id_treeNode[pre[preStart]] = curRoot;

        int rootInInorder=inStart;

        while (rootInInorder < inEnd && pre[preStart] != in[rootInInorder]) rootInInorder++;

        curRoot->lc=buildTree(pre,in,preStart+1,rootInInorder-inStart+preStart,inStart,rootInInorder-1);
        curRoot->rc=buildTree(pre,in,rootInInorder-inStart+preStart+1,preEnd,rootInInorder+1,inEnd);
        if(curRoot->lc)
        {
            curRoot->lc->set_par(curRoot);
            cout<<"I am "<<curRoot->lc->id<<" my par is "<<curRoot->id<<endl;
        }
        if(curRoot->rc)
        {
            curRoot->rc->set_par(curRoot);
            cout<<"I am "<<curRoot->rc->id<<" my par is "<<curRoot->id<<endl;
        }
        return curRoot;
}

void Router::init() {
    // read input and build data structure
    _RUNDIR = "../run_tmp/" + setting.get_case_name() + "/";
    _DRAWDIR = "../draw/" + setting.get_case_name() + "/";
    string cmd_clean = "rm -rf " + _RUNDIR + "; mkdir -p " + _RUNDIR+"; rm -rf " + _DRAWDIR + "; mkdir -p " + _DRAWDIR;
    system(cmd_clean.c_str());
    readInput();
}
void Router::readInput() {
    vector<sink> sinks=parse(setting.input_file_name);

    int tap_id=0;
    for(sink sink:sinks)
    {
        taps.emplace_back(tap_id, sink.x, sink.y,sink.layer,sink.capacitance);
        tap_id++;
    }
    cout.setf(ios::fixed,ios::floatfield);
    for (auto& tap : taps) {
        cout << tap << endl;
    }
    cout << padding << "Finish Reading input" << padding << endl;
}

// Recursive Geometric Matching
void Router::RGM() {
    // TODO
}

// nearest neighbor topology
void Router::NS() {
    // TODO
}

// Agglomerative  Hierarchical Clustering
void Router::HC() {
    // // 1. Create XY array
    // // vector<vector<double>> points;
    // string points_str;
    // points_str += "[";
    // for (int i = 0; i < taps.size(); i++) {
    //     points_str += taps[i].str_xy();
    //     if (i != taps.size() - 1) {
    //         points_str += ",";
    //     }
    // }
    // points_str += "]";
    // // cout <<"Points: " << points_str << endl;

    // // 2. Create clusterizer
    // using namespace alglib;
    // clusterizerstate s;
    // ahcreport rep;
    // real_2d_array xy(points_str.c_str());
    // // int linkage_type = SINGLE_LINKAGE;
    // int linkage_type = COMPLETE_LINKAGE;
    // int metric = L1;

    // clusterizercreate(s);
    // clusterizersetpoints(s, xy, metric);  // manhattan dist
    // clusterizersetahcalgo(s, linkage_type);
    // clusterizerrunahc(s, rep);
    // printf("%s\n", rep.z.tostring().c_str());

    // // 3. Construct binary tree topology
    // //topo = make_shared<TreeTopology>(rep.z);
    // //topo->init(xy.rows(), 2 * xy.rows() - 1);

    // cout << padding << "Finish Topology generation" << padding << endl;
}



Segment Router::TRRintersect(TRR& trr1, TRR& trr2) {
    // get four edges
    // cout << "Merging: " << trr1 << " and " << trr2 << endl;
    vector<GridPoint> trr1_boundary_grid;
    vector<GridPoint> trr2_boundary_grid;
    vector<Segment> trr1_Sides;
    vector<Segment> trr2_Sides;
    assert(trr1.radius!=0||trr2.radius!=0);
    //if there is one trr's radius = 0
    if(trr1.radius==0||trr2.radius==0)
    {
        if(trr1.radius==0)
        {

            if (trr2.core.slope() > 0) {
                trr2_boundary_grid.emplace_back(trr2.core.p1.x, trr2.core.p1.y - trr2.radius);
                trr2_boundary_grid.emplace_back(trr2.core.p2.x + trr2.radius, trr2.core.p2.y);
                trr2_boundary_grid.emplace_back(trr2.core.p2.x, trr2.core.p2.y + trr2.radius);
                trr2_boundary_grid.emplace_back(trr2.core.p1.x - trr2.radius, trr2.core.p1.y);  // clock-wise
            } else if (trr2.core.slope() < 0) {
                trr2_boundary_grid.emplace_back(trr2.core.p1.x + trr2.radius, trr2.core.p1.y);
                trr2_boundary_grid.emplace_back(trr2.core.p2.x, trr2.core.p2.y + trr2.radius);
                trr2_boundary_grid.emplace_back(trr2.core.p2.x - trr2.radius, trr2.core.p2.y);
                trr2_boundary_grid.emplace_back(trr2.core.p1.x, trr2.core.p1.y - trr2.radius);  // clock-wise
            } else {                                                                            // leaf node
                trr2_boundary_grid.emplace_back(trr2.core.p1.x, trr2.core.p1.y - trr2.radius);
                trr2_boundary_grid.emplace_back(trr2.core.p1.x + trr2.radius, trr2.core.p1.y);
                trr2_boundary_grid.emplace_back(trr2.core.p1.x, trr2.core.p1.y + trr2.radius);
                trr2_boundary_grid.emplace_back(trr2.core.p1.x - trr2.radius, trr2.core.p1.y);  // clock-wise
            }

            for (int i = 0; i < 3; i++) {
                trr2_Sides.emplace_back(trr2_boundary_grid[i], trr2_boundary_grid[i + 1]);
            }
            trr2_Sides.emplace_back(trr2_boundary_grid[3], trr2_boundary_grid[0]);

            for (auto& seg2 : trr2_Sides) {
                //cout<<"seg1: "<<seg1<<"seg2: "<<seg2<<endl;
                Segment seg = trr1.core.intersect(seg2);//! seg should be a single point in most cases
                // ? could there be 2 intersection points for core and trr sides?
                if (seg.id == 0) {
                    return seg;
                }
                else if(seg.id==-2)
                {
                    if(trr2.insideTRR(trr1.core.p1))
                    {
                        return Segment(seg.p1,trr1.core.p1);
                    }
                    else if(trr2.insideTRR(trr1.core.p2))
                    {
                        return Segment(seg.p1,trr1.core.p2);
                    }
                }

            }

        }
        else{
            if (trr1.core.slope() > 0) {
                trr1_boundary_grid.emplace_back(trr1.core.p1.x, trr1.core.p1.y - trr1.radius);
                trr1_boundary_grid.emplace_back(trr1.core.p2.x + trr1.radius, trr1.core.p2.y);
                trr1_boundary_grid.emplace_back(trr1.core.p2.x, trr1.core.p2.y + trr1.radius);
                trr1_boundary_grid.emplace_back(trr1.core.p1.x - trr1.radius, trr1.core.p1.y);  // clock-wise
            } else if (trr1.core.slope() < 0) {
                trr1_boundary_grid.emplace_back(trr1.core.p1.x + trr1.radius, trr1.core.p1.y);
                trr1_boundary_grid.emplace_back(trr1.core.p2.x, trr1.core.p2.y + trr1.radius);
                trr1_boundary_grid.emplace_back(trr1.core.p2.x - trr1.radius, trr1.core.p2.y);
                trr1_boundary_grid.emplace_back(trr1.core.p1.x, trr1.core.p1.y - trr1.radius);  // clock-wise
            } else {                                                                            // leaf node
                trr1_boundary_grid.emplace_back(trr1.core.p1.x, trr1.core.p1.y - trr1.radius);
                trr1_boundary_grid.emplace_back(trr1.core.p1.x + trr1.radius, trr1.core.p1.y);
                trr1_boundary_grid.emplace_back(trr1.core.p1.x, trr1.core.p1.y + trr1.radius);
                trr1_boundary_grid.emplace_back(trr1.core.p1.x - trr1.radius, trr1.core.p1.y);  // clock-wise
            }
            
            for (int i = 0; i < 3; i++) {
                trr1_Sides.emplace_back(trr1_boundary_grid[i], trr1_boundary_grid[i + 1]);
            }
            trr1_Sides.emplace_back(trr1_boundary_grid[3], trr1_boundary_grid[0]);

            for (auto& seg1 : trr1_Sides) {
                //cout<<"seg1: "<<seg1<<"seg2: "<<seg2<<endl;
                Segment seg = trr2.core.intersect(seg1);//! seg should be a single point in most cases
                // ? could there be 2 intersection points for core and trr sides?
                if (seg.id == 0) {
                    return seg;
                }
                else if(seg.id==-2)
                {
                    if(trr1.insideTRR(trr2.core.p1))
                    {
                        return Segment(seg.p1,trr2.core.p1);
                    }
                    else if(trr1.insideTRR(trr2.core.p2))
                    {
                        return Segment(seg.p1,trr2.core.p2);
                    }
                }

            }

        }
        
        Segment ret;
        ret.id = -1;
        return ret;
    }
    // if both trr's radius > 0
    if (trr1.core.slope() > 0) {
        trr1_boundary_grid.emplace_back(trr1.core.p1.x, trr1.core.p1.y - trr1.radius);
        trr1_boundary_grid.emplace_back(trr1.core.p2.x + trr1.radius, trr1.core.p2.y);
        trr1_boundary_grid.emplace_back(trr1.core.p2.x, trr1.core.p2.y + trr1.radius);
        trr1_boundary_grid.emplace_back(trr1.core.p1.x - trr1.radius, trr1.core.p1.y);  // clock-wise
    } else if (trr1.core.slope() < 0) {
        trr1_boundary_grid.emplace_back(trr1.core.p1.x + trr1.radius, trr1.core.p1.y);
        trr1_boundary_grid.emplace_back(trr1.core.p2.x, trr1.core.p2.y + trr1.radius);
        trr1_boundary_grid.emplace_back(trr1.core.p2.x - trr1.radius, trr1.core.p2.y);
        trr1_boundary_grid.emplace_back(trr1.core.p1.x, trr1.core.p1.y - trr1.radius);  // clock-wise
    } else {                                                                            // leaf node
        trr1_boundary_grid.emplace_back(trr1.core.p1.x, trr1.core.p1.y - trr1.radius);
        trr1_boundary_grid.emplace_back(trr1.core.p1.x + trr1.radius, trr1.core.p1.y);
        trr1_boundary_grid.emplace_back(trr1.core.p1.x, trr1.core.p1.y + trr1.radius);
        trr1_boundary_grid.emplace_back(trr1.core.p1.x - trr1.radius, trr1.core.p1.y);  // clock-wise
    }

    if (trr2.core.slope() > 0) {
        trr2_boundary_grid.emplace_back(trr2.core.p1.x, trr2.core.p1.y - trr2.radius);
        trr2_boundary_grid.emplace_back(trr2.core.p2.x + trr2.radius, trr2.core.p2.y);
        trr2_boundary_grid.emplace_back(trr2.core.p2.x, trr2.core.p2.y + trr2.radius);
        trr2_boundary_grid.emplace_back(trr2.core.p1.x - trr2.radius, trr2.core.p1.y);  // clock-wise
    } else if (trr2.core.slope() < 0) {
        trr2_boundary_grid.emplace_back(trr2.core.p1.x + trr2.radius, trr2.core.p1.y);
        trr2_boundary_grid.emplace_back(trr2.core.p2.x, trr2.core.p2.y + trr2.radius);
        trr2_boundary_grid.emplace_back(trr2.core.p2.x - trr2.radius, trr2.core.p2.y);
        trr2_boundary_grid.emplace_back(trr2.core.p1.x, trr2.core.p1.y - trr2.radius);  // clock-wise
    } else {                                                                            // leaf node
        trr2_boundary_grid.emplace_back(trr2.core.p1.x, trr2.core.p1.y - trr2.radius);
        trr2_boundary_grid.emplace_back(trr2.core.p1.x + trr2.radius, trr2.core.p1.y);
        trr2_boundary_grid.emplace_back(trr2.core.p1.x, trr2.core.p1.y + trr2.radius);
        trr2_boundary_grid.emplace_back(trr2.core.p1.x - trr2.radius, trr2.core.p1.y);  // clock-wise
    }

    for (int i = 0; i < 3; i++) {
        trr1_Sides.emplace_back(trr1_boundary_grid[i], trr1_boundary_grid[i + 1]);
        trr2_Sides.emplace_back(trr2_boundary_grid[i], trr2_boundary_grid[i + 1]);
    }
    trr1_Sides.emplace_back(trr1_boundary_grid[3], trr1_boundary_grid[0]);
    trr2_Sides.emplace_back(trr2_boundary_grid[3], trr2_boundary_grid[0]);

    // cout << "Print trr1's sides" << endl;
    // for (auto& seg1 : trr1_Sides) {
    //     cout << seg1 << endl;
    // }

    // cout << "Print trr2's sides" << endl;

    // for (auto& seg2 : trr2_Sides) {
    //     cout << seg2 << endl;
    // }

    // for 4*4 check intersect
    for (auto& seg1 : trr1_Sides) {
        for (auto& seg2 : trr2_Sides) {
            //cout<<"seg1: "<<seg1<<"seg2: "<<seg2<<endl;
            Segment seg = seg1.intersect(seg2);
            if (seg.id == 0||seg.id==-2) {
                return seg;
            }
        }
    }
    cout << "Cannot find intersection between two TRRs" << endl;
    Segment ret;
    draw_TRR_pair(trr1,trr2);
    for (auto& seg1 : trr1_Sides) {
        for (auto& seg2 : trr2_Sides) {
            cout<<"seg1: "<<seg1<<"seg2: "<<seg2<<endl;
            Segment seg = seg1.intersect(seg2);
        }
    }
    ret.id = -1;
    return ret;
}
// Deferred-Merge Embedding
void Router::DME() {
    // Segment seg1(GridPoint(1.0,3.0),GridPoint(2.0,4.0));

    // Segment seg2(GridPoint(-0.5,6.5),GridPoint(6,0));

    // cout << seg1.intersect(seg2) << endl;
    // exit(1);
    cout << BLUE << "[Router]" << RESET<<" begin DME\n";

    vertexMS.resize(topo->size);
    vertexTRR.resize(topo->size);
    vertexDistE.resize(topo->size);

    bool RLC=true;

    // 1. Build Tree of Segments (bottom up)
    std::function<void(shared_ptr<TreeNode>)> postOrderTraversal = [&](shared_ptr<TreeNode> curNode) {
        int curId = curNode->id;
        if (curNode->lc != NULL && curNode->rc != NULL) {//!的确不会有只有一个子节点的中间节点
            postOrderTraversal(curNode->lc);
            postOrderTraversal(curNode->rc);

            // create merging segment for curNode
            //auto& ms_a = vertexMS[curNode->lc->id];
            //auto& ms_b = vertexMS[curNode->rc->id];

            auto ms_a = curNode->lc->trr.core;
            auto ms_b = curNode->rc->trr.core;
            // get |e_a|, |e_b|
            double d = min(L1Dist(ms_a.p1, ms_b.p1), L1Dist(ms_a.p1, ms_b.p2));
            d = min(d, L1Dist(ms_a.p2, ms_b.p1));
            d = min(d, L1Dist(ms_a.p2, ms_b.p2));  // but why need to calc 2*2 possiblity?
            
            // double e_a_dist = (ms_b.delay - ms_a.delay + d) / 2;
            // double e_b_dist = (ms_a.delay - ms_b.delay + d) / 2;
            double e_a_dist;
            double e_b_dist;
            //! ea for lc and eb for rc
            if(delay_model==LINEAR_DELAY)// linear delay model
            {
                e_a_dist = (curNode->rc->delay - curNode->lc->delay + d) / 2;
                e_b_dist = (curNode->lc->delay - curNode->rc->delay + d) / 2;
                if (e_a_dist < 0 || e_b_dist < 0) {
                    cout << "Skew too large" << endl;//!
                    exit(1);
                }
            } else if(delay_model==ELMORE_DELAY)// elmore delay(3d)
            {
                double x=calc_x_RC(curNode->lc,curNode->rc,curNode,d);
                if(0 <= x && x <= d){
                    e_a_dist = x;
                    e_b_dist = d - x;
                }
                else if(x < 0){//! 这里是否有假设 a b的相对位置关系？？
                    e_b_dist = calc_L2_RC(curNode->lc,curNode->rc,curNode, 0);
                    //!assert(e_b_dist > d);
                    e_a_dist = 0;
                }
                else if(x > d){
                    e_a_dist = calc_L2_RC(curNode->lc,curNode->rc,curNode, 1);
                    //!assert(e_a_dist > d);
                    e_b_dist = 0;
                }
            }

            RLC_calculation(curNode,curNode->lc,curNode->rc,e_a_dist,e_b_dist);

            // todo : this delay should be changed
            // ms_v.delay = e_a_dist + ms_a.delay; 
            // e_a+ms_a is supposed to be equal to e_b+ms_b
            // todo update treenode delay and capacitance, segment should be a member of tree node, and TRR should be a member of segment, what about rebuild the code structure?
            
            update_merge_Capacitance(curNode,curNode->lc, curNode->rc,e_a_dist,e_b_dist);//? there is a same function in RLC_calculation

            update_merge_Delay(curNode,curNode->lc, curNode->rc,e_a_dist,e_b_dist);

            vertexDistE[curNode->lc->id] = e_a_dist;
            vertexDistE[curNode->rc->id] = e_b_dist;

            // get trr_a, trr_b
            //TRR trr_a(ms_a, e_a_dist);
            //TRR trr_b(ms_b, e_b_dist);
            //vertexTRR[curNode->lc->id] = trr_a;
            //vertexTRR[curNode->rc->id] = trr_b;

            curNode->lc->trr.radius=e_a_dist;
            curNode->rc->trr.radius=e_b_dist;

            // intersect trr_a, trr_b to get ms_v
            Segment ms_v = TRRintersect(curNode->lc->trr, curNode->rc->trr);
            // cout << "Merging result: " << ms_v << endl;
            if (ms_v.id == -1) {
                cout << "Merge failure" << endl;
                exit(1);
            }

            //! new
            //curNode->load_capacitance=curNode->lc->load_capacitance+curNode->rc->load_capacitance+c_v*d;//todo: a+b+wire capacitance, and consider snaking?
            //? not d but max(ea,eb)?
            //! new
            vertexMS[curId] = ms_v;

            curNode->trr.core=ms_v;
            // cout << "Delay diff " << e_a_dist + ms_a.delay - (e_b_dist + ms_b.delay) << endl;
        } else {
            // Create ms for leaf node
            
            vertexMS[curId] = Segment(taps[curId], taps[curId]);
            //vertexMS[curId] = Segment(taps[curId]);
            
            curNode->trr.core = Segment(taps[curId]);
            //cout<<"leaf node id: "<<curId<<endl;
        }
    };
    postOrderTraversal(topo->root);
    cout  << "Finish bottom-up process"  << endl;
    draw_bottom_up();
    //exit(0);

    // 2. Find Exact Placement(top down)
    pl.resize(topo->size);
    sol.resize(topo->size);
    //auto& rootMS = vertexMS[topo->root->id];
    auto& rootMS=topo->root->trr.core;
    std::function<void(shared_ptr<TreeNode>)> preOrderTraversal = [&](shared_ptr<TreeNode> curNode) {
        int curId = curNode->id;

        if (curNode->lc != NULL && curNode->rc != NULL) {
            // handle curNode
            if (curNode == topo->root) {
                GridPoint tmp;
                // tmp.x = (rootMS.p1.x + rootMS.p2.x) /2;
                // tmp.y = (rootMS.p1.y + rootMS.p2.y) /2;
                clockSource = rootMS.p1;
                pl[curId] = rootMS.p1;

                //  clockSource = tmp;
                // pl[curId] = tmp;
            } else {
                auto& par = curNode->par;
                int parId = par->id;
                //auto& trr_par = vertexTRR[parId];
                TRR trr_par;
                trr_par.core = Segment(pl[parId], pl[parId]);
                
                trr_par.radius = vertexDistE[curId];
                trr_par.radius=curNode->trr.radius;
                assert(vertexDistE[curId]==curNode->trr.radius);
                //! vertexDistE[curId] should equal to curNode->trr.radius here

                // cout <<std::fixed<< "Before merge: the value for trr_par is" << setprecision(2) << trr_par << endl;
                // if(trr_par.radius == 122663.50){
                //     cout << 3 << endl;
                // }
                //Segment merged = trr_par.intersect(vertexMS[curId]);
                Segment merged = trr_par.intersect(curNode->trr.core);
               
                // if(merged.isLeaf() == false){    
                //     cout << trr_par << " intersecting "<< vertexMS[curId] <<  endl;
                //     cout << " Not leaf" <<endl;
                //     cout << merged << endl;
                // }
                if (merged.id == -1) {
                    draw_TRR_pair(trr_par,TRR(curNode->trr.core,0));
                    cout << "TRR-MS merging failed" << endl;
                    exit(1);
                }
                pl[curId] = merged.p1;
            }

            cout << "Steiner Point " << curId << " located at " << pl[curId] << endl;
            preOrderTraversal(curNode->lc);
            preOrderTraversal(curNode->rc);
        } else {
            // sinks
            //pl[curId] = vertexMS[curId].p1;
            pl[curId]=curNode->trr.core.p1;
            return;
        }
    };
    preOrderTraversal(topo->root);
    cout  << "Finish top-down process"  << endl;

    cout << padding << "Finished DME" << padding << endl;
}

void Router::route() {
    //HC();  // try hierarchical clustering
    DME();
}
bool db_equal(double a, double b) { return abs(a - b) < eps; }
void Router::buildSolution() {
    // preorder traversal to buil grsteiner structure
    std::function<void(shared_ptr<TreeNode>)> preOrderTraversal = [&](shared_ptr<TreeNode> curNode) {
        int curId = curNode->id;
        if (curNode->lc != NULL && curNode->rc != NULL) {
            // handle curNode
            shared_ptr<GrSteiner>& curSteiner = sol[curId];
            auto& lc = curNode->lc;
            auto& rc = curNode->rc;
            shared_ptr<GrSteiner> lcSteiner = make_shared<GrSteiner>(pl[lc->id]);
            shared_ptr<GrSteiner> rcSteiner = make_shared<GrSteiner>(pl[rc->id]);

            // Connect lc
            if (db_equal(curSteiner->x, lcSteiner->x) || db_equal(curSteiner->y, lcSteiner->y)) {
                lcSteiner->set_par(curSteiner);
            } else {  // Use L-shape
                shared_ptr<GrSteiner> middle = make_shared<GrSteiner>(GridPoint(curSteiner->x, lcSteiner->y));
                lcSteiner->set_par(middle);
                middle->set_par(curSteiner);
            }
            if (db_equal(curSteiner->x, rcSteiner->x) || db_equal(curSteiner->y, rcSteiner->y)) {
                rcSteiner->set_par(curSteiner);
            } else {  // Use L-shape
                shared_ptr<GrSteiner> middle = make_shared<GrSteiner>(GridPoint(curSteiner->x, rcSteiner->y));
                rcSteiner->set_par(middle);
                middle->set_par(curSteiner);
            }
            sol[lc->id] = lcSteiner;
            sol[rc->id] = rcSteiner;
            preOrderTraversal(lc);
            preOrderTraversal(rc);
        } else {
            // sinks
            // pl[curId] = vertexMS[curId].p1;
            return;
        }
    };
    sol[topo->root->id] = make_shared<GrSteiner>(pl[topo->root->id]);
    preOrderTraversal(topo->root);
}


// void Router::reportTotalWL() {


//     std::function<void(shared_ptr<GrSteiner>, double&)> preOrderTraversal = [&](shared_ptr<GrSteiner> curNode,double& wl) {

//         if (curNode->lc != NULL && curNode->rc != NULL) {
        
//         }
//         if (curNode->par == NULL) {  // reached source
//             return;
//         }
//         auto &nxtNode = curNode->par;
//         wl += L1Dist(*curNode,*nxtNode);
//         preOrderTraversal(nxtNode, wl);
//     };
//     double total_wl = 0;
    
    
//     // check wirelength
//     cout << padding << "Finish Write Result" << padding << endl;
// }


void Router::writeSolution() {
    ofstream fout(setting.output_file_name);
    if (fout.fail()) {
        cout << "Fail to open file:" << setting.output_file_name << endl;
        exit(1);
    } else {
        cout << padding << "Successfully open input:" << setting.output_file_name << padding << endl;
    }
    double total_wl = 0;
    std::unordered_set< PointPair, PointPairHasher> calculated_edges;

    std::function<void(shared_ptr<GrSteiner>, double& wl)> traceToSource = [&](shared_ptr<GrSteiner> curNode,
                                                                               double& wl) {
        fout << *curNode << " ";
        if (curNode->par == NULL) {  // reached source
            return;
        }
        auto &nxtNode = curNode->par;
        wl += L1Dist(*curNode,*nxtNode);
        PointPair tmp(curNode->x,curNode->y,nxtNode->x,nxtNode->y);
        if(calculated_edges.find(tmp) == calculated_edges.end()){
            calculated_edges.insert(tmp);
            total_wl += L1Dist(*curNode,*nxtNode);
        }
        traceToSource(nxtNode, wl);
    };

    
    vector<double> wirelenghs(taps.size(), 0);
    fout << std::fixed << setprecision(2) << clockSource << endl;
    for (int tapId = 0; tapId < taps.size(); tapId++) {
        fout << tapId << " ";
        traceToSource(sol[tapId], wirelenghs[tapId]);
        cout << std::fixed << setprecision(2) << "WL for tap" << tapId << ": " << wirelenghs[tapId] << endl;
        fout << endl;
    }
    // check wirelength
    cout << "Total Wirelength: " << total_wl << endl;
    cout << padding << "Finish Write Result" << padding << endl;
}

void Router::buildTopology()
{
    system(("rm -rf "+_RUNDIR+"/preOrder.txt").c_str());
    system(("rm -rf "+_RUNDIR+"/inOrder.txt").c_str());
    // ex: ./bin/ntuplace-r -aux ./run_tmp/die0/die0.aux -out ./run_tmp/die0 > ./run_tmp/die0-ntuplace.log
    string cmd = "../bin/ZST_DME " + setting.input_file_name + " " + _RUNDIR+" > "+_RUNDIR+setting.get_case_name()+".log";

    cout << "topology-cmd: " << cmd << "\n";
    
    system(cmd.c_str());
    //exit(0);
    vector<int> preOrderId;
    vector<int> inOrderId;
    vector<pair<int,int>> IdAndLayer;

    ifstream preOrder(_RUNDIR+"preOrder.txt");
    if (preOrder.fail()) {
        cout << "Fail to open file:" << setting.preOrderfile << endl;
        exit(1);
    } else {
        cout << padding << "Successfully open input:" << setting.preOrderfile << padding << endl;
    }

    string preline;
    while (getline(preOrder, preline)) {
        istringstream iss(preline);
        //cout<<preline;
        int Id;
        iss>>Id;
        preOrderId.push_back(Id);
    }
    //cout<<preOrderId.size();
    ifstream inOrder(_RUNDIR+"inOrder.txt");
    if (inOrder.fail()) {
        cout << "Fail to open file:" << setting.inOrderfile << endl;
        exit(1);
    } else {
        cout << padding << "Successfully open input:" << setting.inOrderfile << padding << endl;
    }

    string line;
    while (getline(inOrder, line)) {
        istringstream in(line);
        //cout<<line;
        int Id;
        in>>Id;
        //cout<<"?";
        inOrderId.push_back(Id);
    }
    
    ifstream layer(_RUNDIR+"preOrder.txt");
    if (layer.fail()) {
        cout << "Fail to open file:" << setting.layerfile << endl;
        exit(1);
    } else {
        cout << padding << "Successfully open input:" << setting.layerfile << padding << endl;
    }

    string layerline;
    while (getline(layer, line)) {
        istringstream in(line);
        //cout<<line<<endl;
        int Id;
        int layer;
        in>>Id>>layer;
        //cout<<"?";
        IdAndLayer.push_back(make_pair(Id,layer));
    }
    
    assert(preOrderId.size()==inOrderId.size());
    topo=make_shared<TreeTopology>();
    assert(topo);
    topo->inittree(taps.size(),preOrderId.size(),preOrderId,inOrderId);
    topo->treeLayerCal();
    topo->layerassignment(IdAndLayer);
    //exit(0);
}

void Router::setdelay_model(int delaymodel)
{   
    delay_model=delaymodel;
}

void Router::draw_bottom_up()
{
    string outFile = _DRAWDIR+ setting.get_case_name() + "_bottom_up.plt";
    ofstream outfile( outFile.c_str() , ios::out );

    outfile << " " << endl;
    outfile << "set terminal png size 4000,4000" << endl;
    outfile << "set output " << "\"" << _DRAWDIR << setting.get_case_name()+"_bottom_up"<<".png\"" << endl;
    // outfile << "set multiplot layout 1, 2" << endl;
    outfile << "set size ratio -1" << endl;
    outfile << "set nokey" << endl << endl;

    // for(int i=0; i<cell_list_top.size(); i++){
    //     outfile << "set label " << i + 2 << " \"" << cell_list_top[i]->get_name() << "\" at " << cell_list_top[i]->get_posX() + cell_list_top[i]->get_width() / 2 << "," << cell_list_top[i]->get_posY() + cell_list_top[i]->get_height() / 2 << " center front" << endl;
    // }
    // outfile << "set xrange [0:" << _pChip->get_width() << "]" << endl;
    // outfile << "set yrange [0:" << _pChip->get_height() << "]" << endl;
    // outfile << "plot[:][:] '-' w l lt 3 lw 2, '-' with filledcurves closed fc \"grey90\" fs border lc \"red\", '-' with filledcurves closed fc \"yellow\" fs border lc \"black\", '-' w l lt 1" << endl << endl;

    outfile << "plot[:][:]  '-' w l lt 3 lw 2, '-' with filledcurves closed fc \"grey90\" fs border lc \"red\", '-' w l lt 1" << endl << endl;
    
    outfile << "# TRR" << endl;
    std::function<void(shared_ptr<TreeNode>)> postOrderTraversal = [&](shared_ptr<TreeNode> curNode) {
        int curId = curNode->id;
        if (curNode->lc != NULL && curNode->rc != NULL) {
            postOrderTraversal(curNode->lc);
            postOrderTraversal(curNode->rc);
            curNode->trr.draw_TRR(outfile);
            return;
        } else {
            curNode->trr.draw_TRR(outfile);
        }
    };
    postOrderTraversal(topo->root);
    outfile << "EOF" << endl;

    outfile << "# TRR cores" << endl;
    std::function<void(shared_ptr<TreeNode>)> postOrderTraversal_core = [&](shared_ptr<TreeNode> curNode) {
        int curId = curNode->id;
        if (curNode->lc != NULL && curNode->rc != NULL) {
            postOrderTraversal_core(curNode->lc);
            postOrderTraversal_core(curNode->rc);
            curNode->trr.draw_core(outfile);
            return;
        } else {
            curNode->trr.draw_core(outfile);
        }
    };
    postOrderTraversal_core(topo->root);
    outfile << "EOF" << endl;
    

    // outfile << "pause -1 'Press any key to close.'" << endl;
    outfile.close();

    system(("gnuplot " + outFile).c_str());

    cout << BLUE << "[Router]" << RESET << " - Visualize the bottom_up graph in \'" << outFile << "\'.\n";
}

void Router::draw_solution()
{
    string outFile = _DRAWDIR+ setting.get_case_name() + "_solution.plt";
    ofstream outfile( outFile.c_str() , ios::out );

    outfile << " " << endl;
    outfile << "set terminal png size 4000,4000" << endl;
    outfile << "set output " << "\"" << _DRAWDIR << setting.get_case_name()+"_solution"<<".png\"" << endl;
    // outfile << "set multiplot layout 1, 2" << endl;
    outfile << "set size ratio -1" << endl;
    outfile << "set nokey" << endl << endl;

    // for(int i=0; i<cell_list_top.size(); i++){
    //     outfile << "set label " << i + 2 << " \"" << cell_list_top[i]->get_name() << "\" at " << cell_list_top[i]->get_posX() + cell_list_top[i]->get_width() / 2 << "," << cell_list_top[i]->get_posY() + cell_list_top[i]->get_height() / 2 << " center front" << endl;
    // }
    // outfile << "set xrange [0:" << _pChip->get_width() << "]" << endl;
    // outfile << "set yrange [0:" << _pChip->get_height() << "]" << endl;
    // outfile << "plot[:][:] '-' w l lt 3 lw 2, '-' with filledcurves closed fc \"grey90\" fs border lc \"red\", '-' with filledcurves closed fc \"yellow\" fs border lc \"black\", '-' w l lt 1" << endl << endl;

    outfile << "plot[:][:]  '-' w l lt 3 lw 2, '-' with filledcurves closed fc \"grey90\" fs border lc \"red\", '-' w l lt 1" << endl << endl;
    
    outfile << "# TREE" << endl;
    std::function<void(shared_ptr<GrSteiner>)> traceToSource = [&](shared_ptr<GrSteiner> curNode) {
        if (curNode->par == NULL) {  // reached source
            return;
        }
        auto &nxtNode = curNode->par;
        plotLinePLT(outfile,curNode->x,curNode->y,nxtNode->x,nxtNode->y);
        traceToSource(nxtNode);
    };
    for (int tapId = 0; tapId < taps.size(); tapId++) {
        traceToSource(sol[tapId]);
    }
    outfile << "EOF" << endl;
    
    // outfile << "pause -1 'Press any key to close.'" << endl;
    outfile.close();

    system(("gnuplot " + outFile).c_str());

    cout << BLUE << "[Router]" << RESET << " - Visualize the bottom_up graph in \'" << outFile << "\'.\n";
}

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

void update_merge_Capacitance(shared_ptr<TreeNode> nodeMerge, shared_ptr<TreeNode> nodeLeft, shared_ptr<TreeNode> nodeRight, double ea, double eb){
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

void RLC_calculation(shared_ptr<TreeNode> nodeMerge, shared_ptr<TreeNode> nodeLeft, shared_ptr<TreeNode> nodeRight, double& ea, double& eb){
    float t_a, t_b;
    //! ea for left and eb for right
    t_a = calc_delay_RLC(nodeMerge, nodeLeft, ea);
    t_b = calc_delay_RLC(nodeMerge, nodeRight, eb);
    //需要考虑extra wirelength存在的情况。但是尽量不要引入额外的线长。
    while(fabs(t_a + nodeLeft->delay - (t_b + nodeRight->delay)) >= 1){//? >=1?
        if(t_a + nodeLeft->delay > t_b + nodeRight->delay)// ? does this gurantee that ea>eb?
            modify_coord_by_L1(ea, eb, nodeLeft, nodeRight, skewModifyStep);
        else
            modify_coord_by_L2(ea, eb, nodeLeft, nodeRight, skewModifyStep);// ? does this gurantee that eb>ea?

        t_a = calc_delay_RLC(nodeMerge, nodeLeft, ea);
        t_b = calc_delay_RLC(nodeMerge, nodeRight, eb);
    }
    //printf("原本的 delay: left: %f, right: %f\n", nodeLeft->delay, nodeRight->delay);
    printf("left: %f, right: %f\n t_a: %f, t_b: %f\n total left: %f, total right: %f\n", nodeLeft->delay, nodeRight->delay, t_a, t_b, t_a+nodeLeft->delay, t_b + nodeRight->delay);
    t_a = calc_delay_RLC(nodeMerge, nodeLeft, ea);
    t_b = calc_delay_RLC(nodeMerge, nodeRight, eb);
    printf("RC delay: %f, %f\n", t_a, t_b);

    update_merge_Capacitance(nodeMerge, nodeLeft, nodeRight, ea, eb);
    //update_merge_Delay(nodeMerge, nodeLeft, nodeRight, ea, eb);
    if(t_a + nodeLeft->delay > t_b + nodeRight->delay)
        nodeMerge->delay = t_a + nodeLeft->delay;
    else
        nodeMerge->delay = t_b + nodeRight->delay;
}

//delay_a > delay_b, ea -= x/2, eb+= x/2
void modify_coord_by_L1(double& ea, double& eb, shared_ptr<TreeNode> nodeLeft, shared_ptr<TreeNode> nodeRight, float x){
    double min_manhattan=min_manhattan_dist(nodeLeft,nodeRight); 

    assert(ea!=0||eb!=0);
    
    if(ea!=0 && eb!=0)    // if no node absolutely has extra wirelength(else, one of ea and eb should be 0) one special case: ea/eb=0 and eb/ea=d
    {   if(ea>=x/2&&eb<=min_manhattan-x/2)
        {
            ea -= x/2;
            eb += x/2;
        }
        else if(ea<x/2)// then eb must > min_manhattan-x/2, in this case, either decrease step or allow extra wirelength
        {// here we allow extra wirelength
            double delta=x-ea;
            ea=0;
            eb+=delta;
        }

    }
    else if(ea==0)
    {
        cout<<"should eb has extra wirelength?\n";
    }
    else if(eb==0)
    {
        if((ea-min_manhattan)>x)
        {
            ea-=x;
        }
        else
        {
            double delta=x-(ea-min_manhattan);
            ea-=x;
            eb+=delta;
        }
    }
    // if one node already has extra wirelength(must be ea??? or add assert?)
}

//delay_a < delay_b, ea += x/2, eb -= x/2
void modify_coord_by_L2(double& ea, double& eb, shared_ptr<TreeNode> nodeLeft, shared_ptr<TreeNode> nodeRight, float x){
    assert(nodeLeft!=NULL||nodeRight!=NULL);
    double min_manhattan=min_manhattan_dist(nodeLeft,nodeRight); 

    assert(ea!=0||eb!=0);

    if(ea!=0 && eb!=0)    // if no node absolutely has extra wirelength(else, one of ea and eb should be 0) one special case: ea/eb=0 and eb/ea=d
    {   if(eb>=x/2&&ea<=min_manhattan-x/2)
        {
            eb -= x/2;
            ea += x/2;
        }
        else if(eb<x/2)// then eb must > min_manhattan-x/2, in this case, either decrease step or allow extra wirelength
        {// here we allow extra wirelength
            double delta=x-eb;
            eb=0;
            ea+=delta;
        }

    }
    else if(eb==0)
    {
        cout<<"should ea has extra wirelength?\n";
    }
    else if(ea==0)
    {
        if((eb-min_manhattan)>x)
        {
            eb-=x;
        }
        else
        {
            double delta=x-(eb-min_manhattan);
            eb-=x;
            ea+=delta;
        }
    }
    // if one node already has extra wirelength(must be ea??? or add assert?)
}

double L1Dist(GridPoint p1, GridPoint p2) { return abs(p1.x - p2.x) + abs(p1.y - p2.y); }

double min_manhattan_dist(shared_ptr<TreeNode> nodeLeft, shared_ptr<TreeNode> nodeRight)
{
    auto ms_a = nodeLeft->trr.core;
    auto ms_b = nodeRight->trr.core;
    // get |e_a|, |e_b|
    double d = min(L1Dist(ms_a.p1, ms_b.p1), L1Dist(ms_a.p1, ms_b.p2));
    d = min(d, L1Dist(ms_a.p2, ms_b.p1));
    d = min(d, L1Dist(ms_a.p2, ms_b.p2));  // but why need to calc 2*2 possiblity?
    assert(d>0);
    return d;
}

void TRR::draw_TRR(ofstream& stream)
{
    vector<GridPoint> trr1_boundary_grid;
    if (core.slope() > 0) {
        trr1_boundary_grid.emplace_back(core.p1.x, core.p1.y - radius);
        trr1_boundary_grid.emplace_back(core.p2.x + radius, core.p2.y);
        trr1_boundary_grid.emplace_back(core.p2.x, core.p2.y + radius);
        trr1_boundary_grid.emplace_back(core.p1.x - radius, core.p1.y);  // clock-wise
    } else if (core.slope() < 0) {
        trr1_boundary_grid.emplace_back(core.p1.x + radius, core.p1.y);
        trr1_boundary_grid.emplace_back(core.p2.x, core.p2.y + radius);
        trr1_boundary_grid.emplace_back(core.p2.x - radius, core.p2.y);
        trr1_boundary_grid.emplace_back(core.p1.x, core.p1.y - radius);  // clock-wise
    } else {                                                                            // leaf node
        trr1_boundary_grid.emplace_back(core.p1.x, core.p1.y - radius);
        trr1_boundary_grid.emplace_back(core.p1.x + radius, core.p1.y);
        trr1_boundary_grid.emplace_back(core.p1.x, core.p1.y + radius);
        trr1_boundary_grid.emplace_back(core.p1.x - radius, core.p1.y);  // clock-wise
    }

    plotBoxPLT(stream,
    trr1_boundary_grid[0].x,trr1_boundary_grid[0].y,
    trr1_boundary_grid[1].x,trr1_boundary_grid[1].y,
    trr1_boundary_grid[2].x,trr1_boundary_grid[2].y,
    trr1_boundary_grid[3].x,trr1_boundary_grid[3].y    
    );
}

void TRR::draw_core(ofstream& stream)
{
    plotLinePLT(stream,core.p1.x,
    core.p1.y,
    core.p2.x,
    core.p2.y
    );
}

bool TRR::insideTRR(GridPoint point)
{
    vector<GridPoint> trr1_boundary_grid;
    vector<Segment> trr1_Sides;
    if(db_equal(radius,0))// when radius of TRR==0
    {
        Segment pointSeg=Segment(point,point);
        Segment seg=core.intersect(pointSeg);
        return seg.id==-2;
    }

    if (core.slope() > 0) {
        trr1_boundary_grid.emplace_back(core.p1.x, core.p1.y - radius);
        trr1_boundary_grid.emplace_back(core.p2.x + radius, core.p2.y);
        trr1_boundary_grid.emplace_back(core.p2.x, core.p2.y + radius);
        trr1_boundary_grid.emplace_back(core.p1.x - radius, core.p1.y);  // clock-wise
    } else if (core.slope() < 0) {
        trr1_boundary_grid.emplace_back(core.p1.x + radius, core.p1.y);
        trr1_boundary_grid.emplace_back(core.p2.x, core.p2.y + radius);
        trr1_boundary_grid.emplace_back(core.p2.x - radius, core.p2.y);
        trr1_boundary_grid.emplace_back(core.p1.x, core.p1.y - radius);  // clock-wise
    } else {  // core is leaf node
        trr1_boundary_grid.emplace_back(core.p1.x, core.p1.y - radius);
        trr1_boundary_grid.emplace_back(core.p1.x + radius, core.p1.y);
        trr1_boundary_grid.emplace_back(core.p1.x, core.p1.y + radius);
        trr1_boundary_grid.emplace_back(core.p1.x - radius, core.p1.y);  // clock-wise
    }

    
    for (int i = 0; i < 3; i++) {
        trr1_Sides.emplace_back(trr1_boundary_grid[i], trr1_boundary_grid[i + 1]);
    }
    trr1_Sides.emplace_back(trr1_boundary_grid[3], trr1_boundary_grid[0]);

    vector<double> interceps;// use interceps to determine if a point is in a TRR, just like linear programming, but here all slopes of constraints are 1 or -1

    for(Segment side:trr1_Sides)
    {
        if(side.slope()>0)
        {
            interceps.emplace_back(side.p1.y-side.p1.x);
        }
        else if(side.slope()<0)
        {
            interceps.emplace_back(side.p1.y+side.p1.x);
        }
    }
    sort(interceps.begin(),interceps.end());

    return (point.x+point.y)>=interceps[2]&&(point.x+point.y)<=interceps[3]&&(point.y-point.x)>=interceps[0]&&(point.y-point.x)<=interceps[1];

}

void Router::draw_TRR_pair(TRR trr1,TRR trr2)
{
    string outFile = _DRAWDIR+ "TRR_PAIR"+ ".plt";
    ofstream outfile( outFile.c_str() , ios::out );

    outfile << " " << endl;
    outfile << "set terminal png size 4000,4000" << endl;
    outfile << "set output " << "\"" << _DRAWDIR << "TRR_PAIR"<<".png\"" << endl;
    // outfile << "set multiplot layout 1, 2" << endl;
    outfile << "set size ratio -1" << endl;
    outfile << "set nokey" << endl << endl;

    // for(int i=0; i<cell_list_top.size(); i++){
    //     outfile << "set label " << i + 2 << " \"" << cell_list_top[i]->get_name() << "\" at " << cell_list_top[i]->get_posX() + cell_list_top[i]->get_width() / 2 << "," << cell_list_top[i]->get_posY() + cell_list_top[i]->get_height() / 2 << " center front" << endl;
    // }
    // outfile << "set xrange [0:" << _pChip->get_width() << "]" << endl;
    // outfile << "set yrange [0:" << _pChip->get_height() << "]" << endl;
    // outfile << "plot[:][:] '-' w l lt 3 lw 2, '-' with filledcurves closed fc \"grey90\" fs border lc \"red\", '-' with filledcurves closed fc \"yellow\" fs border lc \"black\", '-' w l lt 1" << endl << endl;

    outfile << "plot[:][:]  '-' w l lt 3 lw 2, '-' with filledcurves closed fc \"grey90\" fs border lc \"red\", '-' w l lt 1" << endl << endl;
    
    outfile << "# TRR" << endl;
    trr1.draw_TRR(outfile);
    trr2.draw_TRR(outfile);
    outfile << "EOF" << endl;

    outfile << "# TRR cores" << endl;
    trr1.draw_core(outfile);
    trr2.draw_core(outfile);
    outfile << "EOF" << endl;
    

    // outfile << "pause -1 'Press any key to close.'" << endl;
    outfile.close();

    system(("gnuplot " + outFile).c_str());

    cout << BLUE << "[Router]" << RESET << " - Visualize the TRR pair in \'" << outFile << "\'.\n";
}