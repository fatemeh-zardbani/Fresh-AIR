#ifndef RTREE_723OCT21_H
#define RTREE_723OCT21_H

// copied from 132Sep stochastic version
// adding minimum threshold
// and no shuffle just move to the end

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <cstdlib>

#include <algorithm>
#include <functional>
#include <vector>

// TODO delete later
#include <iostream>

using namespace std;

// FATEMEH
#include <utility>
#include <queue>
#include <string>
#include <chrono>
#include <string.h>
// FATEMEH OUT
// 2023
#include <valarray>
// #include <bits/stdc++.h>
#include "valarray"
// 2023

#define ASSERT assert // RTree uses ASSERT( condition )
#ifndef Min
#define Min std::min
#endif //Min
#ifndef Max
#define Max std::max
#endif //Max

#ifdef timebreak
    clock_t scan_pending_time = 0.0;
    clock_t delete_pending_time = 0.0;
    clock_t search_time = 0.0;
    clock_t add_ltas_time = 0.0;
    clock_t ripple_time = 0.0;

    // in search
    clock_t tba_assignment_time = 0.0;
    clock_t cracking_time = 0.0;
    clock_t traversedown_time = 0.0;

    // in ripple
    clock_t shuffle_time = 0.0;
#endif 

//
// RTree.h
//

#define RTREE_TEMPLATE template<class DATATYPE, int TMAXNODES, int TMINNODES, int TMAXDATAPOINTS, int TMINDATAPOINTS>
#define RTREE_QUAL RTree<DATATYPE, TMAXNODES, TMINNODES, TMAXDATAPOINTS, TMINDATAPOINTS>

#define RTREE_DONT_USE_MEMPOOLS // This version does not contain a fixed memory allocator, fill in lines with EXAMPLE to implement one.
//#define RTREE_USE_SPHERICAL_VOLUME // Better split classification, may be slower on some systems

// FATEMEH DEFINE
//#define DATA_COUNT 19349214
//#define DATA_COUNT 2000000
// #define DATA_COUNT 70380191
//#define DATA_COUNT 64000000
//#define DATA_COUNT 153206662
//#define DATA_COUNT 75000000
//#define DATA_COUNT 2000000
//#define DATA_COUNT 7550261
//#define DATA_COUNT 70000v14
#ifdef data_size6M
    #define DATA_COUNT 6000000
#elif data_size6_4M
    #define DATA_COUNT 6400000
#elif data_size6_8M
    #define DATA_COUNT 6800000
#elif data_size7_2M
    #define DATA_COUNT 7200000
#elif data_size7_6M
    #define DATA_COUNT 7600000
#elif data_size7_92M
    #define DATA_COUNT 79200000
#elif data_size500K
    #define DATA_COUNT 500000
#elif data_size14_25M
    #define DATA_COUNT 14250000
#elif data_size15_2M
    #define DATA_COUNT 15200000
#elif data_size16_15M
    #define DATA_COUNT 16150000
#elif data_size17_1M
    #define DATA_COUNT 17100000
#elif data_size18_05M
    #define DATA_COUNT 18050000
#elif data_size18_81M
    #define DATA_COUNT 18810000
#endif  
// #define CRACK_SWITCH_SIZE_THRESHOLD 8000
#define INITIAL_HOLES 0

#define NUMDIMS 2

#ifdef dh0
    #define DEFAULT_HOLES 0
#elif dh64
    #define DEFAULT_HOLES 64 
#endif

#ifdef stats
    int count_total_holes = INITIAL_HOLES;
    int count_internal_nodes = 0;
    int count_regular_leaves = 0;
    int count_irregular_leaves = 0;
#endif

//#define DATA_COUNT 10000
//#define NUMDIMS 100

static float m_data_arr_mins[DATA_COUNT * 5][NUMDIMS];
static float m_data_arr_maxes[DATA_COUNT * 5][NUMDIMS];
static int m_data_arr_ids[DATA_COUNT * 5];

// 2023
static float m_pending_insertions_mins[DATA_COUNT][NUMDIMS];
static float m_pending_insertions_maxes[DATA_COUNT][NUMDIMS];
static int m_pending_insertions_ids[DATA_COUNT];

static float tba_mins[DATA_COUNT][NUMDIMS];
static float tba_maxes[DATA_COUNT][NUMDIMS];
static int tba_ids[DATA_COUNT];
// 2023


//static int m_data_arr_datainfo[DATA_COUNT];

//float **m_data_arr_mins = new float*[DATA_COUNT];
//float **m_data_arr_maxes = new float*[DATA_COUNT];
//int m_data_arr_ids[DATA_COUNT];


//float** m_data_arr_mins = new float*[DATA_COUNT];
//float** m_data_arr_maxes = new float*[DATA_COUNT];
//int* m_data_arr_ids = new int[DATA_COUNT];

// FATEMEH OUT


int partition(int *mapping, int low, int high, int offset)
{
    // Choosing the pivot
    int folan;
    int pivot = mapping[high - offset];
 
    // Index of smaller element and indicates
    // the right position of pivot found so far
    int i = (low - 1);
 
    for (int j = low; j <= high - 1; j++) {
 
        // If current element is smaller than the pivot
        if (mapping[j - offset] < pivot) {
 
            // Increment index of smaller element
            i++;
            for(folan = 0; folan< NUMDIMS; folan++){
                swap(tba_mins[i][folan], tba_mins[j][folan]);
                swap(tba_maxes[i][folan], tba_maxes[j][folan]);
            }
            swap(tba_ids[i], tba_ids[j]);
            swap(mapping[i-offset],mapping[j-offset]);
        }
    }
    for(folan = 0; folan< NUMDIMS; folan++){
        swap(tba_mins[(i + 1)][folan], tba_mins[high][folan]);
        swap(tba_maxes[(i + 1)][folan], tba_maxes[high][folan]);
    }
    swap(tba_ids[i + 1], tba_ids[high]);
    swap(mapping[i + 1 - offset], mapping[high - offset]);
    return (i + 1);
}


void quickSort(int *mapping, int low, int high, int offset)
{
    if (low < high) {
 
        // pi is partitioning index, arr[p]
        // is now at right place
        int pi = partition(mapping, low, high, offset);
 
        // Separately sort elements before
        // partition and after partition
        quickSort(mapping, low, pi - 1, offset);
        quickSort(mapping, pi + 1, high, offset);
    }
}



// FATEMEH NOTES:
// QueryAdaptive: never cretaes the Q_rect
// QueryAdaptive_v2: Always creates the q_rect, adds to NN
// QueryAdaptive_v3: sometimes creates the q_rect, adds to NN
// QueryAdaptive_v4: always, just inserts

// TODO there is something wrong with the RemoveBtch function in the reinsert part, have to fix that
// Only qa uses the 2dc_v5, the rest are wrong

// 2dc_v4, multiple passes, but ok
// 2dc_v5, one pass all done!

/// \class RTree
/// Implementation of RTree, a multidimensional bounding rectangle tree.
/// Example usage: For a 3-dimensional tree use RTree<Object*, float, 3> myTree;
///
/// This modified, templated C++ version by Greg Douglas at Auran (http://www.auran.com)
///
/// DATATYPE Referenced data, should be int, void*, obj* etc. no larger than sizeof<void*> and simple type
/// float Type of element such as int or float
/// NUMDIMS Number of dimensions such as 2 or 3
/// float Type of element that allows fractional and large values such as float or double, for use in volume calcs
///
/// NOTES: Inserting and removing data requires the knowledge of its constant Minimal Bounding Rectangle.
///        This version uses new/delete for nodes, I recommend using a fixed size allocator for efficiency.
///        Instead of using a callback function for returned results, I recommend and efficient pre-sized, grow-only memory
///        array similar to MFC CArray or STL Vector for returning search query result.
///
template<class DATATYPE,
        int TMAXNODES = 8, int TMINNODES = TMAXNODES / 2, int TMAXDATAPOINTS = 256, int TMINDATAPOINTS = TMAXDATAPOINTS / 2>
                class RTree {
                protected:

                    struct Node;  // Fwd decl.  Used by other internal structs and iterator

                public:

                    // These constant must be declared after Branch and before Node struct
                    // Stuck up here for MSVC 6 compiler.  NSVC .NET 2003 is much happier.
                    enum {
                        MAXNODES = TMAXNODES,                         ///< Max elements in node
                        MINNODES = TMINNODES,                         ///< Min elements in node
                        MAXDATAPOINTS = TMAXDATAPOINTS,
                        MINDATAPOINTS = TMINDATAPOINTS
                    };



                public:

                    RTree();

                    RTree(const RTree &other);

                    RTree(const float data_set_min[][NUMDIMS], const float data_set_max[][NUMDIMS], int data_set_ids[],
                          int data_set_size);

                    // initialize rtree from data_file
                    // also do it with static memory
                    RTree(std::string data_file_name);

                    virtual ~RTree();

                    /// Insert entry
                    /// \param a_min Min of bounding rect
                    /// \param a_max Max of bounding rect
                    /// \param a_dataId Positive Id of data.  Maybe zero, but negative numbers not allowed.
                    void Insert_old(const float a_min[NUMDIMS], const float a_max[NUMDIMS], const DATATYPE &a_dataId);

                    /// \param a_min Min of bounding rect
                    /// \param a_max Max of bounding rect
                    /// \param a_dataId Positive Id of data.  Maybe zero, but negative numbers not allowed.
                    void Insert(const float a_min[NUMDIMS], const float a_max[NUMDIMS], const DATATYPE &a_dataId);

                    /// Remove entry
                    /// \param a_min Min of bounding rect
                    /// \param a_max Max of bounding rect
                    /// \param a_dataId Positive Id of data.  Maybe zero, but negative numbers not allowed.
                    void Remove(const float a_min[NUMDIMS], const float a_max[NUMDIMS], const DATATYPE &a_dataId);

                    /// Find all within search rectangle
                    /// \param a_min Min of search bounding rect
                    /// \param a_max Max of search bounding rect
                    /// \param a_searchResult Search result array.  Caller should set grow size. Function will reset, not append to array.
                    /// \param a_resultCallback Callback function to return result.  Callback should return 'true' to continue searching
                    /// \param a_context User context to pass as parameter to a_resultCallback
                    /// \return Returns the number of entries found
                    int Search(const float a_min[NUMDIMS], const float a_max[NUMDIMS]);

                    // FATEMEH
                    int CountDataPoints();
                    int SumDataPointIDs();
                    int CountNodes();
                    int TreeHeight();
                    std::pair<int, int> CountLeaves();
                    // FATEMEH OUT

                    int QueryAdaptive(const float a_min[NUMDIMS], const float a_max[NUMDIMS]);

                    // FATEMEH TEST
                    int QueryAdaptive_v2(const float a_min[NUMDIMS], const float a_max[NUMDIMS]);
                    int QueryAdaptive_v3(const float a_min[NUMDIMS], const float a_max[NUMDIMS]);
                    int QueryAdaptive_v4(const float a_min[NUMDIMS], const float a_max[NUMDIMS]);
                    // for search_9 and 2dc8
                    int QueryAdaptive_v5(const float a_min[NUMDIMS], const float a_max[NUMDIMS]);
                    // search_10 and 2dc9
                    int QueryAdaptive_v6(const float a_min[NUMDIMS], const float a_max[NUMDIMS]);
                    // search_11 and choose which 2dc(8 or 9)
                    int QueryAdaptive_v7(const float a_min[NUMDIMS], const float a_max[NUMDIMS]);
                    // search_12 and choose which 2dc8, but slightly different conditions
                    int QueryAdaptive_v8(const float a_min[NUMDIMS], const float a_max[NUMDIMS]);
                    // search_14 and choose which 2dc10
                    int QueryAdaptive_v9(const float a_min[NUMDIMS], const float a_max[NUMDIMS]);
                    // search_15 and choose which 2dc8 for 3d data
                    int QueryAdaptive_v10(const float a_min[NUMDIMS], const float a_max[NUMDIMS]);
                    // for the _v8 that resturns results
                    int QueryAdaptive_v11(const float a_min[NUMDIMS], const float a_max[NUMDIMS]);
                    // FATEMEH TEST OUT

                    // ELEGANT STUFF
                    int QueryAdaptive_v12(const float a_min[NUMDIMS], const float a_max[NUMDIMS]);
                    // elegant v2. with search_20, add_ltas_v2
                    int QueryAdaptive_v13(const float a_min[NUMDIMS], const float a_max[NUMDIMS]);

                    // 2023
                    // with the to_be_inserted stuff
                    int QueryAdaptive_v14(const float a_min[NUMDIMS], const float a_max[NUMDIMS]);

                    void getLeafArea(string filename);
                    // 2023

                    /// Remove all entries from tree
                    void RemoveAll();

                    /// Count the data elements in this container.  This is slow as no internal counter is maintained.
                    int Count();



                    /// Create 2d-crack defined in adaptive tree creation, temporarily public
                    //    void TwoDCrack(Node* node_parent_of_start, int parent_of_start_branch_index, Node* start, const float a_min[NUMDIMS], const float a_max[NUMDIMS]);

                    // tester
                    bool SearchNNtester(const float a_min[NUMDIMS], const float a_max[NUMDIMS]);

                    // print out the tree structure
                    void printTree(string file_name);

                    // print for debug
                    void printTreeDebug(string file_name);

                    void printLeafSizes(string file_name);

                    // 2023
                    void printLeafLinkedList();
                    void printIDs();
                    int getPendingCount(){return m_pending_insertions_count;}
                    int getRightestRight(){return m_youngest_leaf->m_R;}
                    void printLeafBoundsOfData(int data_index);
                    bool printLeafBoundsOfDataRec(Node *a_node, int data_index);
                    void findAllLeavesWithDataInThem(int data_index);
                    void lookForEmptyNodes();
                    // 2023

                    // testing, should not be public
                    //    bool Insert_anylevel(const Branch &a_branch, Node *start, int a_level);



                protected:

                    /// Minimal bounding rectangle (n-dimensional)
                    struct Rect {
                        Rect(){}
                        Rect(float mins[NUMDIMS], float maxes[NUMDIMS]){
                            for(int i=0; i< NUMDIMS; i++) {
                                m_min[i] = mins[i];
                                m_max[i] = maxes[i];
                            }
                        }
                        float m_min[NUMDIMS];                      ///< Min dimensions of bounding box
                        float m_max[NUMDIMS];                      ///< Max dimensions of bounding box
                    };

                    /// May be data or may be another subtree
                    /// The parents level determines this.
                    /// If the parents level is 0, then this is data
                    struct Branch {
                        Rect m_rect;                                  ///< Bounds
                        Node *m_child;                                ///< Child node
                        DATATYPE m_data;                              ///< Data Id
                        Branch(){}
                        Branch(Node *the_node, Rect the_rect, DATATYPE the_data){
                            m_child = the_node;
                            m_data = the_data;
                            for(int i=0; i< NUMDIMS; i++) {
                                m_rect.m_min[i] = the_rect.m_min[i];
                                m_rect.m_max[i] = the_rect.m_max[i];
                            }
                        }
                        bool operator<(const Branch& a) const{
                            // cout << "calling < operaor in branch" << endl;
                            return m_child->m_L < a.m_child->m_L;
                        }
                    };

                    /// Node for each branch level
                    struct Node {
                        bool IsInternalNode() { return (m_level > 0); } // Not a leaf, but a internal node
                        bool IsLeaf() { return (m_level == 0); } // A leaf, contains data
                        // 2023
                        // we always use the isReg function to check regularity
                        // to make sure the holes are always correct, let's just check here
                        // bool isRegular() { return m_regular; } // leaf is regular is m_regular is true, irregular ow
                        // we also only ask leaves
                        bool isRegular() { return ((m_R - m_L - m_holes) <= MAXDATAPOINTS); } // leaf is regular is m_regular is true, irregular ow

                        // for debug
                        std::string name;
                        int m_count;                                  ///< Count
                        int m_level;                                  ///< Leaf is zero, others positive
                        // 2023:  SHOULD NOT USE THIS M_REGULAR THING.
                        // IT IS CORRECT I THINK
                        // BUT USE THE ISREGULAR FUNCTION!!
                        bool m_regular;                               ///< Some leaves can be irregular. Leaf is irregular if this is false, regular ow
                        //        std::vector<Branch> m_branch;                 ///< Branch
                        Branch m_branch[MAXNODES];                    ///< Branch
//                        Branch m_data_points[MAXDATAPOINTS];            ///< data points, couldnt keep them in m_branch since the limits differ
                        int m_L;                                        ///< left side of data_points interval in m_data_arr
                        int m_R;                                        ///< right side of data_points interval in m_data_arr
                        bool m_transferred;                             ///< whether data has been copied to m_data_points, if true has been copied, if false data is in m_data_arr
                        Node* m_parent;
                        int m_id;

                        // MARCH 2023
                        // pointer to leaf on the right side
                        Node* m_younger_brother;
                        // pointer to leaf on the left side
                        Node* m_older_brother;
                        // number of empty spaces at the front of this leaf's L.
                        int m_holes;
                        // MARCH 2023


                        Node(){}
                    };

                    /// A link list of nodes for reinsertion after a delete operation
                    struct ListNode {
                        ListNode *m_next;                             ///< Next in list
                        Node *m_node;                                 ///< Node
                    };

                    /// Variables for finding a split partition
                    struct PartitionVars {
                        enum {
                            NOT_TAKEN = -1
                        }; // indicates that position

                        int m_partition[MAXDATAPOINTS + 1];
                        //        std::vector<int> m_partition;
                        int m_total;
                        int m_minFill;
                        int m_count[2];
                        Rect m_cover[2];
                        float m_area[2];

                        Branch m_branchBuf[MAXDATAPOINTS + 1];
                        //        std::vector<Branch> m_branchBuf;
                        int m_branchCount;
                        Rect m_coverSplit;
                        float m_coverSplitArea;
                    };

                    // FATEMEH TYPEDEF
                    typedef std::pair<bool,bool> bools;
                    typedef std::pair<std::vector<Branch>, std::vector<int>> tbd_vectors;
                    typedef std::pair<int, int> position;
                    typedef std::pair<std::vector<Branch>, std::pair<position, Rect>> fragment;
                    typedef std::vector<fragment> fragments;
                    typedef std::pair<Node*, std::pair<tbd_vectors, fragments>> gathered_node_info;
                    //    typedef std::pair<Node*, std::pair<tbd_vectors, fragment[]>> gathered_node_info;

                    struct Fragment{
                        int position_row;
                        int position_colomn;
                        std::vector<Branch> data_points;
                        Rect cover;
                    };

                    struct NodeInfo{
                        Node* this_node;
                        std::vector<Branch> tbd_branches;
                        std::vector<int> tbd_indexes;
                        Fragment fragments[9];
                    };

                    struct NodeCracks{
                        Node* this_node;

                        //        int cracks[2*NUMDIMS + 3 - 1];
                        //        int crack_covers[2*NUMDIMS + 3][2 * NUMDIMS];

                        int cracks[2*NUMDIMS + 1 - 1];
                        float crack_covers_mins[2*NUMDIMS + 1][NUMDIMS];
                        float crack_covers_maxes[2*NUMDIMS + 1][NUMDIMS];
                        bool three_or_five; // true = 3, false = 5

                    };



                    // fatemeh's doing
                    struct CrackVars {
                        //        std::vector<Branch> data_points; // data points in irregular leaves, that intersect with the query rect
                        //        std::vector<Node*> data_point_parents; // each irregular leaf the data was originally in
                        // FATEMEH NEW
                        std::vector<gathered_node_info> gathered_nodes;

                        std::vector<NodeCracks> gathered_node_cracks;
                        // FATEMEH NEW OUT
                        //        int irregulars_encoutered_count = 0; // count of irregular rectangles encountered while searching
                        //        Node *last_irregular; // last irregular leaf intersecting with the query rect
                        //        Node *last_irregular_parent = NULL; // parent of last irregular leaf intersecting with the query rect
                        //        Node *potential_parent = NULL; // bad coding
                        //        int potential_branch_index; // bad coding
                        //        int parent_of_last_irregular_branch_index; // branch index of last irregular in parent of last irregular leaf intersecting with the query rect
                        bool encountered_regular = false;
                        //        CrackVars() : irregulars_encoutered_count(0) {}


                        ~CrackVars(){
                            std::vector<gathered_node_info>().swap(gathered_nodes);
                        }
                    };

                    // ELEGANT SETUP
                    struct LTA{
                        Node* this_leaf;
                        vector<Branch> pieces_created;
                        LTA(Node* folan, vector<Branch> filan){
                            this_leaf = folan;
                            pieces_created = filan;
                        }
                    };

                    typedef vector<LTA> LeavesToAdd;

                    // 2023

                    struct LTA_v2{
                        Node* this_leaf;
                        int Ls[2*NUMDIMS + 1 + 1];
                        int Rs[2*NUMDIMS + 1 + 1];
                        int holes[2*NUMDIMS + 1 + 1];
                        float crack_covers_mins[2*NUMDIMS + 1 + 1][NUMDIMS];
                        float crack_covers_maxes[2*NUMDIMS + 1 + 1][NUMDIMS];
                        int how_many_created = 0;
                        // need for later
                        int qp_L;
                        // vector<Branch> qp_tba;
                        // float* tba_mins;
                        // float* tba_maxes;
                        // int* tba_ids;
                        int tba_start;
                        int tba_count;

                    };
                    // 2023

                    typedef vector<LTA_v2> LeavesToAdd_v2;


                    // END ELEGANT

                    // 2023
                    struct overlapping_leaf_tbas{
                        Node *this_leaf;
                        // vector<Branch> this_tba;
                        // float* this_tba_mins;
                        // float* this_tba_maxes;
                        // int* this_tba_ids;
                        int this_tba_start;
                        int this_tba_count;
                        overlapping_leaf_tbas(Node* l, int tba_start, int tba_count){
                            this_leaf = l;
                            // this_tba = tba;
                            // this_tba_mins = tba_mins;
                            // this_tba_maxes = tba_maxes;
                            // this_tba_ids = tba_ids;
                            this_tba_start = tba_start;
                            this_tba_count = tba_count;
                        }
                    };
                    // 2023

                    struct NodeCracks_9{
                        Node* this_node;
                        int cracks[8];
                        float crack_covers_mins[9][NUMDIMS];
                        float crack_covers_maxes[9][NUMDIMS];


                    };

                    struct CrackVars_9{
                        std::vector<NodeCracks_9> gathered_node_cracks;
                        bool encountered_regular = false;
                    };

                    // FATEMEH NN
                    struct CompareNode {
                        int dist(Rect q_rect, Rect r) {
                            // turn the query rectangle into a point
                            int q_point[NUMDIMS];

                            for (int i = 0; i < NUMDIMS; i++) {
                                q_point[i] = (q_rect.m_min[i] + q_rect.m_max[i]) / 2;
                            }
                            // calculate the MINDIST of the q_point and r
                            // as defined by Nick Roussopoulos Stephen Kelley Frederic Vincent
                            int distance = 0;
                            for (int i = 0; i < NUMDIMS; i++) {
                                if (q_point[i] < r.m_min[i]) { distance += ((q_point[i] - r.m_min[i]) * (q_point[i] - r.m_min[i])); }
                                else if (q_point[i] > r.m_max[i]) {
                                    distance += ((q_point[i] - r.m_max[i]) * (q_point[i] - r.m_max[i]));
                                }
                            }

                            return distance;
                        }

                        bool operator()(const pair<Node *, Rect> &lhs, const pair<Node *, Rect> &rhs) {
                            Rect query_rect = lhs.second;
                            //            printf("calling NC in CompareNode operator 1\n");
                            int dist1 = dist(query_rect, RTREE_QUAL::NodeCover(lhs.first));
                            //            printf("calling NC in CompareNode operator 2\n");
                            int dist2 = dist(query_rect, RTREE_QUAL::NodeCover(rhs.first));
                            return dist1 > dist2;
                        }

                    };

                    struct BranchList{
                        std::vector<Branch> branch_list;
                    };


                    // FATEMEH OUT

                    Node *AllocNode();

                    void FreeNode(Node *a_node);

                    void InitNode(Node *a_node);

                    void InitRect(Rect *a_rect);

                    bool InsertRectRec(const Branch &a_branch, Node *a_node, Node **a_newNode, int a_level);

                    bool InsertRectRec(const Branch &a_branch, Node *a_node, Node **a_newNode, int a_level, Rect &a_coverRect, Rect &a_newCoverRect);


                    bool InsertRect(const Branch &a_branch, Node **a_root, int a_level);

                    // FATEMEH
                    bool InsertRectRecReg(const Branch &a_branch, Node *a_node, Node **a_newNode, int a_level);

                    bool InsertRectReg(const Branch &a_branch, Node **a_root, int a_level);


                    bools InsertRectOnlyRegRec(const Branch &a_branch, Node *a_node, Node **a_newNode, int a_level);

                    bools InsertRectOnlyRegRec(const Branch &a_branch, Node *a_node, Node **a_newNode, int a_level, Rect &a_coverRect, Rect &a_newCoverRect);



                    // returns whether we inserted or not
                    bools InsertRectOnlyReg(const Branch &a_branch, Node **a_root, int a_level);

                    // FATEMEH OUT

                    // FATEMEH
                    void getLeafArea(ofstream &file, Node *a_node);
                    void Insert_anylevel(const Branch &a_branch, Node *start, int a_level);
                    int getNodeBranchIndex(Node* a_node);
                    // FATEMEH OUT

                    Rect NodeCover(Node *a_node);

                    //11 oct
                    Rect LeafCover(int L, int R);

                    // FATEMEH DEBUG
                    static Rect FatemehCover(std::vector<Branch> m_branch);

                    bool AddBranch(const Branch *a_branch, Node *a_node, Node **a_newNode);


                    bool AddBranch(const Branch *a_branch, Node *a_node, Node **a_newNode, Rect &a_coverRect, Rect &a_newCoverRect);

                    void DisconnectBranch(Node *a_node, int a_index);

                    int PickBranch(const Rect *a_rect, Node *a_node);

                    // FATEMEH
                    int PickBranchReg(const Rect *a_rect, Node *a_node);
                    // FATEMEH OUT

                    static Rect CombineRect(const Rect *a_rectA, const Rect *a_rectB);

                    void SplitNode(Node *a_node, const Branch *a_branch, Node **a_newNode);

                    void SplitNode(Node *a_node, const Branch *a_branch, Node **a_newNode, Rect &a_coverRect, Rect &a_newCoverRect);

                    float RectSphericalVolume(Rect *a_rect);

                    float RectVolume(Rect *a_rect);

                    float CalcRectVolume(Rect *a_rect);

                    void GetBranches(Node *a_node, const Branch *a_branch, PartitionVars *a_parVars);

                    void ChoosePartition(PartitionVars *a_parVars, int a_minFill);

                    void LoadNodes(Node *a_nodeA, Node *a_nodeB, PartitionVars *a_parVars);

                    void InitParVars(PartitionVars *a_parVars, int a_maxRects, int a_minFill);

                    void PickSeeds(PartitionVars *a_parVars);

                    void Classify(int a_index, int a_group, PartitionVars *a_parVars);

                    bool RemoveRect(Rect *a_rect, const DATATYPE &a_id, Node **a_root);

                    bool RemoveBatchFromNode(Node *a_node, std::vector<Branch> a_ids, std::vector<int> a_indexes, Node **a_root);

                    bool RemoveRectRec(Rect *a_rect, const DATATYPE &a_id, Node *a_node, ListNode **a_listNode);



                    bools RemoveBatchRec(Node *a_from_node, std::vector<Branch> a_ids, std::vector<int> a_indexes, Node *a_node, ListNode **a_listNode);

                    ListNode *AllocListNode();

                    void FreeListNode(ListNode *a_listNode);

                    bool Overlap(Rect *a_rectA, Rect *a_rectB) const;

                    bool Overlap(Rect *a_rectA, const float a_min[NUMDIMS], const float a_max[NUMDIMS]) const;

                    void ReInsert(Node *a_node, ListNode **a_listNode);

                    bool Search(Node *a_node, Rect *a_rect, int &a_foundCount, CrackVars *crack_vars);

                    // FATEMEH
                    bool Search_2(Node *a_node, Rect *a_rect, int &a_foundCount, CrackVars *crack_vars);

                    bool Search_old(Node *a_node, Rect *a_rect, int &a_foundCount, CrackVars *crack_vars);
                    bool Search_3(Node *a_node, Rect *a_rect, int &a_foundCount, CrackVars *crack_vars);

                    bool Search_4(Node *a_node, Rect *a_rect, int &a_foundCount, CrackVars *crack_vars);

                    bool Search_5(Node *a_node, Rect *a_rect, int &a_foundCount, CrackVars *crack_vars);
                    bool Search_6(Node *a_node, Rect *a_rect, int &a_foundCount, CrackVars *crack_vars);
                    bool Search_7(Node *a_node, Rect *a_rect, int &a_foundCount, CrackVars *crack_vars);

                    // for 3d
                    bool Search_8(Node *a_node, Rect *a_rect, int &a_foundCount, CrackVars *crack_vars);

                    // search 6 without scracks and with inline node-cover, 5 cracks
                    bool Search_9(Node *a_node, Rect *a_rect, int &a_foundCount, CrackVars *crack_vars);
                    // 3 cracks
                    bool Search_10(Node *a_node, Rect *a_rect, int &a_foundCount, CrackVars *crack_vars);
                    // sometimes 3 sometimes 5 cracks
                    bool Search_11(Node *a_node, Rect *a_rect, int &a_foundCount, CrackVars *crack_vars);
                    // the new cracking strategy: the one similar to kd-tree TODO: the size checks are incomplete so it's not okay!!, DONE
                    bool Search_12(Node *a_node, Rect *a_rect, int &a_foundCount, CrackVars *crack_vars);

                    // same as 12, just no checking the min condition
                    bool Search_122(Node *a_node, Rect *a_rect, int &a_foundCount, CrackVars *crack_vars);
                    // same as 12, just using the weird heuristic to figure out what is wrong with edges performance
                    bool Search_123(Node *a_node, Rect *a_rect, int &a_foundCount, CrackVars *crack_vars);

                    // baseline, 9 way cracking in static setting
                    bool Search_14(Node *a_node, Rect *a_rect, int &a_foundCount, CrackVars_9 *crack_vars);

                    // for 3d, 7 ind cracks
                    bool Search_15(Node *a_node, Rect *a_rect, int &a_foundCount, CrackVars *crack_vars);

                    // for 3d, exactly 15, just no checking if the pieces are too small
                    bool Search_152(Node *a_node, Rect *a_rect, int &a_foundCount, CrackVars *crack_vars);



                    //exactly search_12, but returns the results
                    bool Search_16(Node *a_node, Rect *a_rect, vector<DATATYPE> &results, CrackVars *crack_vars);


                    // ELEGANT SETUP
                    // for elegant version, this one just finds the irreg, calls do_2dc, which returns some branches
                    // the branches are then added to the tree in add_leaves_to_add
                    bool Search_17(Node *a_node, Rect *a_rect, int &a_foundCount, LeavesToAdd *all_lta);

                    int Elegant_2dc(Node *a_node, Rect *a_rect, LeavesToAdd *all_lta);

                    // 2023
                    // void TraverseDownTreeTilLeaf(Node* a_node, vector<Branch> *tba, vector<overlapping_leaf_tbas> *overlapping_leaves);
                    // void TraverseDownTreeTilLeaf(Node* a_node, float* tba_mins, float* tba_maxes, int* tba_ids, int tba_count, vector<overlapping_leaf_tbas> *overlapping_leaves);
                    void TraverseDownTreeTilLeaf(Node* a_node, int tba_start, int tba_count, vector<overlapping_leaf_tbas> *overlapping_leaves);

                    // 2023


                    bool Search_18(Node *a_node, Rect *a_rect, int &a_foundCount, LeavesToAdd *all_lta);
                    int Elegant_2dc_v2(Node *a_node, Rect *a_rect, LeavesToAdd *all_lta);
                    // exactly the same as 18, just not calling elegant_2dc_v2, instead copy the code in
                    bool Search_19(Node *a_node, Rect *a_rect, int &a_foundCount, LeavesToAdd *all_lta);
                    // same as 19, but with no min_l check
                    bool Search_192(Node *a_node, Rect *a_rect, int &a_foundCount, LeavesToAdd *all_lta);
                    // same as 192, just one this_piece and not a piece_for_crack and piece_for_choice
                    bool Search_1925(Node *a_node, Rect *a_rect, int &a_foundCount, LeavesToAdd *all_lta);
                    // same as 19 --> wml, just no special treatment
                    bool Search_193(Node *a_node, Rect *a_rect, int &a_foundCount, LeavesToAdd *all_lta);

                    // wml, no special treatment, but with the lta_v2, works with add_ltas_v2, based off 19
                    bool Search_20(Node *a_node, Rect *a_rect, int &a_foundCount, LeavesToAdd_v2 *all_lta);
                    // woml, ns, lta_v2, add_ltas_v2, based off 1925
                    bool Search_2025(Node *a_node, Rect *a_rect, int &a_foundCount, LeavesToAdd_v2 *all_lta);
                    // same as 2025, but with weird choice
                    bool Search_2026(Node *a_node, Rect *a_rect, int &a_foundCount, LeavesToAdd_v2 *all_lta);

                    // 2023
                    // for the pending_insertions things
                    // bool Search_2027(Node *a_node, Rect *a_rect, int &a_foundCount, LeavesToAdd_v2 *all_lta, vector<Branch> *tba, vector<overlapping_leaf_tbas> *overlapping_leaves);
                    // bool Search_2027(Node *a_node, Rect *a_rect, int &a_foundCount, LeavesToAdd_v2 *all_lta, float* tba_mins, float* tba_maxes, int* tba_ids, int tba_count, vector<overlapping_leaf_tbas> *overlapping_leaves);
                    bool Search_2027(Node *a_node, Rect *a_rect, int &a_foundCount, LeavesToAdd_v2 *all_lta, int tba_start, int tba_count, vector<overlapping_leaf_tbas> *overlapping_leaves);



                    void Add_ltas(LeavesToAdd *all_lta);
                    // sort the branches based on their R value

                    static bool compareByL(const Branch a, const Branch b);

                    // to work with lta_v2
                    void Add_ltas_v2(LeavesToAdd_v2 *all_lta);
                    // 2023
                    void Add_ltas_v3(LeavesToAdd_v2 *all_lta, vector<overlapping_leaf_tbas> *overlapping_leaves);
                    // the rippling things
                    void ripple(vector<overlapping_leaf_tbas> *overlapping_leaves);
                    // new psuedo code, ripple v4 in overleaf
                    
                    void deleteEmptyNodesUp(Node* start);
                    bool is_leaf_in_OL(Node *a_leaf, vector<overlapping_leaf_tbas> *overlapping_leaves);
                    // void shuffle_right(vector<overlapping_leaf_tbas> *overlapping_leaves, Node *this_leaf, vector<Branch> this_tba, int k);
                    // void shuffle_left(vector<overlapping_leaf_tbas> *overlapping_leaves, Node *this_leaf, vector<Branch> this_tba, int k);
                    
                    void shuffle_right(vector<overlapping_leaf_tbas> *overlapping_leaves, Node *this_leaf, int this_tba_start, int k);
                    void shuffle_left(vector<overlapping_leaf_tbas> *overlapping_leaves, Node *this_leaf, int this_tba_start, int k);
                    void ripple_v2(vector<overlapping_leaf_tbas> *overlapping_leaves);

                    void ripple_v3(vector<overlapping_leaf_tbas> *overlapping_leaves);






                    // FATEMEH OUT


                    // FATEMEH
                    bool CountDataPoints(Node *a_node, int &a_foundCount);
                    bool SumDataPointIDs(Node *a_node, int &a_foundCount);
                    bool CountNodes(Node *a_node, int &a_foundCount);
                    void PrintLeafSizesRec(Node *a_node, ofstream &myfile);
                    void lookForEmptyNodesRec(Node* a_node);
                    void findNodeCount(int node_id, Node* this_node);
                    // FATEMEH OUT

                    bool SearchNN(Rect query_rect, Node **result);

                    void RemoveAllRec(Node *a_node);

                    void Reset();

                    void CountRec(Node *a_node, int &a_count);

                    void CopyBranch(Branch &current, const Branch &other);


                    /// Create 2d-crack defined in adaptive tree creation, temporarily public, let's see how this one works
                    void TwoDCrack_v2(std::vector<Branch> data_points, Node *node_parent_of_start, int parent_of_start_branch_index,
                                      Node *start, const float a_min[NUMDIMS], const float a_max[NUMDIMS], BranchList* reinsert_list);

                    void TwoDCrack_v3(std::vector<Branch> data_points, Node *node_parent_of_start, int parent_of_start_branch_index,
                                      Node *start, const float a_min[NUMDIMS], const float a_max[NUMDIMS], BranchList* reinsert_list);

                    // FATEMEH NEW
                    // returns whether or not one of the leaves created was regular
                    bool TwoDCrack_v4(Node *start, const float a_min[NUMDIMS], const float a_max[NUMDIMS], BranchList* reinsert_list);
                    bool TwoDCrack_v5(Node *start, fragments frags, const float a_min[NUMDIMS], const float a_max[NUMDIMS], BranchList* reinsert_list);
                    bool TwoDCrack_v6(Node *start, fragment frags[], const float a_min[NUMDIMS], const float a_max[NUMDIMS], BranchList* reinsert_list);
                    bool TwoDCrack_v7(Node *start, int crack_indexes[], const float a_min[NUMDIMS], const float a_max[NUMDIMS]);

                    // for search_9, without scrack and nodecover in line
                    bool TwoDCrack_v8(Node *start, int crack_indexes[2*NUMDIMS], float cover_mins[2*NUMDIMS + 1][NUMDIMS], float cover_maxes[2*NUMDIMS + 1][NUMDIMS], const float a_min[NUMDIMS], const float a_max[NUMDIMS]);
                    // for search_10, same as search_9, just 3-way crack
                    bool TwoDCrack_v9(Node *start, int crack_indexes[2*NUMDIMS], float cover_mins[2*NUMDIMS + 1][NUMDIMS], float cover_maxes[2*NUMDIMS + 1][NUMDIMS], const float a_min[NUMDIMS], const float a_max[NUMDIMS]);

                    //9 way crack search_14
                    bool TwoDCrack_v10(Node *start, int crack_indexes[8], float cover_mins[9][NUMDIMS], float cover_maxes[9][NUMDIMS], const float a_min[NUMDIMS], const float a_max[NUMDIMS]);


                    // FATEMEH NEW OUT

                    // FATEMEH
                    int ChooseAxis(Rect node_cover, Rect query_cover);
                    int ChooseAxisAgain(Rect node_cover, Rect query_cover, int ignore_this_dim);

                    // returns 2*axis + {0 if min or 1 if max}, for all choices ignore = -3
                    int ChooseAxisLargestSideMid(Rect node_cover, Rect query_cover);

                    // figuring out what heuristic edges results was using, so this is ChooseAxis and then choosing min or max
                    int ChooseAxis_weird(Rect node_cover, Rect query_cover);



                    // returns the index where the crack is at
                    int CrackOnAxisCompMin(int L, int R, const float& crack_value, int axis);
                    int CrackOnAxisCompMax(int L, int R, const float& crack_value, int axis);

                    // for in line node cover
                    int CrackOnAxisCompMin_v2(int L, int R, const float& crack_value, int axis, float *cover_mins_left, float *cover_maxes_left, float *cover_mins_right, float *cover_maxes_right);
                    int CrackOnAxisCompMax_v2( int L, int R, const float& crack_value, int axis, float *cover_mins_left, float *cover_maxes_left, float *cover_mins_right, float *cover_maxes_right);

                    // just different way of returning the right and left rects
                    int CrackOnAxisCompMin_v3( int L, int R, const float& crack_value, int axis, Rect *left_rect, Rect *right_rect);
                    int CrackOnAxisCompMax_v3( int L, int R, const float& crack_value, int axis, Rect *left_rect, Rect *right_rect);
                    // just calls partition_v4
                    int CrackOnAxisCompMin_v4( int L, int R, const float& crack_value, int axis, Rect *left_rect, Rect *right_rect);
                    int CrackOnAxisCompMax_v4( int L, int R, const float& crack_value, int axis, Rect *left_rect, Rect *right_rect);
                    // just calls partition v45
                    int CrackOnAxisCompMin_v45( int L, int R, const float& crack_value, int axis, Rect *left_rect, Rect *right_rect);
                    int CrackOnAxisCompMax_v45( int L, int R, const float& crack_value, int axis, Rect *left_rect, Rect *right_rect);
                    // just calls partition v46
                    int CrackOnAxisCompMin_v46( int L, int R, const float& crack_value, int axis, Rect *left_rect, Rect *right_rect);
                    int CrackOnAxisCompMax_v46( int L, int R, const float& crack_value, int axis, Rect *left_rect, Rect *right_rect);


                    void swap(float &a, float &b);
                    //    void swap(Branch* a, Branch *b);
                    void swap_index(int i, int j);
                    int MyPartition(int L, int R, float pivot_value, int axis, bool min_or_max);

                    // for in line node cover
                    int MyPartition_v2(int L, int R, float pivot_value, int axis, bool min_or_max, float *cover_mins_left, float *cover_maxes_left, float *cover_mins_right, float *cover_maxes_right);
                    //different way of passing the variables
                    int MyPartition_v3(int L, int R, float pivot_value, int axis, bool min_or_max, Rect *left_rect, Rect *right_rect);
                    // imitating ckd code, less swaps
                    int MyPartition_v4(int L, int R, float pivot_value, int axis, bool min_or_max, Rect *left_rect, Rect *right_rect);
                    // same as v4, just assigns the x1_rect for less indirection
                    int MyPartition_v45(int L, int R, float pivot_value, int axis, bool min_or_max, Rect *left_rect, Rect *right_rect);
                    // for testing the partition theory, completely wrong!!!!!!!!!!!!!! does not update left and right rects, so completely wrong tree creation
                    int MyPartition_v46(int L, int R, float pivot_value, int axis, bool min_or_max, Rect *left_rect, Rect *right_rect);


                    // FATEMEH OUT


                    // print tree structure
                    void printTreeRec(Node* a_node, ofstream &myfile);

                    // print for debug
                    void printTreeDebugRec(Node* a_node, ofstream &myfile);




                    //    void TwoDCrack(Node* start, const float a_min[NUMDIMS], const float a_max[NUMDIMS]);

                    Node *m_root;                                    ///< Root of tree
                    float m_unitSphereVolume;                 ///< Unit sphere constant for required number of dimensions
                    int last_id;
//                    Branch **m_data_arr = new Branch*[DATA_COUNT];    // static memory allocation for data
                    // MARCH 2023
                    // the first leaf, the left most leaf
                    Node* m_oldest_leaf;
                    // last leaf, right most leaf
                    Node* m_youngest_leaf;
                    // MARCH 2023


                    BranchList m_pendings;
                    Rect root_cover;

                    // 2023
                    int m_pending_insertions_count;
                    // 2023

                };



RTREE_TEMPLATE
RTREE_QUAL::RTree() {
    ASSERT(MAXNODES > MINNODES);
    ASSERT(MINNODES > 0);

    // Precomputed volumes of the unit spheres for the first few dimensions
    const float UNIT_SPHERE_VOLUMES[] = {
            0.000000f, 2.000000f, 3.141593f, // Dimension  0,1,2
            4.188790f, 4.934802f, 5.263789f, // Dimension  3,4,5
            5.167713f, 4.724766f, 4.058712f, // Dimension  6,7,8
            3.298509f, 2.550164f, 1.884104f, // Dimension  9,10,11
            1.335263f, 0.910629f, 0.599265f, // Dimension  12,13,14
            0.381443f, 0.235331f, 0.140981f, // Dimension  15,16,17
            0.082146f, 0.046622f, 0.025807f, // Dimension  18,19,20
    };

    m_root = AllocNode();
    m_root->m_level = 0;
    m_root->m_regular = false;
    m_root->m_parent = NULL;
    last_id = 0;
    m_root->m_id = last_id;

    // MARCH 2023
    m_youngest_leaf = m_root;
    m_oldest_leaf = m_root;
    m_root->m_older_brother = NULL;
    m_root->m_younger_brother = NULL;
    // MARCH 2023

    last_id++;
    m_unitSphereVolume = (float) UNIT_SPHERE_VOLUMES[NUMDIMS];
}

RTREE_TEMPLATE
RTREE_QUAL::RTree(const float data_set_min[][NUMDIMS], const float data_set_max[][NUMDIMS], int data_set_ids[],
                  int data_set_size) {
    ASSERT(MAXNODES > MINNODES);
    ASSERT(MINNODES > 0);

    // Precomputed volumes of the unit spheres for the first few dimensions
    const float UNIT_SPHERE_VOLUMES[] = {
            0.000000f, 2.000000f, 3.141593f, // Dimension  0,1,2
            4.188790f, 4.934802f, 5.263789f, // Dimension  3,4,5
            5.167713f, 4.724766f, 4.058712f, // Dimension  6,7,8
            3.298509f, 2.550164f, 1.884104f, // Dimension  9,10,11
            1.335263f, 0.910629f, 0.599265f, // Dimension  12,13,14
            0.381443f, 0.235331f, 0.140981f, // Dimension  15,16,17
            0.082146f, 0.046622f, 0.025807f, // Dimension  18,19,20
    };

    m_root = AllocNode();
    m_root->m_level = 0;
    m_root->m_regular = false;
    m_root->m_parent = NULL;
    last_id = 0;
    m_root->m_id = last_id;
    last_id ++;
    m_unitSphereVolume = (float) UNIT_SPHERE_VOLUMES[NUMDIMS];

    for (int i = 0; i < data_set_size; i++) {
        Branch branch;
        branch.m_data = data_set_ids[i];
        branch.m_child = NULL;

        for (int axis = 0; axis < NUMDIMS; ++axis) {
            branch.m_rect.m_min[axis] = data_set_min[i][axis];
            branch.m_rect.m_max[axis] = data_set_min[i][axis];
        }
        // the root is irregular and no new_nodes are needed
        AddBranch(&branch, m_root, NULL);
    }

    // MARCH 2023
    m_youngest_leaf = m_root;
    m_oldest_leaf = m_root;
    m_root->m_older_brother = NULL;
    m_root->m_younger_brother = NULL;
    // MARCH 2023
}

RTREE_TEMPLATE
RTREE_QUAL::RTree(std::string data_file_name){
    // srand(time(NULL));
    std::ifstream data_file(data_file_name.c_str());
    if(data_file.is_open()){cout << data_file_name << endl;}
    // int data_count = std::count(std::istreambuf_iterator<char>(data_file),
    //                             std::istreambuf_iterator<char>(), '\n');
    // cout << "FOR FATEMEH " << data_count << endl;
    // if(data_count < DATA_COUNT){
    //     cout << "NOT ENOUGH DATA!!\n";
    //     exit(3);
    // }
    cout << "I am fatemeh\n";

    data_file.clear();
    data_file.seekg(0, ios::beg);

    // float min_x, min_y, max_x, max_y;
    float min[NUMDIMS]; float max[NUMDIMS];

    for(int i = INITIAL_HOLES; i < DATA_COUNT+INITIAL_HOLES; i++){
//        Branch branch;
//        m_data_arr[i] = new Branch;

        for(int j = 0; j < NUMDIMS; j++){
            data_file >> min[j];
        }
        for(int j = 0; j < NUMDIMS; j++){
            data_file >> max[j];
        }

        // data_file >> min_x >> min_y >> max_x >> max_y;
        // min[0] = min_x; min[1] = min_y; max[0] = max_x; max[1] = max_y;

//        m_data_arr[i]->m_data = i;
        m_data_arr_ids[i] = i - INITIAL_HOLES;
//        m_data_arr[i]->m_child = NULL;
//        m_data_arr_mins[i] = new float[NUMDIMS];
//        m_data_arr_maxes[i] = new float[NUMDIMS];
//        m_data_arr_mins[i] = new float[NUMDIMS];
//        m_data_arr_maxes[i] = new float[NUMDIMS];
        for (int axis = 0; axis < NUMDIMS; ++axis) {
//            m_data_arr[i]->m_rect.m_min[axis] = min[axis];
            m_data_arr_mins[i][axis] = min[axis];
//            m_data_arr[i]->m_rect.m_max[axis] = max[axis];
            m_data_arr_maxes[i][axis] = max[axis];
        }

//        m_data_arr_datainfo[i] = i;

//        m_data_arr[i] = &branch;
        //        cout << "in for testing " << m_data_arr[i].m_rect.m_min[1] << " " << min_y << endl;
    }

    data_file.close();

    m_root = AllocNode();
    m_root->m_level = 0;
    m_root->m_regular = false;
    m_root->m_parent = NULL;
    last_id = 0;
    m_root->m_id = last_id;

    last_id++;

    m_root->m_transferred = false;
    m_root->m_L = 0;
    m_root->m_R = DATA_COUNT + INITIAL_HOLES;
    m_root->m_holes = INITIAL_HOLES;
    m_root->m_count = m_root->m_R - m_root->m_L - m_root->m_holes;
    ASSERT(m_root->m_count == DATA_COUNT);

    root_cover = NodeCover(m_root);

    // MARCH 2023
    m_youngest_leaf = m_root;
    m_oldest_leaf = m_root;
    m_root->m_older_brother = NULL;
    m_root->m_younger_brother = NULL;
    m_pending_insertions_count = 0;
    // MARCH 2023

    #ifdef stats
        count_irregular_leaves++;
        // auto t = CountLeaves();
        // cout << "in RTRee added one leaf " << count_leaves << " " << t.first + t.second << endl;
    #endif

    //    cout << "testing " << m_data_arr[0].m_rect.m_min[1] << endl;

}




RTREE_TEMPLATE
RTREE_QUAL::~RTree() {
    Reset(); // Free, or reset node memory
}


RTREE_TEMPLATE
void RTREE_QUAL::Insert_old(const float a_min[NUMDIMS], const float a_max[NUMDIMS], const DATATYPE &a_dataId) {


    Branch branch;
    branch.m_data = a_dataId;
    branch.m_child = NULL;

    for (int axis = 0; axis < NUMDIMS; ++axis) {
        branch.m_rect.m_min[axis] = a_min[axis];
        branch.m_rect.m_max[axis] = a_max[axis];
    }

    InsertRect(branch, &m_root, 0);
}

// 2023

// just add it to the pendings
RTREE_TEMPLATE
void RTREE_QUAL::Insert(const float a_min[NUMDIMS], const float a_max[NUMDIMS], const DATATYPE &a_dataId) {

    // cout << "inserting " << a_dataId << " into loc " << m_pending_insertions_count << endl;

    for (int axis = 0; axis < NUMDIMS; ++axis) {
        m_pending_insertions_mins[m_pending_insertions_count][axis] = a_min[axis];
        m_pending_insertions_maxes[m_pending_insertions_count][axis] = a_max[axis];
    }

    m_pending_insertions_ids[m_pending_insertions_count] = a_dataId;

    m_pending_insertions_count ++;

}
// 2023


RTREE_TEMPLATE
void RTREE_QUAL::Remove(const float a_min[NUMDIMS], const float a_max[NUMDIMS], const DATATYPE &a_dataId) {


    Rect rect;

    for (int axis = 0; axis < NUMDIMS; ++axis) {
        rect.m_min[axis] = a_min[axis];
        rect.m_max[axis] = a_max[axis];
    }

    RemoveRect(&rect, a_dataId, &m_root);
}


RTREE_TEMPLATE
int RTREE_QUAL::Search(const float a_min[NUMDIMS], const float a_max[NUMDIMS]) {


    Rect rect;

    for (int axis = 0; axis < NUMDIMS; ++axis) {
        rect.m_min[axis] = a_min[axis];
        rect.m_max[axis] = a_max[axis];
    }

    // NOTE: May want to return search result another way, perhaps returning the number of found elements here.

    int foundCount = 0;
    CrackVars cv;
    Search(m_root, &rect, foundCount, &cv);

    // TODO
    // now that cv is full, we have to take action
    // if count == 1 and no_regulars: 2dCrack
    // if no regular and more than one irr: put it all in one query_rect
    // if some regular, some irr: do insert for each point again:)
    // TODO somehow, our insert has to only insert into regulars now

    return foundCount;
}


RTREE_TEMPLATE
int RTREE_QUAL::CountDataPoints() {
    int foundCount = 0;
    CountDataPoints(m_root, foundCount);
    return foundCount;
}


RTREE_TEMPLATE
int RTREE_QUAL::CountNodes() {
    int foundCount = 0;
    CountNodes(m_root, foundCount);
    return foundCount;
}

RTREE_TEMPLATE
int RTREE_QUAL::TreeHeight() {
    return (m_root->m_level + 1);
}

RTREE_TEMPLATE
int RTREE_QUAL::SumDataPointIDs() {
    int foundCount = 0;
    SumDataPointIDs(m_root, foundCount);
    return foundCount;
}



// FATEMEH TEST
RTREE_TEMPLATE
int RTREE_QUAL::QueryAdaptive(const float *a_min, const float *a_max) {

    //    printf("in query_adapt\n");
    Rect rect;

    for (int axis = 0; axis < NUMDIMS; ++axis) {
        rect.m_min[axis] = a_min[axis];
        rect.m_max[axis] = a_max[axis];
    }

    // NOTE: May want to return search result another way, perhaps returning the number of found elements here.

    int foundCount = 0;
    CrackVars cv;
    //    printf("before calling search in q_adapt\n");
    //    auto start_time = std::chrono::steady_clock::now();
    //    printTreeDebug("../tree_right_before_search.txt");
    Search_old(m_root, &rect, foundCount, &cv);

    //    std::chrono::duration<double> search_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);


    // TODO
    // now that cv is full, we have to take action
    // if count == 1 and no_regulars: 2dCrack
    // if no regular and more than one irr: put it all in one query_rect
    // if some regular, some irr: do insert for each point again:)
    // TODO somehow, our insert has to only insert into regulars now

    // FATEMEH NEW
    int gathered_nodes_size = cv.gathered_nodes.size();
    if (gathered_nodes_size == 0 && !cv.encountered_regular) {
        //        printf("the query_rect overlaps nothing\n");
        return foundCount;
    } else if (gathered_nodes_size == 0 && cv.encountered_regular) {
        //        printf("only regulars, result is in the nowhere\n"); //TODO DO WE CARE?
        return foundCount;
    }
    // the case in between does not happen beacuase we ony add to gathered_node if something was found
    else if (gathered_nodes_size == 1 && !cv.encountered_regular) {
        // only one irreg
        //        printf("IN 2DCRACK\n");

        BranchList reinsert_list;

        //        bool is_regular_in_neighborhood = TwoDCrack_v4(cv.gathered_nodes[0].first, a_min, a_max, &reinsert_list);
        //        start_time = std::chrono::steady_clock::now();

        //        bool is_regular_in_neighborhood = TwoDCrack_v6(cv.gathered_nodes[0].first, (cv.gathered_nodes[0].second).second, a_min, a_max, &reinsert_list);
        bool is_regular_in_neighborhood = TwoDCrack_v4(cv.gathered_nodes[0].first, a_min, a_max, &reinsert_list);


        //        std::chrono::duration<double> two_dc_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);

        //        start_time = std::chrono::steady_clock::now();


        for(int i= 0; i < reinsert_list.branch_list.size(); i++){
            if(is_regular_in_neighborhood){
                InsertRectReg(reinsert_list.branch_list[i], &m_root, 0);
            }
            else{
                InsertRect(reinsert_list.branch_list[i], &m_root, 0);
            }
        }

        //        std::chrono::duration<double> reinsertion_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);

        //        cout << "TIMES search: " << search_time.count() << " 2dc: " << two_dc_time.count() << " reinsertion: " << reinsertion_time.count() << endl;


        return foundCount;
    }
    // the case in between will not happen for the same reason as before
    else if ((gathered_nodes_size >= 1 && cv.encountered_regular) || gathered_nodes_size > 1) {
        // just need to do 2dc on every node in gathered_nodes
        // no need to remove anything

        //        printf("REINSERTING INTO REGULAR BRANCHES\n");

        bool is_regular_in_neighborhood = false; bool tmp = false;
        BranchList reinsert_list;

        for(int i=0; i < cv.gathered_nodes.size(); i++){
            tmp = TwoDCrack_v4(cv.gathered_nodes[i].first, a_min, a_max, &reinsert_list);
            //            tmp = TwoDCrack_v6(cv.gathered_nodes[i].first, (cv.gathered_nodes[i].second).second, a_min, a_max, &reinsert_list);
            is_regular_in_neighborhood = tmp || is_regular_in_neighborhood;
        }

        // reinsert orphans

        for(int i= 0; i < reinsert_list.branch_list.size(); i++){
            if(is_regular_in_neighborhood){
                InsertRectReg(reinsert_list.branch_list[i], &m_root, 0);
            }
            else{
                InsertRect(reinsert_list.branch_list[i], &m_root, 0);
            }
        }
        return foundCount;
    }

    // FATEMEH NEW OUT

    return foundCount;
}

// FATEMEH TEST OUT



// FATEMEH TEST
RTREE_TEMPLATE
int RTREE_QUAL::QueryAdaptive_v5(const float *a_min, const float *a_max) {

    //printf("in query_adapt\n");
    Rect rect;

    for (int axis = 0; axis < NUMDIMS; ++axis) {
        rect.m_min[axis] = a_min[axis];
        rect.m_max[axis] = a_max[axis];
    }

    // NOTE: May want to return search result another way, perhaps returning the number of found elements here.

    int foundCount = 0;
    CrackVars cv;

    Search_9(m_root, &rect, foundCount, &cv);
    //    Search_12(m_root, &rect, foundCount, &cv);


    //    std::chrono::duration<double> search_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);


    // TODO
    // now that cv is full, we have to take action
    // if count == 1 and no_regulars: 2dCrack
    // if no regular and more than one irr: put it all in one query_rect
    // if some regular, some irr: do insert for each point again:)
    // TODO somehow, our insert has to only insert into regulars now

    // FATEMEH NEW
    int gathered_nodes_size = cv.gathered_node_cracks.size();
    if (gathered_nodes_size == 0 && !cv.encountered_regular) {
        //        printf("the query_rect overlaps nothing\n");
        return foundCount;
    } else if (gathered_nodes_size == 0 && cv.encountered_regular) {
        //        printf("only regulars, result is in the nowhere\n"); //TODO DO WE CARE?
        return foundCount;
    }
    // the case in between does not happen beacuase we ony add to gathered_node if something was found
    else if (gathered_nodes_size == 1 && !cv.encountered_regular) {
        // only one irreg
        //        printf("IN 2DCRACK\n");

        BranchList reinsert_list;

        //        bool is_regular_in_neighborhood = TwoDCrack_v4(cv.gathered_nodes[0].first, a_min, a_max, &reinsert_list);
        //        start_time = std::chrono::steady_clock::now();

        //        bool is_regular_in_neighborhood = TwoDCrack_v6(cv.gathered_nodes[0].first, (cv.gathered_nodes[0].second).second, a_min, a_max, &reinsert_list);
        bool is_regular_in_neighborhood = TwoDCrack_v8(cv.gathered_node_cracks[0].this_node, cv.gathered_node_cracks[0].cracks, cv.gathered_node_cracks[0].crack_covers_mins, cv.gathered_node_cracks[0].crack_covers_maxes, a_min, a_max);


        //        std::chrono::duration<double> two_dc_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);

        //        start_time = std::chrono::steady_clock::now();

        //        for(int i= 0; i < reinsert_list.branch_list.size(); i++){
        //            if(is_regular_in_neighborhood){
        //                InsertRectReg(reinsert_list.branch_list[i], &m_root, 0);
        //            }
        //            else{
        //                InsertRect(reinsert_list.branch_list[i], &m_root, 0);
        //            }
        //        }

        //        std::chrono::duration<double> reinsertion_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);

        //        cout << "TIMES search: " << search_time.count() << " 2dc: " << two_dc_time.count() << " reinsertion: " << reinsertion_time.count() << endl;


        return foundCount;
    }
    // the case in between will not happen for the same reason as before
    else if ((gathered_nodes_size >= 1 && cv.encountered_regular) || gathered_nodes_size > 1) {
        // just need to do 2dc on every node in gathered_nodes
        // no need to remove anything

        //        cout << "MULTIPLE IRREGS " << cv.gathered_node_cracks.size() << endl;

        bool is_regular_in_neighborhood = false; bool tmp = false;
        BranchList reinsert_list;

        for(int i=0; i < cv.gathered_node_cracks.size(); i++){
            //            cout << "one irregular at a time "<< endl;
            //            tmp = TwoDCrack_v4(cv.gathered_nodes[i].first, a_min, a_max, &reinsert_list);
            //            tmp = TwoDCrack_v6(cv.gathered_nodes[i].first, (cv.gathered_nodes[i].second).second, a_min, a_max, &reinsert_list);
            tmp = TwoDCrack_v8(cv.gathered_node_cracks[i].this_node, cv.gathered_node_cracks[i].cracks,cv.gathered_node_cracks[i].crack_covers_mins, cv.gathered_node_cracks[i].crack_covers_maxes, a_min, a_max);
            is_regular_in_neighborhood = tmp || is_regular_in_neighborhood;
        }

        // reinsert orphans

        //        for(int i= 0; i < reinsert_list.branch_list.size(); i++){
        //            if(is_regular_in_neighborhood){
        //                InsertRectReg(reinsert_list.branch_list[i], &m_root, 0);
        //            }
        //            else{
        //                InsertRect(reinsert_list.branch_list[i], &m_root, 0);
        //            }
        //        }
        return foundCount;
    }

    // FATEMEH NEW OUT

    return foundCount;
}




// FATEMEH TEST
RTREE_TEMPLATE
int RTREE_QUAL::QueryAdaptive_v6(const float *a_min, const float *a_max) {

    Rect rect;

    for (int axis = 0; axis < NUMDIMS; ++axis) {
        rect.m_min[axis] = a_min[axis];
        rect.m_max[axis] = a_max[axis];
    }

    int foundCount = 0;
    CrackVars cv;

    Search_10(m_root, &rect, foundCount, &cv);



    // FATEMEH NEW
    int gathered_nodes_size = cv.gathered_node_cracks.size();
    if (gathered_nodes_size == 0 && !cv.encountered_regular) {
        return foundCount;
    } else if (gathered_nodes_size == 0 && cv.encountered_regular) {
        return foundCount;
    }
    // the case in between does not happen beacuase we ony add to gathered_node if something was found
    else if (gathered_nodes_size == 1 && !cv.encountered_regular) {
        // only one irreg

        bool is_regular_in_neighborhood = TwoDCrack_v9(cv.gathered_node_cracks[0].this_node, cv.gathered_node_cracks[0].cracks, cv.gathered_node_cracks[0].crack_covers_mins, cv.gathered_node_cracks[0].crack_covers_maxes, a_min, a_max);

        return foundCount;
    }
    // the case in between will not happen for the same reason as before
    else if ((gathered_nodes_size >= 1 && cv.encountered_regular) || gathered_nodes_size > 1) {
        // just need to do 2dc on every node in gathered_nodes
        // no need to remove anything


        bool is_regular_in_neighborhood = false; bool tmp = false;
        BranchList reinsert_list;

        for(int i=0; i < cv.gathered_node_cracks.size(); i++){
            tmp = TwoDCrack_v9(cv.gathered_node_cracks[i].this_node, cv.gathered_node_cracks[i].cracks,cv.gathered_node_cracks[i].crack_covers_mins, cv.gathered_node_cracks[i].crack_covers_maxes, a_min, a_max);
            is_regular_in_neighborhood = tmp || is_regular_in_neighborhood;
        }

        return foundCount;
    }

    // FATEMEH NEW OUT

    return foundCount;
}


// FATEMEH TEST
RTREE_TEMPLATE
int RTREE_QUAL::QueryAdaptive_v7(const float *a_min, const float *a_max) {

    Rect rect;

    for (int axis = 0; axis < NUMDIMS; ++axis) {
        rect.m_min[axis] = a_min[axis];
        rect.m_max[axis] = a_max[axis];
    }

    int foundCount = 0;
    CrackVars cv;

    Search_11(m_root, &rect, foundCount, &cv);



    // FATEMEH NEW
    int gathered_nodes_size = cv.gathered_node_cracks.size();
    if (gathered_nodes_size == 0 && !cv.encountered_regular) {
        return foundCount;
    } else if (gathered_nodes_size == 0 && cv.encountered_regular) {
        return foundCount;
    }
    // the case in between does not happen beacuase we ony add to gathered_node if something was found
    else if (gathered_nodes_size == 1 && !cv.encountered_regular) {
        // only one irreg
        if(cv.gathered_node_cracks[0].three_or_five) {
            // three cracks
            TwoDCrack_v9(cv.gathered_node_cracks[0].this_node, cv.gathered_node_cracks[0].cracks,
                         cv.gathered_node_cracks[0].crack_covers_mins, cv.gathered_node_cracks[0].crack_covers_maxes,
                         a_min, a_max);
        } else {
            // five cracks
            TwoDCrack_v8(cv.gathered_node_cracks[0].this_node,
                         cv.gathered_node_cracks[0].cracks,
                         cv.gathered_node_cracks[0].crack_covers_mins,
                         cv.gathered_node_cracks[0].crack_covers_maxes, a_min, a_max);
        }
        return foundCount;
    }
    // the case in between will not happen for the same reason as before
    else if ((gathered_nodes_size >= 1 && cv.encountered_regular) || gathered_nodes_size > 1) {
        // just need to do 2dc on every node in gathered_nodes
        // no need to remove anything


        bool is_regular_in_neighborhood = false; bool tmp = false;
        BranchList reinsert_list;

        for(int i=0; i < cv.gathered_node_cracks.size(); i++){
            if(cv.gathered_node_cracks[i].three_or_five) {
                // three cracks
                tmp = TwoDCrack_v9(cv.gathered_node_cracks[i].this_node, cv.gathered_node_cracks[i].cracks,
                                   cv.gathered_node_cracks[i].crack_covers_mins,
                                   cv.gathered_node_cracks[i].crack_covers_maxes, a_min, a_max);

            } else{
                // five cracks
                tmp = TwoDCrack_v8(cv.gathered_node_cracks[i].this_node, cv.gathered_node_cracks[i].cracks,
                                   cv.gathered_node_cracks[i].crack_covers_mins,
                                   cv.gathered_node_cracks[i].crack_covers_maxes, a_min, a_max);

            }
            is_regular_in_neighborhood = tmp || is_regular_in_neighborhood;
        }

        return foundCount;
    }

    // FATEMEH NEW OUT

    return foundCount;
}




// FATEMEH TEST
RTREE_TEMPLATE
int RTREE_QUAL::QueryAdaptive_v8(const float *a_min, const float *a_max) {

    //printf("in query_adapt\n");
    Rect rect;

    for (int axis = 0; axis < NUMDIMS; ++axis) {
        rect.m_min[axis] = a_min[axis];
        rect.m_max[axis] = a_max[axis];
    }

    // NOTE: May want to return search result another way, perhaps returning the number of found elements here.

    int foundCount = 0;
    CrackVars cv;

    Search_12(m_root, &rect, foundCount, &cv);


    //    std::chrono::duration<double> search_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);


    // TODO
    // now that cv is full, we have to take action
    // if count == 1 and no_regulars: 2dCrack
    // if no regular and more than one irr: put it all in one query_rect
    // if some regular, some irr: do insert for each point again:)
    // TODO somehow, our insert has to only insert into regulars now

    // FATEMEH NEW
    int gathered_nodes_size = cv.gathered_node_cracks.size();

    bool is_regular_in_neighborhood = false; bool tmp = false;

    for(int i=0; i < cv.gathered_node_cracks.size(); i++){
        tmp = TwoDCrack_v8(cv.gathered_node_cracks[i].this_node, cv.gathered_node_cracks[i].cracks,cv.gathered_node_cracks[i].crack_covers_mins, cv.gathered_node_cracks[i].crack_covers_maxes, a_min, a_max);
        is_regular_in_neighborhood = tmp || is_regular_in_neighborhood;
    }


    return foundCount;


    // FATEMEH NEW OUT

}



// FATEMEH TEST
RTREE_TEMPLATE
int RTREE_QUAL::QueryAdaptive_v9(const float *a_min, const float *a_max) {

    //printf("in query_adapt\n");
    Rect rect;

    for (int axis = 0; axis < NUMDIMS; ++axis) {
        rect.m_min[axis] = a_min[axis];
        rect.m_max[axis] = a_max[axis];
    }

    // NOTE: May want to return search result another way, perhaps returning the number of found elements here.

    int foundCount = 0;
    CrackVars_9 cv;

    Search_14(m_root, &rect, foundCount, &cv);


    //    std::chrono::duration<double> search_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);


    // TODO
    // now that cv is full, we have to take action
    // if count == 1 and no_regulars: 2dCrack
    // if no regular and more than one irr: put it all in one query_rect
    // if some regular, some irr: do insert for each point again:)
    // TODO somehow, our insert has to only insert into regulars now

    // FATEMEH NEW
    int gathered_nodes_size = cv.gathered_node_cracks.size();

    //    bool is_regular_in_neighborhood = false;
    bool tmp = false;

    for(int i=0; i < cv.gathered_node_cracks.size(); i++){
        tmp = TwoDCrack_v10(cv.gathered_node_cracks[i].this_node, cv.gathered_node_cracks[i].cracks,cv.gathered_node_cracks[i].crack_covers_mins, cv.gathered_node_cracks[i].crack_covers_maxes, a_min, a_max);
        //        is_regular_in_neighborhood = tmp || is_regular_in_neighborhood;
    }


    return foundCount;


    // FATEMEH NEW OUT

}


RTREE_TEMPLATE
int RTREE_QUAL::QueryAdaptive_v10(const float *a_min, const float *a_max) {

    //printf("in query_adapt\n");
    Rect rect;

    for (int axis = 0; axis < NUMDIMS; ++axis) {
        rect.m_min[axis] = a_min[axis];
        rect.m_max[axis] = a_max[axis];
    }

    // NOTE: May want to return search result another way, perhaps returning the number of found elements here.

    int foundCount = 0;
    CrackVars cv;

    Search_15(m_root, &rect, foundCount, &cv);


    //    std::chrono::duration<double> search_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);


    // TODO
    // now that cv is full, we have to take action
    // if count == 1 and no_regulars: 2dCrack
    // if no regular and more than one irr: put it all in one query_rect
    // if some regular, some irr: do insert for each point again:)
    // TODO somehow, our insert has to only insert into regulars now

    // FATEMEH NEW
    int gathered_nodes_size = cv.gathered_node_cracks.size();

    //    bool is_regular_in_neighborhood = false;
    bool tmp = false;

    for(int i=0; i < cv.gathered_node_cracks.size(); i++){
        tmp = TwoDCrack_v8(cv.gathered_node_cracks[i].this_node, cv.gathered_node_cracks[i].cracks,cv.gathered_node_cracks[i].crack_covers_mins, cv.gathered_node_cracks[i].crack_covers_maxes, a_min, a_max);
        //        is_regular_in_neighborhood = tmp || is_regular_in_neighborhood;
    }
    return foundCount;

    // FATEMEH NEW OUT

}



// FATEMEH TEST
RTREE_TEMPLATE
int RTREE_QUAL::QueryAdaptive_v11(const float *a_min, const float *a_max) {

    //printf("in query_adapt\n");
    Rect rect;

    for (int axis = 0; axis < NUMDIMS; ++axis) {
        rect.m_min[axis] = a_min[axis];
        rect.m_max[axis] = a_max[axis];
    }

    // NOTE: May want to return search result another way, perhaps returning the number of found elements here.

    int foundCount = 0;
    vector<DATATYPE> results;
    CrackVars cv;

    Search_16(m_root, &rect, results, &cv);


    //    std::chrono::duration<double> search_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);


    // TODO
    // now that cv is full, we have to take action
    // if count == 1 and no_regulars: 2dCrack
    // if no regular and more than one irr: put it all in one query_rect
    // if some regular, some irr: do insert for each point again:)
    // TODO somehow, our insert has to only insert into regulars now

    // FATEMEH NEW
    int gathered_nodes_size = cv.gathered_node_cracks.size();

    bool is_regular_in_neighborhood = false; bool tmp = false;

    for(int i=0; i < cv.gathered_node_cracks.size(); i++){
        tmp = TwoDCrack_v8(cv.gathered_node_cracks[i].this_node, cv.gathered_node_cracks[i].cracks,cv.gathered_node_cracks[i].crack_covers_mins, cv.gathered_node_cracks[i].crack_covers_maxes, a_min, a_max);
        is_regular_in_neighborhood = tmp || is_regular_in_neighborhood;
    }


    return results.size();


    // FATEMEH NEW OUT

}

RTREE_TEMPLATE
int RTREE_QUAL::QueryAdaptive_v12(const float *a_min, const float *a_max) {

    Rect rect;

    for (int axis = 0; axis < NUMDIMS; ++axis) {
        rect.m_min[axis] = a_min[axis];
        rect.m_max[axis] = a_max[axis];
    }


    int foundCount = 0;
    LeavesToAdd all_lta;

    Search_19(m_root, &rect, foundCount, &all_lta);
    Add_ltas(&all_lta);

    return foundCount;

}


RTREE_TEMPLATE
int RTREE_QUAL::QueryAdaptive_v13(const float *a_min, const float *a_max) {

    Rect rect;

    for (int axis = 0; axis < NUMDIMS; ++axis) {
        rect.m_min[axis] = a_min[axis];
        rect.m_max[axis] = a_max[axis];
    }


    int foundCount = 0;
    LeavesToAdd_v2 all_lta;

    cout << "in v13, before search..." << endl;

    Search_2026(m_root, &rect, foundCount, &all_lta);
    cout << "in v13, before add_ltas... " << endl;
    Add_ltas_v2(&all_lta);
    cout << "in v13, after add_ltas." << endl;
    return foundCount;

}


RTREE_TEMPLATE
int RTREE_QUAL::QueryAdaptive_v14(const float *a_min, const float *a_max) {

    #ifdef timebreak
        scan_pending_time = 0.0;
        delete_pending_time = 0.0;
        search_time = 0.0;
        add_ltas_time = 0.0;
        ripple_time = 0.0;

        // in search
        tba_assignment_time = 0.0;
        cracking_time = 0.0;
        traversedown_time = 0.0;

        // in ripple
        shuffle_time = 0.0;
    #endif 

    // cout << "before before item 16763" << endl;
    // if(m_pending_insertions_ids[15] == 16763){
    //     cout << "bounds after replacement shit " << m_pending_insertions_mins[15][0] << " " << m_pending_insertions_maxes[15][0] << " , " << m_pending_insertions_mins[15][1] << " " << m_pending_insertions_maxes[15][1] << endl;
    // }

    Rect rect;

    for (int axis = 0; axis < NUMDIMS; ++axis) {
        rect.m_min[axis] = a_min[axis];
        rect.m_max[axis] = a_max[axis];
    }

    // TODO: look through the pending insertions
    // if it overlaps query, add it to a vector tba and give to search
    // TODO: make new search function that takes the tba vector and distributes it around

    // old
    // vector<Branch> tba;

    // vector<int> tba_indexes;

    // for(int i = 0; i < m_pending_insertions_count; i++){
    //     if(Overlap(&rect, m_pending_insertions_mins[i], m_pending_insertions_maxes[i])){
    //         tba_indexes.push_back(i);

    //         Branch branch;
    //         branch.m_data = m_pending_insertions_ids[i];
    //         branch.m_child = NULL;

    //         for (int axis = 0; axis < NUMDIMS; ++axis) {
    //             branch.m_rect.m_min[axis] = m_pending_insertions_mins[i][axis];
    //             branch.m_rect.m_max[axis] = m_pending_insertions_maxes[i][axis];
    //         }

    //         tba.push_back(branch);

    //         // TODO: I need to adjust the mbr of the root with the new tbis
    //         for (int axis = 0; axis < NUMDIMS; ++axis) {
    //             if(root_cover.m_min[axis] > m_pending_insertions_mins[i][axis])
    //                 {root_cover.m_min[axis] = m_pending_insertions_mins[i][axis];}
    //             if(root_cover.m_max[axis] < m_pending_insertions_maxes[i][axis])
    //                 {root_cover.m_max[axis] = m_pending_insertions_maxes[i][axis];}
    //         }

    //     }
    // } 
    // old

    // new
    // float tba_mins[m_pending_insertions_count][NUMDIMS];
    // float tba_maxes[m_pending_insertions_count][NUMDIMS];
    // int tba_ids[m_pending_insertions_count];

    // float* tba_mins;
    // tba_mins = (float*)malloc(m_pending_insertions_count * NUMDIMS* sizeof(float));
    // float* tba_maxes;
    // tba_maxes = (float*)malloc(m_pending_insertions_count *NUMDIMS* sizeof(float));

    // int *tba_ids;
    // tba_ids = (int*)malloc(m_pending_insertions_count * sizeof(int));


    // int tba_indexes[m_pending_insertions_count];
    int * tba_indexes = new int[m_pending_insertions_count];
    int tba_count = 0;
    #ifdef timebreak
        auto start_time = clock();
    #endif
    for(int i = 0; i < m_pending_insertions_count; i++){
        if(Overlap(&rect, m_pending_insertions_mins[i], m_pending_insertions_maxes[i])){

            tba_indexes[tba_count] = i;
            tba_ids[tba_count] = m_pending_insertions_ids[i];
            // if(tba_ids[tba_count] == 16763){
            //     cout << "in pending 16763 " << m_pending_insertions_mins[i][0] << " " << m_pending_insertions_maxes[i][0] <<  ", " << m_pending_insertions_mins[i][1] << " " << m_pending_insertions_maxes[i][1] << endl;
            // }

            for (int axis = 0; axis < NUMDIMS; ++axis) {
                tba_mins[tba_count][axis] = m_pending_insertions_mins[i][axis];
                tba_maxes[tba_count][axis] = m_pending_insertions_maxes[i][axis];
            }

            tba_count++;

            // TODO: I need to adjust the mbr of the root with the new tbis
            for (int axis = 0; axis < NUMDIMS; ++axis) {
                if(root_cover.m_min[axis] > m_pending_insertions_mins[i][axis])
                    {root_cover.m_min[axis] = m_pending_insertions_mins[i][axis];}
                if(root_cover.m_max[axis] < m_pending_insertions_maxes[i][axis])
                    {root_cover.m_max[axis] = m_pending_insertions_maxes[i][axis];}
            }
        }
    }

    #ifdef timebreak
        scan_pending_time = clock() - start_time;
    #endif


    // new


    // TODOM DONE
    // now i have the index of the guys that will be inserted this time
    // I have to remove them
    // replace them from the back

    // cout << "pending before removing stuff &&&&&&&&&&&" << endl;
    // for(int i = 0; i < m_pending_insertions_count; i++){
    //     cout << m_pending_insertions_ids[i] << endl;
    // }
    // cout << "&&&&&&&&&&" << endl;

    // for(int i = 0; i < tba.size(); i++){
    //     // replace pendings[tba_index[i]] with pendings[count - tba.size + i]

    //     for (int axis = 0; axis < NUMDIMS; ++axis) {
    //         m_pending_insertions_mins[tba_indexes[i]][axis] = m_pending_insertions_mins[m_pending_insertions_count - tba.size() + i][axis];
    //         m_pending_insertions_maxes[tba_indexes[i]][axis] = m_pending_insertions_maxes[m_pending_insertions_count - tba.size() + i][axis];
    //     }

    //     m_pending_insertions_ids[tba_indexes[i]] = m_pending_insertions_ids[m_pending_insertions_count - tba.size() + i];
    // }

    // old
    // for(int i = tba.size() - 1; i > -1; i--){
    //     if(tba_indexes[i] != m_pending_insertions_count - 1){
    //         for (int axis = 0; axis < NUMDIMS; ++axis) {
    //             m_pending_insertions_mins[tba_indexes[i]][axis] = m_pending_insertions_mins[m_pending_insertions_count -1][axis];
    //             m_pending_insertions_maxes[tba_indexes[i]][axis] = m_pending_insertions_maxes[m_pending_insertions_count - 1][axis];
    //         }
    //         m_pending_insertions_ids[tba_indexes[i]] = m_pending_insertions_ids[m_pending_insertions_count - 1];

    //         m_pending_insertions_count--;
    //     }
    //     else{
    //         m_pending_insertions_count--;
    //     }
    // }
    // old
    // cout << "before item 16763" << endl;
    // if(m_pending_insertions_ids[15] == 16763){
    //     cout << "bounds after replacement shit " << m_pending_insertions_mins[15][0] << " " << m_pending_insertions_maxes[15][0] << " , " << m_pending_insertions_mins[15][1] << " " << m_pending_insertions_maxes[15][1] << endl;
    // }

    #ifdef timebreak
        start_time = clock();
    #endif

    for(int i = tba_count -1; i > -1; i--){
        if(tba_indexes[i] != m_pending_insertions_count - 1){
            for (int axis = 0; axis < NUMDIMS; ++axis) {
                m_pending_insertions_mins[tba_indexes[i]][axis] = m_pending_insertions_mins[m_pending_insertions_count -1][axis];
                m_pending_insertions_maxes[tba_indexes[i]][axis] = m_pending_insertions_maxes[m_pending_insertions_count - 1][axis];
            }
            m_pending_insertions_ids[tba_indexes[i]] = m_pending_insertions_ids[m_pending_insertions_count - 1];

            m_pending_insertions_count--;
        }
        else{
            m_pending_insertions_count--;
        }
    }

    #ifdef timebreak
        delete_pending_time = clock() - start_time;
    #endif
    

    // cout << "pending after removing stuff ^^^^^^^^^^^^^^^^^^^" << endl;
    // for(int i = 0; i < m_pending_insertions_count; i++){
    //     cout << m_pending_insertions_ids[i] << endl;
    // }
    // cout << "^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;


    
    // findAllLeavesWithDataInThem(45780);

    // TODO, DONE
    // I also need to -- the count of pendings

    // cout << "pending_insert_count before update " << m_pending_insertions_count << " tba.size() " << tba.size() << endl;

    // m_pending_insertions_count -= tba.size();

    // cout << "pending insert count after update " << m_pending_insertions_count << endl;
    // cout << "***************** " << tba_count << endl;
    // if(tba_count > 0){
    //     for(int ii = 0; ii < tba_count; ii++){
    //         // cout << tba_ids[ii] << " " << tba_mins[ii][0] << " " << tba_maxes[ii][0] << " , " << tba_mins[ii][1] << " " << tba_maxes[ii][1] << endl;
    //         cout << tba_ids[ii] << endl;
    //     }
    // }
    // cout << "*****************" << endl;

    // cout << "item 16763" << endl;
    // if(m_pending_insertions_ids[15] == 16763){
    //     cout << "bounds after replacement shit " << m_pending_insertions_mins[15][0] << " " << m_pending_insertions_maxes[15][0] << " , " << m_pending_insertions_mins[15][1] << " " << m_pending_insertions_maxes[15][1] << endl;
    // }

    


    int foundCount = 0;
    // foundCount += tba.size();
    foundCount += tba_count;
    LeavesToAdd_v2 all_lta;

    vector<overlapping_leaf_tbas> ol;

    // findNodeCount(15518, m_root);

    // cout << "before search, ol.size() " << ol.size() << endl;
    // cout << "search..." << endl;
    // Search_2027(m_root, &rect, foundCount, &all_lta, &tba, &ol);
    #ifdef timebreak
        start_time = clock();
    #endif
    Search_2027(m_root, &rect, foundCount, &all_lta, 0, tba_count, &ol);
    #ifdef timebreak
        search_time = clock() - start_time;
    #endif
    // cout << "after search, before add_ltas, ol.size() " << ol.size() << " ol[0].this_leaf.id " << ol[0].this_leaf->m_id << endl;
    // cout << "after search, before add_ltas, ol.size() " << ol.size() << endl;
    // if(ol.size() > 0 ){
    //     cout << "ol[0].this_leaf.id " << ol[0].this_leaf->m_id << " ol[0].this_tba.size() " << ol[0].this_tba.size() << endl;
    // }

    // findAllLeavesWithDataInThem(45780);
    // findNodeCount(15518, m_root);



    // cout << "add ltas..." << endl;
    #ifdef timebreak
        start_time = clock();
    #endif
    Add_ltas_v3(&all_lta, &ol);
    #ifdef timebreak
        add_ltas_time = clock() - start_time;
    #endif
    // cout << "after add_ltas, ol.size() " << ol.size()  << endl;
    // if(ol.size() > 0 ){
    //     cout << "ol[0].this_leaf.id " << ol[0].this_leaf->m_id << " ol[0].this_tba.size() " << ol[0].this_tba.size() << endl;
    // }
    // findAllLeavesWithDataInThem(45780);
    // findNodeCount(15518, m_root);



    // cout << "tree leaves: " << endl;
    // printLeafLinkedList();
    // printIDs();

    // ripple
    // cout << "before ripple, ol.size() " << ol.size() << endl;
    // if(ol.size() > 0 ){
    //     cout << "ol[0].this_leaf.id " << ol[0].this_leaf->m_id << " ol[0].this_tba.size() " << ol[0].this_tba.size() << endl;
    // }
    // for(int ii = 0; ii < ol.size(); ii++){
    //     cout << "leaf_id, tba size: " << ol[ii].this_leaf->m_id << " " <<  ol[ii].this_tba.size() << endl;
    // }
    // cout << "ripple ..." << endl;
    #ifdef timebreak
        start_time = clock();
    #endif
    ripple_v3(&ol);
    #ifdef timebreak
        ripple_time = clock() - start_time;
        cout <<  (double)search_time/CLOCKS_PER_SEC << " " << (double)tba_assignment_time/CLOCKS_PER_SEC << " " << (double)cracking_time/CLOCKS_PER_SEC << " " << (double)traversedown_time/CLOCKS_PER_SEC << " ";
        cout << (double)add_ltas_time/CLOCKS_PER_SEC << " " << (double)ripple_time/CLOCKS_PER_SEC << " " << (double)shuffle_time/CLOCKS_PER_SEC << " ";
        cout << (double)scan_pending_time/CLOCKS_PER_SEC << " " << (double)delete_pending_time/CLOCKS_PER_SEC << endl;
    #endif
    // cout << "After ripple " << endl;
    // findNodeCount(15518, m_root);


    // findAllLeavesWithDataInThem(45780);

    // cout << "Before free in v14\n";
    // for(int i =0; i < m_pending_insertions_count; i++){
    //     free(tba_mins[i]);
    //     free(tba_maxes[i]);
    // }
    // free(tba_mins);
    // free(tba_maxes);
    // free(tba_ids);
    // cout << "After free in v14\n";

    //  cout << "end of query item 16763" << endl;
    // if(m_pending_insertions_ids[15] == 16763){
    //     cout << "bounds after replacement shit " << m_pending_insertions_mins[15][0] << " " << m_pending_insertions_maxes[15][0] << " , " << m_pending_insertions_mins[15][1] << " " << m_pending_insertions_maxes[15][1] << endl;
    // }

    return foundCount;

}




RTREE_TEMPLATE
int RTREE_QUAL::Count() {
    int count = 0;
    CountRec(m_root, count);

    return count;
}


RTREE_TEMPLATE
void RTREE_QUAL::CountRec(Node *a_node, int &a_count) {
    if (a_node->IsInternalNode())  // not a leaf node
        {
        for (int index = 0; index < a_node->m_count; ++index) {
            CountRec(a_node->m_branch[index].m_child, a_count);
        }
        //        for (typename std::vector<Branch>::iterator it = a_node->m_branch.begin(); it != a_node->m_branch.end(); it++) {
        //            CountRec(it->m_child, a_count);
        //        }
        } else // A leaf node
        {
        a_count += a_node->m_count;
        }
}




RTREE_TEMPLATE
void RTREE_QUAL::RemoveAll() {
    // Delete all existing nodes
    Reset();

    m_root = AllocNode();
    m_root->m_level = 0;
}


RTREE_TEMPLATE
void RTREE_QUAL::Reset() {
#ifdef RTREE_DONT_USE_MEMPOOLS
    // Delete all existing nodes
    RemoveAllRec(m_root);
#else // RTREE_DONT_USE_MEMPOOLS
    // Just reset memory pools.  We are not using complex types
    // EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}


RTREE_TEMPLATE
void RTREE_QUAL::RemoveAllRec(Node *a_node) {

    //    printf("HEEEEEY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    //            ASSERT(a_node);
    //            ASSERT(a_node->m_level >= 0);
    // I think we have to make this in decreasing order and it shoudl be fine, but Im testing for now
    if (a_node->IsInternalNode()) // This is an internal node in the tree
        {
        //        printf("calling NC in RempveAllRec\n");
        // Rect for_me = NodeCover(a_node);
        for (int index = 0; index < a_node->m_count; ++index) {
            //            printf("removing child %d with level: %d from node: %d, %d, %d, %d\n", index, a_node->m_branch[index].m_child->m_level, for_me.m_min[0], for_me.m_min[1], for_me.m_max[0], for_me.m_max[1]);
            RemoveAllRec(a_node->m_branch[index].m_child);
        }
        }
    FreeNode(a_node);
}


RTREE_TEMPLATE
typename RTREE_QUAL::Node *RTREE_QUAL::AllocNode() {
    Node *newNode;
#ifdef RTREE_DONT_USE_MEMPOOLS
    newNode = new Node;
#else // RTREE_DONT_USE_MEMPOOLS
    // EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
InitNode(newNode);
return newNode;
}


RTREE_TEMPLATE
void RTREE_QUAL::FreeNode(Node *a_node) {
    //            ASSERT(a_node);

#ifdef RTREE_DONT_USE_MEMPOOLS
delete a_node;
#else // RTREE_DONT_USE_MEMPOOLS
// EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}


// Allocate space for a node in the list used in DeletRect to
// store Nodes that are too empty.
RTREE_TEMPLATE
typename RTREE_QUAL::ListNode *RTREE_QUAL::AllocListNode() {
#ifdef RTREE_DONT_USE_MEMPOOLS
    return new ListNode;
#else // RTREE_DONT_USE_MEMPOOLS
    // EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}


RTREE_TEMPLATE
void RTREE_QUAL::FreeListNode(ListNode *a_listNode) {
#ifdef RTREE_DONT_USE_MEMPOOLS
    delete a_listNode;
#else // RTREE_DONT_USE_MEMPOOLS
    // EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}


RTREE_TEMPLATE
void RTREE_QUAL::InitNode(Node *a_node) {
    a_node->m_count = 0;
    a_node->m_level = -1;
    a_node->m_regular = true;
    a_node->name = "init";
    a_node->m_id = last_id;
    a_node->m_transferred = false;
    last_id++;

    // 2023
    a_node->m_younger_brother = NULL;
    a_node->m_older_brother = NULL;
    a_node->m_holes = 0;
    // 2023
}


RTREE_TEMPLATE
void RTREE_QUAL::InitRect(Rect *a_rect) {
    for (int index = 0; index < NUMDIMS; ++index) {
        a_rect->m_min[index] = (float) 0;
        a_rect->m_max[index] = (float) 0;
    }
}


// Inserts a new data rectangle into the index structure.
// Recursively descends tree, propagates splits back up.
// Returns 0 if node was not split.  Old node updated.
// If node was split, returns 1 and sets the pointer pointed to by
// new_node to point to the new node.  Old node updated to become one of two.
// The level argument specifies the number of steps up from the leaf
// level to insert; e.g. a data rectangle goes in at level = 0.
RTREE_TEMPLATE
bool RTREE_QUAL::InsertRectRec(const Branch &a_branch, Node *a_node, Node **a_newNode, int a_level) {
    //            ASSERT(a_node && a_newNode);
    //            ASSERT(a_level >= 0 && a_level <= a_node->m_level);

    // recurse until we reach the correct level for the new record. data records
    // will always be called with a_level == 0 (leaf)
    if (a_node->m_level > a_level) {
        // Still above level for insertion, go down tree recursively
        Node *otherNode;

        // find the optimal branch for this record
        int index = PickBranch(&a_branch.m_rect, a_node);

        // recursively insert this record into the picked branch
        bool childWasSplit = InsertRectRec(a_branch, a_node->m_branch[index].m_child, &otherNode, a_level);

        if (!childWasSplit) {
            // Child was not split. Merge the bounding box of the new record with the
            // existing bounding box
            a_node->m_branch[index].m_rect = CombineRect(&a_branch.m_rect, &(a_node->m_branch[index].m_rect));
            return false;
        } else {
            // Child was split. The old branches are now re-partitioned to two nodes
            // so we have to re-calculate the bounding boxes of each node
            //            printf("calling NC in InsertRectRec\n");
            a_node->m_branch[index].m_rect = NodeCover(a_node->m_branch[index].m_child);
            Branch branch;
            branch.m_child = otherNode;
            //            printf("calling NC in InsertRectRec again\n");
            branch.m_rect = NodeCover(otherNode);

            // The old node is already a child of a_node. Now add the newly-created
            // node to a_node as well. a_node might be split because of that.
            return AddBranch(&branch, a_node, a_newNode);
        }
    }
    else if (a_node->m_level == a_level) {
        // We have reached level for insertion. Add rect, split if necessary
        //        cout << "levels were equal calling addbranch\n";
        return AddBranch(&a_branch, a_node, a_newNode);
    } else {
        // Should never occur
        //                ASSERT(0);
        return false;
    }
}


// Inserts a new data rectangle into the index structure.
// Recursively descends tree, propagates splits back up.
// Returns 0 if node was not split.  Old node updated.
// If node was split, returns 1 and sets the pointer pointed to by
// new_node to point to the new node.  Old node updated to become one of two.
// The level argument specifies the number of steps up from the leaf
// level to insert; e.g. a data rectangle goes in at level = 0.
RTREE_TEMPLATE
bool RTREE_QUAL::InsertRectRec(const Branch &a_branch, Node *a_node, Node **a_newNode, int a_level, Rect &a_coverRect, Rect &a_newCoverRect) {
    //            ASSERT(a_node && a_newNode);
    //            ASSERT(a_level >= 0 && a_level <= a_node->m_level);

    // recurse until we reach the correct level for the new record. data records
    // will always be called with a_level == 0 (leaf)
    if (a_node->m_level > a_level) {
        // Still above level for insertion, go down tree recursively
        Node *otherNode;

        // find the optimal branch for this record
        int index = PickBranch(&a_branch.m_rect, a_node);

        // recursively insert this record into the picked branch
        bool childWasSplit = InsertRectRec(a_branch, a_node->m_branch[index].m_child, &otherNode, a_level, a_coverRect, a_newCoverRect);

        if (!childWasSplit) {
            // Child was not split. Merge the bounding box of the new record with the
            // existing bounding box
            a_node->m_branch[index].m_rect = CombineRect(&a_branch.m_rect, &(a_node->m_branch[index].m_rect));
            return false;
        } else {
            // Child was split. The old branches are now re-partitioned to two nodes
            // so we have to re-calculate the bounding boxes of each node
            //            printf("calling NC in InsertRectRec\n");
            //            a_node->m_branch[index].m_rect = NodeCover(a_node->m_branch[index].m_child);
            a_node->m_branch[index].m_rect = Rect(a_coverRect.m_max, a_coverRect.m_max);
            Branch branch;
            branch.m_child = otherNode;
            //            printf("calling NC in InsertRectRec again\n");
            //            branch.m_rect = NodeCover(otherNode);
            branch.m_rect = Rect(a_newCoverRect.m_min, a_newCoverRect.m_max);

            // The old node is already a child of a_node. Now add the newly-created
            // node to a_node as well. a_node might be split because of that.
            return AddBranch(&branch, a_node, a_newNode, a_coverRect, a_newCoverRect);
        }
    }
    else if (a_node->m_level == a_level) {
        // We have reached level for insertion. Add rect, split if necessary
        //        cout << "levels were equal calling addbranch\n";
        return AddBranch(&a_branch, a_node, a_newNode, a_coverRect, a_newCoverRect);
    } else {
        // Should never occur
        //                ASSERT(0);
        return false;
    }
}



// Insert a data rectangle into an index structure.
// InsertRect provides for splitting the root;
// returns 1 if root was split, 0 if it was not.
// The level argument specifies the number of steps up from the leaf
// level to insert; e.g. a data rectangle goes in at level = 0.
// InsertRect2 does the recursion.
//
RTREE_TEMPLATE
bool RTREE_QUAL::InsertRect(const Branch &a_branch, Node **a_root, int a_level) {
    //            ASSERT(a_root);
    //            ASSERT(a_level >= 0 && a_level <= (*a_root)->m_level);


    Node *newNode;

    //    cout << "in insrertrect\n";

    Rect rect1, rect2;

    if (InsertRectRec(a_branch, *a_root, &newNode, a_level, rect1, rect2))  // Root split
        {
        // Grow tree taller and new root
        Node *newRoot = AllocNode();
        newRoot->m_level = (*a_root)->m_level + 1;
        newRoot->m_parent = NULL;
        // DEBUG
        newRoot->name = "new root";

        Branch branch;

        // add old root node as a child of the new root
        //        printf("calling NC in InsertRect\n");

        //        branch.m_rect = NodeCover(*a_root);
        branch.m_rect = Rect(rect1.m_min, rect1.m_max);
        branch.m_child = *a_root;
        AddBranch(&branch, newRoot, NULL);

        // add the split node as a child of the new root
        //        printf("calling NC in InsertRectRec again\n");

        //        branch.m_rect = NodeCover(newNode);
        branch.m_rect = Rect(rect2.m_min, rect2.m_max);
        branch.m_child = newNode;
        AddBranch(&branch, newRoot, NULL);

        // set the new root as the root node
        *a_root = newRoot;

        return true;
        }

    return false;
}


// returns <was_split, inserted_or_not>
RTREE_TEMPLATE
typename RTREE_QUAL::bools RTREE_QUAL::InsertRectOnlyRegRec(const Branch &a_branch, Node *a_node, Node **a_newNode, int a_level) {
    //            ASSERT(a_node && a_newNode);
    //            ASSERT(a_level >= 0 && a_level <= a_node->m_level);

    // recurse until we reach the correct level for the new record. data records
    // will always be called with a_level == 0 (leaf)
    if (a_node->m_level > a_level) {
        // Still above level for insertion, go down tree recursively
        Node *otherNode;

        // find the optimal branch for this record
        int index = PickBranch(&a_branch.m_rect, a_node);

        // recursively insert this record into the picked branch
        bools childWasSplit = InsertRectOnlyRegRec(a_branch, a_node->m_branch[index].m_child, &otherNode, a_level);

        if (!childWasSplit.first) {
            // Child was not split. Merge the bounding box of the new record with the
            // existing bounding box
            a_node->m_branch[index].m_rect = CombineRect(&a_branch.m_rect, &(a_node->m_branch[index].m_rect));
            return std::make_pair(false, childWasSplit.second);
        } else {
            // Child was split. The old branches are now re-partitioned to two nodes
            // so we have to re-calculate the bounding boxes of each node
            //            printf("calling NC in InsertRectRec\n");
            a_node->m_branch[index].m_rect = NodeCover(a_node->m_branch[index].m_child);
            Branch branch;
            branch.m_child = otherNode;
            //            printf("calling NC in InsertRectRec again\n");
            branch.m_rect = NodeCover(otherNode);

            // The old node is already a child of a_node. Now add the newly-created
            // node to a_node as well. a_node might be split because of that.
            return std::make_pair(AddBranch(&branch, a_node, a_newNode), childWasSplit.second);
        }
    }
    else if (a_node->m_level == a_level) {
        // We have reached level for insertion. Add rect, split if necessary
        //        cout << "levels were equal calling addbranch\n";

        if(a_level == 0 && a_node->isRegular()) {
            return std::make_pair(AddBranch(&a_branch, a_node, a_newNode), true);
        } else{
            return std::make_pair(false, false);
        }
    } else {
        // Should never occur
        //                ASSERT(0);
        return std::make_pair(false, false);
    }
}


// returns <was_split, inserted_or_not>
RTREE_TEMPLATE
typename RTREE_QUAL::bools RTREE_QUAL::InsertRectOnlyRegRec(const Branch &a_branch, Node *a_node, Node **a_newNode, int a_level, Rect &a_coverRect, Rect &a_newCoverRect) {
    //            ASSERT(a_node && a_newNode);
    //            ASSERT(a_level >= 0 && a_level <= a_node->m_level);

    // recurse until we reach the correct level for the new record. data records
    // will always be called with a_level == 0 (leaf)
    if (a_node->m_level > a_level) {
        // Still above level for insertion, go down tree recursively
        Node *otherNode;

        // find the optimal branch for this record
        int index = PickBranch(&a_branch.m_rect, a_node);

        // recursively insert this record into the picked branch
        bools childWasSplit = InsertRectOnlyRegRec(a_branch, a_node->m_branch[index].m_child, &otherNode, a_level, a_coverRect, a_newCoverRect);

        if (!childWasSplit.first) {
            // Child was not split. Merge the bounding box of the new record with the
            // existing bounding box
            a_node->m_branch[index].m_rect = CombineRect(&a_branch.m_rect, &(a_node->m_branch[index].m_rect));
            return std::make_pair(false, childWasSplit.second);
        } else {
            // Child was split. The old branches are now re-partitioned to two nodes
            // so we have to re-calculate the bounding boxes of each node
            //            printf("calling NC in InsertRectRec\n");
            //            a_node->m_branch[index].m_rect = NodeCover(a_node->m_branch[index].m_child);
            a_node->m_branch[index].m_rect = Rect(a_coverRect.m_min, a_coverRect.m_max);
            Branch branch;
            branch.m_child = otherNode;
            //            printf("calling NC in InsertRectRec again\n");
            branch.m_rect = NodeCover(otherNode);
            branch.m_rect = Rect(a_newCoverRect.m_min, a_newCoverRect.m_max);

            // The old node is already a child of a_node. Now add the newly-created
            // node to a_node as well. a_node might be split because of that.
            return std::make_pair(AddBranch(&branch, a_node, a_newNode, a_coverRect, a_newCoverRect), childWasSplit.second);
        }
    }
    else if (a_node->m_level == a_level) {
        // We have reached level for insertion. Add rect, split if necessary
        //        cout << "levels were equal calling addbranch\n";

        if(a_level == 0 && a_node->isRegular()) {
            return std::make_pair(AddBranch(&a_branch, a_node, a_newNode, a_coverRect, a_newCoverRect), true);
        } else{
            return std::make_pair(false, false);
        }
    } else {
        // Should never occur
        //                ASSERT(0);
        return std::make_pair(false, false);
    }
}



RTREE_TEMPLATE
typename RTREE_QUAL::bools RTREE_QUAL::InsertRectOnlyReg(const Branch &a_branch, Node **a_root, int a_level) {
    //            ASSERT(a_root);
    //            ASSERT(a_level >= 0 && a_level <= (*a_root)->m_level);


    Node *newNode;

    //    cout << "in insert rect\n";
    Rect rect1, rect2;
    bools root_was_split_and_found = InsertRectOnlyRegRec(a_branch, *a_root, &newNode, a_level, rect1, rect2);
    if (root_was_split_and_found.first)  // Root split
        {
        // Grow tree taller and new root
        Node *newRoot = AllocNode();
        newRoot->m_level = (*a_root)->m_level + 1;
        newRoot->m_parent = NULL;
        // DEBUG
        newRoot->name = "new root";

        Branch branch;

        // add old root node as a child of the new root
        //        printf("calling NC in InsertRect\n");
        //        branch.m_rect = NodeCover(*a_root);
        branch.m_rect = Rect(rect1.m_min, rect1.m_max);
        branch.m_child = *a_root;
        AddBranch(&branch, newRoot, NULL);

        // add the split node as a child of the new root
        //        printf("calling NC in InsertRectRec again\n");
        //        branch.m_rect = NodeCover(newNode);
        branch.m_rect = Rect(rect2.m_min, rect2.m_max);
        branch.m_child = newNode;
        AddBranch(&branch, newRoot, NULL);

        // set the new root as the root node
        *a_root = newRoot;

        return std::make_pair(true, root_was_split_and_found.second);
        }

    return std::make_pair(false, root_was_split_and_found.second);
}


// Find the smallest rectangle that includes all rectangles in branches of a node.
RTREE_TEMPLATE
typename RTREE_QUAL::Rect RTREE_QUAL::NodeCover(Node *a_node) {
    //            ASSERT(a_node);
    Rect rect;
    Rect alaki_rect;
    if( !a_node->IsLeaf() || a_node->m_transferred) {
        //        cout << "NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO\n";
        rect = a_node->m_branch[0].m_rect;
        for (int index = 1; index < a_node->m_count; ++index) {
            rect = CombineRect(&rect, &(a_node->m_branch[index].m_rect));
        }
    } else{
//        rect = m_data_arr[a_node->m_L]->m_rect;
        rect = Rect(m_data_arr_mins[a_node->m_L], m_data_arr_maxes[a_node->m_L]);
        for (int index = a_node->m_L + 1; index < a_node->m_R; ++index) {
//            rect = CombineRect(&rect, &(m_data_arr[index]->m_rect));
            alaki_rect = Rect(m_data_arr_mins[index], m_data_arr_maxes[index]);
            rect = CombineRect(&rect, &(alaki_rect));
        }
    }
    //    cout << "\n" << "here!!! " << rect.m_min[0] << rect.m_min[1] << ", " << rect.m_max[0] << rect.m_max[1] << "\n";
    return rect;
}


// 11 oct
RTREE_TEMPLATE
typename RTREE_QUAL::Rect RTREE_QUAL::LeafCover(int L, int R){
    Rect rect = Rect(m_data_arr_mins[L], m_data_arr_maxes[L]);
    for (int index = L + 1; index < R; ++index) {
        for(int i = 0; i < NUMDIMS; i++){
            if(rect.m_min[i] > m_data_arr_mins[index][i]) rect.m_min[i] = m_data_arr_mins[index][i];
            if(rect.m_max[i] < m_data_arr_maxes[index][i]) rect.m_max[i] = m_data_arr_maxes[index][i];
        }

    }

    // debug
    //    for(int i = 0; i < NUMDIMS; i++){
    //        if(rect.m_min[i] > rect.m_max[i]) {cout << "INJA GHALATE; CHE GOHI KHORDI AKHE. " << i << " " << rect.m_min[i] << " " << rect.m_max[i] << endl;}
    //    }

    return rect;
}

// FATEMEH DEBUG
// Find the smallest rectangle that includes all rectangles in branches of a node.
RTREE_TEMPLATE
typename RTREE_QUAL::Rect RTREE_QUAL::FatemehCover(std::vector<Branch> m_branch) {
    if(m_branch.size() == 0){return Rect();}

    Rect rect = m_branch[0].m_rect;
    for (int index = 1; index < m_branch.size(); ++index) {
        rect = CombineRect(&rect, &(m_branch[index].m_rect));
    }
    //    cout << "\n" << "here!!! " << rect.m_min[0] << rect.m_min[1] << ", " << rect.m_max[0] << rect.m_max[1] << "\n";
    return rect;
}

// FATEMEH OUT


// Add a branch to a node.  Split the node if necessary.
// Returns 0 if node not split.  Old node updated.
// Returns 1 if node split, sets *new_node to address of new node.
// Old node updated, becomes one of two.
// RTREE_TEMPLATE
// bool RTREE_QUAL::AddBranch(const Branch *a_branch, Node *a_node, Node **a_newNode) {
//     //            ASSERT(a_branch);
//     //            ASSERT(a_node);
//     // regularity check
//     if (a_node->isRegular()) {
//         int max_children;
//         if(a_node->IsLeaf()) {max_children = MAXDATAPOINTS;} else {max_children = MAXNODES;}

// //        if (a_node->IsLeaf() && !a_node->m_transferred){
// //            //                    ASSERT(a_node->isRegular());
// //            // transfer the stuff into array
// //            for(int i=a_node->m_L; i < a_node->m_R; i++){
// //                a_node->m_data_points[(i - a_node->m_L)] = m_data_arr[i];
// //            }
// //            a_node->m_transferred = true;
// //        }

//         if (a_node->m_count < max_children)  // Split won't be necessary
//             {
//             //            printf("simple regular AddBranch\n");
//             if(a_node->IsLeaf()){
// //                a_node->m_data_points[a_node->m_count] = *a_branch;
//                 cout << "INSERTING INTO A LEAF! THIS IS NOT SUPPORTED NOW!!!! \n";
//             } else {
//                 a_node->m_branch[a_node->m_count] = *a_branch;
//             }
//             if (a_node->m_level > 0) (a_branch->m_child)->m_parent = a_node;
//             //            a_node->m_branch.push_back(*a_branch);
//             //            printf("tracking line 1\n");
//             ++a_node->m_count;
//             //            printf("tracking line 2\n");
//             return false;
//             //            printf("for fatemeh: %d, %d, %d, %d\n", a_branch->m_rect.m_min[0],a_branch->m_rect.m_min[1],a_branch->m_rect.m_max[0],a_branch->m_rect.m_max[1]);

//             } else {
//             //                    ASSERT(a_newNode);
//             // only called on regular nodes!:)
//             //            printf("!!!!!!!!!!!!! had to split !!!!!!!!!!!!!!!!!!!!!!\n");
//             SplitNode(a_node, a_branch, a_newNode);
//             return true;
//         }
//     } else {
//         // is irregular
//         // should only be called on irr leaves from the 2dcrack func

//         //must be a leaf
//         //        printf("adding branch to irregular leaf %d, %d, %d, %d\n", a_branch->m_rect.m_min[0], a_branch->m_rect.m_min[1], a_branch->m_rect.m_max[0], a_branch->m_rect.m_max[1]);
//         //                ASSERT(a_node->IsLeaf());
//         //        a_node->m_branch.push_back(*a_branch);
//         //        ++a_node->m_count;
//         cout << "REMEMBER YOU HAVE NOT FIXED THIS!! ADDBRANCH FOR IRREG LEAF\n";
//         cout << "NEVER MIND SHOULD NEVER HAVE HAPPENED!\n";
//         return false;
//     }
// }


// 2023
// new formatted addBranch
RTREE_TEMPLATE
bool RTREE_QUAL::AddBranch(const Branch *a_branch, Node *a_node, Node **a_newNode) {
    //            ASSERT(a_branch);
    //            ASSERT(a_node);
    #ifdef stats
    if(a_node->m_level > 1)
        {count_internal_nodes++;}
    else{
        if(a_branch->m_child->isRegular())
            count_regular_leaves++;
        else
            count_irregular_leaves++;
        // cout << "in addbrach 1 added one leaf " << count_leaves <<  endl;
    }
    #endif

    if(a_node->IsInternalNode()){
        // regular internal node
        if (a_node->m_count < MAXNODES) // Split won't be necessary
        {
            a_node->m_branch[a_node->m_count] = *a_branch;
            (a_branch->m_child)->m_parent = a_node;
            ++a_node->m_count;
            return false;
        }
        else{
            // only called on regular nodes!:)
            SplitNode(a_node, a_branch, a_newNode);
            return true;
        }
    }
    else
    {
        // it's a leaf, so we have to check regularity
        if (a_node->isRegular()) 
        {
            cout << "INSERTING INTO REGUALR LEAF, NOT SUPPORTED LIKE THIS, SHOULD NOT HAPPEN!!!\n";
            return false;
        }
        else
        {
            cout << "INSERTING INTO IRREGUALR LEAF, NOT SUPPORTED LIKE THIS, SHOULD NOT HAPPEN!!!\n";
            return false;
        }
    }
    
}





// Add a branch to a node.  Split the node if necessary.
// Returns 0 if node not split.  Old node updated.
// Returns 1 if node split, sets *new_node to address of new node.
// Old node updated, becomes one of two.
// RTREE_TEMPLATE
// bool RTREE_QUAL::AddBranch(const Branch *a_branch, Node *a_node, Node **a_newNode, Rect &a_coverRect, Rect &a_newCoverRect) {
//     //            ASSERT(a_branch);
//     //            ASSERT(a_node);
//     // regularity check
//     if (a_node->isRegular()) {
//         int max_children;
//         if(a_node->IsLeaf()) {max_children = MAXDATAPOINTS;} else {max_children = MAXNODES;}

// //        if (a_node->IsLeaf() && !a_node->m_transferred){
// //            //                    ASSERT(a_node->isRegular());
// //            // transfer the stuff into array
// //            for(int i=a_node->m_L; i < a_node->m_R; i++){
// //                a_node->m_data_points[(i - a_node->m_L)] = m_data_arr[i];
// //            }
// //            a_node->m_transferred = true;
// //        }

//         if (a_node->m_count < max_children)  // Split won't be necessary
//             {
//             //            printf("simple regular AddBranch\n");
//             if(a_node->IsLeaf()){
// //                a_node->m_data_points[a_node->m_count] = *a_branch;
//                 cout << "INSERTING INTO LEAF, NOT SUPPORTED!!!\n";
//             } else {
//                 a_node->m_branch[a_node->m_count] = *a_branch;
//             }
//             if (a_node->m_level > 0) (a_branch->m_child)->m_parent = a_node;
//             //            a_node->m_branch.push_back(*a_branch);
//             //            printf("tracking line 1\n");
//             ++a_node->m_count;
//             //            printf("tracking line 2\n");
//             return false;
//             //            printf("for fatemeh: %d, %d, %d, %d\n", a_branch->m_rect.m_min[0],a_branch->m_rect.m_min[1],a_branch->m_rect.m_max[0],a_branch->m_rect.m_max[1]);

//             } else {
//             //                    ASSERT(a_newNode);
//             // only called on regular nodes!:)
//             //            printf("!!!!!!!!!!!!! had to split !!!!!!!!!!!!!!!!!!!!!!\n");
//             SplitNode(a_node, a_branch, a_newNode, a_coverRect, a_newCoverRect);
//             return true;
//         }
//     } else {
//         // is irregular
//         // should only be called on irr leaves from the 2dcrack func

//         //must be a leaf
//         //        printf("adding branch to irregular leaf %d, %d, %d, %d\n", a_branch->m_rect.m_min[0], a_branch->m_rect.m_min[1], a_branch->m_rect.m_max[0], a_branch->m_rect.m_max[1]);
//         //                ASSERT(a_node->IsLeaf());
//         //        a_node->m_branch.push_back(*a_branch);
//         //        ++a_node->m_count;
//         cout << "REMEMBER YOU HAVE NOT FIXED THIS!! ADDBRACNH FOR IRREG LEAF\n";
//         cout << "NEVER MIND SHOULD NEVER HAVE HAPPENED!\n";
//         return false;
//     }
// }


// 2023
// new formatted addBranch
RTREE_TEMPLATE
bool RTREE_QUAL::AddBranch(const Branch *a_branch, Node *a_node, Node **a_newNode, Rect &a_coverRect, Rect &a_newCoverRect) {
    //            ASSERT(a_branch);
    //            ASSERT(a_node);
    #ifdef stats
    if(a_node->m_level > 1)
        {count_internal_nodes++;}
    else{
        if(a_branch->m_child->isRegular())
            count_regular_leaves++;
        else
            count_irregular_leaves++;
        // cout << "in addbrach 1 added one leaf " << count_leaves <<  endl;
    }
    #endif
    if(a_node->IsInternalNode()){
        // regular internal node
        if (a_node->m_count < MAXNODES) // Split won't be necessary
        {
            a_node->m_branch[a_node->m_count] = *a_branch;
            (a_branch->m_child)->m_parent = a_node;
            ++a_node->m_count;
            return false;
        }
        else{
            // only called on regular nodes!:)
            SplitNode(a_node, a_branch, a_newNode, a_coverRect, a_newCoverRect);
            return true;
        }
    }
    else
    {
        // it's a leaf, so we have to check regularity
        if (a_node->isRegular()) 
        {
            cout << "INSERTING INTO REGUALR LEAF, NOT SUPPORTED LIKE THIS, SHOULD NOT HAPPEN!!!\n";
            return false;
        }
        else
        {
            cout << "INSERTING INTO IRREGUALR LEAF, NOT SUPPORTED LIKE THIS, SHOULD NOT HAPPEN!!!\n";
            return false;
        }
    }
    
}


// Disconnect a dependent node.
// Caller must return (or stop using iteration index) after this as count has changed
RTREE_TEMPLATE
void RTREE_QUAL::DisconnectBranch(Node *a_node, int a_index) {
    //            ASSERT(a_node && (a_index >= 0) && (a_index < MAXNODES));
    //            ASSERT(a_node && (a_index >= 0));
    //            ASSERT(a_node);
    //            ASSERT(a_index >=0);
    //            ASSERT(a_node->m_count > 0);

    #ifdef stats
        if(a_node->m_level > 1) count_internal_nodes--;
        else
        { 
            if((a_node->m_branch[a_index]).m_child->isRegular())
                count_regular_leaves--;
            else
                count_irregular_leaves--;
            // auto t = CountLeaves();
            // cout << "in Disconnect removed one leaf " << count_leaves<< " " << t.first + t.second << endl;
        }
    #endif

    // MARCH 2023
    if(a_node->m_level == 1){
        // cout << "in disconnect, before oldest check" << endl;
        // let's fix the linked list things first
        if((a_node->m_branch[a_index]).m_child == m_oldest_leaf){
            // is head
            // cout << "was head" << endl;
            m_oldest_leaf = ((a_node->m_branch[a_index]).m_child)->m_younger_brother;
        }
        else{
            // cout << "was not head" << endl;
            (((a_node->m_branch[a_index]).m_child)->m_older_brother)->m_younger_brother = ((a_node->m_branch[a_index]).m_child)->m_younger_brother;
        }
        // cout << "before youngest check" << endl;
        if((a_node->m_branch[a_index]).m_child == m_youngest_leaf){
            // cout << "was tail" << endl;
            m_youngest_leaf = ((a_node->m_branch[a_index]).m_child)->m_older_brother;
        }
        else{
            // cout << "was not tail" << endl;
            (((a_node->m_branch[a_index]).m_child)->m_younger_brother)->m_older_brother = ((a_node->m_branch[a_index]).m_child)->m_older_brother;
        }
    }
    // MARCH 2023
    // cout << "after linkedlist fixes" << endl;

    // cout << "a_node.count" << a_node->m_count << " swapping " << a_index << " and " << a_node->m_count - 1 << endl;

    // Remove element by swapping with the last element to prevent gaps in array
    a_node->m_branch[a_index] = a_node->m_branch[a_node->m_count - 1];
    //    if(a_node->m_branch[a_index].m_child != NULL){
    //        std::vector<Branch>().swap((a_node->m_branch[a_index].m_child)->m_branch);
    //    }
    //    a_node->m_branch[a_index] = a_node->m_branch.back();
    //    a_node->m_branch.pop_back();
    a_node->m_count -= 1;
    // cout << "new count " << a_node->m_count << endl;
    //    if(a_node->m_count == 0) {printf("count zero in disB \n");}
    int max_children;
    if(a_node->IsLeaf()) { max_children = MAXDATAPOINTS;} else {max_children = MAXNODES;}
    if (a_node->m_count <= max_children) a_node->m_regular = true;


}




// Pick a branch.  Pick the one that will need the smallest increase
// in area to accomodate the new rectangle.  This will result in the
// least total area for the covering rectangles in the current node.
// In case of a tie, pick the one which was smaller before, to get
// the best resolution when searching.
RTREE_TEMPLATE
int RTREE_QUAL::PickBranch(const Rect *a_rect, Node *a_node) {
    //            ASSERT(a_rect && a_node);

    bool firstTime = true;
    float increase;
    float bestIncr = (float) -1;
    float area;
    float bestArea;
    int best = 0;
    Rect tempRect;
    // cout << "**********************" << endl;
    // cout << "in pick branch, a_rect: (" << a_rect->m_min[0] << ", " << a_rect->m_min[1] << ") (" << a_rect->m_max[0] << ", " << a_rect->m_max[1] << ")" << endl; 
    for (int index = 0; index < a_node->m_count; ++index) {
        
        Rect *curRect = &a_node->m_branch[index].m_rect;
        // cout << "in loop," << index << " id: " << a_node->m_branch[index].m_child->m_id << " cureRect: (" << curRect->m_min[0] << ", " << curRect->m_min[1] << ") (" << curRect->m_max[0] << ", " << curRect->m_max[1] << ")" << endl;
        area = CalcRectVolume(curRect);
        tempRect = CombineRect(a_rect, curRect);
        increase = CalcRectVolume(&tempRect) - area;
        // cout << "increase: " << increase<< endl;
        if ((increase < bestIncr) || firstTime) {
            best = index;
            bestArea = area;
            bestIncr = increase;
            firstTime = false;
        } else if ((increase == bestIncr) && (area < bestArea)) {
            best = index;
            bestArea = area;
            bestIncr = increase;
        }
        // cout << "best: " << best << endl;
    }
    // cout << "**********************" << endl;
    return best;
}

RTREE_TEMPLATE
int RTREE_QUAL::PickBranchReg(const Rect *a_rect, Node *a_node) {
    //            ASSERT(a_rect && a_node);

    bool firstTime = true;
    float increase;
    float bestIncr = (float) -1;
    float area;
    float bestArea;
    int best = 0;
    Rect tempRect;

    int index = 0;
    for(int index=0; index < a_node->m_count; ++index)
    {

        //    for (typename std::vector<Branch>::iterator it = a_node->m_branch.begin(); it != a_node->m_branch.end(); it++) {
        if (a_node->m_branch[index].m_child->isRegular()) {
            Rect* curRect = &(a_node->m_branch[index].m_rect);
            area = CalcRectVolume(curRect);
            tempRect = CombineRect(a_rect, curRect);
            increase = CalcRectVolume(&tempRect) - area;
            if ((increase < bestIncr) || firstTime) {
                best = index;
                bestArea = area;
                bestIncr = increase;
                firstTime = false;
            } else if ((increase == bestIncr) && (area < bestArea)) {
                best = index;
                bestArea = area;
                bestIncr = increase;
            }
        }
        index++;
    }
    return best;
}


// Combine two rectangles into larger one containing both
RTREE_TEMPLATE
typename RTREE_QUAL::Rect RTREE_QUAL::CombineRect(const Rect *a_rectA, const Rect *a_rectB) {
    //            ASSERT(a_rectA && a_rectB);

    Rect newRect;

    for (int index = 0; index < NUMDIMS; ++index) {
        //        printf("IDKY %d, %d \n", a_rectA->m_min[index], a_rectB->m_min[index]);
        newRect.m_min[index] = Min(a_rectA->m_min[index], a_rectB->m_min[index]);
        newRect.m_max[index] = Max(a_rectA->m_max[index], a_rectB->m_max[index]);
    }

    return newRect;
}



// Split a node.
// Divides the nodes branches and the extra one between two nodes.
// Old node is one of the new ones, and one really new one is created.
// Tries more than one method for choosing a partition, uses best result.
RTREE_TEMPLATE
void RTREE_QUAL::SplitNode(Node *a_node, const Branch *a_branch, Node **a_newNode) {
    //            ASSERT(a_node);
    //            ASSERT(a_branch);

    // Could just use local here, but member or external is faster since it is reused
    PartitionVars localVars;
    PartitionVars *parVars = &localVars;

    // Load all the branches into a buffer, initialize old node
    GetBranches(a_node, a_branch, parVars);

    // Find partition
    int min_fill;
    if(a_node->IsLeaf()) {min_fill = MINDATAPOINTS;} else {min_fill = MINNODES;}
    ChoosePartition(parVars, min_fill);

    if(parVars->m_count[0] == 1 || parVars->m_count[1] == 1){
        printf("WHAT THE FUDGE\n");
    }

    // Create a new node to hold (about) half of the branches
    *a_newNode = AllocNode();
    (*a_newNode)->m_level = a_node->m_level;
    (*a_newNode)->name = "from split node";

    // Put branches from buffer into 2 nodes according to the chosen partition
    a_node->m_count = 0;
    //    a_node->m_branch.clear();
    LoadNodes(a_node, *a_newNode, parVars);

    //            ASSERT((a_node->m_count + (*a_newNode)->m_count) == parVars->m_total);
}


RTREE_TEMPLATE
void RTREE_QUAL::SplitNode(Node *a_node, const Branch *a_branch, Node **a_newNode, Rect &a_coverRect, Rect &a_newCoverRect) {
    //            ASSERT(a_node);
    //            ASSERT(a_branch);

    // Could just use local here, but member or external is faster since it is reused
    PartitionVars localVars;
    PartitionVars *parVars = &localVars;

    // Load all the branches into a buffer, initialize old node
    GetBranches(a_node, a_branch, parVars);

    // Find partition
    int min_fill;
    if(a_node->IsLeaf()) {min_fill = MINDATAPOINTS;} else {min_fill = MINNODES;}
    ChoosePartition(parVars, min_fill);

    if(parVars->m_count[0] == 1 || parVars->m_count[1] == 1){
        printf("WHAT THE FUDGE\n");
    }

    // Create a new node to hold (about) half of the branches
    *a_newNode = AllocNode();
    (*a_newNode)->m_level = a_node->m_level;
    (*a_newNode)->name = "from split node";

    // Put branches from buffer into 2 nodes according to the chosen partition
    a_node->m_count = 0;
    //    a_node->m_branch.clear();
    LoadNodes(a_node, *a_newNode, parVars);
    a_coverRect = (parVars->m_cover[0]);
    a_newCoverRect = (parVars->m_cover[1]);

    //            ASSERT((a_node->m_count + (*a_newNode)->m_count) == parVars->m_total);
}


// Calculate the n-dimensional volume of a rectangle
RTREE_TEMPLATE
float RTREE_QUAL::RectVolume(Rect *a_rect) {
    //            ASSERT(a_rect);

    float volume = (float) (a_rect->m_max[0] - a_rect->m_min[0]);
    float len;

    for (int index = 1; index < NUMDIMS; ++index) {
        //        if(volume < 0) {
        //            cout << index << " volume went below zero here\n";
        //            cout << "the rect was: " << a_rect->m_min[0] << " " << a_rect->m_min[1] << " " << a_rect->m_max[0] << " " << a_rect->m_max[1] << endl;
        //        }

        len = a_rect->m_max[index] - a_rect->m_min[index];
        //        if(len < 0) {cout << "len below zero " << a_rect->m_max[index] << " " << a_rect->m_min[index] << " " << len << endl;}
        //                ASSERT(len >= (float) 0);
        //        cout << index << " len here: " << len << endl;
        volume *= (len);
    }
        if(volume < 0) {
            cout << " after for volume went below zero here\n";
            cout << "the rect was: " << a_rect->m_min[0] << " " << a_rect->m_min[1] << " " << a_rect->m_min[2] << " " << a_rect->m_max[0] << " " << a_rect->m_max[1] << " " << a_rect->m_max[2] << endl;
        }

    ASSERT(volume >= (float) 0);

    return volume;
}


// The exact volume of the bounding sphere for the given Rect
RTREE_TEMPLATE
float RTREE_QUAL::RectSphericalVolume(Rect *a_rect) {
    //            ASSERT(a_rect);

    float sumOfSquares = (float) 0;
    float radius;

    for (int index = 0; index < NUMDIMS; ++index) {
        float halfExtent = ((float) a_rect->m_max[index] - (float) a_rect->m_min[index]) * 0.5f;
        sumOfSquares += halfExtent * halfExtent;
    }

    radius = (float) sqrt(sumOfSquares);

    // Pow maybe slow, so test for common dims like 2,3 and just use x*x, x*x*x.
    if (NUMDIMS == 3) {
        return (radius * radius * radius * m_unitSphereVolume);
    } else if (NUMDIMS == 2) {
        return (radius * radius * m_unitSphereVolume);
    } else {
        return (float) (pow(radius, NUMDIMS) * m_unitSphereVolume);
    }
}


// Use one of the methods to calculate retangle volume
RTREE_TEMPLATE
float RTREE_QUAL::CalcRectVolume(Rect *a_rect) {
#ifdef RTREE_USE_SPHERICAL_VOLUME
    return RectSphericalVolume(a_rect); // Slower but helps certain merge cases
#else // RTREE_USE_SPHERICAL_VOLUME
return RectVolume(a_rect); // Faster but can cause poor merges
#endif // RTREE_USE_SPHERICAL_VOLUME
}


// Load branch buffer with branches from full node plus the extra branch.
RTREE_TEMPLATE
void RTREE_QUAL::GetBranches(Node *a_node, const Branch *a_branch, PartitionVars *a_parVars) {
    //            ASSERT(a_node);
    //            ASSERT(a_branch);

    //            ASSERT(a_node->m_count == MAXNODES);
    // TODO technically this should work, seeing as partitioning is just for regular nodes, but make sure
    // Load the branch buffer
    //    for (int index = 0; index < MAXNODES; ++index) {
    //        a_parVars->m_branchBuf[index] = a_node->m_branch[index];
    //    }
    //    a_parVars->m_branchBuf[MAXNODES] = *a_branch;
    int max_children;
    if(a_node->IsLeaf()){max_children = MAXDATAPOINTS;} else {max_children = MAXNODES;}
    for (int index = 0; index < max_children; ++index) {
        //        a_parVars->m_branchBuf.push_back(a_node->m_branch[index]);
        if(a_node->IsLeaf()){
//            a_parVars->m_branchBuf[index] = a_node->m_data_points[index];
// this never happens i think
// but I have to include it any way
//            a_parVars->m_branchBuf[index] = *m_data_arr[a_node->m_L + index];
            a_parVars->m_branchBuf[index] = Branch(NULL, Rect(m_data_arr_mins[a_node->m_L + index], m_data_arr_maxes[a_node->m_L + index]), m_data_arr_ids[a_node->m_L + index]);
        } else {
            a_parVars->m_branchBuf[index] = a_node->m_branch[index];
        }
    }
    //    a_parVars->m_branchBuf.push_back(*a_branch);
    a_parVars->m_branchBuf[max_children] = *a_branch;
    a_parVars->m_branchCount = max_children + 1;

    // Calculate rect containing all in the set
    a_parVars->m_coverSplit = a_parVars->m_branchBuf[0].m_rect;
    for (int index = 1; index < max_children + 1; ++index) {
        a_parVars->m_coverSplit = CombineRect(&a_parVars->m_coverSplit, &a_parVars->m_branchBuf[index].m_rect);
    }
    //    cout << "calling crv in get branches\n";
    a_parVars->m_coverSplitArea = CalcRectVolume(&a_parVars->m_coverSplit);
}


// Method #0 for choosing a partition:
// As the seeds for the two groups, pick the two rects that would waste the
// most area if covered by a single rectangle, i.e. evidently the worst pair
// to have in the same group.
// Of the remaining, one at a time is chosen to be put in one of the two groups.
// The one chosen is the one with the greatest difference in area expansion
// depending on which group - the rect most strongly attracted to one group
// and repelled from the other.
// If one group gets too full (more would force other group to violate min
// fill requirement) then other group gets the rest.
// These last are the ones that can go in either group most easily.
RTREE_TEMPLATE
void RTREE_QUAL::ChoosePartition(PartitionVars *a_parVars, int a_minFill) {
    //            ASSERT(a_parVars);

    float biggestDiff;
    int group, chosen = 0, betterGroup = 0;

    InitParVars(a_parVars, a_parVars->m_branchCount, a_minFill);
    PickSeeds(a_parVars);

    //    cout << "in CHOOSE PARTITION m_total " << a_parVars->m_total << endl;

    while (((a_parVars->m_count[0] + a_parVars->m_count[1]) < a_parVars->m_total)
    && (a_parVars->m_count[0] < (a_parVars->m_total - a_parVars->m_minFill))
    && (a_parVars->m_count[1] < (a_parVars->m_total - a_parVars->m_minFill))) {
        biggestDiff = (float) -1;
        for (int index = 0; index < a_parVars->m_total; ++index) {
            if (PartitionVars::NOT_TAKEN == a_parVars->m_partition[index]) {
                Rect *curRect = &a_parVars->m_branchBuf[index].m_rect;
                Rect rect0 = CombineRect(curRect, &a_parVars->m_cover[0]);
                Rect rect1 = CombineRect(curRect, &a_parVars->m_cover[1]);
                //                cout << "calling crv in choose partition\n";
                float growth0 = CalcRectVolume(&rect0) - a_parVars->m_area[0];
                float growth1 = CalcRectVolume(&rect1) - a_parVars->m_area[1];
                float diff = growth1 - growth0;
                if (diff >= 0) {
                    group = 0;
                } else {
                    group = 1;
                    diff = -diff;
                }

                if (diff > biggestDiff) {
                    biggestDiff = diff;
                    chosen = index;
                    betterGroup = group;
                } else if ((diff == biggestDiff) && (a_parVars->m_count[group] < a_parVars->m_count[betterGroup])) {
                    chosen = index;
                    betterGroup = group;
                }
            }
        }
        //        cout << "CALLING CLASSIFY IN CHOOSEPARTITION FOR " << chosen << endl;
        Classify(chosen, betterGroup, a_parVars);
    }

    // If one group too full, put remaining rects in the other
    if ((a_parVars->m_count[0] + a_parVars->m_count[1]) < a_parVars->m_total) {
        if (a_parVars->m_count[0] >= a_parVars->m_total - a_parVars->m_minFill) {
            group = 1;
        } else {
            group = 0;
        }
        for (int index = 0; index < a_parVars->m_total; ++index) {
            if (PartitionVars::NOT_TAKEN == a_parVars->m_partition[index]) {
                //                cout << "CALLING CLASSIFY IN CHOOSEPARTITION FOR THE REST " << index << endl;

                Classify(index, group, a_parVars);
            }
        }
    }

    //            ASSERT((a_parVars->m_count[0] + a_parVars->m_count[1]) == a_parVars->m_total);
    //            ASSERT((a_parVars->m_count[0] >= a_parVars->m_minFill) &&
    //                   (a_parVars->m_count[1] >= a_parVars->m_minFill));
}


// Copy branches from the buffer into two nodes according to the partition.
RTREE_TEMPLATE
void RTREE_QUAL::LoadNodes(Node *a_nodeA, Node *a_nodeB, PartitionVars *a_parVars) {
    //            ASSERT(a_nodeA);
    //            ASSERT(a_nodeB);
    //            ASSERT(a_parVars);

    for (int index = 0; index < a_parVars->m_total; ++index) {
        //                ASSERT(a_parVars->m_partition[index] == 0 || a_parVars->m_partition[index] == 1);

        int targetNodeIndex = a_parVars->m_partition[index];
        Node *targetNodes[] = {a_nodeA, a_nodeB};

        // It is assured that AddBranch here will not cause a node split.
        bool nodeWasSplit = AddBranch(&a_parVars->m_branchBuf[index], targetNodes[targetNodeIndex], NULL);
        //                ASSERT(!nodeWasSplit);
         #ifdef stats
            if(a_nodeA->m_level == 1){
                if(a_parVars->m_branchBuf[index].m_child->isRegular())
                    count_regular_leaves--;
                else
                    count_irregular_leaves--;
            } 
        #endif
    }
     #ifdef stats
        if(a_nodeA->m_level > 1) 
            count_internal_nodes -= a_parVars->m_total;
        // else {
        //     count_leaves -= a_parVars->m_total;
        // }
    #endif
}


// Initialize a PartitionVars structure.
RTREE_TEMPLATE
void RTREE_QUAL::InitParVars(PartitionVars *a_parVars, int a_maxRects, int a_minFill) {
    //            ASSERT(a_parVars);

    a_parVars->m_count[0] = a_parVars->m_count[1] = 0;
    a_parVars->m_area[0] = a_parVars->m_area[1] = (float) 0;
    a_parVars->m_total = a_maxRects;
    a_parVars->m_minFill = a_minFill;
    for (int index = 0; index < a_maxRects; ++index) {
        a_parVars->m_partition[index] = PartitionVars::NOT_TAKEN;
    }
    //    for (int index = 0; index < a_maxRects; ++index) {
    //        a_parVars->m_partition.push_back(PartitionVars::NOT_TAKEN);
    //    }
}


RTREE_TEMPLATE
void RTREE_QUAL::PickSeeds(PartitionVars *a_parVars) {
    int seed0 = 0, seed1 = 0;
    float worst, waste;
    float area[a_parVars->m_total + 1];

    for (int index = 0; index < a_parVars->m_total; ++index) {
        //        cout << "calling crv in pick seeds\n";
        area[index] = CalcRectVolume(&a_parVars->m_branchBuf[index].m_rect);
    }


    worst = -1 * a_parVars->m_coverSplitArea - 1;
    //    cout << "in pickseeds before for " << a_parVars->m_total << " " << worst << " " << a_parVars->m_coverSplitArea << endl;

    for (int indexA = 0; indexA < a_parVars->m_total - 1; ++indexA) {
        for (int indexB = indexA + 1; indexB < a_parVars->m_total; ++indexB) {
            Rect oneRect = CombineRect(&a_parVars->m_branchBuf[indexA].m_rect, &a_parVars->m_branchBuf[indexB].m_rect);
            //            cout << "calling crv in pick seeds\n";
            waste = CalcRectVolume(&oneRect) - area[indexA] - area[indexB];
            //            cout << "waste: " << waste << endl;
            if (waste > worst) {
                worst = waste;
                seed0 = indexA;
                seed1 = indexB;
                //                cout << "in pickseeds, setting seeds: " << seed0 << " " << seed1 << endl;
            }
        }
    }
    //    cout << "CALLING CLASSIFY IN PICKSEEDS  " << a_parVars->m_total << " " << seed0 << " " << seed1 << endl;
    Classify(seed0, 0, a_parVars);
    Classify(seed1, 1, a_parVars);
}


// Put a branch in one of the groups.
RTREE_TEMPLATE
void RTREE_QUAL::Classify(int a_index, int a_group, PartitionVars *a_parVars) {
    //    cout << "CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << a_index << " " << a_parVars->m_partition[a_index] << endl;
    //            ASSERT(a_parVars);
    //            ASSERT(PartitionVars::NOT_TAKEN == a_parVars->m_partition[a_index]);

    a_parVars->m_partition[a_index] = a_group;

    // Calculate combined rect
    if (a_parVars->m_count[a_group] == 0) {
        a_parVars->m_cover[a_group] = a_parVars->m_branchBuf[a_index].m_rect;
    } else {
        a_parVars->m_cover[a_group] = CombineRect(&a_parVars->m_branchBuf[a_index].m_rect,
                                                  &a_parVars->m_cover[a_group]);
    }

    // Calculate volume of combined rect
    //    cout << "calling crv in classify\n";
    a_parVars->m_area[a_group] = CalcRectVolume(&a_parVars->m_cover[a_group]);

    ++a_parVars->m_count[a_group];
}




// Delete a data rectangle from an index structure.
// Pass in a pointer to a Rect, the tid of the record, ptr to ptr to root node.
// Returns 1 if record not found, 0 if success.
// RemoveRect provides for eliminating the root.
RTREE_TEMPLATE
bool RTREE_QUAL::RemoveRect(Rect *a_rect, const DATATYPE &a_id, Node **a_root) {
    //            ASSERT(a_rect && a_root);
    //            ASSERT(*a_root);

    ListNode *reInsertList = NULL;

    if (!RemoveRectRec(a_rect, a_id, *a_root, &reInsertList)) {
        // Found and deleted a data item
        // Reinsert any branches from eliminated nodes
        //        printf("in RemoveRect: found and deleted a data item\n");
        while (reInsertList) {

            Node *tempNode = reInsertList->m_node;
            //            printf("in RemoveRect: reinserting\n");

            for (int index = 0; index < tempNode->m_count; ++index) {
                // TODO go over this code. should I use (tempNode->m_level - 1)?
                InsertRect(tempNode->m_branch[index],
                           a_root,
                           tempNode->m_level);
            }

            ListNode *remLNode = reInsertList;
            reInsertList = reInsertList->m_next;

            FreeNode(remLNode->m_node);
            FreeListNode(remLNode);
        }

        // Check for redundant root (not leaf, 1 child) and eliminate TODO replace
        // if with while? In case there is a whole branch of redundant roots...
        if ((*a_root)->m_count == 1 && (*a_root)->IsInternalNode()) {
            Node *tempNode = (*a_root)->m_branch[0].m_child;

            //                    ASSERT(tempNode);
            FreeNode(*a_root);
            *a_root = tempNode;
        }
        return false;
    } else {
        return true;
    }
}



// removes all datapoint in a_ids from the node a_from_node
// return whether or not a_from_node was deleted
RTREE_TEMPLATE
bool RTREE_QUAL::RemoveBatchFromNode(Node *a_node, std::vector<Branch> a_ids, std::vector<int> a_indexes, Node **a_root){
    //            ASSERT(a_node && a_root);
    //            ASSERT(*a_root);

    ListNode *reInsertList = NULL;
    bools res = RemoveBatchRec(a_node, a_ids, a_indexes, *a_root, &reInsertList);
    if (!res.first) {
        // Found and deleted a data item
        // Reinsert any branches from eliminated nodes
        //        printf("in RemoveRect: found and deleted a data item\n");
        while (reInsertList) {

            Node *tempNode = reInsertList->m_node;
            //            printf("in RemoveRect: reinserting\n");

            for (int index = 0; index < tempNode->m_count; ++index) {
                // TODO go over this code. should I use (tempNode->m_level - 1)?
                printf("FALALALA %d\n", tempNode->m_level);
                InsertRect(tempNode->m_branch[index],
                           a_root,
                           (tempNode->m_level - 1));
            }

            ListNode *remLNode = reInsertList;
            reInsertList = reInsertList->m_next;

            FreeNode(remLNode->m_node);
            FreeListNode(remLNode);
        }

        // Check for redundant root (not leaf, 1 child) and eliminate TODO replace
        // if with while? In case there is a whole branch of redundant roots...
        if ((*a_root)->m_count == 1 && (*a_root)->IsInternalNode()) {
            //            cout << "NO MAN'S LAND...\n";
            Node *tempNode = (*a_root)->m_branch[0].m_child;
            tempNode->m_parent = NULL;

            //                    ASSERT(tempNode);
            FreeNode(*a_root);
            *a_root = tempNode;
        }
        //        return false;
        return res.second;
    } else {
        //        return true;
        return res.second;
    }
}



// Delete a rectangle from non-root part of an index structure.
// Called by RemoveRect.  Descends tree recursively,
// merges branches on the way back up.
// Returns 1 if record not found, 0 if success.
RTREE_TEMPLATE
bool RTREE_QUAL::RemoveRectRec(Rect *a_rect, const DATATYPE &a_id, Node *a_node, ListNode **a_listNode) {
    //            ASSERT(a_rect && a_node && a_listNode);
    //            ASSERT(a_node->m_level >= 0);

    if (a_node->IsInternalNode())  // not a leaf node
        {
        for (int index = 0; index < a_node->m_count; ++index) {
            //            cout  << "remove a_rect: " << a_rect->m_min[0] << " " << a_rect->m_min[1] << " " << a_rect->m_max[0] << " " << a_rect->m_max[1] << "\n";
            if (Overlap(a_rect, &(a_node->m_branch[index].m_rect))) {
                if (!RemoveRectRec(a_rect, a_id, a_node->m_branch[index].m_child, a_listNode)) {
                    //                    printf("returned from deleting from child, child's count is: %d  hey \n", a_node->m_branch[index].m_child->m_count);
                    int min_children;
                    if(a_node->m_branch[index].m_child->IsLeaf()) {min_children = MINDATAPOINTS;} else {min_children = MINNODES;}
                    if (a_node->m_branch[index].m_child->m_count >= min_children) {
                        // child removed, just resize parent rect
                        //                        printf("calling NC in RemoveRectRec\n");
                        a_node->m_branch[index].m_rect = NodeCover(a_node->m_branch[index].m_child);
                    } else {
                        // child removed, not enough entries in node, eliminate node
                        //                        printf("cascade delete\n");
                        ReInsert(a_node->m_branch[index].m_child, a_listNode);
                        DisconnectBranch(a_node, index); // Must return after this call as count has changed
                    }
                    return false;
                }
            }
        }
        return true;
        } else // A leaf node
        {
        //        printf("deleting from leaf, count is: %d  hey \n", a_node->m_count);
        //        if(a_node->m_count == 1){
        //            printf("here here \n");
        //        }
        for (int index = 0; index < a_node->m_count; ++index) {
            if (a_node->m_branch[index].m_data == a_id) {
                DisconnectBranch(a_node, index); // Must return after this call as count has changed
                return false;
            }
        }
        return true;
        }
}


// removes all datapoint in a_ids from the node a_from_node
// returns a pair of bools
// <!whether we found a_from_node, whether a_from_node was deleted>
RTREE_TEMPLATE
typename RTREE_QUAL::bools RTREE_QUAL::RemoveBatchRec(Node *a_from_node, std::vector<Branch> a_ids, std::vector<int> a_indexes, Node *a_node, ListNode **a_listNode){
    //            ASSERT(a_from_node && a_node && a_listNode);
    //            ASSERT(a_node->m_level >= 0);
    //    auto start_time = std::chrono::steady_clock::now();
    Rect rect = NodeCover(a_from_node);
    //    std::chrono::duration<double> res = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);
    //    cout << res.count() << " nodecover in remove batch" << endl;
    if (a_node->IsInternalNode())  // not a leaf node
        {
        for (int index = 0; index < a_node->m_count; ++index) {
            //            cout  << "remove a_rect: " << a_rect->m_min[0] << " " << a_rect->m_min[1] << " " << a_rect->m_max[0] << " " << a_rect->m_max[1] << "\n";
            if (Overlap(&rect, &(a_node->m_branch[index].m_rect))) {
                bools res = RemoveBatchRec(a_from_node, a_ids, a_indexes, a_node->m_branch[index].m_child, a_listNode);
                if (!res.first) {
                    //                    printf("returned from deleting from child, child's count is: %d  hey \n", a_node->m_branch[index].m_child->m_count);
                    int min_children;
                    if(a_node->m_branch[index].m_child->IsLeaf()) {min_children = MINDATAPOINTS;} else {min_children = MINNODES;}
                    if (a_node->m_branch[index].m_child->m_count >= min_children) {
                        // child removed, just resize parent rect
                        //                        printf("calling NC in RemoveRectRec\n");
                        a_node->m_branch[index].m_rect = NodeCover(a_node->m_branch[index].m_child);
                        return res;
                    } else {
                        // child removed, not enough entries in node, eliminate node
                        //                        printf("cascade delete\n");

                        ReInsert(a_node->m_branch[index].m_child, a_listNode);
                        DisconnectBranch(a_node, index); // Must return after this call as count has changed


                        // FATEMEH TIP:
                        // i think we should check if a_from_node is a child of this guy, I mean if we returned from him
                        // and then also fix the return values to return the fact that the a_from_leaf was destroyed
                        // or we could pass around a bool

                        if (a_node->m_level == 1){
                            // FATEMEH TIP
                            // this could be the a_from_node's parent
                            // actually this definitely is the parent
                            // if it returned false, and is level 1, it must be a_from_node's father
                            return std::make_pair(false, true);
                        }
                        return res;
                    }
                    //                    return false;
                }
            }
        }
        return std::make_pair(true, false);
        } else // A leaf node
        {
        //        printf("deleting from leaf, count is: %d  hey \n", a_node->m_count);
        //        if(a_node->m_count == 1){
        //            printf("here here \n");
        //        }
        //  TODO have to check that a_node and a_from_node are exactly the same, not just overlap DONE not tested
        Rect a_rect = NodeCover(a_from_node); Rect b_rect = NodeCover(a_node);
        if(a_rect.m_min[0] == b_rect.m_min[0] &&
        a_rect.m_min[1] == b_rect.m_min[1] &&
        a_rect.m_max[0] == b_rect.m_max[0] &&
        a_rect.m_max[1] == b_rect.m_max[1]) {
            // TODO check that all have been found and deleted
            //            start_time = std::chrono::steady_clock::now();
            //            for (int id_index = 0; id_index < a_ids.size(); id_index++) {
            //                for (int index = 0; index < a_node->m_count; ++index) {
            //                    if (a_node->m_branch[index].m_data == a_ids[id_index]) {
            //                        DisconnectBranch(a_node, index); // Must break after this call as count has changed
            //                        break;
            //                    }
            //                }
            //            }

            for(int i = a_indexes.size() - 1; i > -1; i--){
                //                ASSERT(a_node->m_branch[a_indexes[i]].m_data == a_ids[i].m_data);
                DisconnectBranch(a_node, a_indexes[i]);

            }
            //            cout << "back here \n";
            //            res = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);
            //            cout << res.count() << " remove for in batchremove" << endl;
            return std::make_pair(false, false); // found the node
        }
        return std::make_pair(true, false);
        }
}


// Decide whether two rectangles overlap.
RTREE_TEMPLATE
bool RTREE_QUAL::Overlap(Rect *a_rectA, Rect *a_rectB) const {
    //            ASSERT(a_rectA && a_rectB);

    for (int index = 0; index < NUMDIMS; ++index) {
        if (a_rectA->m_min[index] > a_rectB->m_max[index] ||
        a_rectB->m_min[index] > a_rectA->m_max[index]) {
            return false;
        }
    }
    return true;
}

// Decide whether two rectangles overlap.
RTREE_TEMPLATE
bool RTREE_QUAL::Overlap(Rect *a_rectA, const float a_min[NUMDIMS], const float a_max[NUMDIMS]) const {
    //            ASSERT(a_rectA && a_rectB);

    for (int index = 0; index < NUMDIMS; ++index) {
        if (a_rectA->m_min[index] > a_max[index] ||
        a_min[index] > a_rectA->m_max[index]) {
            return false;
        }
    }
    return true;
}


// Add a node to the reinsertion list.  All its branches will later
// be reinserted into the index structure.
RTREE_TEMPLATE
void RTREE_QUAL::ReInsert(Node *a_node, ListNode **a_listNode) {
    ListNode *newListNode;

    newListNode = AllocListNode();
    newListNode->m_node = a_node;
    newListNode->m_next = *a_listNode;
    *a_listNode = newListNode;
}


// Search in an index tree or subtree for all data retangles that overlap the argument rectangle.
RTREE_TEMPLATE
bool RTREE_QUAL::Search(Node *a_node, Rect *a_rect, int &a_foundCount, CrackVars *crackVars) {
    //            ASSERT(a_node);
    //            ASSERT(a_node->m_level >= 0);
    //            ASSERT(a_rect);
    //            ASSERT(a_node->m_count == a_node->m_branch.size());
    //    printf("in search\n");
    if (a_node->IsInternalNode()) {
        //        printf("was inetrnal node\n");
        // This is an internal node in the tree
        for (int index = 0; index < a_node->m_count; ++index) {
            // TODO I think here we have to compare with the NodeCover of the child, since the calculation might be old
            // this can happen because of the shrinks in the search
            //            if (Overlap(a_rect, &a_node->m_branch[index].m_rect))
            //            printf("calling NC in Search\n");
            //            if(a_node->m_branch[index].m_child->m_branch == NULL) {printf("HEY THERE !!\n");}
            //            printf("I THINK IT'S HERE\n");
            //            if(a_node->m_branch[index].m_child->m_branch.empty()){printf("HERE HERE HERE count = %d \n", a_node->m_branch[index].m_child->m_count);}
            //            a_node->m_branch[index].m_rect = NodeCover(a_node->m_branch[index].m_child);
            //            cout  << "a_rect: " << a_rect->m_min[0] << " " << a_rect->m_min[1] << " " << a_rect->m_max[0] << " " << a_rect->m_max[1] << "\n";
            if (Overlap(a_rect, &a_node->m_branch[index].m_rect)) {
                //                printf("overlapped with this guy at pos: %d, %d, %d, %d\n", a_node->m_branch[index].m_rect.m_min[0], a_node->m_branch[index].m_rect.m_min[1], a_node->m_branch[index].m_rect.m_max[0], a_node->m_branch[index].m_rect.m_max[1]);
                // TODO update cv, with last_encountered parent node and branch number, DONE
                //                crackVars->potential_parent = a_node;
                //                crackVars->potential_branch_index = index;
                //                if(a_node->m_branch[index].m_child->m_count > 100) printf("this could be fucked up %d hey \n", a_node->m_branch[index].m_child->m_count);
                if (!Search(a_node->m_branch[index].m_child, a_rect, a_foundCount, crackVars)) {
                    // The callback indicated to stop searching
                    return false;
                }
            }
        }
    } else {
        //        printf("was not internal node\n");
        // This is a leaf node
        //        if (!a_node->isRegular()) {
        ////            printf("was not regular node\n");
        ////            crackVars->parent_of_last_irregular_branch_index = crackVars->potential_branch_index;
        ////            crackVars->last_irregular_parent = crackVars->potential_parent;
        ////            crackVars->last_irregular = a_node;
        ////            crackVars->irregulars_encoutered_count = crackVars->irregulars_encoutered_count + 1;
        //        } else {
        //            crackVars->encountered_regular = true;
        //        }

        if(a_node->isRegular()){
            crackVars->encountered_regular = true;
        }

        std::vector<int> to_be_deleted_indexes;
        std::vector<Branch> to_be_deleted_ids;
        // TODO
        // we have to create the fragments here
        // one fragment is the query frag which is in the tbd vectors
        // the rest we have to create here
        fragment fs[9];
        int count_of_peices = 9;
        //        cout << "creating fragments: " << fs.size() << endl;
        // TODO make this reasonable, this is probably the worst way to do this:/

        fs[0] = (std::make_pair(std::vector<Branch>(), std::make_pair(std::make_pair(0, 0), Rect())));
        fs[1] = (std::make_pair(std::vector<Branch>(), std::make_pair(std::make_pair(0, 1), Rect())));
        fs[2] = (std::make_pair(std::vector<Branch>(), std::make_pair(std::make_pair(0, 2), Rect())));
        fs[3] = (std::make_pair(std::vector<Branch>(), std::make_pair(std::make_pair(1, 0), Rect())));
        fs[4] = (std::make_pair(std::vector<Branch>(), std::make_pair(std::make_pair(1, 1), Rect())));
        fs[5] = (std::make_pair(std::vector<Branch>(), std::make_pair(std::make_pair(1, 2), Rect())));
        fs[6] = (std::make_pair(std::vector<Branch>(), std::make_pair(std::make_pair(2, 0), Rect())));
        fs[7] = (std::make_pair(std::vector<Branch>(), std::make_pair(std::make_pair(2, 1), Rect())));
        fs[8] = (std::make_pair(std::vector<Branch>(), std::make_pair(std::make_pair(2, 2), Rect())));
        //        auto start_time = std::chrono::steady_clock::now();
        //        std::chrono::duration<double> insertion_dur = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);
        //        auto insertion_time = insertion_dur.count();
        //
        //        std::chrono::duration<double> cover_dur = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);
        //        auto cover_time = cover_dur.count();
        //        cout << "creating fragments: " << fs.size() << endl;
        // TODO make this reasonable, this is probably the worst way to do this:/
        //        auto start_time = std::chrono::steady_clock::now();
        if(a_node->isRegular()) {
            for (int index = 0; index < a_node->m_count; ++index) {
                //            printf("looking through the children\n");
                // TODO I still think that even here we should look at the NodeCover instead of the saved Rect DONE, NOT NEEDED
                // I was wrong, the children here are data points, so the rect is correct no need to update to the NodeCover
                //            cout  << "a_rect: " << a_rect->m_min[0] << " " << a_rect->m_min[1] << " " << a_rect->m_max[0] << " " << a_rect->m_max[1] << "\n";
                if (Overlap(a_rect, &a_node->m_branch[index].m_rect)) {
                    //                printf("overlapped with this guy at pos: %d, %d, %d, %d its regularity: %d\n", a_node->m_branch[index].m_rect.m_min[0], a_node->m_branch[index].m_rect.m_min[1], a_node->m_branch[index].m_rect.m_max[0], a_node->m_branch[index].m_rect.m_max[1], a_node->m_regular);
                    DATATYPE &id = a_node->m_branch[index].m_data;
                    ++a_foundCount;

                }
            }
        } else{
            for (int index = 0; index < a_node->m_count; ++index) {
                //            printf("looking through the children\n");
                // TODO I still think that even here we should look at the NodeCover instead of the saved Rect DONE, NOT NEEDED
                // I was wrong, the children here are data points, so the rect is correct no need to update to the NodeCover
                //            cout  << "a_rect: " << a_rect->m_min[0] << " " << a_rect->m_min[1] << " " << a_rect->m_max[0] << " " << a_rect->m_max[1] << "\n";
                if (Overlap(a_rect, &a_node->m_branch[index].m_rect)) {
                    DATATYPE &id = a_node->m_branch[index].m_data;
                    ++a_foundCount;
                    //TODO, DONE
                    // if it is irregular, update crack_vars
                    // and remove this branch from node

                    to_be_deleted_indexes.push_back(index);
                    to_be_deleted_ids.push_back(a_node->m_branch[index]);
                    fs[4].first.push_back(a_node->m_branch[index]);
                    if(fs[4].first.size() == 1){
                        // first one
                        fs[4].second.second = a_node->m_branch[index].m_rect;
                        //                        cout << "THE NEXT TWO LINES SHOULD BE IDENTICAL\n";
                        //                        cout << fs[4].second.second.m_min[0] << " " << fs[4].second.second.m_min[1] << " " << fs[4].second.second.m_max[0] << " " << fs[4].second.second.m_max[1] << endl;
                        //                        cout << a_node->m_branch[index].m_rect.m_min[0] << " " << a_node->m_branch[index].m_rect.m_min[1] << " " << a_node->m_branch[index].m_rect.m_max[0] << " " << a_node->m_branch[index].m_rect.m_max[1] << endl;

                    } else{
                        // combine them
                        fs[4].second.second = CombineRect(&((fs[4].second).second), &(a_node->m_branch[index].m_rect));
                    }
                } else {
                    // we have to fill out frags here
                    if(a_rect->m_max[1] < a_node->m_branch[index].m_rect.m_min[1]){
                        // three boxes: 00, 01, 02
                        if(a_rect->m_min[0] > a_node->m_branch[index].m_rect.m_max[0]){
                            //                            start_time = std::chrono::steady_clock::now();
                            fs[0].first.push_back(a_node->m_branch[index]);
                            //                            insertion_dur = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);
                            //                            insertion_time += insertion_dur.count();
                            //                            start_time = std::chrono::steady_clock::now();
                            if(fs[0].first.size() == 1){
                                // first one
                                fs[0].second.second = a_node->m_branch[index].m_rect;
                            } else{
                                // combine them
                                fs[0].second.second = CombineRect(&((fs[0].second).second), &(a_node->m_branch[index].m_rect));
                            }
                            //                            cover_dur = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);
                            //                            cover_time += cover_dur.count();
                        }
                        else if(a_rect->m_max[0] > a_node->m_branch[index].m_rect.m_max[0]){
                            fs[1].first.push_back(a_node->m_branch[index]);
                            if(fs[1].first.size() == 1){
                                // first one
                                fs[1].second.second = a_node->m_branch[index].m_rect;
                            } else{
                                // combine them
                                fs[1].second.second = CombineRect(&((fs[1].second).second), &(a_node->m_branch[index].m_rect));
                            }
                        }
                        else{
                            fs[2].first.push_back(a_node->m_branch[index]);
                            if(fs[2].first.size() == 1){
                                // first one
                                fs[2].second.second = a_node->m_branch[index].m_rect;
                            } else{
                                // combine them
                                fs[2].second.second = CombineRect(&((fs[2].second).second), &(a_node->m_branch[index].m_rect));
                            }
                        }
                    }
                    else if(a_rect->m_min[1] < a_node->m_branch[index].m_rect.m_min[1]){
                        // three boxes: 10, 11, 12
                        if(a_rect->m_min[0] > a_node->m_branch[index].m_rect.m_max[0]){
                            fs[3].first.push_back(a_node->m_branch[index]);
                            if(fs[3].first.size() == 1){
                                // first one
                                fs[3].second.second = a_node->m_branch[index].m_rect;
                            } else{
                                // combine them
                                fs[3].second.second = CombineRect(&((fs[3].second).second), &(a_node->m_branch[index].m_rect));
                            }
                        }
                        else if(a_rect->m_max[0] > a_node->m_branch[index].m_rect.m_max[0]){
                            // this should be empty as it is the query bounds!
                            // wel not empty, it should just be filled out in the if statement and not in else
                            fs[4].first.push_back(a_node->m_branch[index]);
                            if(fs[4].first.size() == 1){
                                // first one
                                fs[4].second.second = a_node->m_branch[index].m_rect;
                            } else{
                                // combine them
                                fs[4].second.second = CombineRect(&((fs[4].second).second), &(a_node->m_branch[index].m_rect));
                            }
                        }
                        else{
                            fs[5].first.push_back(a_node->m_branch[index]);
                            if(fs[5].first.size() == 1){
                                // first one
                                fs[5].second.second = a_node->m_branch[index].m_rect;
                            } else{
                                // combine them
                                fs[5].second.second = CombineRect(&((fs[5].second).second), &(a_node->m_branch[index].m_rect));
                            }
                        }
                    }
                    else{
                        // three boxes: 20, 21, 22
                        if(a_rect->m_min[0] > a_node->m_branch[index].m_rect.m_max[0]){
                            fs[6].first.push_back(a_node->m_branch[index]);
                            if(fs[6].first.size() == 1){
                                // first one
                                fs[6].second.second = a_node->m_branch[index].m_rect;
                            } else{
                                // combine them
                                fs[6].second.second = CombineRect(&((fs[6].second).second), &(a_node->m_branch[index].m_rect));
                            }
                        }
                        else if(a_rect->m_max[0] > a_node->m_branch[index].m_rect.m_max[0]){
                            fs[7].first.push_back(a_node->m_branch[index]);
                            if(fs[7].first.size() == 1){
                                // first one
                                fs[7].second.second = a_node->m_branch[index].m_rect;
                            } else{
                                // combine them
                                fs[7].second.second = CombineRect(&((fs[7].second).second), &(a_node->m_branch[index].m_rect));
                            }
                        }
                        else{
                            fs[8].first.push_back(a_node->m_branch[index]);
                            if(fs[8].first.size() == 1){
                                // first one
                                fs[8].second.second = a_node->m_branch[index].m_rect;
                            } else{
                                // combine them
                                fs[8].second.second = CombineRect(&((fs[8].second).second), &(a_node->m_branch[index].m_rect));
                            }
                        }
                    }
                }

            }
        }
        tbd_vectors folan;
        if(!a_node->isRegular() &&  !to_be_deleted_indexes.empty()){
            //            ASSERT(to_be_deleted_indexes.size() == to_be_deleted_ids.size());
            folan = std::make_pair(to_be_deleted_ids, to_be_deleted_indexes);
            //            cout << "creating fragments after: " << fs.size() << endl;
            crackVars->gathered_nodes.push_back(std::make_pair(a_node, std::make_pair(folan, fs)));
        }

        //        cout << "IN SEARCH, insertion times: " << insertion_time << " cover_time: " << cover_time <<  endl;

        //        std::chrono::duration<double> res = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);
        //        cout << res.count() << " search through data points in leaf" << endl;
        //        printf("after searching through the children deleting the t_b_d list\n");
        //         TODO this is wrong, the indexes change after any disconnect:(( FIXED
        // but since they are sorted, we can do a trick
        //        int cnt = 0; int correct_index = 0;
        //        cout << "BLOB" << a_node->m_count << "\n";
        //        start_time = std::chrono::steady_clock::now();

        // FATEMEH
        // this should not be done here, it should be done after the search is complete
        //        if(!a_node->isRegular()){
        //            cout << "ROOT ID BEFORE: " << m_root->m_id << "\n";
        //            bool a_node_was_deleted = RemoveBatchFromNode(a_node, to_be_deleted_ids, to_be_deleted_indexes, &m_root);
        //
        //            cout << "ROOT ID AFTER: " << m_root->m_id << "\n";
        //            if (a_node_was_deleted){
        //                int count = to_be_deleted_ids.size();
        //                crackVars->data_point_parents.erase(crackVars->data_point_parents.end() - count, crackVars->data_point_parents.end());
        //                for (int i = 0; i < count; i++){
        //                    crackVars->data_point_parents.push_back(NULL);
        //                }
        //                ASSERT(crackVars->data_points.size() == crackVars->data_point_parents.size());
        //            }
        //        }
        // FATEMEH OUT


        //        res = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);
        //        cout << res.count() << " removing data points from irregulars!" << endl;
        //        printf("HEY %d hey\n", a_node->m_count);
        //        for (int i = (to_be_deleted_indexes.size() - 1); i > -1; i--) {
        //            cout << "in deleting list in search" << a_node->m_branch[to_be_deleted_indexes[i]].m_rect.m_min[0] << "\n";
        //            RemoveRect(&(a_node->m_branch[to_be_deleted_indexes[i]].m_rect),
        //                       a_node->m_branch[to_be_deleted_indexes[i]].m_data, &m_root);
        //        }

    }

    return true; // Continue searching
}



// Search in an index tree or subtree for all data retangles that overlap the argument rectangle.
RTREE_TEMPLATE
bool RTREE_QUAL::Search_2(Node *a_node, Rect *a_rect, int &a_foundCount, CrackVars *crackVars) {
    //            ASSERT(a_node);
    //            ASSERT(a_node->m_level >= 0);
    //            ASSERT(a_rect);
    //            ASSERT(a_node->m_count == a_node->m_branch.size());
    printf("in search\n");
    if (a_node->IsInternalNode()) {
        //        printf("was inetrnal node\n");
        // This is an internal node in the tree
        for (int index = 0; index < a_node->m_count; ++index) {
            // TODO I think here we have to compare with the NodeCover of the child, since the calculation might be old
            // this can happen because of the shrinks in the search
            //            if (Overlap(a_rect, &a_node->m_branch[index].m_rect))
            //            printf("calling NC in Search\n");
            //            if(a_node->m_branch[index].m_child->m_branch == NULL) {printf("HEY THERE !!\n");}
            //            printf("I THINK IT'S HERE\n");
            //            if(a_node->m_branch[index].m_child->m_branch.empty()){printf("HERE HERE HERE count = %d \n", a_node->m_branch[index].m_child->m_count);}
            //            a_node->m_branch[index].m_rect = NodeCover(a_node->m_branch[index].m_child);
            //            cout  << "a_rect: " << a_rect->m_min[0] << " " << a_rect->m_min[1] << " " << a_rect->m_max[0] << " " << a_rect->m_max[1] << "\n";
            if (Overlap(a_rect, &a_node->m_branch[index].m_rect)) {
                //                printf("overlapped with this guy at pos: %d, %d, %d, %d\n", a_node->m_branch[index].m_rect.m_min[0], a_node->m_branch[index].m_rect.m_min[1], a_node->m_branch[index].m_rect.m_max[0], a_node->m_branch[index].m_rect.m_max[1]);
                // TODO update cv, with last_encountered parent node and branch number, DONE
                //                crackVars->potential_parent = a_node;
                //                crackVars->potential_branch_index = index;
                //                if(a_node->m_branch[index].m_child->m_count > 100) printf("this could be fucked up %d hey \n", a_node->m_branch[index].m_child->m_count);
                if (!Search_2(a_node->m_branch[index].m_child, a_rect, a_foundCount, crackVars)) {
                    // The callback indicated to stop searching
                    return false;
                }
            }
        }
    }
    else {
        //        printf("was not internal node\n");
        // This is a leaf node
        //        if (!a_node->isRegular()) {
        ////            printf("was not regular node\n");
        ////            crackVars->parent_of_last_irregular_branch_index = crackVars->potential_branch_index;
        ////            crackVars->last_irregular_parent = crackVars->potential_parent;
        ////            crackVars->last_irregular = a_node;
        ////            crackVars->irregulars_encoutered_count = crackVars->irregulars_encoutered_count + 1;
        //        } else {
        //            crackVars->encountered_regular = true;
        //        }

        if(a_node->isRegular()){
            crackVars->encountered_regular = true;
        }

        std::vector<int> to_be_deleted_indexes;
        std::vector<Branch> to_be_deleted_ids;
        // TODO
        // we have to create the fragments here
        // one fragment is the query frag which is in the tbd vectors
        // the rest we have to create here
        fragment fs[9];
        int count_of_peices = 9;
        //        cout << "creating fragments: " << fs.size() << endl;
        // TODO make this reasonable, this is probably the worst way to do this:/

        fs[0] = (std::make_pair(std::vector<Branch>(), std::make_pair(std::make_pair(0, 0), Rect())));
        fs[1] = (std::make_pair(std::vector<Branch>(), std::make_pair(std::make_pair(0, 1), Rect())));
        fs[2] = (std::make_pair(std::vector<Branch>(), std::make_pair(std::make_pair(0, 2), Rect())));
        fs[3] = (std::make_pair(std::vector<Branch>(), std::make_pair(std::make_pair(1, 0), Rect())));
        fs[4] = (std::make_pair(std::vector<Branch>(), std::make_pair(std::make_pair(1, 1), Rect())));
        fs[5] = (std::make_pair(std::vector<Branch>(), std::make_pair(std::make_pair(1, 2), Rect())));
        fs[6] = (std::make_pair(std::vector<Branch>(), std::make_pair(std::make_pair(2, 0), Rect())));
        fs[7] = (std::make_pair(std::vector<Branch>(), std::make_pair(std::make_pair(2, 1), Rect())));
        fs[8] = (std::make_pair(std::vector<Branch>(), std::make_pair(std::make_pair(2, 2), Rect())));
        //        auto start_time = std::chrono::steady_clock::now();
        //        std::chrono::duration<double> insertion_dur = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);
        //        auto insertion_time = insertion_dur.count();
        //
        //        std::chrono::duration<double> cover_dur = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);
        //        auto cover_time = cover_dur.count();
        //        cout << "creating fragments: " << fs.size() << endl;
        // TODO make this reasonable, this is probably the worst way to do this:/
        //        auto start_time = std::chrono::steady_clock::now();
        if(a_node->isRegular()) {
            for (int index = 0; index < a_node->m_count; ++index) {
                //            printf("looking through the children\n");
                // TODO I still think that even here we should look at the NodeCover instead of the saved Rect DONE, NOT NEEDED
                // I was wrong, the children here are data points, so the rect is correct no need to update to the NodeCover
                //            cout  << "a_rect: " << a_rect->m_min[0] << " " << a_rect->m_min[1] << " " << a_rect->m_max[0] << " " << a_rect->m_max[1] << "\n";
                if (Overlap(a_rect, &a_node->m_branch[index].m_rect)) {
                    //                printf("overlapped with this guy at pos: %d, %d, %d, %d its regularity: %d\n", a_node->m_branch[index].m_rect.m_min[0], a_node->m_branch[index].m_rect.m_min[1], a_node->m_branch[index].m_rect.m_max[0], a_node->m_branch[index].m_rect.m_max[1], a_node->m_regular);
                    DATATYPE &id = a_node->m_branch[index].m_data;
                    ++a_foundCount;

                }
            }
        } else{
            for (int index = 0; index < a_node->m_count; ++index) {
                //            printf("looking through the children\n");
                // TODO I still think that even here we should look at the NodeCover instead of the saved Rect DONE, NOT NEEDED
                // I was wrong, the children here are data points, so the rect is correct no need to update to the NodeCover
                //            cout  << "a_rect: " << a_rect->m_min[0] << " " << a_rect->m_min[1] << " " << a_rect->m_max[0] << " " << a_rect->m_max[1] << "\n";
                if (Overlap(a_rect, &a_node->m_branch[index].m_rect)) {
                    DATATYPE &id = a_node->m_branch[index].m_data;
                    ++a_foundCount;
                    //TODO, DONE
                    // if it is irregular, update crack_vars
                    // and remove this branch from node

                    to_be_deleted_indexes.push_back(index);
                    to_be_deleted_ids.push_back(a_node->m_branch[index]);
                    fs[4].first.push_back(a_node->m_branch[index]);
                    if(fs[4].first.size() == 1){
                        // first one
                        fs[4].second.second = a_node->m_branch[index].m_rect;
                        //                        cout << "THE NEXT TWO LINES SHOULD BE IDENTICAL\n";
                        //                        cout << fs[4].second.second.m_min[0] << " " << fs[4].second.second.m_min[1] << " " << fs[4].second.second.m_max[0] << " " << fs[4].second.second.m_max[1] << endl;
                        //                        cout << a_node->m_branch[index].m_rect.m_min[0] << " " << a_node->m_branch[index].m_rect.m_min[1] << " " << a_node->m_branch[index].m_rect.m_max[0] << " " << a_node->m_branch[index].m_rect.m_max[1] << endl;

                    } else{
                        // combine them
                        fs[4].second.second = CombineRect(&((fs[4].second).second), &(a_node->m_branch[index].m_rect));
                    }
                } else {
                    // we have to fill out frags here
                    if(a_rect->m_max[1] < a_node->m_branch[index].m_rect.m_min[1]){
                        // three boxes: 00, 01, 02
                        //                        if(a_rect->m_min[0] > a_node->m_branch[index].m_rect.m_max[0]){
                        //                            start_time = std::chrono::steady_clock::now();
                        fs[0].first.push_back(a_node->m_branch[index]);
                        //                            insertion_dur = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);
                        //                            insertion_time += insertion_dur.count();
                        //                            start_time = std::chrono::steady_clock::now();
                        if(fs[0].first.size() == 1){
                            // first one
                            fs[0].second.second = a_node->m_branch[index].m_rect;
                        } else{
                            // combine them
                            fs[0].second.second = CombineRect(&((fs[0].second).second), &(a_node->m_branch[index].m_rect));
                        }
                        //                            cover_dur = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);
                        //                            cover_time += cover_dur.count();
                        //                        }
                        //                        else if(a_rect->m_max[0] > a_node->m_branch[index].m_rect.m_max[0]){
                        //                            fs[1].first.push_back(a_node->m_branch[index]);
                        //                            if(fs[1].first.size() == 1){
                        //                                // first one
                        //                                fs[1].second.second = a_node->m_branch[index].m_rect;
                        //                            } else{
                        //                                // combine them
                        //                                fs[1].second.second = CombineRect(&((fs[1].second).second), &(a_node->m_branch[index].m_rect));
                        //                            }
                        //                        }
                        //                        else{
                        //                            fs[2].first.push_back(a_node->m_branch[index]);
                        //                            if(fs[2].first.size() == 1){
                        //                                // first one
                        //                                fs[2].second.second = a_node->m_branch[index].m_rect;
                        //                            } else{
                        //                                // combine them
                        //                                fs[2].second.second = CombineRect(&((fs[2].second).second), &(a_node->m_branch[index].m_rect));
                        //                            }
                        //                        }
                    }
                    else if(a_rect->m_min[1] < a_node->m_branch[index].m_rect.m_min[1]){
                        // three boxes: 10, 11, 12
                        if(a_rect->m_min[0] > a_node->m_branch[index].m_rect.m_max[0]){
                            fs[3].first.push_back(a_node->m_branch[index]);
                            if(fs[3].first.size() == 1){
                                // first one
                                fs[3].second.second = a_node->m_branch[index].m_rect;
                            } else{
                                // combine them
                                fs[3].second.second = CombineRect(&((fs[3].second).second), &(a_node->m_branch[index].m_rect));
                            }
                        }
                        //                        else if(a_rect->m_max[0] > a_node->m_branch[index].m_rect.m_max[0]){
                        //                            // this should be empty as it is the query bounds!
                        //                            // wel not empty, it should just be filled out in the if statement and not in else
                        //                            fs[4].first.push_back(a_node->m_branch[index]);
                        //                            if(fs[4].first.size() == 1){
                        //                                // first one
                        //                                fs[4].second.second = a_node->m_branch[index].m_rect;
                        //                            } else{
                        //                                // combine them
                        //                                fs[4].second.second = CombineRect(&((fs[4].second).second), &(a_node->m_branch[index].m_rect));
                        //                            }
                        //                        }
                        else{
                            fs[5].first.push_back(a_node->m_branch[index]);
                            if(fs[5].first.size() == 1){
                                // first one
                                fs[5].second.second = a_node->m_branch[index].m_rect;
                            } else{
                                // combine them
                                fs[5].second.second = CombineRect(&((fs[5].second).second), &(a_node->m_branch[index].m_rect));
                            }
                        }
                    }
                    else{
                        // three boxes: 20, 21, 22
                        //                        if(a_rect->m_min[0] > a_node->m_branch[index].m_rect.m_max[0]){
                        fs[6].first.push_back(a_node->m_branch[index]);
                        if(fs[6].first.size() == 1){
                            // first one
                            fs[6].second.second = a_node->m_branch[index].m_rect;
                        } else{
                            // combine them
                            fs[6].second.second = CombineRect(&((fs[6].second).second), &(a_node->m_branch[index].m_rect));
                        }
                        //                        }
                        //                        else if(a_rect->m_max[0] > a_node->m_branch[index].m_rect.m_max[0]){
                        //                            fs[7].first.push_back(a_node->m_branch[index]);
                        //                            if(fs[7].first.size() == 1){
                        //                                // first one
                        //                                fs[7].second.second = a_node->m_branch[index].m_rect;
                        //                            } else{
                        //                                // combine them
                        //                                fs[7].second.second = CombineRect(&((fs[7].second).second), &(a_node->m_branch[index].m_rect));
                        //                            }
                        //                        }
                        //                        else{
                        //                            fs[8].first.push_back(a_node->m_branch[index]);
                        //                            if(fs[8].first.size() == 1){
                        //                                // first one
                        //                                fs[8].second.second = a_node->m_branch[index].m_rect;
                        //                            } else{
                        //                                // combine them
                        //                                fs[8].second.second = CombineRect(&((fs[8].second).second), &(a_node->m_branch[index].m_rect));
                        //                            }
                        //                        }
                    }
                }

            }
        }
        tbd_vectors folan;
        if(!a_node->isRegular() &&  !to_be_deleted_indexes.empty()){
            //                    ASSERT(to_be_deleted_indexes.size() == to_be_deleted_ids.size());
            folan = std::make_pair(to_be_deleted_ids, to_be_deleted_indexes);
            //            cout << "creating fragments after: " << fs.size() << endl;
            //cout << "adding to cv\n";
            crackVars->gathered_nodes.push_back(std::make_pair(a_node, std::make_pair(folan, fs)));
        }

        //        cout << "IN SEARCH, insertion times: " << insertion_time << " cover_time: " << cover_time <<  endl;

        //        std::chrono::duration<double> res = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);
        //        cout << res.count() << " search through data points in leaf" << endl;
        //        printf("after searching through the children deleting the t_b_d list\n");
        //         TODO this is wrong, the indexes change after any disconnect:(( FIXED
        // but since they are sorted, we can do a trick
        //        int cnt = 0; int correct_index = 0;
        //        cout << "BLOB" << a_node->m_count << "\n";
        //        start_time = std::chrono::steady_clock::now();

        // FATEMEH
        // this should not be done here, it should be done after the search is complete
        //        if(!a_node->isRegular()){
        //            cout << "ROOT ID BEFORE: " << m_root->m_id << "\n";
        //            bool a_node_was_deleted = RemoveBatchFromNode(a_node, to_be_deleted_ids, to_be_deleted_indexes, &m_root);
        //
        //            cout << "ROOT ID AFTER: " << m_root->m_id << "\n";
        //            if (a_node_was_deleted){
        //                int count = to_be_deleted_ids.size();
        //                crackVars->data_point_parents.erase(crackVars->data_point_parents.end() - count, crackVars->data_point_parents.end());
        //                for (int i = 0; i < count; i++){
        //                    crackVars->data_point_parents.push_back(NULL);
        //                }
        //                ASSERT(crackVars->data_points.size() == crackVars->data_point_parents.size());
        //            }
        //        }
        // FATEMEH OUT


        //        res = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - start_time);
        //        cout << res.count() << " removing data points from irregulars!" << endl;
        //        printf("HEY %d hey\n", a_node->m_count);
        //        for (int i = (to_be_deleted_indexes.size() - 1); i > -1; i--) {
        //            cout << "in deleting list in search" << a_node->m_branch[to_be_deleted_indexes[i]].m_rect.m_min[0] << "\n";
        //            RemoveRect(&(a_node->m_branch[to_be_deleted_indexes[i]].m_rect),
        //                       a_node->m_branch[to_be_deleted_indexes[i]].m_data, &m_root);
        //        }

    }

    return true; // Continue searching
}



// Search in an index tree or subtree for all data retangles that overlap the argument rectangle.
RTREE_TEMPLATE
bool RTREE_QUAL::Search_old(Node *a_node, Rect *a_rect, int &a_foundCount, CrackVars *crackVars) {
    //            ASSERT(a_node);
    //            ASSERT(a_node->m_level >= 0);
    //            ASSERT(a_rect);
    //            ASSERT(a_node->m_count == a_node->m_branch.size());
    //    printf("in search\n");
    if (a_node->IsInternalNode()) {
        //        printf("was inetrnal node\n");
        // This is an internal node in the tree
        for (int index = 0; index < a_node->m_count; ++index) {
            if (Overlap(a_rect, &a_node->m_branch[index].m_rect)) {
                if (!Search_old(a_node->m_branch[index].m_child, a_rect, a_foundCount, crackVars)) {
                    // The callback indicated to stop searching
                    return false;
                }
            }
        }
    }
    else {
        if(a_node->isRegular()){
            crackVars->encountered_regular = true;
        }

        std::vector<int> to_be_deleted_indexes;
        std::vector<Branch> to_be_deleted_ids;

        if(a_node->isRegular()) {
            for (int index = 0; index < a_node->m_count; ++index) {
                if (Overlap(a_rect, &a_node->m_branch[index].m_rect)) {
                    DATATYPE &id = a_node->m_branch[index].m_data;
                    ++a_foundCount;
                }
            }
        } else{
            for (int index = 0; index < a_node->m_count; ++index) {
                if (Overlap(a_rect, &a_node->m_branch[index].m_rect)) {
                    DATATYPE &id = a_node->m_branch[index].m_data;
                    ++a_foundCount;
                    //TODO, DONE
                    // if it is irregular, update crack_vars
                    // and remove this branch from node

                    to_be_deleted_indexes.push_back(index);
                    to_be_deleted_ids.push_back(a_node->m_branch[index]);

                }
            }

        }
        fragments fs;
        tbd_vectors folan;
        if(!a_node->isRegular() &&  !to_be_deleted_indexes.empty()){
            //                    ASSERT(to_be_deleted_indexes.size() == to_be_deleted_ids.size());
            folan = std::make_pair(to_be_deleted_ids, to_be_deleted_indexes);
            //            cout << "creating fragments after: " << fs.size() << endl;
            crackVars->gathered_nodes.push_back(std::make_pair(a_node, std::make_pair(folan, fs)));
        }


    }

    return true; // Continue searching
}




RTREE_TEMPLATE
bool RTREE_QUAL::Search_2026(Node *a_node, Rect *a_rect, int &a_foundCount, LeavesToAdd_v2 *all_lta) {
    //            ASSERT(a_node);
    //            ASSERT(a_node->m_level >= 0);
    //            ASSERT(a_rect);
    //            ASSERT(a_node->m_count == a_node->m_branch.size());
    //        cout << "looked_through_a_node" << endl;

    if (a_node->IsInternalNode()) {
        // This is an internal node in the tree
        for (int index = 0; index < a_node->m_count; ++index) {
            if (Overlap(a_rect, &a_node->m_branch[index].m_rect)) {
                if (!Search_2026(a_node->m_branch[index].m_child, a_rect, a_foundCount, all_lta)) {
                    // The callback indicated to stop searching
                    return false;
                }
            }
        }
    }
    else {
        if(a_node->isRegular()) {
            for (int index = a_node->m_L; index < a_node->m_R; ++index) {
                if (Overlap(a_rect, m_data_arr_mins[index], m_data_arr_maxes[index])) {
                    ++a_foundCount;
                }
            }

        }
        else{
            // call elegant_2dc
            //            cout << "found count before " << a_foundCount << endl;
            //            a_foundCount += Elegant_2dc_v2(a_node, a_rect, all_lta);
            //            cout << "found count after " << a_foundCount << endl;

            Rect node_cover;
            if(a_node->m_parent != NULL) {
                node_cover = ((a_node->m_parent)->m_branch[getNodeBranchIndex(a_node)]).m_rect;
            } else{
                node_cover = root_cover;
            }


            Rect this_piece_rect = Rect(node_cover.m_min, node_cover.m_max);
            int this_piece_L = a_node->m_L;
            int this_piece_R = a_node->m_R;


            int choice, mm, chosen_axis;

            LTA_v2 this_lta;
            this_lta.how_many_created =0;
            this_lta.this_leaf = a_node;
            Rect potential_query_piece = Rect(node_cover.m_min, node_cover.m_max); 
            Rect other_piece;


            int crack_index, filan, ashghal;
            int sum_of_choices = 0; int sum_of_mms = 0;

            // 18NOV
            int largest_piece_index = 0;
            int stochastic_crack_axis = 0; float stochastic_crack_pivot;
            float X, X2, X3;
            int stochastic_crack_index;

            for(int i =0; i < 2*NUMDIMS; i++){
                //            for(int i =0; i < 4; i++){

                for(filan=0; filan < NUMDIMS; filan++) {
                    other_piece.m_min[filan] = potential_query_piece.m_min[filan];
                    other_piece.m_max[filan] = potential_query_piece.m_max[filan];
                }



                if(i == ((2 * NUMDIMS) - 1)){
//                    cout << "last choice\n";
                    choice = (((2*NUMDIMS) - 1) * NUMDIMS) - sum_of_choices;
                }
                else {
                    choice = ChooseAxisLargestSideMid(this_piece_rect, *a_rect);
                }
                // 11oct 16:26
                if(choice < 0) break;
                //
                if(sum_of_mms >= NUMDIMS) mm = 1;
                else if (sum_of_mms <= (-1*NUMDIMS)) mm = 0;
                else {
                    mm = choice % 2;
                    if(mm == 0) sum_of_mms++;
                    else sum_of_mms--;
                }

//                cout << "sum of choices "<< sum_of_choices << " choice " << choice << endl;
                chosen_axis = choice / 2;
//                cout << "chosen axis " << chosen_axis << endl;

//                ASSERT(chosen_axis >= 0);
//                ASSERT(chosen_axis < NUMDIMS);

                if(chosen_axis >= NUMDIMS || chosen_axis < 0) break;

                sum_of_choices += choice;

//                potential_query_piece.m_min[chosen_axis] = 21474836;
//                potential_query_piece.m_max[chosen_axis] = -21474836;
                other_piece.m_min[chosen_axis] = 21474836;
                other_piece.m_max[chosen_axis] = -21474836;

                if(mm == 0) {
                    potential_query_piece.m_min[chosen_axis] = 21474836;
                    crack_index = CrackOnAxisCompMax_v45( this_piece_L, this_piece_R, a_rect->m_min[chosen_axis],
                                                        chosen_axis, &other_piece, &potential_query_piece);
                }
                else{
                    potential_query_piece.m_max[chosen_axis] = -21474836;
                    crack_index = CrackOnAxisCompMin_v45( this_piece_L, this_piece_R, a_rect->m_max[chosen_axis],
                                                        chosen_axis, &potential_query_piece, &other_piece);
                }

//                cout << "INJA DARIM CHECK MIKONIM " << potential_query_piece.m_min[chosen_axis] << " " << potential_query_piece.m_max[chosen_axis] << endl;
//                if(potential_query_piece.m_min[chosen_axis] == 21474836 || potential_query_piece.m_max[chosen_axis] == -21474836) {
//                    cout << "naaaa nabayd in etefagh mioftad " << chosen_axis << " mm " << mm << " L " << this_piece_L << " crackindex " << crack_index << " R " << this_piece_R << endl;
//                }
//                if(other_piece.m_min[chosen_axis] > other_piece.m_max[chosen_axis] && other_piece.m_min[chosen_axis] != 21474836) {
//                    cout << "naaaaaaa nabayd in etefagh mioftad " << chosen_axis << " mm " << mm << " L " << this_piece_L << " crackindex " << crack_index << " R " << this_piece_R << endl;
//                    cout << "values " << other_piece.m_min[chosen_axis] << " " << other_piece.m_max[chosen_axis] << endl;
//                }


                if((mm == 0 && crack_index != this_piece_L) ||
                (mm == 1 && this_piece_R != crack_index)){

                    if(mm == 0){

                        this_lta.Ls[this_lta.how_many_created] = this_piece_L;
                        this_lta.Rs[this_lta.how_many_created] = crack_index;
                        // set up ptcrack
                        this_piece_L = crack_index;
                    }
                    else{
                        this_lta.Ls[this_lta.how_many_created] = crack_index;
                        this_lta.Rs[this_lta.how_many_created] = this_piece_R;
                        // set up ptcrack
                        this_piece_R = crack_index;
                    }
                    this_piece_rect = Rect(potential_query_piece.m_min, potential_query_piece.m_max);


                    for(ashghal = 0; ashghal < NUMDIMS; ashghal ++){
//                        if(other_piece.m_min[ashghal]  > other_piece.m_max[ashghal]){
//                            cout << "IN SEARCH2026, ADDING OTHERPIECE, MIN>MAX dim " <<  ashghal << " chosen_axis " << chosen_axis  << " mm " << mm << " L " <<  this_lta.Ls[this_lta.how_many_created] << " R " << this_lta.Rs[this_lta.how_many_created] << endl;
//                            cout << "values  min " << other_piece.m_min[ashghal] <<  " max " << other_piece.m_max[ashghal] << endl;
//                        }

                        this_lta.crack_covers_mins[this_lta.how_many_created][ashghal] = other_piece.m_min[ashghal];
                        this_lta.crack_covers_maxes[this_lta.how_many_created][ashghal] = other_piece.m_max[ashghal];
                    }

                    this_lta.how_many_created++;
                }
                else{
                    // 14 oct
                    // baghiash bi ahamiate chon inja mitoonim break konim
                    // baresh dashtam dobare Panos movafegh nabood
//                    break;

                    // in baraye ine k oon right va left cover k partition midim ye meghdar avaliyeE dashte bashan
                    // in joori b 2 tashoon potential query midim

                    // FK KONAM ALREADY MOSAVIAN VALI MOTMAEN NISTAM
                    // vali hala k in taghyir ro dadim doroste k ghablesh ba meghdar ghabli set kardim
                    potential_query_piece = Rect(this_piece_rect.m_min, this_piece_rect.m_max);
                    // HALA BARAYE CHOICE BADI THIS_PIECE_RECT RO AVAZ MIKONIM

                    if( mm == 0 ) {this_piece_rect.m_min[chosen_axis] = a_rect->m_min[chosen_axis]; }
                    else {this_piece_rect.m_max[chosen_axis] = a_rect->m_max[chosen_axis];}

                }



                if((this_piece_R - this_piece_L) <= MAXDATAPOINTS) {
                    break;
                }

            }

            // scan last piece
            // age mikhay result baghean adadesh dorost bashe bayad injoori scan koni
            // vali injoori ye adad k hamishe bozorgtar mosavi javab vagheEye b dast miad
            // DEBUG
            for(int folan = this_piece_L; folan < this_piece_R; folan++){
            
                if (Overlap(a_rect, m_data_arr_mins[folan], m_data_arr_maxes[folan])) {
                    ++a_foundCount;
                }
            }
            // a_foundCount += (this_piece_R - this_piece_L);
            // DEBUG

            if(this_piece_L != this_piece_R) {
                this_lta.Ls[this_lta.how_many_created] = this_piece_L;
                this_lta.Rs[this_lta.how_many_created] = this_piece_R;
                for(ashghal = 0; ashghal < NUMDIMS; ashghal ++){

//                    if(this_piece_rect.m_min[ashghal]  > this_piece_rect.m_max[ashghal]){
//                        cout << "IN SEARCH2026, ADDING LAST PIECE, MIN>MAX dim " << ashghal  << " mm " << mm << " L " <<  this_lta.Ls[this_lta.how_many_created] << " R " << this_lta.Rs[this_lta.how_many_created] << endl;
//                        cout << "values  min " << this_piece_rect.m_min[ashghal] <<  " max " << this_piece_rect.m_max[ashghal] << endl;
//                    }

                    this_lta.crack_covers_mins[this_lta.how_many_created][ashghal] = this_piece_rect.m_min[ashghal];
                    this_lta.crack_covers_maxes[this_lta.how_many_created][ashghal] = this_piece_rect.m_max[ashghal];
                }
                this_lta.how_many_created++;
            }
            // aval bayad bozorgtarin tikke ro peyda konim
            largest_piece_index = 0;
//            cout << "0: L " << this_lta.Ls[0] << " R " << this_lta.Rs[0] << endl;
            for(ashghal = 1; ashghal < this_lta.how_many_created; ashghal++){
//                cout << ashghal << ": L " << this_lta.Ls[ashghal] << " R " << this_lta.Rs[ashghal] << endl;
                if((this_lta.Rs[ashghal] - this_lta.Ls[ashghal]) > (this_lta.Rs[largest_piece_index] - this_lta.Ls[largest_piece_index])) largest_piece_index = ashghal;
            }
//            cout << "the chosen largest piece: " << largest_piece_index << endl;
            // hala bayad in tike ro beshkanim roo y jaye random
            stochastic_crack_axis = 0;
            // faghat bar 2d emtehan mikonim pas hamin basse
            // if((this_lta.crack_covers_maxes[largest_piece_index][1] - this_lta.crack_covers_mins[largest_piece_index][1])  > (this_lta.crack_covers_maxes[largest_piece_index][0] - this_lta.crack_covers_mins[largest_piece_index][0]) ){
            //     stochastic_crack_axis = 1;_arr_mins[this_lta.Ls[largest_piece_index] + rand() % (this_lta.Rs[largest_piece_index] - this_lta.Ls[largest_piece_index])]
            // }
            for(ashghal = 1; ashghal < NUMDIMS; ashghal++){
                if((this_lta.crack_covers_maxes[largest_piece_index][ashghal] - this_lta.crack_covers_mins[largest_piece_index][ashghal])  > (this_lta.crack_covers_maxes[largest_piece_index][stochastic_crack_axis] - this_lta.crack_covers_mins[largest_piece_index][stochastic_crack_axis]) ){
                    stochastic_crack_axis = ashghal;
                }
            }

            X = m_data_arr_mins[this_lta.Ls[largest_piece_index] + rand() % (this_lta.Rs[largest_piece_index] - this_lta.Ls[largest_piece_index])][stochastic_crack_axis];    // split in random position
            X2 = m_data_arr_mins[this_lta.Ls[largest_piece_index] + rand() % (this_lta.Rs[largest_piece_index] - this_lta.Ls[largest_piece_index])][stochastic_crack_axis];    // split in random position
            X3 = m_data_arr_mins[this_lta.Ls[largest_piece_index] + rand() % (this_lta.Rs[largest_piece_index] - this_lta.Ls[largest_piece_index])][stochastic_crack_axis];    // split in random position

           // make X the median of the samples
           if (X2 > X) swap(X, X2);
           if (X > X3) swap(X, X3);
           if (X2 > X) swap(X, X2);
           stochastic_crack_pivot = X;

           // hala bayad vaghean crack konim
           Rect l_rect = Rect(this_lta.crack_covers_mins[largest_piece_index], this_lta.crack_covers_maxes[largest_piece_index]);
           Rect r_rect = Rect(this_lta.crack_covers_mins[largest_piece_index], this_lta.crack_covers_maxes[largest_piece_index]);
           l_rect.m_min[stochastic_crack_axis] = 21474836;
           l_rect.m_max[stochastic_crack_axis] = -21474836;
           r_rect.m_min[stochastic_crack_axis] = 21474836;

           stochastic_crack_index = CrackOnAxisCompMax_v45( this_lta.Ls[largest_piece_index], this_lta.Rs[largest_piece_index], stochastic_crack_pivot,
                                                            stochastic_crack_axis, &l_rect, &r_rect);
            if(stochastic_crack_index != this_lta.Ls[largest_piece_index] && stochastic_crack_index != this_lta.Rs[largest_piece_index]){
//               cout << "creating the pieces\n";

               // right ro ezafe kon
               this_lta.Ls[this_lta.how_many_created] = stochastic_crack_index;
               this_lta.Rs[this_lta.how_many_created] = this_lta.Rs[largest_piece_index];
               for(ashghal = 0; ashghal < NUMDIMS; ashghal ++){
                   this_lta.crack_covers_mins[this_lta.how_many_created][ashghal] = r_rect.m_min[ashghal];
                   this_lta.crack_covers_maxes[this_lta.how_many_created][ashghal] = r_rect.m_max[ashghal];
               }
               this_lta.how_many_created++;


               // left ro bezar sar e jay e ooni k shekoondim
               this_lta.Rs[largest_piece_index] = stochastic_crack_index;
               for(ashghal = 0; ashghal < NUMDIMS; ashghal ++){
                   this_lta.crack_covers_mins[largest_piece_index][ashghal] = l_rect.m_min[ashghal];
                   this_lta.crack_covers_maxes[largest_piece_index][ashghal] = l_rect.m_max[ashghal];
               }

           }


            all_lta->push_back(this_lta);
        }
    }
    return true; // Continue searching
}


RTREE_TEMPLATE
void RTREE_QUAL::TraverseDownTreeTilLeaf(Node *a_node, int tba_start, int tba_count, vector<overlapping_leaf_tbas> *overlapping_leaves){
    // cout << "in traverse down till leaf" << endl;
    
    if(!a_node->IsInternalNode()){
        // cout << "was leaf" << endl;

        // for(int i = 0; i < tba->size(); i++){
        //     if((tba->at(i)).m_data == 104841){
        //         cout << "was assigned to leaf " << a_node->m_id << " in tdttl" << endl;
        //     }
        // }

        // is leaf, search is over
        // overlapping_leaves->push_back(overlapping_leaf_tbas(a_node, *tba));
        // if(a_node->m_id == 2026){
        //     cout << "in traverse, putting 2026 in ol, start: " << tba_start << endl;
        // }
        overlapping_leaves->push_back(overlapping_leaf_tbas(a_node, tba_start, tba_count));
        // cout << "after push back" << endl;
        if(a_node->m_parent != NULL) {
            // cout << "was not root tbacnt " << tba_count << endl;
            int a_node_branch_index = getNodeBranchIndex(a_node);
            Rect trash_rect;
            for(int i =tba_start; i < tba_start + tba_count; i++){
                trash_rect = (Rect(tba_mins[i], tba_maxes[i]));
                (a_node->m_parent)->m_branch[a_node_branch_index].m_rect = CombineRect(&trash_rect, &((a_node->m_parent)->m_branch[a_node_branch_index].m_rect));
            }
        }
        // else{
        //     for(int i =0; i < tba->size(); i++) root_cover = CombineRect(&root_cover, &tba->at(i).m_rect);
        // }
        
        // for(int j = 0; j < tba_count; j++)
        // {
        //     free(tba_mins[j]);
        //     free(tba_maxes[j]);
        // }
        // free(tba_ids);
        // free(tba_mins);
        // free(tba_maxes);
        


        return;
    }
    // cout << "was inetranl " << a_node->m_level << endl;
    // is internal
    // vector<Branch> tba_to_children[a_node->m_count];
    // int a_node_branch_index;
    // if(a_node->m_parent != NULL) a_node_branch_index = getNodeBranchIndex(a_node);
    // int branch_choice;
    // for(int i =0; i < tba->size(); i++){
    //     branch_choice = PickBranch(&tba->at(i).m_rect, a_node);
    //     tba_to_children[branch_choice].push_back(tba->at(i));
    //     // LET'S EXTEND THE MBR OF THE INNER NODES HERE
    //     // this node's mbr would get extended
    //     // i'll do it stupidly for now
    //     // TODO there must be a better way
    //     if(a_node->m_parent != NULL) {
    //         (a_node->m_parent)->m_branch[a_node_branch_index].m_rect = CombineRect(&tba->at(i).m_rect, &((a_node->m_parent)->m_branch[a_node_branch_index].m_rect));
    //     }
    //     else root_cover = CombineRect(&root_cover, &tba->at(i).m_rect);
    // }

    int *tba_to_children_counts = (int *)calloc(a_node->m_count, sizeof(int));
    int *tba_map = (int *)calloc(tba_count,sizeof(int));


    // cout << "after mallocs tba_count " << tba_count << endl;

    int a_node_branch_index;
    if(a_node->m_parent != NULL) a_node_branch_index = getNodeBranchIndex(a_node);
    // int branch_choice;
    Rect this_rect;
    for(int i = tba_start; i < tba_start + tba_count; i++){
        this_rect = (Rect(tba_mins[i], tba_maxes[i]));
        tba_map[i-tba_start] = PickBranch(&this_rect, a_node);

        // tba_to_children_ids[branch_choice][tba_to_children_counts[branch_choice]] = tba_ids[i];
        // for(int folan = 0; folan < NUMDIMS; folan++){
        //     tba_to_children_mins[branch_choice][tba_to_children_counts[branch_choice]][folan] = tba_mins[i][folan];
        //     tba_to_children_maxes[branch_choice][tba_to_children_counts[branch_choice]][folan] = tba_maxes[i][folan];
        // }
        // memcpy(tba_to_children_mins[branch_choice][tba_to_children_counts[branch_choice]], tba_mins[i], NUMDIMS*sizeof(float));
        // memcpy(tba_to_children_maxes[branch_choice][tba_to_children_counts[branch_choice]], tba_maxes[i], NUMDIMS*sizeof(float));


        tba_to_children_counts[tba_map[i-tba_start]]++;


        // cout << "RECT is "  << this_rect.m_min[0] << " " << this_rect.m_min[1] << " " << this_rect.m_max[0] << " " << this_rect.m_max[1] << endl;
        
        if(a_node->m_parent != NULL) {
            (a_node->m_parent)->m_branch[a_node_branch_index].m_rect = CombineRect(&this_rect, &((a_node->m_parent)->m_branch[a_node_branch_index].m_rect));
        }
        // else root_cover = CombineRect(&root_cover, &(Rect(tba_mins[i], tba_maxes[i])));
    }




    // accumulated of counts
    // sort mins, maxes, and ids based on map

    // cout << "before sort in traverse" << endl;
    // for(int i = 0; i < tba_count; i++)
    //     cout << tba_ids[i] << " ";
    // cout <<"\n";
    
    quickSort(tba_map, tba_start,tba_start + tba_count-1, tba_start);
    
    
    int *tba_to_children_acc_counts = (int *)calloc(a_node->m_count, sizeof(int));
    tba_to_children_acc_counts[0] = 0;
    for(int i = 1; i < a_node->m_count ; i++)
        tba_to_children_acc_counts[i] = tba_to_children_acc_counts[i-1] + tba_to_children_counts[i-1];
    // cout << "after sort in traverse" << endl;
    // for(int i = 0; i < tba_count; i++)
    //     cout << tba_ids[i] << " ";
    // cout <<"\n";
    // for(int i = 0; i < a_node->m_count; i++)
    //     cout << tba_to_children_counts[i] << " ";
    // cout <<"\n";

    //  for(int i = 0; i < a_node->m_count; i++)
    //     cout << tba_to_children_acc_counts[i] << " ";
    // cout <<"\n";
        




    for(int i =0; i < a_node->m_count; i++){
        if(tba_to_children_counts[i] > 0){
            // cout << "traverse recur call, tba start " << tba_start << " " << tba_to_children_acc_counts[i] << endl;
            ASSERT(0 <= tba_start + tba_to_children_acc_counts[i] && tba_start + tba_to_children_acc_counts[i] < DATA_COUNT);
            TraverseDownTreeTilLeaf(a_node->m_branch[i].m_child, tba_start + tba_to_children_acc_counts[i], tba_to_children_counts[i], overlapping_leaves);
        }
    }


    // for(int j = 0; j < tba_count; j++)
    // {
    //     free(tba_mins[j]);
    //     free(tba_maxes[j]);
    // }
    // free(tba_ids);
    // free(tba_mins);
    // free(tba_maxes);

    // free(tba_to_children_counts);
    // if(a_node->m_level != 1){
    //     for(int i = 0; i < a_node->m_count;i++)
    //     {
    //         for(int j = 0; j < tba_count; j++)
    //         {
    //             free(tba_to_children_mins[i][j]);
    //             free(tba_to_children_maxes[i][j]);
    //         }
    //         free(tba_to_children_ids[i]);
    //         free(tba_to_children_mins[i]);
    //         free(tba_to_children_maxes[i]);
    //     }
    // }
    // cout << "Before free in traverse\n";
    free(tba_to_children_counts);
    free(tba_map);
    free(tba_to_children_acc_counts);
    // cout << "After free in traverse\n";
}



RTREE_TEMPLATE
bool RTREE_QUAL::Search_2027(Node *a_node, Rect *a_rect, int &a_foundCount, LeavesToAdd_v2 *all_lta, int tba_start, int tba_count, vector<overlapping_leaf_tbas> *overlapping_leaves) {
    //            ASSERT(a_node);
    //            ASSERT(a_node->m_level >= 0);
    //            ASSERT(a_rect);
    //            ASSERT(a_node->m_count == a_node->m_branch.size());
    //        cout << "looked_through_a_node" << endl;

    if (a_node->IsInternalNode()) {
        // cout << "was internal level: " << a_node->m_level << endl;
        // This is an internal node in the tree
        // TODO:
        // i think I should assign the tba to each of the children
        // then I also need to put that as the input of that search recurse
        // then we also have to expand the mbr of the children here

        // very very bad code
        if(a_node->m_count > 0){
            // old
            // vector<Branch> tba_to_children[a_node->m_count];
            // int a_node_branch_index;
            // if(a_node->m_parent != NULL) a_node_branch_index = getNodeBranchIndex(a_node);
            // int branch_choice;
            // for(int i =0; i < tba->size(); i++){
            //     branch_choice = PickBranch(&tba->at(i).m_rect, a_node);
            //     tba_to_children[branch_choice].push_back(tba->at(i));
            //     // LET'S EXTEND THE MBR OF THE INNER NODES HERE
            //     // this node's mbr would get extended
            //     // i'll do it stupidly for now
            //     // TODO there must be a better way
            //     if(a_node->m_parent != NULL) {
            //         (a_node->m_parent)->m_branch[a_node_branch_index].m_rect = CombineRect(&tba->at(i).m_rect, &((a_node->m_parent)->m_branch[a_node_branch_index].m_rect));
            //     }
            //     else root_cover = CombineRect(&root_cover, &tba->at(i).m_rect);
            // }
            // old


            // new
            // float tba_to_children_mins[a_node->m_count][tba_count][NUMDIMS];
            // float tba_to_children_maxes[a_node->m_count][tba_count][NUMDIMS];
            // int tba_to_children_ids[a_node->m_count][tba_count];

            // cout << "internal node search, before mallocs" << endl;

            #ifdef timebreak
                clock_t start_time = clock();
            #endif

            // float** tba_to_children_mins[a_node->m_count];
            // float** tba_to_children_maxes[a_node->m_count];
            // int* tba_to_children_ids[a_node->m_count];
            // for(int i =0; i < a_node->m_count; i++){
            //     tba_to_children_ids[i] = (int*)malloc(tba_count*sizeof(int));
            //     tba_to_children_mins[i] = (float**)malloc(tba_count*sizeof(float*));
            //     tba_to_children_maxes[i] = (float**)malloc(tba_count*sizeof(float*));
            //     for(int j = 0; j < tba_count; j++){
            //         tba_to_children_mins[i][j] = (float*)malloc(NUMDIMS*sizeof(float));
            //         tba_to_children_maxes[i][j] = (float*)malloc(NUMDIMS*sizeof(float));

            //     }
            // }

            

            // int tba_to_children_counts[a_node->m_count];
            int *tba_to_children_counts = (int *)calloc(a_node->m_count, sizeof(int));
            int *tba_map = (int *)calloc(tba_count,sizeof(int));


            // cout << "after mallocs tba_count " << tba_count << endl;

            int a_node_branch_index;
            if(a_node->m_parent != NULL) a_node_branch_index = getNodeBranchIndex(a_node);
            // int branch_choice;
            Rect this_rect;
            for(int i = tba_start; i < tba_start + tba_count; i++){
                // this_rect = (Rect(tba_mins+i*NUMDIMS, tba_maxes+i*NUMDIMS));
                this_rect = (Rect(tba_mins[i], tba_maxes[i]));

                tba_map[i-tba_start] = PickBranch(&this_rect, a_node);

                // tba_to_children_ids[branch_choice][tba_to_children_counts[branch_choice]] = tba_ids[i];
                // for(int folan = 0; folan < NUMDIMS; folan++){
                //     tba_to_children_mins[branch_choice][tba_to_children_counts[branch_choice]][folan] = tba_mins[i][folan];
                //     tba_to_children_maxes[branch_choice][tba_to_children_counts[branch_choice]][folan] = tba_maxes[i][folan];
                // }
                // memcpy(tba_to_children_mins[branch_choice][tba_to_children_counts[branch_choice]], tba_mins[i], NUMDIMS*sizeof(float));
                // memcpy(tba_to_children_maxes[branch_choice][tba_to_children_counts[branch_choice]], tba_maxes[i], NUMDIMS*sizeof(float));


                tba_to_children_counts[tba_map[i-tba_start]]++;


                // cout << "RECT is "  << this_rect.m_min[0] << " " << this_rect.m_min[1] << " " << this_rect.m_max[0] << " " << this_rect.m_max[1] << endl;
                
                if(a_node->m_parent != NULL) {
                    (a_node->m_parent)->m_branch[a_node_branch_index].m_rect = CombineRect(&this_rect, &((a_node->m_parent)->m_branch[a_node_branch_index].m_rect));
                }
                // else root_cover = CombineRect(&root_cover, &(Rect(tba_mins[i], tba_maxes[i])));
            }

            #ifdef timebreak
                tba_assignment_time += (clock() - start_time); 
            #endif

            // accumulated of counts
            // sort mins, maxes, and ids based on map
            // cout << "before sort in 2027" << endl;
            // for(int i = 0; i < tba_count; i++)
            //     cout << tba_ids[i] << " ";
            // cout <<"\n";   
            quickSort(tba_map, tba_start, tba_start + tba_count - 1, tba_start);
            // cout << "after sort in 2027" << endl;
            // for(int i = 0; i < tba_count; i++)
            //     cout << tba_ids[i] << " ";
            // cout <<"\n";
            // for(int i = 0; i < a_node->m_count; i++)
            //     cout << tba_to_children_counts[i] << " ";
            // cout <<"\n";
            int *tba_to_children_acc_counts = (int *)calloc(a_node->m_count, sizeof(int));
            tba_to_children_acc_counts[0] = 0;
            for(int i = 1; i < a_node->m_count; i++)
                tba_to_children_acc_counts[i] = tba_to_children_acc_counts[i-1] + tba_to_children_counts[i-1];
            // for(int i = 0; i < a_node->m_count; i++)
            //     cout << tba_to_children_acc_counts[i] << " ";
            // cout <<"\n";



            // new


            // int check_sum = 0;
            // for(int i = 0; i < a_node->m_count; i++){
            //     check_sum += tba_to_children[i].size();
            // }
            // if(check_sum != tba->size()){
            //     cout << "TERRIBLE JOB PEOPLE. WHAT IS HAPPENING HERE... " << check_sum << " " << tba->size() << endl;
            // }

            // for (int index = 0; index < a_node->m_count; ++index) {
            //     cout << "^^^^^^^^^^^ tba_children[" << index << "].size()  " << tba_to_children[index].size() << endl;
            // }
            for (int index = 0; index < a_node->m_count; ++index) {
                // cout << "___________ tba_children[" << index << "].size()  " << tba_to_children[index].size() << endl;
                // TODO: then we also have to expand the mbr of the children here
                // I will do it after it returns so it doesn't mess up the overlap condition check

                // if((a_node->m_branch[index].m_child)->m_id == 59078) {
                //     cout << "in parent of 59078, checking overlap..." << endl;
                //     cout << "its rect is: " << (a_node->m_branch[index].m_rect).m_min[0] << " " << (a_node->m_branch[index].m_rect).m_max[0] << " ";
                //     cout << (a_node->m_branch[index].m_rect).m_min[1] << " " << (a_node->m_branch[index].m_rect).m_max[1] << endl;
                // }

                if (Overlap(a_rect, &a_node->m_branch[index].m_rect)) {
                    // if((a_node->m_branch[index].m_child)->m_id == 1575) cout << "it overlapped, rec call" << endl;
                    // cout << "overlapped internaL node" << endl;
                    ASSERT(0 <= tba_start + tba_to_children_acc_counts[index] && tba_start + tba_to_children_acc_counts[index] < DATA_COUNT);

                    if (!Search_2027(a_node->m_branch[index].m_child, a_rect, a_foundCount, all_lta, tba_start + tba_to_children_acc_counts[index], tba_to_children_counts[index], overlapping_leaves)) {
                        // The callback indicated to stop searching
                        return false;
                    }
                }
                else{
                    // even if it doesn't overlap, we still have to add it to the overlapping_leaves to add the tba
                    
                    // TODO: fixt this disaster!!!!!
                    // cout << "did not overlap, tba_children[" << index << "].size() " << tba_to_children[index].size() << endl;
                    // if(tba_to_children[index].size() > 0) {
                    if(tba_to_children_counts[index] > 0) {
                        // cout << "adding non-overlapping child to overlapping_leaves " << endl;
                        // ASSERT(a_node->m_level == 1);
                        // cout << "I just realised that it may not be a leaf, and then I don't know what to do...." << endl;
                        // overlapping_leaves->push_back(overlapping_leaf_tbas(a_node->m_branch[index].m_child, (tba_to_children[index])));
                        // cout << "after adding to ol, ol[-1].tba.size() " << (overlapping_leaves->back()).this_tba.size() << endl;
                        // THIS IS VERY BAD
                        // HAVE TO FIGURE THIS OUT AND MAKE SURE IT DOESN'T HAPPEN
                        // BANDAID SOLUTION FOR NOW:
                        // EXPAND THE MBRS HERE TOO...

                        // for(int ii = 0; ii < tba_to_children[index].size(); ii++){
                        //     a_node->m_branch[index].m_rect = CombineRect(&(a_node->m_branch[index].m_rect), &(tba_to_children[index][ii].m_rect));
                        // }
                        // cout << "before traverse call " << endl;
                        #ifdef timebreak
                            start_time = clock();
                        #endif
                        // cout << "calling traverse with start: " << tba_start << " " << tba_to_children_acc_counts[index] << endl;
                        ASSERT(0 <= tba_start + tba_to_children_acc_counts[index] && tba_start + tba_to_children_acc_counts[index] < DATA_COUNT);

                        TraverseDownTreeTilLeaf(a_node->m_branch[index].m_child, tba_start + tba_to_children_acc_counts[index], tba_to_children_counts[index], overlapping_leaves);
                        #ifdef timebreak
                        traversedown_time += (clock() - start_time);
                        #endif
                        // cout << "after traverse call " << endl;
                    }
                }
                // cout << "after all" << endl;
                // expand the mbrs
                // for(int ii=0; ii < tba_to_children[index].size(); ii++){
                //     a_node->m_branch[index].m_rect = CombineRect(&tba_to_children[index][ii].m_rect, &(a_node->m_branch[index].m_rect));
                // }
                // THIS WOULD NOT WORK
                // THE LEAVES ALSO RETURN HERE
                // THEY DON'T EVEN BECAUSE THERE IS THE RETURN
                // BUT EVEN WITHOUT THAT A LEAF THAT HAS GOTTEN CRACKED COULD RETURN HERE, WE WOULD EXTEND ITS MBR FOR NO REASON

            }


            // FREE THE ARRAYS
            // if(a_node->m_level != 1){
            //     for(int i = 0; i< a_node->m_count; i++){
            //         for(int j = 0; j < tba_count; j++){
            //             free(tba_to_children_mins[i][j]);
            //             free(tba_to_children_maxes[i][j]);
            //         }
            //         free(tba_to_children_mins[i]);
            //         free(tba_to_children_maxes[i]);
            //         free(tba_to_children_ids[i]);
            //     }
            // }
            // the way we did in traverse down
            // if(a_node->m_parent != NULL){
            //     for(int ii = 0; ii < tba_count; ii++){
            //         free(tba_mins[ii]);
            //         free(tba_maxes[ii]);
            //     }
            //     free(tba_ids);
            //     free(tba_mins);
            //     free(tba_maxes);
            // }

            // cout << "Before free in search\n";
            free(tba_to_children_counts);
            free(tba_map);
            free(tba_to_children_acc_counts);
            // cout << "After free in search\n";

        }
        else{
            cout << "BAD BAD PEOPLE. BAD BAD CODE. FIX ME!! tba.size() " << tba_count << " " << a_node->m_id << " count " << a_node->m_count << endl;
        }
    }
    else {
        // cout << "©" << endl;
        // if(a_node->m_id == 39024) cout << "checking out leaf 39024 L, h, R " << a_node->m_L << " " << a_node->m_holes << " " << a_node->m_R << endl;
        if(a_node->isRegular()) {
            // cout <<"it was regular" << endl;
            int a_node_branch_index;
            if(a_node->m_parent != NULL) a_node_branch_index = getNodeBranchIndex(a_node);   
            // cout << "#############" << endl;   
            for (int index = a_node->m_L + a_node->m_holes; index < a_node->m_R; index++) {
                // if(m_data_arr_ids[index] == 24115){
                //         cout << "checking this guy out ebfore overlap check" << endl;
                // }
                if (Overlap(a_rect, m_data_arr_mins[index], m_data_arr_maxes[index])) {
                    // if(m_data_arr_ids[index] == 24115){
                    //     cout << "checkign this guy out" << endl;
                    //     cout << m_data_arr_mins[index][0] << " " << m_data_arr_mins[index][1] << " " << m_data_arr_maxes[index][0] << " " << m_data_arr_maxes[index][1] << endl;
                    // }
                    // cout << m_data_arr_ids[index] << endl;
                    a_foundCount++;
                }
                // else{
                //     cout << "did not overlap " << m_data_arr_ids[index] << " index: " << index << endl;
                //     cout << "shape: " ;
                //     for(int ii = 0; ii < NUMDIMS; ii++) cout << m_data_arr_mins[index][ii] << " " << m_data_arr_maxes[index][ii] << ", ";
                //     cout << endl;
                // }
            }
            // cout << "#############" << endl;
            // 2023
            // for(int ii =0; ii < tba->size(); ii++){
            //     if(Overlap(a_rect, &tba->at(ii).m_rect)){
            //         ++a_foundCount;
            //     }
            // }
            // they all shoudl overlap
            // that is the whole reason they eneded up in the tba
            // so I'll just do it in QueryAdaptive when I find them...
            // for(int i = 0; i < tba->size(); i++){
            //     if((tba->at(i)).m_data == 104841){
            //         cout << "was assigned to regular leaf " << a_node->m_id << endl;
            //     }
            // }
            // if(a_node->m_id == 2026 || a_node->m_id == 3762){
            //     cout << "pushing " << a_node->m_id << " into ol with start " << tba_start << endl;
            // }
            overlapping_leaves->push_back(overlapping_leaf_tbas(a_node, tba_start, tba_count));
            if(tba_count > 0) {
                // cout << "adding regular leaf to ol" << endl;
                // overlapping_leaves->push_back(overlapping_leaf_tbas(a_node, *tba));
                // and let's extend the MBR
                // TODO: I am not sure if this assignment thing will actually copy the rects ro not....
                if(a_node->m_parent != NULL) {
                    Rect new_rect = (a_node->m_parent)->m_branch[a_node_branch_index].m_rect;
                    Rect trash_rect;
                    for(int ii = tba_start; ii < tba_start + tba_count; ii++){
                        // if(tba->at(ii).m_data == 100352){
                        //     cout << "inserting 100351 into regular leaf " << a_node->m_id << endl;
                        // }
                        trash_rect = (Rect(tba_mins[ii], tba_maxes[ii]));
                        new_rect = CombineRect(&trash_rect, &new_rect);
                    }
                    // cout << "new cover rect: "<< endl;
                    // for(int ii = 0 ;ii < NUMDIMS; ii++) cout << new_rect.m_min[ii] << " " << new_rect.m_max[ii] << " , ";
                    // cout << endl;
                    ((a_node->m_parent)->m_branch[a_node_branch_index]).m_rect = new_rect;
                    // cout << "new cover rect in tree: "<< endl;
                    // for(int ii = 0 ;ii < NUMDIMS; ii++) cout << (a_node->m_parent)->m_branch[a_node_branch_index].m_rect.m_min[ii] << " " << (a_node->m_parent)->m_branch[a_node_branch_index].m_rect.m_max[ii] << " , ";
                    // cout << endl;

                }
                // else{
                //     // it's the root that is a leaf still...
                //     Rect new_rect = root_cover;
                //     for(int ii =0; ii < tba->size(); ii++){
                //         new_rect = CombineRect(&tba->at(ii).m_rect, &new_rect);
                //     }
                //     root_cover = new_rect;
                // }
            }
            // 2023

        }
        else{

            #ifdef timebreak
                auto start_time = clock();
            #endif

            // if(a_node->m_id == 39024){
            //     cout << "cracking leaf " << a_node->m_id << "!!!!!! L, h, R " << a_node->m_L << " " <<  a_node->m_holes <<  " " << a_node->m_R << endl;
            // }

            // for(int i = 0; i < tba->size(); i++){
            //     if((tba->at(i)).m_data == 104841){
            //         cout << "was assigned to irregular leaf " << a_node->m_id << " so must be inone fo the children " << endl;
            //     }
            // }
            // if(a_node->m_id == 627){
            //     cout << "CRACKING 627!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
            // }
            // call elegant_2dc
            //            cout << "found count before " << a_foundCount << endl;
            //            a_foundCount += Elegant_2dc_v2(a_node, a_rect, all_lta);
            //            cout << "found count after " << a_foundCount << endl;

            // srand(time(NULL));

            Rect node_cover;
            if(a_node->m_parent != NULL) {
                node_cover = ((a_node->m_parent)->m_branch[getNodeBranchIndex(a_node)]).m_rect;
            } else{
                node_cover = root_cover;
            }


            Rect this_piece_rect = Rect(node_cover.m_min, node_cover.m_max);
            // 2023
            // int this_piece_L = a_node->m_L;
            int this_piece_L = a_node->m_L;
            int this_piece_R = a_node->m_R;

            int this_piece_holes = a_node->m_holes;


            int choice, mm, chosen_axis;

            LTA_v2 this_lta;
            this_lta.how_many_created =0;
            this_lta.this_leaf = a_node;
            Rect potential_query_piece = Rect(node_cover.m_min, node_cover.m_max); 
            Rect other_piece;


            int crack_index, filan, ashghal;
            int sum_of_choices = 0; int sum_of_mms = 0;

            // 18NOV
            int largest_piece_index = 0;
            int stochastic_crack_axis = 0; float stochastic_crack_pivot;
            float X, X2, X3;
            int stochastic_crack_index;

            // 14 sep
            float old_pqr_chosen_axis[2]; // 0: min, 1: max

            for(int i =0; i < 2*NUMDIMS; i++){
                // if(a_node->m_id == 39024){
                //     cout << "in crack loop tpl, tpr, tph " << this_piece_L << " " << this_piece_R << " " << this_piece_holes << endl;
                // }
                //            for(int i =0; i < 4; i++){

                for(filan=0; filan < NUMDIMS; filan++) {
                    other_piece.m_min[filan] = potential_query_piece.m_min[filan];
                    other_piece.m_max[filan] = potential_query_piece.m_max[filan];
                }



                if(i == ((2 * NUMDIMS) - 1)){
//                    cout << "last choice\n";
                    choice = (((2*NUMDIMS) - 1) * NUMDIMS) - sum_of_choices;
                }
                else {
                    choice = ChooseAxisLargestSideMid(this_piece_rect, *a_rect);
                }
                // 11oct 16:26
                // if(a_node->m_id == 39024) cout << "cp1 " << endl;
                if(choice < 0) break;
                //
                if(sum_of_mms >= NUMDIMS) mm = 1;
                else if (sum_of_mms <= (-1*NUMDIMS)) mm = 0;
                else {
                    mm = choice % 2;
                    if(mm == 0) sum_of_mms++;
                    else sum_of_mms--;
                }

//                cout << "sum of choices "<< sum_of_choices << " choice " << choice << endl;
                chosen_axis = choice / 2;
//                cout << "chosen axis " << chosen_axis << endl;

//                ASSERT(chosen_axis >= 0);
//                ASSERT(chosen_axis < NUMDIMS);
                // if(a_node->m_id == 39024) cout << "cp 2" << endl;
                if(chosen_axis >= NUMDIMS || chosen_axis < 0) break;

                sum_of_choices += choice;

//                potential_query_piece.m_min[chosen_axis] = 21474836;
//                potential_query_piece.m_max[chosen_axis] = -21474836;
                other_piece.m_min[chosen_axis] = 21474836;
                other_piece.m_max[chosen_axis] = -21474836;

                // 14 Sep
                old_pqr_chosen_axis[0] = potential_query_piece.m_min[chosen_axis];
                old_pqr_chosen_axis[1] = potential_query_piece.m_max[chosen_axis];

                // if(a_node->m_id == 39024) cout << "cp 3 mm " << mm << " caxis: " << chosen_axis << endl;
 
                if(mm == 0) {
                    // if(a_node->m_id == 39024) cout << "pivot: " << a_rect->m_min[chosen_axis] << endl;
                    potential_query_piece.m_min[chosen_axis] = 21474836;
                    crack_index = CrackOnAxisCompMax_v45( this_piece_L + this_piece_holes, this_piece_R, a_rect->m_min[chosen_axis],
                                                        chosen_axis, &other_piece, &potential_query_piece);
                }
                else{
                    // if(a_node->m_id == 39024) cout << "pivot: " << a_rect->m_max[chosen_axis] << endl;
                    potential_query_piece.m_max[chosen_axis] = -21474836;
                    crack_index = CrackOnAxisCompMin_v45( this_piece_L + this_piece_holes, this_piece_R, a_rect->m_max[chosen_axis],
                                                        chosen_axis, &potential_query_piece, &other_piece);
                }

//                cout << "INJA DARIM CHECK MIKONIM " << potential_query_piece.m_min[chosen_axis] << " " << potential_query_piece.m_max[chosen_axis] << endl;
//                if(potential_query_piece.m_min[chosen_axis] == 21474836 || potential_query_piece.m_max[chosen_axis] == -21474836) {
//                    cout << "naaaa nabayd in etefagh mioftad " << chosen_axis << " mm " << mm << " L " << this_piece_L << " crackindex " << crack_index << " R " << this_piece_R << endl;
//                }
//                if(other_piece.m_min[chosen_axis] > other_piece.m_max[chosen_axis] && other_piece.m_min[chosen_axis] != 21474836) {
//                    cout << "naaaaaaa nabayd in etefagh mioftad " << chosen_axis << " mm " << mm << " L " << this_piece_L << " crackindex " << crack_index << " R " << this_piece_R << endl;
//                    cout << "values " << other_piece.m_min[chosen_axis] << " " << other_piece.m_max[chosen_axis] << endl;
//                }

                // if(a_node->m_id == 39024) cout << "cp 4" << endl;
                if((mm == 0 && (crack_index - this_piece_L - this_piece_holes) >= MINDATAPOINTS) ||
                (mm == 1 && (this_piece_R - crack_index) >= MINDATAPOINTS)){

                    // distribute holes
                    // if(a_node->m_id == 39024) {
                    //     cout << "after min check if cracking 39024 at " << crack_index << endl;
                    // }


                    if(mm == 0){
                        
                        int half_the_holes = (int)(this_piece_holes / 2);
                        this_lta.Ls[this_lta.how_many_created] = this_piece_L;

                        if(this_piece_R != crack_index){
                            if(crack_index - this_piece_holes - this_piece_L >= half_the_holes)
                            {
                                // this_lta.Ls[this_lta.how_many_created] = this_piece_L;
                                // this_lta.Rs[this_lta.how_many_created] = crack_index;

                                
                                this_lta.holes[this_lta.how_many_created] = this_piece_holes - half_the_holes;

                                // moving h/2 from the end of this piece to the start, to leave the rest of the holes to be distributed again
                                memcpy(m_data_arr_mins[this_piece_L + this_piece_holes - half_the_holes], m_data_arr_mins[crack_index - half_the_holes], NUMDIMS*half_the_holes*sizeof(float));
                                memcpy(m_data_arr_maxes[this_piece_L + this_piece_holes - half_the_holes], m_data_arr_maxes[crack_index - half_the_holes], NUMDIMS*half_the_holes*sizeof(float));
                                memcpy(m_data_arr_ids + this_piece_L + this_piece_holes - half_the_holes, m_data_arr_ids + crack_index - half_the_holes, half_the_holes*sizeof(int));

                                this_lta.Rs[this_lta.how_many_created] = crack_index - half_the_holes;

                                // if(a_node->m_id == 1052) {
                                //     cout << "it was mm = 0, copying from " <<  crack_index - half_the_holes << " to " << this_piece_L + this_piece_holes - half_the_holes << " this many " << half_the_holes << endl;
                                // }

                                // set up ptcrack
                                this_piece_L = crack_index - half_the_holes;
                                this_piece_holes = half_the_holes;
                            }
                            else
                            {
                                // this_lta.Ls[this_lta.how_many_created] = this_piece_L;

                                this_lta.holes[this_lta.how_many_created] = this_piece_holes - half_the_holes;
                                int how_many_to_move = crack_index - this_piece_holes - this_piece_L;

                                memcpy(m_data_arr_mins[this_piece_L + this_piece_holes - half_the_holes], m_data_arr_mins[crack_index - how_many_to_move], NUMDIMS*how_many_to_move*sizeof(float));
                                memcpy(m_data_arr_maxes[this_piece_L + this_piece_holes - half_the_holes], m_data_arr_maxes[crack_index - how_many_to_move], NUMDIMS*how_many_to_move*sizeof(float));
                                memcpy(m_data_arr_ids + this_piece_L + this_piece_holes - half_the_holes, m_data_arr_ids + crack_index - how_many_to_move, how_many_to_move*sizeof(int));

                                this_lta.Rs[this_lta.how_many_created] = crack_index - half_the_holes;
                                this_piece_L = crack_index - half_the_holes;
                                this_piece_holes = half_the_holes;

                            }

                        }else {
                            this_lta.Rs[this_lta.how_many_created] = crack_index;
                            this_lta.holes[this_lta.how_many_created] = this_piece_holes;

                            this_piece_holes = 0; 
                            this_piece_L = crack_index;
                        }

                    }
                    else{
                        // this_lta.Ls[this_lta.how_many_created] = crack_index;
                        int half_the_holes = (int)(this_piece_holes / 2);
                        this_lta.Rs[this_lta.how_many_created] = this_piece_R;

                        if(crack_index != this_piece_L + this_piece_holes){
                            if(crack_index - this_piece_holes - this_piece_L >= half_the_holes)
                            {
                                // this_lta.Rs[this_lta.how_many_created] = this_piece_R;

                                

                                this_lta.holes[this_lta.how_many_created] = half_the_holes;

                                // moving h/2 from the end of this piece to the start, to leave the rest of the holes to be distributed again
                                memcpy(m_data_arr_mins[this_piece_L + this_piece_holes - half_the_holes], m_data_arr_mins[crack_index - half_the_holes], NUMDIMS*half_the_holes*sizeof(float));
                                memcpy(m_data_arr_maxes[this_piece_L + this_piece_holes - half_the_holes], m_data_arr_maxes[crack_index - half_the_holes], NUMDIMS*half_the_holes*sizeof(float));
                                memcpy(m_data_arr_ids + this_piece_L + this_piece_holes - half_the_holes, m_data_arr_ids + crack_index - half_the_holes, half_the_holes*sizeof(int));

                                this_lta.Ls[this_lta.how_many_created] = crack_index - half_the_holes;

                                // if(a_node->m_id == 1052) {
                                //     cout << "case 1 it was mm = 1, copying from " <<  crack_index - half_the_holes << " to " << this_piece_L + this_piece_holes - half_the_holes << " this many " << half_the_holes << endl;
                                //     cout << "copied " << m_data_arr_ids[this_piece_L + this_piece_holes - half_the_holes] << endl; 
                                // }

                                // set up ptcrack
                                this_piece_R = crack_index - half_the_holes;
                                this_piece_holes -= half_the_holes;
                            }
                            else
                            {
                                // this_lta.Ls[this_lta.how_many_created] = crack_index;
                                // this_lta.Rs[this_lta.how_many_created] = this_piece_R;

                                // int half_the_holes = (int)(this_piece_holes / 2);

                                this_lta.holes[this_lta.how_many_created] = this_piece_holes - half_the_holes;
                                int how_many_to_move = crack_index - this_piece_holes - this_piece_L;

                                // moving h/2 from the end of this piece to the start, to leave the rest of the holes to be distributed again
                                memcpy(m_data_arr_mins[this_piece_L + half_the_holes], m_data_arr_mins[crack_index - how_many_to_move], NUMDIMS*how_many_to_move*sizeof(float));
                                memcpy(m_data_arr_maxes[this_piece_L + half_the_holes], m_data_arr_maxes[crack_index - how_many_to_move], NUMDIMS*how_many_to_move*sizeof(float));
                                memcpy(m_data_arr_ids + this_piece_L + half_the_holes, m_data_arr_ids + crack_index - how_many_to_move, how_many_to_move*sizeof(int));

                                this_lta.Ls[this_lta.how_many_created] = this_piece_L + half_the_holes + how_many_to_move;

                                // if(a_node->m_id == 1052) {
                                //     cout << "case 2 it was mm = 1, copying from " <<  crack_index - half_the_holes << " to " << this_piece_L + this_piece_holes - half_the_holes << " this many " << how_many_to_move << endl;
                                //     cout << "copied " << m_data_arr_ids[this_piece_L + this_piece_holes - half_the_holes] << endl; 
                                // }

                                // set up ptcrack
                                this_piece_R = this_piece_L + half_the_holes + how_many_to_move;
                                this_piece_holes = half_the_holes;
                            }
                        } else {
                            ASSERT(crack_index - this_piece_L == this_piece_holes);
                            this_lta.Ls[this_lta.how_many_created] = this_piece_L;
                            this_lta.holes[this_lta.how_many_created] = this_piece_holes;

                            this_piece_R = crack_index;
                        }
                    }
                    this_piece_rect = Rect(potential_query_piece.m_min, potential_query_piece.m_max);

                    // 2023
                    // dummy
                    // this_lta.holes[this_lta.how_many_created] = 0;


                    for(ashghal = 0; ashghal < NUMDIMS; ashghal ++){
//                        if(other_piece.m_min[ashghal]  > other_piece.m_max[ashghal]){
//                            cout << "IN SEARCH2026, ADDING OTHERPIECE, MIN>MAX dim " <<  ashghal << " chosen_axis " << chosen_axis  << " mm " << mm << " L " <<  this_lta.Ls[this_lta.how_many_created] << " R " << this_lta.Rs[this_lta.how_many_created] << endl;
//                            cout << "values  min " << other_piece.m_min[ashghal] <<  " max " << other_piece.m_max[ashghal] << endl;
//                        }

                        this_lta.crack_covers_mins[this_lta.how_many_created][ashghal] = other_piece.m_min[ashghal];
                        this_lta.crack_covers_maxes[this_lta.how_many_created][ashghal] = other_piece.m_max[ashghal];
                    }

                    // if(a_node->m_id == 39024) {
                    //     cout << "cover added to lta after update" << endl;
                    //     cout << this_lta.crack_covers_mins[this_lta.how_many_created][0] << " " << this_lta.crack_covers_maxes[this_lta.how_many_created][0] << " ";
                    //     cout << this_lta.crack_covers_mins[this_lta.how_many_created][1] << " " << this_lta.crack_covers_maxes[this_lta.how_many_created][1] << endl;


                    //     cout << "tpr:" << endl;
                    //     cout << this_piece_rect.m_min[0] << " " << this_piece_rect.m_max[0] << " " << this_piece_rect.m_min[1] << " " << this_piece_rect.m_max[1] << endl;
                    
                    //     cout << "pqr:" << endl;
                    //     cout << potential_query_piece.m_min[0] << " " << potential_query_piece.m_max[0] << " " << potential_query_piece.m_min[1] << " " << potential_query_piece.m_max[1] << endl;
                    
                    // }

                    this_lta.how_many_created++;
                }
                else{
                    // 14 oct
                    // baghiash bi ahamiate chon inja mitoonim break konim
                    // baresh dashtam dobare Panos movafegh nabood
//                    break;

                    // in baraye ine k oon right va left cover k partition midim ye meghdar avaliyeE dashte bashan
                    // in joori b 2 tashoon potential query midim

                    // FK KONAM ALREADY MOSAVIAN VALI MOTMAEN NISTAM
                    // vali hala k in taghyir ro dadim doroste k ghablesh ba meghdar ghabli set kardim
                    // potential_query_piece = Rect(this_piece_rect.m_min, this_piece_rect.m_max);
                    // HALA BARAYE CHOICE BADI THIS_PIECE_RECT RO AVAZ MIKONIM

                    // if(a_node->m_id == 39024){
                    //     cout << "in case 123 before update" << endl;
                
                    //     cout << "tpr:" << endl;
                    //     cout << this_piece_rect.m_min[0] << " " << this_piece_rect.m_max[0] << " " << this_piece_rect.m_min[1] << " " << this_piece_rect.m_max[1] << endl;
                    
                    //     cout << "pqr:" << endl;
                    //     cout << potential_query_piece.m_min[0] << " " << potential_query_piece.m_max[0] << " " << potential_query_piece.m_min[1] << " " << potential_query_piece.m_max[1] << endl;
                    
                    // }

                    if( mm == 0 ) {
                        this_piece_rect.m_min[chosen_axis] = a_rect->m_min[chosen_axis];
                        potential_query_piece.m_min[chosen_axis] = old_pqr_chosen_axis[0]; 
                        }
                    else {
                        this_piece_rect.m_max[chosen_axis] = a_rect->m_max[chosen_axis];
                        potential_query_piece.m_max[chosen_axis] =old_pqr_chosen_axis[1];
                        }

                    // if(a_node->m_id == 39024){
                    //     cout << "in case 123 after update" << endl;
                
                    //     cout << "tpr:" << endl;
                    //     cout << this_piece_rect.m_min[0] << " " << this_piece_rect.m_max[0] << " " << this_piece_rect.m_min[1] << " " << this_piece_rect.m_max[1] << endl;
                    
                    //     cout << "pqr:" << endl;
                    //     cout << potential_query_piece.m_min[0] << " " << potential_query_piece.m_max[0] << " " << potential_query_piece.m_min[1] << " " << potential_query_piece.m_max[1] << endl;
                    
                    // }

                }



                if((this_piece_R - this_piece_L - this_piece_holes) <= MAXDATAPOINTS) {
                    break;
                }

            }

            // scan last piece
            // age mikhay result baghean adadesh dorost bashe bayad injoori scan koni
            // vali injoori ye adad k hamishe bozorgtar mosavi javab vagheEye b dast miad
            // DEBUG
            // cout << "$$$$$$$$$$" << endl;
            // cout << "from " << this_piece_L + this_piece_holes << " to " << this_piece_R << endl;
            for(int folan = this_piece_L + this_piece_holes; folan < this_piece_R; folan++){
                // if( m_data_arr_ids[folan] == 110652)
                // {
                //     cout << "Item 110652 --> " << m_data_arr_mins[folan][0] << " " << m_data_arr_mins[folan][1] << " " << m_data_arr_maxes[folan][0] << " " << m_data_arr_maxes[folan][1] << endl;
                // }
                if (Overlap(a_rect, m_data_arr_mins[folan], m_data_arr_maxes[folan])) {
                    // if(m_data_arr_ids[folan] == 101221)
                    //     cout << m_data_arr_mins[folan][0] << " " << m_data_arr_mins[folan][1] << " " << m_data_arr_maxes[folan][0] << " " << m_data_arr_maxes[folan][1] << endl;
                    // cout << m_data_arr_ids[folan] << endl;
                    a_foundCount++;
                }
            }
            // cout << "$$$$$$$$$$" << endl;
            // a_foundCount += (this_piece_R - this_piece_L);
            // DEBUG

            // 2023
            // let's also set the qp_stuff
            this_lta.qp_L = this_piece_L;
            // cout << "in search, adding the tba.size() " << tba->size() << endl;
            // if(tba->size() > 0)  this_lta.qp_tba = *tba;
            // else this_lta.qp_tba = NULL;
            // cout << "Before " << this_lta.qp_tba.size() << endl;
            // this_lta.qp_tba = *tba;
            // this_lta.tba_mins = tba_mins;
            // this_lta.tba_maxes = tba_maxes;
            // this_lta.tba_ids = tba_ids;
            this_lta.tba_start = tba_start;
            this_lta.tba_count = tba_count;
            // cout << "After " << this_lta.qp_tba.size() << endl;

            // cout << "in search, adding the qp_tba " << this_lta.qp_tba->size() << " tba.size() " << tba->size() << endl;

            if(this_piece_L + this_piece_holes != this_piece_R) {

                // if(a_node->m_id == 39024) cout << "making last qiery piece" << endl;
                this_lta.Ls[this_lta.how_many_created] = this_piece_L;
                this_lta.Rs[this_lta.how_many_created] = this_piece_R;

                // 2023
                // dummy
                this_lta.holes[this_lta.how_many_created] = this_piece_holes;

                for(ashghal = 0; ashghal < NUMDIMS; ashghal ++){

//                    if(this_piece_rect.m_min[ashghal]  > this_piece_rect.m_max[ashghal]){
//                        cout << "IN SEARCH2026, ADDING LAST PIECE, MIN>MAX dim " << ashghal  << " mm " << mm << " L " <<  this_lta.Ls[this_lta.how_many_created] << " R " << this_lta.Rs[this_lta.how_many_created] << endl;
//                        cout << "values  min " << this_piece_rect.m_min[ashghal] <<  " max " << this_piece_rect.m_max[ashghal] << endl;
//                    }

                    this_lta.crack_covers_mins[this_lta.how_many_created][ashghal] = potential_query_piece.m_min[ashghal];
                    this_lta.crack_covers_maxes[this_lta.how_many_created][ashghal] = potential_query_piece.m_max[ashghal];
                }

                // if(a_node->m_id == 39024) {
                //     cout << "cover added to lta after update" << endl;
                //     cout << this_lta.crack_covers_mins[this_lta.how_many_created][0] << " " << this_lta.crack_covers_maxes[this_lta.how_many_created][0] << " ";
                //     cout << this_lta.crack_covers_mins[this_lta.how_many_created][1] << " " << this_lta.crack_covers_maxes[this_lta.how_many_created][1] << endl;

                // }
                 // TODO: we should also adjust the mbr of the qpeice
                 // done in add_ltas
                // for(int ii = 0; ii < tba_count; ii++){
                //     for(ashghal = 0; ashghal < NUMDIMS; ashghal ++){
                //         if(this_lta.crack_covers_mins[this_lta.how_many_created][ashghal] > tba_mins[ii][ashghal])
                //         {
                //             this_lta.crack_covers_mins[this_lta.how_many_created][ashghal] = tba_mins[ii][ashghal];
                //         }
                //         if(this_lta.crack_covers_maxes[this_lta.how_many_created][ashghal] < tba_maxes[ii][ashghal])
                //         {
                //             this_lta.crack_covers_maxes[this_lta.how_many_created][ashghal] = tba_maxes[ii][ashghal];
                //         }
                //     }
                // }   
                this_lta.how_many_created++;
            }
            else{
                this_lta.qp_L = this_lta.Ls[this_lta.how_many_created - 1];
            }

            
           

            // TODO: now I have to figure out what to do if the query piece is the one being stochastically cracked
            // I think for now, I will just not do it if it is the biggest piece...

            // aval bayad bozorgtarin tikke ro peyda konim
            /// START POF STOCH
            largest_piece_index = 0;
//            cout << "0: L " << this_lta.Ls[0] << " R " << this_lta.Rs[0] << endl;
            for(ashghal = 1; ashghal < this_lta.how_many_created; ashghal++){
//                cout << ashghal << ": L " << this_lta.Ls[ashghal] << " R " << this_lta.Rs[ashghal] << endl;
                if((this_lta.Rs[ashghal] - this_lta.Ls[ashghal]) > (this_lta.Rs[largest_piece_index] - this_lta.Ls[largest_piece_index])) largest_piece_index = ashghal;
            }

            // 2023
            if(this_lta.Ls[largest_piece_index] != this_lta.qp_L){

                // if(a_node->m_id == 39024) cout << "making doing stoch crack" << endl;
            //    cout << "the chosen largest piece: " << largest_piece_index << endl;
            //    cout << "L, R, h: " << this_lta.Ls[largest_piece_index] << " " << this_lta.Rs[largest_piece_index] << " " << this_lta.holes[largest_piece_index] << endl;
                // hala bayad in tike ro beshkanim roo y jaye random
                stochastic_crack_axis = 0;
                // faghat bar 2d emtehan mikonim pas hamin basse
                // if((this_lta.crack_covers_maxes[largest_piece_index][1] - this_lta.crack_covers_mins[largest_piece_index][1])  > (this_lta.crack_covers_maxes[largest_piece_index][0] - this_lta.crack_covers_mins[largest_piece_index][0]) ){
                //     stochastic_crack_axis = 1;_arr_mins[this_lta.Ls[largest_piece_index] + rand() % (this_lta.Rs[largest_piece_index] - this_lta.Ls[largest_piece_index])]
                // }
                for(ashghal = 1; ashghal < NUMDIMS; ashghal++){
                    if((this_lta.crack_covers_maxes[largest_piece_index][ashghal] - this_lta.crack_covers_mins[largest_piece_index][ashghal])  > (this_lta.crack_covers_maxes[largest_piece_index][stochastic_crack_axis] - this_lta.crack_covers_mins[largest_piece_index][stochastic_crack_axis]) ){
                        stochastic_crack_axis = ashghal;
                    }
                }

                X = m_data_arr_mins[this_lta.Ls[largest_piece_index] + rand() % (this_lta.Rs[largest_piece_index] - this_lta.Ls[largest_piece_index])][stochastic_crack_axis];    // split in random position
                X2 = m_data_arr_mins[this_lta.Ls[largest_piece_index] + rand() % (this_lta.Rs[largest_piece_index] - this_lta.Ls[largest_piece_index])][stochastic_crack_axis];    // split in random position
                X3 = m_data_arr_mins[this_lta.Ls[largest_piece_index] + rand() % (this_lta.Rs[largest_piece_index] - this_lta.Ls[largest_piece_index])][stochastic_crack_axis];    // split in random position

                // make X the median of the samples
                if (X2 > X) swap(X, X2);
                if (X > X3) swap(X, X3);
                if (X2 > X) swap(X, X2);
                stochastic_crack_pivot = X;

                // hala bayad vaghean crack konim
                Rect l_rect = Rect(this_lta.crack_covers_mins[largest_piece_index], this_lta.crack_covers_maxes[largest_piece_index]);
                Rect r_rect = Rect(this_lta.crack_covers_mins[largest_piece_index], this_lta.crack_covers_maxes[largest_piece_index]);
                l_rect.m_min[stochastic_crack_axis] = 21474836;
                l_rect.m_max[stochastic_crack_axis] = -21474836;
                r_rect.m_min[stochastic_crack_axis] = 21474836;

                stochastic_crack_index = CrackOnAxisCompMax_v45( this_lta.Ls[largest_piece_index] + this_lta.holes[largest_piece_index], this_lta.Rs[largest_piece_index], stochastic_crack_pivot,
                                                                stochastic_crack_axis, &l_rect, &r_rect);
                // cout << "stochastic crack index: " << stochastic_crack_index << endl;
                if((stochastic_crack_index - this_lta.Ls[largest_piece_index] - this_lta.holes[largest_piece_index]) >= MINDATAPOINTS 
                        && (stochastic_crack_index - this_lta.Rs[largest_piece_index] ) >= MINDATAPOINTS){

                    // if(a_node->m_id == 39024) cout << "making doing stoch crack item" << endl;
                //   cout << "creating the pieces\n";
                    // distribute the holes
                    int half_the_holes = (int)(this_lta.holes[largest_piece_index] / 2);
                    // this_lta.Ls[this_lta.how_many_created] = this_piece_L;

                    if(stochastic_crack_index - this_lta.holes[largest_piece_index] - this_lta.Ls[largest_piece_index] >= half_the_holes)
                    {
                        // this_lta.Ls[this_lta.how_many_created] = this_piece_L;
                        // this_lta.Rs[this_lta.how_many_created] = crack_index;
                        // cout << "case 1" << endl;


                        this_lta.Ls[this_lta.how_many_created] = stochastic_crack_index - half_the_holes;
                        this_lta.Rs[this_lta.how_many_created] = this_lta.Rs[largest_piece_index];
                        this_lta.holes[this_lta.how_many_created] = half_the_holes;
                        for(ashghal = 0; ashghal < NUMDIMS; ashghal ++){
                            this_lta.crack_covers_mins[this_lta.how_many_created][ashghal] = r_rect.m_min[ashghal];
                            this_lta.crack_covers_maxes[this_lta.how_many_created][ashghal] = r_rect.m_max[ashghal];
                        }
                        // cout << "new piece, L, R, h: " << this_lta.Ls[this_lta.how_many_created] << " " << this_lta.Rs[this_lta.how_many_created] << " " << this_lta.holes[this_lta.how_many_created] << endl;

                        this_lta.how_many_created++;

                        
                        this_lta.holes[largest_piece_index] = this_lta.holes[largest_piece_index] - half_the_holes;


                        // cout << "copying from " << stochastic_crack_index - half_the_holes << " to " << this_lta.Ls[largest_piece_index] +  this_lta.holes[largest_piece_index] << " this many " << half_the_holes << endl;
                        // moving h/2 from the end of this piece to the start, to leave the rest of the holes to be distributed again
                        memcpy(m_data_arr_mins[this_lta.Ls[largest_piece_index] +  this_lta.holes[largest_piece_index]], m_data_arr_mins[stochastic_crack_index - half_the_holes], NUMDIMS*half_the_holes*sizeof(float));
                        memcpy(m_data_arr_maxes[this_lta.Ls[largest_piece_index] +  this_lta.holes[largest_piece_index]], m_data_arr_maxes[stochastic_crack_index - half_the_holes], NUMDIMS*half_the_holes*sizeof(float));
                        memcpy(m_data_arr_ids + this_lta.Ls[largest_piece_index] +  this_lta.holes[largest_piece_index], m_data_arr_ids + stochastic_crack_index - half_the_holes, half_the_holes*sizeof(int));

                        this_lta.Rs[largest_piece_index] = stochastic_crack_index - half_the_holes;

                        for(ashghal = 0; ashghal < NUMDIMS; ashghal ++){
                            this_lta.crack_covers_mins[largest_piece_index][ashghal] = l_rect.m_min[ashghal];
                            this_lta.crack_covers_maxes[largest_piece_index][ashghal] = l_rect.m_max[ashghal];
                        }

                        // cout << "old piece, L, R, h: " << this_lta.Ls[largest_piece_index] << " " << this_lta.Rs[largest_piece_index] << " " << this_lta.holes[largest_piece_index] << endl;


                        // if(a_node->m_id == 16502) {
                        //     cout << "it was mm = 0, copying from " <<  crack_index - half_the_holes << " to " << this_piece_L + this_piece_holes - half_the_holes << " this many " << half_the_holes << endl;
                        // }

                        // set up ptcrack
                        // this_piece_L = crack_index - half_the_holes;
                        // this_piece_holes = half_the_holes;

                       
                    }
                    else
                    {

                        // cout << "case 2" << endl;

                        this_lta.Ls[this_lta.how_many_created] = stochastic_crack_index - half_the_holes;
                        this_lta.Rs[this_lta.how_many_created] = this_lta.Rs[largest_piece_index];
                        this_lta.holes[this_lta.how_many_created] = half_the_holes;
                        for(ashghal = 0; ashghal < NUMDIMS; ashghal ++){
                            this_lta.crack_covers_mins[this_lta.how_many_created][ashghal] = r_rect.m_min[ashghal];
                            this_lta.crack_covers_maxes[this_lta.how_many_created][ashghal] = r_rect.m_max[ashghal];
                        }
                        // cout << "new piece, L, R, h: " << this_lta.Ls[this_lta.how_many_created] << " " << this_lta.Rs[this_lta.how_many_created] << " " << this_lta.holes[this_lta.how_many_created] << endl;

                        this_lta.how_many_created++;



                        // this_lta.Ls[this_lta.how_many_created] = this_piece_L;
                        int old_holes = this_lta.holes[largest_piece_index];

                        int how_many_to_move = stochastic_crack_index - this_lta.holes[largest_piece_index] - this_lta.Ls[largest_piece_index];
                        this_lta.holes[largest_piece_index] = this_lta.holes[largest_piece_index] - half_the_holes;

                        // cout << "copying from " << stochastic_crack_index - how_many_to_move << " to " << this_lta.Ls[largest_piece_index] +  this_lta.holes[largest_piece_index] << " this many " << how_many_to_move << endl;

                        memcpy(m_data_arr_mins[this_lta.Ls[largest_piece_index] + this_lta.holes[largest_piece_index]], m_data_arr_mins[stochastic_crack_index - how_many_to_move], NUMDIMS*how_many_to_move*sizeof(float));
                        memcpy(m_data_arr_maxes[this_lta.Ls[largest_piece_index] + this_lta.holes[largest_piece_index]], m_data_arr_maxes[stochastic_crack_index - how_many_to_move], NUMDIMS*how_many_to_move*sizeof(float));
                        memcpy(m_data_arr_ids + this_lta.Ls[largest_piece_index] + this_lta.holes[largest_piece_index], m_data_arr_ids + stochastic_crack_index - how_many_to_move, how_many_to_move*sizeof(int));

                        this_lta.Rs[largest_piece_index] = stochastic_crack_index - half_the_holes;
                        // cout << "old piece, L, R, h: " << this_lta.Ls[largest_piece_index] << " " << this_lta.Rs[largest_piece_index] << " " << this_lta.holes[largest_piece_index] << endl;

                        
                        // this_piece_L = crack_index - half_the_holes;
                        // this_piece_holes = half_the_holes;

                    }

                    
                    // end of distribute


                    // // right ro ezafe kon
                    // this_lta.Ls[this_lta.how_many_created] = stochastic_crack_index;
                    // this_lta.Rs[this_lta.how_many_created] = this_lta.Rs[largest_piece_index];
                    // for(ashghal = 0; ashghal < NUMDIMS; ashghal ++){
                    //     this_lta.crack_covers_mins[this_lta.how_many_created][ashghal] = r_rect.m_min[ashghal];
                    //     this_lta.crack_covers_maxes[this_lta.how_many_created][ashghal] = r_rect.m_max[ashghal];
                    // }
                    // this_lta.how_many_created++;


                    // // left ro bezar sar e jay e ooni k shekoondim
                    // this_lta.Rs[largest_piece_index] = stochastic_crack_index;
                    // for(ashghal = 0; ashghal < NUMDIMS; ashghal ++){
                    //     this_lta.crack_covers_mins[largest_piece_index][ashghal] = l_rect.m_min[ashghal];
                    //     this_lta.crack_covers_maxes[largest_piece_index][ashghal] = l_rect.m_max[ashghal];
                    // }

                }

           }
            // cout << "end of stch" << endl;
            // END OF STOCH 
            all_lta->push_back(this_lta);
            // cout << "after puch back" << endl;
            #ifdef timebreak
                cracking_time += (clock() - start_time);
            #endif

        }
    }
    return true; // Continue searching
}




RTREE_TEMPLATE
bool RTREE_QUAL::compareByL(const Branch a, const Branch b){
    if(a.m_child->m_L < b.m_child->m_L) return true;
    return false;
}


RTREE_TEMPLATE
void RTREE_QUAL::Add_ltas_v2(LeavesToAdd_v2 *all_lta){

    for(int i = 0; i < all_lta->size(); i++){

        // cout << "in first for loop, " << i << endl;

        vector<Branch> tba((all_lta->at(i)).how_many_created);
        
        cout  << "how_many_cre " << ((all_lta->at(i)).how_many_created) << endl;
        cout << "tba size: " << tba.size() << endl;
        cout  << "how_many_cre " << ((all_lta->at(i)).how_many_created) << endl;
        Node* this_leaf_older_brother; Node* this_leaf_younger_brother;
        if((all_lta->at(i).this_leaf)->m_parent != NULL){
            Node* parent_of_start = (all_lta->at(i).this_leaf)->m_parent;
            int branch_index = getNodeBranchIndex((all_lta->at(i).this_leaf));
            this_leaf_older_brother = (all_lta->at(i).this_leaf)->m_older_brother;
            this_leaf_younger_brother = (all_lta->at(i).this_leaf)->m_younger_brother;
            cout << "before disconnect branch..." << endl;
            DisconnectBranch(parent_of_start, branch_index);
            cout << "after disconnect branch, before second for loop..." << endl;
            for(int j = 0 ; j < ((all_lta->at(i)).how_many_created); j++){
                cout << "in second loop " << j << endl;
                Branch this_branch;
                this_branch.m_child = AllocNode();
                this_branch.m_child->m_L = (all_lta->at(i)).Ls[j];
                this_branch.m_child->m_R = (all_lta->at(i)).Rs[j];
//                this_branch.m_rect = Rect((all_lta->at(i)).crack_covers_mins[j], (all_lta->at(i)).crack_covers_maxes[j]);

                this_branch.m_child->m_level = 0;
                // 2023
                // in the naive solution to crackign with holes, at this point, ther are no holes, so this coutn would be correct, and the assignment of regularity would also be correct
                this_branch.m_child->m_count = this_branch.m_child->m_R - this_branch.m_child->m_L;
//                if(this_branch.m_child->m_count > MAXDATAPOINTS) this_branch.m_child->m_regular = false;
//                else this_branch.m_child->m_regular = true;

                if(this_branch.m_child->m_count > MAXDATAPOINTS) {
                    this_branch.m_child->m_regular = false;
                    this_branch.m_rect = Rect((all_lta->at(i)).crack_covers_mins[j], (all_lta->at(i)).crack_covers_maxes[j]);
                }
                else {
                    this_branch.m_child->m_regular = true;
                    this_branch.m_rect = LeafCover(this_branch.m_child->m_L, this_branch.m_child->m_R);
                }
//                for(int ashghal = 0; ashghal < NUMDIMS; ashghal++){
//                    if(this_branch.m_rect.m_min[ashghal] > this_branch.m_rect.m_max[ashghal]){
//                        cout << "TOOYE ADDLTA, MIN>MAX " << ashghal << " min " << this_branch.m_rect.m_min[ashghal] << " max " << this_branch.m_rect.m_max[ashghal] << endl;
//                    }
//                }

                tba[j] = this_branch;
                // tba.push_back(this_branch);
                // Insert_anylevel(this_branch, parent_of_start, 1);
            }
            cout << "end of second loop, before sort..." << endl;
            cout << "tba[0, 1].child.L " << tba[0].m_child->m_L << " " << tba[0].m_child->m_R << ", " << tba[1].m_child->m_L << " " << tba[1].m_child->m_R  << endl;

            // TODO: sort them
            // nor dure if ascending or desc
            std::sort(tba.begin(), tba.end());
            // std::sort(tba.begin(), tba.end(), compareByL);
            cout << "after sort, before insertions..." << endl;
            cout << "tba[0, 1].child.L " << tba[0].m_child->m_L << " " << tba[0].m_child->m_R << ", " << tba[1].m_child->m_L << " " << tba[1].m_child->m_R  << endl;

            // TODO: then do another loop and insert them
            // before inserting also add the older and younger sisters now that they are ordered
            for(int j = 0 ; j < ((all_lta->at(i)).how_many_created); j++){
                // TODO: WE NEED TO KEEP THE NEXT AND PREV OF THE ONE WE JUST DOSCONNECTED TO KNOW WHERE TO PUT THE FIRST ONE...
                tba.at(j).m_child->m_younger_brother = this_leaf_younger_brother;
                tba.at(j).m_child->m_older_brother = this_leaf_older_brother;

                if(this_leaf_older_brother == NULL){
                    // is head
                    m_oldest_leaf = tba.at(j).m_child;
                }
                else{
                    this_leaf_older_brother->m_younger_brother = tba.at(j).m_child;
                }

                if(this_leaf_younger_brother == NULL){
                    // is tail
                    m_youngest_leaf = tba.at(j).m_child;
                }
                else{
                    this_leaf_younger_brother->m_older_brother = tba.at(j).m_child;
                }


                Insert_anylevel(tba.at(j), parent_of_start, 1);

                // also, if it is the qpiece, then add it to the overlapping_leaves

                // for the next round:
                this_leaf_older_brother = tba.at(j).m_child;

                // TODO: this is stupid, I keep re-writing the older and younger brothers so this is bad. I will have to optimise it
                // i think it's correct, but not smart.
            }
        }
        else{
            cout << "in add_ltas, leaf was root..." << endl;
            m_root->m_count = 0;
            m_root->m_level++;
            m_root->m_regular = true;
            m_root->m_parent = NULL;

            this_leaf_older_brother = NULL; this_leaf_younger_brother = NULL;

            for(int j = 0 ; j < ((all_lta->at(i)).how_many_created); j++){
                cout << "was root, in loop... " << j << endl;
                Branch this_branch;
                this_branch.m_child = AllocNode();
                this_branch.m_child->m_L = (all_lta->at(i)).Ls[j];
                this_branch.m_child->m_R = (all_lta->at(i)).Rs[j];
//                this_branch.m_rect = Rect((all_lta->at(i)).crack_covers_mins[j], (all_lta->at(i)).crack_covers_maxes[j]);

                this_branch.m_child->m_level = 0;
                this_branch.m_child->m_count = this_branch.m_child->m_R - this_branch.m_child->m_L;
//                if(this_branch.m_child->m_count > MAXDATAPOINTS) this_branch.m_child->m_regular = false;
//                else this_branch.m_child->m_regular = true;
                if(this_branch.m_child->m_count > MAXDATAPOINTS) {
                    this_branch.m_child->m_regular = false;
                    this_branch.m_rect = Rect((all_lta->at(i)).crack_covers_mins[j], (all_lta->at(i)).crack_covers_maxes[j]);
                }
                else {
                    this_branch.m_child->m_regular = true;
                    this_branch.m_rect = LeafCover(this_branch.m_child->m_L, this_branch.m_child->m_R);
                }
//                for(int ashghal = 0; ashghal < NUMDIMS; ashghal++){
//                    if(this_branch.m_rect.m_min[ashghal] > this_branch.m_rect.m_max[ashghal]){
//                        cout << "TOOYE ADDLTA, MIN>MAX " << ashghal << " min " << this_branch.m_rect.m_min[ashghal] << " max " << this_branch.m_rect.m_max[ashghal] << endl;
//                    }
//                }
                tba[j] = this_branch;
                // tba.push_back(this_branch);
                // InsertRect(this_branch, &m_root, 1);
            }
            cout << "after loop, sorting..." << endl;
            // TODO: sort them
            // nor dure if ascending or desc
            // std::sort(tba.begin(), tba.end(), compareByL);
            cout << "tba size: " << tba.size() << endl;
            cout  << " how_many_cre " << ((all_lta->at(i)).how_many_created) << endl;
            cout << "tba[0, 1].child.L " << tba[0].m_child->m_L << " " << tba[0].m_child->m_R << ", " << tba[1].m_child->m_L << " " << tba[1].m_child->m_R << ", " << tba[2].m_child->m_L << " " << tba[2].m_child->m_R << endl;
            std::sort(tba.begin(), tba.end());
            cout << "after sort..." << endl;
            cout << "tba[0, 1].child.L " << tba[0].m_child->m_L << " " << tba[0].m_child->m_R << ", " << tba[1].m_child->m_L << " " << tba[1].m_child->m_R << ", " << tba[2].m_child->m_L << " " << tba[2].m_child->m_R << endl;
            cout << "sorted, insert loop... " << endl;
            // TODO: then do another loop and insert them
            // before inserting also add the older and younger sisters now that they are ordered
            for(int j = 0 ; j < ((all_lta->at(i)).how_many_created); j++){
                cout << "insisde insert loop, " << j << endl;
                // TODO: WE NEED TO KEEP THE NEXT AND PREV OF THE ONE WE JUST DOSCONNECTED TO KNOW WHERE TO PUT THE FIRST ONE...
                tba.at(j).m_child->m_younger_brother = this_leaf_younger_brother;
                tba.at(j).m_child->m_older_brother = this_leaf_older_brother;

                if(this_leaf_older_brother == NULL){
                    // is head
                    cout << "case 1" << endl;
                    m_oldest_leaf = tba.at(j).m_child;
                }
                else{
                    cout << "case 2" << endl;

                    this_leaf_older_brother->m_younger_brother = tba.at(j).m_child;
                }

                if(this_leaf_younger_brother == NULL){
                    cout << "case 3" << endl;

                    // is tail
                    m_youngest_leaf = tba.at(j).m_child;
                }
                else{
                    cout << "case 4" << endl;

                    this_leaf_younger_brother->m_older_brother = tba.at(j).m_child;
                }

                cout << "set brother pointers, inserting..." << endl;
                InsertRect(tba.at(j), &m_root, 1);
                cout << "inserted, setting new oppa..." << endl;

                // for the next round:
                this_leaf_older_brother = tba.at(j).m_child;
                cout << "time for next loop" << endl;

                // TODO: this is stupid, I keep re-writing the older and younger brothers so this is bad. I will have to optimise it
                // i think it's correct, but not smart.
            }

        }
    }


}


/*
// also, if it is the qpiece, then add it to the overlapping_leaves
                if(all_lta->at(i).qp_L == tba.at(j).m_child->m_L){
                    overlapping_leaves->push_back(overlapping_leaf_tbas(tba.at(j).m_child, all_lta->at(i).qp_tba));
                }


// also, if it is the qpiece, then add it to the overlapping_leaves
                if(all_lta->at(i).qp_L == tba.at(j).m_child->m_L){
                    overlapping_leaves->push_back(overlapping_leaf_tbas(tba.at(j).m_child, all_lta->at(i).qp_tba));
                }
*/




RTREE_TEMPLATE
void RTREE_QUAL::Add_ltas_v3(LeavesToAdd_v2 *all_lta, vector<overlapping_leaf_tbas> *overlapping_leaves){

    // cout << "start of add_lta, ol.size() " << overlapping_leaves->size() << endl;
    // if(overlapping_leaves->size() > 0){
    //     cout << "ol[0].leaf.id " << (overlapping_leaves->at(0)).this_leaf->m_id << " ol[0].tba.size() " << (overlapping_leaves->at(0)).this_tba.size() << endl;
    // } 
    
    Node* this_leaf_older_brother; Node* this_leaf_younger_brother;

    Rect trash_rect;
    // int original_L, original_holes;

    for(int i = 0; i < all_lta->size(); i++){

        // cout << "in first for loop, " << i << endl;

        // if((all_lta->at(i)).how_many_created == 1) {
        //     overlapping_leaves->push_back(overlapping_leaf_tbas(all_lta->at(i).this_leaf, all_lta->at(i).tba_mins, all_lta->at(i).tba_maxes, all_lta->at(i).tba_ids, all_lta->at(i).tba_count));
        //     continue;
        // }

        if((all_lta->at(i)).how_many_created == 1) 
        {
            // if(all_lta->at(i).this_leaf->m_id == 2026){
            //     cout << "putting irreg leaf 2026 into ol woth start: " <<  all_lta->at(i).tba_start << endl;
            // }
            overlapping_leaves->push_back(overlapping_leaf_tbas(all_lta->at(i).this_leaf, all_lta->at(i).tba_start, all_lta->at(i).tba_count));
            Rect new_rect; int branch_index;
            Rect trash_rect;
            if((all_lta->at(i).this_leaf)->m_parent != NULL){
                branch_index = getNodeBranchIndex(all_lta->at(i).this_leaf);
                new_rect = ((all_lta->at(i).this_leaf)->m_parent)->m_branch[branch_index].m_rect;
            
                for(int ii = (all_lta->at(i).tba_start); ii < (all_lta->at(i).tba_start) + (all_lta->at(i).tba_count); ii++){
                        trash_rect = (Rect(tba_mins[ii], tba_maxes[ii]));
                        new_rect = CombineRect(&(new_rect), &trash_rect);
                    }

                ((all_lta->at(i).this_leaf)->m_parent)->m_branch[branch_index].m_rect = new_rect;
            }

            continue;
        }

        vector<Branch> tba((all_lta->at(i)).how_many_created);

        int old_id = (all_lta->at(i).this_leaf)->m_id;
        
        // cout  << "how_many_cre " << ((all_lta->at(i)).how_many_created) << endl;
        // cout << "tba size: " << tba.size() << endl;
        // cout  << "how_many_cre " << ((all_lta->at(i)).how_many_created) << endl;
        
        if((all_lta->at(i).this_leaf)->m_parent != NULL){
            Node* parent_of_start = (all_lta->at(i).this_leaf)->m_parent;
            int branch_index = getNodeBranchIndex((all_lta->at(i).this_leaf));
            this_leaf_older_brother = (all_lta->at(i).this_leaf)->m_older_brother;
            this_leaf_younger_brother = (all_lta->at(i).this_leaf)->m_younger_brother;
            // original_L = (all_lta->at(i).this_leaf)->m_L;
            // original_holes = (all_lta->at(i).this_leaf)->m_holes;
            // cout << "before disconnect branch..." << endl;
            DisconnectBranch(parent_of_start, branch_index);
            FreeNode(all_lta->at(i).this_leaf);
            // cout << "after disconnect branch, before second for loop..." << endl;
            for(int j = 0 ; j < ((all_lta->at(i)).how_many_created); j++){
                // cout << "in second loop " << j << endl;
                Branch this_branch;
                this_branch.m_child = AllocNode();
                this_branch.m_child->m_L = (all_lta->at(i)).Ls[j];
                this_branch.m_child->m_R = (all_lta->at(i)).Rs[j];

                // 2023
                this_branch.m_child->m_holes = (all_lta->at(i)).holes[j];
//                this_branch.m_rect = Rect((all_lta->at(i)).crack_covers_mins[j], (all_lta->at(i)).crack_covers_maxes[j]);

                this_branch.m_child->m_level = 0;
                this_branch.m_child->m_count = this_branch.m_child->m_R - this_branch.m_child->m_L - this_branch.m_child->m_holes;
//                if(this_branch.m_child->m_count > MAXDATAPOINTS) this_branch.m_child->m_regular = false;
//                else this_branch.m_child->m_regular = true;

                if(this_branch.m_child->m_count > MAXDATAPOINTS) {
                    this_branch.m_child->m_regular = false;
                    this_branch.m_rect = Rect((all_lta->at(i)).crack_covers_mins[j], (all_lta->at(i)).crack_covers_maxes[j]);
                }
                else {
                    this_branch.m_child->m_regular = true;
                    this_branch.m_rect = LeafCover(this_branch.m_child->m_L, this_branch.m_child->m_R);
                }

                
//                for(int ashghal = 0; ashghal < NUMDIMS; ashghal++){
//                    if(this_branch.m_rect.m_min[ashghal] > this_branch.m_rect.m_max[ashghal]){
//                        cout << "TOOYE ADDLTA, MIN>MAX " << ashghal << " min " << this_branch.m_rect.m_min[ashghal] << " max " << this_branch.m_rect.m_max[ashghal] << endl;
//                    }
//                }
                // let's also extend the MBR here
                if(all_lta->at(i).qp_L == this_branch.m_child->m_L && (all_lta->at(i).tba_count) > 0){
                    for(int ii = (all_lta->at(i).tba_start); ii < (all_lta->at(i).tba_start) + (all_lta->at(i).tba_count); ii++){
                        // if((all_lta->at(i).qp_tba)[ii].m_data == 100352) {
                            // cout << "inserting " << (all_lta->at(i).qp_tba)[ii].m_data << " into irreg leaf " << this_branch.m_child->m_L << " " << this_branch.m_child->m_R << endl;
                        // }
                        trash_rect = (Rect(tba_mins[ii], tba_maxes[ii]));
                        // cout << "trash_rect is "  << trash_rect.m_min[0] << " " << trash_rect.m_min[1] << " " << trash_rect.m_max[0] << " " << trash_rect.m_max[1] << endl;
                        this_branch.m_rect = CombineRect(&(this_branch.m_rect), &trash_rect);
                    }
                    // cout << "cover after adding the tba:" << endl;
                    // for(int ii = 0; ii < NUMDIMS; ii++) cout << this_branch.m_rect.m_min[ii] << " " << this_branch.m_rect.m_max[ii] << " , ";
                    // cout << endl;
                }

                


                tba[j] = this_branch;
                // tba.push_back(this_branch);
                // Insert_anylevel(this_branch, parent_of_start, 1);
            }
            // cout << "end of second loop, before sort..." << endl;
            // cout << "tba[0, 1].child.L " << tba[0].m_child->m_L << " " << tba[0].m_child->m_R << ", " << tba[1].m_child->m_L << " " << tba[1].m_child->m_R  << endl;

            // TODO: sort them
            // nor dure if ascending or desc
            std::sort(tba.begin(), tba.end());

            // for(int o = 0; o < tba.size(); o++){
            //     if(tba[o].m_child->m_id == 16807){
            //         cout << "in add ltas, created leaf 16807, has " << tba[o].m_child->m_holes << " holes, parent had " << (all_lta->at(i).this_leaf)->m_holes << " holes." << endl;
            //         cout << "parent bounds id L, h, R " << (all_lta->at(i).this_leaf)->m_id << " " << (all_lta->at(i).this_leaf)->m_L << " " << (all_lta->at(i).this_leaf)->m_holes << " " << (all_lta->at(i).this_leaf)->m_R << endl;
            //         cout << " this leaf bounds: " << tba[o].m_child->m_L << " " << tba[o].m_child->m_holes << " " << tba[o].m_child->m_R << endl;
            //         cout << "tba.size() " << tba.size() << " other bounds: " << endl;
            //         // for(int p =0; p < tba.size();p++){
            //         //     cout << tba[p].m_child->m_L << " " << tba[p].m_child->m_holes << " " << tba[p].m_child->m_R << endl;
            //         // }
            //         break;
            //     }
            // }

            // if(old_id == 15372 || old_id == 14089){
            //     cout << "object ids in this leaf: .............." << endl;
            //     for(int o = (all_lta->at(i).this_leaf)->m_L + (all_lta->at(i).this_leaf)->m_holes; o < (all_lta->at(i).this_leaf)->m_R; o++){
            //         cout << m_data_arr_ids[o] << endl;
            //     }
            //     cout << "............" << endl;
            //     for(int o = 0; o < tba.size(); o++){
            //         cout << "in add ltas, cracked leaf leaf 14089, has" << endl;
            //         cout << "parent bounds id L, h, R " << (all_lta->at(i).this_leaf)->m_id << " " << (all_lta->at(i).this_leaf)->m_L << " " << (all_lta->at(i).this_leaf)->m_holes << " " << (all_lta->at(i).this_leaf)->m_R << endl;
            //         cout << " this leaf bounds and id : " << tba[o].m_child->m_id << " " << tba[o].m_child->m_L << " " << tba[o].m_child->m_holes << " " << tba[o].m_child->m_R << endl;
            //         cout << "tba.size() " << tba.size() << " other bounds: " << endl;
            //     }
            // }

            

            //2023
            // cout << "hhhhhhhhhhhhh" << endl;
            // for(int u = 0; u < tba.size(); u++){
            //     cout << tba[u].m_child->m_L << " " << tba[u].m_child->m_R << endl;
            // }
            // cout << "hhhhhhhhhhhhh" << endl;
            //2023
            // std::sort(tba.begin(), tba.end(), compareByL);
            // cout << "after sort, before insertions..." << endl;
            // cout << "tba[0, 1].child.L " << tba[0].m_child->m_L << " " << tba[0].m_child->m_R << ", " << tba[1].m_child->m_L << " " << tba[1].m_child->m_R  << endl;

            // TODO: done, not tested
            // naive solution for holes when cracking.
            // add the holes to the first guy
            // now that the pieces are sorted it will be the first piece.
            // ASSERT(tba[0].m_child->m_L == original_L + original_holes);
            // tba[0].m_child->m_L = original_L;
            // tba[0].m_child->m_holes = original_holes;

            // TODO: then do another loop and insert them
            // before inserting also add the older and younger sisters now that they are ordered
            int tba_size;
            // if((all_lta->at(i).qp_tba) == NULL) tba_size = 0;
            // else tba_size = (all_lta->at(i).qp_tba).size();
            tba_size = (all_lta->at(i).tba_count);
            // cout << "set the tba_size " << tba_size << endl;
            // NEW
            int hmc = ((all_lta->at(i)).how_many_created);
            tba[0].m_child->m_older_brother = this_leaf_older_brother;
            if(this_leaf_older_brother != NULL) this_leaf_older_brother->m_younger_brother = tba[0].m_child;
            else m_oldest_leaf = tba.at(0).m_child;
            tba[hmc - 1].m_child->m_younger_brother = this_leaf_younger_brother;
            if(this_leaf_younger_brother != NULL) this_leaf_younger_brother->m_older_brother = tba[hmc- 1].m_child;
            else m_youngest_leaf = tba[hmc-1].m_child;


            for(int j = 0; j < hmc; j++){
                if(j ==0){
                    tba.at(j).m_child->m_younger_brother = tba.at(j + 1).m_child;
                }
                else if(j == hmc - 1){
                    tba.at(j).m_child->m_older_brother = tba.at(j-1).m_child;
                }
                else{
                    tba.at(j).m_child->m_younger_brother = tba.at(j + 1).m_child;
                    tba.at(j).m_child->m_older_brother = tba.at(j-1).m_child;
                }

                ASSERT(tba.at(j).m_child->m_count > 0);

                Insert_anylevel(tba.at(j), parent_of_start, 1);

                // if( (j == 0 && all_lta->at(i).qp_L == ((tba.at(j)).m_child)->m_L + original_holes && (all_lta->at(i).qp_tba).size() > 0)
                //    || (all_lta->at(i).qp_L == ((tba.at(j)).m_child)->m_L &&  (all_lta->at(i).qp_tba).size() > 0)){
                // if( (j == 0 && all_lta->at(i).qp_L == ((tba.at(j)).m_child)->m_L + original_holes)
                //    || (all_lta->at(i).qp_L == ((tba.at(j)).m_child)->m_L)){
                if((all_lta->at(i).qp_L == ((tba.at(j)).m_child)->m_L)){
                    // cout << "in add_ltas, adding to ol, tba.size() " << (all_lta->at(i).qp_tba).size() << endl;
                    // if(tba.at(j).m_child->m_id == 2026){
                    //     cout << "2 putting irreg leaf 2026 into ol woth start: " <<  all_lta->at(i).tba_start << endl;
                    // }
                    overlapping_leaves->push_back(overlapping_leaf_tbas(tba.at(j).m_child, all_lta->at(i).tba_start, all_lta->at(i).tba_count));
                    // cout << "in the same block, after push back to ol, ol[-1].tba.size() " << ((overlapping_leaves->back()).this_tba).size() << endl;                    
                }
            }

            
            // for(int o = 0; o < tba.size(); o++){
            //     if(tba[o].m_child->m_id == 16807){
            //         // cout << "after insertion, in add ltas, created leaf 16807, has " << tba[o].m_child->m_holes << " holes, parent had " << (all_lta->at(i).this_leaf)->m_holes << " holes." << endl;
            //         if(tba[o].m_child->m_younger_brother != NULL) cout << "right bro id: " << (tba[o].m_child->m_younger_brother)->m_id << " " << (tba[o].m_child->m_younger_brother)->m_L << " " << (tba[o].m_child->m_younger_brother)->m_R <<  endl;
            //         // else cout << "right bro NULL" << endl;
            //         if(tba[o].m_child->m_older_brother != NULL) cout << "left bro id: " << (tba[o].m_child->m_older_brother)->m_id << " " <<  (tba[o].m_child->m_older_brother)->m_L << " " << (tba[o].m_child->m_older_brother)->m_R << endl;
            //         else cout << "left bro NULL" << endl;

            //         // for(int p =0; p < tba.size();p++){
            //         //     cout << tba[p].m_child->m_L << " " << tba[p].m_child->m_holes << " " << tba[p].m_child->m_R << endl;
            //         // }
            //         break;
            //     }
            // }
            // END OF NEW
            // OLD
            // for(int j = 0 ; j < ((all_lta->at(i)).how_many_created); j++){
            //     // TODO: WE NEED TO KEEP THE NEXT AND PREV OF THE ONE WE JUST DOSCONNECTED TO KNOW WHERE TO PUT THE FIRST ONE...
            //     tba.at(j).m_child->m_younger_brother = this_leaf_younger_brother;
            //     tba.at(j).m_child->m_older_brother = this_leaf_older_brother;

            //     if(this_leaf_older_brother == NULL){
            //         // is head
            //         m_oldest_leaf = tba.at(j).m_child;
            //     }
            //     else{
            //         this_leaf_older_brother->m_younger_brother = tba.at(j).m_child;
            //     }

            //     if(this_leaf_younger_brother == NULL){
            //         // is tail
            //         m_youngest_leaf = tba.at(j).m_child;
            //     }
            //     else{
            //         this_leaf_younger_brother->m_older_brother = tba.at(j).m_child;
            //     }


            //     Insert_anylevel(tba.at(j), parent_of_start, 1);

            //     // also, if it is the qpiece, then add it to the overlapping_leaves
            //     // cout << "!!!!!!!!!!!!!before the if statement, tba.size() " << (all_lta->at(i).qp_tba)->size() << " tba_size " << tba_size << endl;
            //     //     if(all_lta->at(i).qp_L == this_branch.m_child->m_L && (all_lta->at(i).qp_tba).size() > 0){

            //     if( (j == 0 && all_lta->at(i).qp_L == ((tba.at(j)).m_child)->m_L + original_holes && (all_lta->at(i).qp_tba).size() > 0)
            //        || (all_lta->at(i).qp_L == ((tba.at(j)).m_child)->m_L &&  (all_lta->at(i).qp_tba).size() > 0)){
            //         // cout << "in add_ltas, adding to ol, tba.size() " << (all_lta->at(i).qp_tba).size() << endl;
            //         overlapping_leaves->push_back(overlapping_leaf_tbas(tba.at(j).m_child, all_lta->at(i).qp_tba));
            //         // cout << "in the same block, after push back to ol, ol[-1].tba.size() " << ((overlapping_leaves->back()).this_tba).size() << endl;                    
            //     }

            //     // for the next round:
            //     this_leaf_older_brother = tba.at(j).m_child;

            //     // TODO: this is stupid, I keep re-writing the older and younger brothers so this is bad. I will have to optimise it
            //     // i think it's correct, but not smart.
            // }
            // END OF OLD


        }
        else{
            // cout << "in add_ltas, leaf was root..." << endl;
            m_root->m_count = 0;
            m_root->m_level++;
            m_root->m_regular = true;
            m_root->m_parent = NULL;

            #ifdef stats
                count_irregular_leaves--;
                count_internal_nodes++;
            #endif

            this_leaf_older_brother = NULL; this_leaf_younger_brother = NULL;

            for(int j = 0 ; j < ((all_lta->at(i)).how_many_created); j++){
                // cout << "was root, in loop... " << j << endl;
                Branch this_branch;
                this_branch.m_child = AllocNode();
                this_branch.m_child->m_L = (all_lta->at(i)).Ls[j];
                this_branch.m_child->m_R = (all_lta->at(i)).Rs[j];
//                this_branch.m_rect = Rect((all_lta->at(i)).crack_covers_mins[j], (all_lta->at(i)).crack_covers_maxes[j]);

                // 2023
                this_branch.m_child->m_holes = (all_lta->at(i)).holes[j];

                this_branch.m_child->m_level = 0;
                this_branch.m_child->m_count = this_branch.m_child->m_R - this_branch.m_child->m_L  - this_branch.m_child->m_holes;
//                if(this_branch.m_child->m_count > MAXDATAPOINTS) this_branch.m_child->m_regular = false;
//                else this_branch.m_child->m_regular = true;

                

                if(this_branch.m_child->m_count > MAXDATAPOINTS) {
                    this_branch.m_child->m_regular = false;
                    this_branch.m_rect = Rect((all_lta->at(i)).crack_covers_mins[j], (all_lta->at(i)).crack_covers_maxes[j]);
                }
                else {
                    this_branch.m_child->m_regular = true;
                    this_branch.m_rect = LeafCover(this_branch.m_child->m_L, this_branch.m_child->m_R);
                }
//                for(int ashghal = 0; ashghal < NUMDIMS; ashghal++){
//                    if(this_branch.m_rect.m_min[ashghal] > this_branch.m_rect.m_max[ashghal]){
//                        cout << "TOOYE ADDLTA, MIN>MAX " << ashghal << " min " << this_branch.m_rect.m_min[ashghal] << " max " << this_branch.m_rect.m_max[ashghal] << endl;
//                    }
//                }

                // let's also extend the MBR here
                if(all_lta->at(i).qp_L == this_branch.m_child->m_L && (all_lta->at(i).tba_count) > 0){
                    Rect trash_rect;
                    for(int ii = (all_lta->at(i).tba_start); ii < (all_lta->at(i).tba_start) + (all_lta->at(i).tba_count); ii++){
                        trash_rect = (Rect(tba_mins[ii], tba_maxes[ii]));
                        this_branch.m_rect = CombineRect(&(this_branch.m_rect), &trash_rect);
                    }
                }
                tba[j] = this_branch;
                // tba.push_back(this_branch);
                // InsertRect(this_branch, &m_root, 1);
            }
            // TODO: sort them
            // nor dure if ascending or desc
            // std::sort(tba.begin(), tba.end(), compareByL);
            std::sort(tba.begin(), tba.end());

            // TODO: naive solution to crackign with holes
            // tba[0].m_child->m_holes = tba[0].m_child->m_L;
            // tba[0].m_child->m_L = 0;

            //2023
            // cout << "hhhhhhhhhhhhh" << endl;
            // for(int u = 0; u < tba.size(); u++){
            //     cout << tba[u].m_child->m_L << " " << tba[u].m_child->m_R << endl;
            // }
            // cout << "hhhhhhhhhhhhh" << endl;
            //2023
            
            // TODO: then do another loop and insert them
            // before inserting also add the older and younger sisters now that they are ordered
            // int tba_size;
            // if((all_lta->at(i).qp_tba) == NULL) tba_size = 0;
            // else tba_size = (all_lta->at(i).qp_tba).size();
            // NEW
            int hmc = ((all_lta->at(i)).how_many_created);
            tba[0].m_child->m_older_brother = this_leaf_older_brother;
            if(this_leaf_older_brother != NULL) this_leaf_older_brother->m_younger_brother = tba[0].m_child;
            else m_oldest_leaf = tba.at(0).m_child;
            tba[hmc - 1].m_child->m_younger_brother = this_leaf_younger_brother;
            if(this_leaf_younger_brother != NULL) this_leaf_younger_brother->m_older_brother = tba[hmc- 1].m_child;
            else m_youngest_leaf = tba[hmc-1].m_child;

            for(int j = 0; j < hmc; j++){
                if(j ==0){
                    tba.at(j).m_child->m_younger_brother = tba.at(j + 1).m_child;
                }
                else if(j == hmc - 1){
                    tba.at(j).m_child->m_older_brother = tba.at(j-1).m_child;
                }
                else{
                    tba.at(j).m_child->m_younger_brother = tba.at(j + 1).m_child;
                    tba.at(j).m_child->m_older_brother = tba.at(j-1).m_child;
                }
                ASSERT(tba.at(j).m_child->m_count > 0);
                InsertRect(tba.at(j), &m_root, 1);

                // if( (j == 0 && all_lta->at(i).qp_L == ((tba.at(j)).m_child)->m_L + original_holes && (all_lta->at(i).qp_tba).size() > 0)
                //    || (all_lta->at(i).qp_L == ((tba.at(j)).m_child)->m_L &&  (all_lta->at(i).qp_tba).size() > 0)){
                // if( (j == 0 && all_lta->at(i).qp_L == ((tba.at(j)).m_child)->m_L + original_holes)
                //    || (all_lta->at(i).qp_L == ((tba.at(j)).m_child)->m_L)){
                if( (all_lta->at(i).qp_L == ((tba.at(j)).m_child)->m_L)){
                    // cout << "in add_ltas, adding to ol, tba.size() " << (all_lta->at(i).qp_tba).size() << endl;
                    overlapping_leaves->push_back(overlapping_leaf_tbas(tba.at(j).m_child, all_lta->at(i).tba_start, all_lta->at(i).tba_count));
                    // cout << "in the same block, after push back to ol, ol[-1].tba.size() " << ((overlapping_leaves->back()).this_tba).size() << endl;                    
                }
            }
            ASSERT(m_root->m_count == hmc);
            // END OF NEW



            // OLD
            // tba_size = (all_lta->at(i).qp_tba).size();
            // for(int j = 0 ; j < ((all_lta->at(i)).how_many_created); j++){
            //     // TODO: WE NEED TO KEEP THE NEXT AND PREV OF THE ONE WE JUST DOSCONNECTED TO KNOW WHERE TO PUT THE FIRST ONE...
            //     tba.at(j).m_child->m_younger_brother = this_leaf_younger_brother;
            //     tba.at(j).m_child->m_older_brother = this_leaf_older_brother;

            //     if(this_leaf_older_brother == NULL){
            //         // is head
            //         // cout << "case 1" << endl;
            //         m_oldest_leaf = tba.at(j).m_child;
            //     }
            //     else{
            //         // cout << "case 2" << endl;

            //         this_leaf_older_brother->m_younger_brother = tba.at(j).m_child;
            //     }

            //     if(this_leaf_younger_brother == NULL){
            //         // cout << "case 3" << endl;

            //         // is tail
            //         m_youngest_leaf = tba.at(j).m_child;
            //     }
            //     else{
            //         // cout << "case 4" << endl;

            //         this_leaf_younger_brother->m_older_brother = tba.at(j).m_child;
            //     }

            //     // cout << "set brother pointers, inserting..." << endl;
            //     InsertRect(tba.at(j), &m_root, 1);
            //     // cout << "inserted, setting new oppa..." << endl;

            //     // also, if it is the qpiece, then add it to the overlapping_leaves
            //     // if(all_lta->at(i).qp_L == tba.at(j).m_child->m_L  &&  tba_size > 0)
            //     if( (j == 0 && all_lta->at(i).qp_L == ((tba.at(j)).m_child)->m_L + original_holes && (all_lta->at(i).qp_tba).size() > 0)
            //        || (all_lta->at(i).qp_L == ((tba.at(j)).m_child)->m_L &&  (all_lta->at(i).qp_tba).size() > 0))
            //     {
            //         overlapping_leaves->push_back(overlapping_leaf_tbas(tba.at(j).m_child, all_lta->at(i).qp_tba));
            //     }

            //     // for the next round:
            //     this_leaf_older_brother = tba.at(j).m_child;
            //     // cout << "time for next loop" << endl;

            //     // TODO: this is stupid, I keep re-writing the older and younger brothers so this is bad. I will have to optimise it
            //     // i think it's correct, but not smart.
            // }
            // END OF OLD

        }
    }

    // cout << "end of add_lta, ol.size() " << overlapping_leaves->size() << endl;
    // if(overlapping_leaves->size() > 0){
    //     cout << "ol[0].leaf.id " << (overlapping_leaves->at(0)).this_leaf->m_id << " ol[0].tba.size() " << (overlapping_leaves->at(0)).this_tba.size() << endl;
    // } 

}



RTREE_TEMPLATE
void RTREE_QUAL::deleteEmptyNodesUp(Node *start){
    // cout << "$$$$$$$$$$" << endl;
    if(start == NULL){
        // cout << "start was null..." << endl;
        return;
    }
    // cout << "in deleteEmptyNodes, start->id " << start->m_id << " count: " << start->m_count << " lvl: " << start->m_level << endl;
    
    if(start->m_parent == NULL){
        // is the root. we don't do anything.
        // cout << "It's the root. It's empty. It really shouldn't be..." << endl;
        return;
    }
    if(start->m_count > 0){
        return;
    }
    
    // it is empty
    Node* parent_of_start = start->m_parent;
    int branch_index = getNodeBranchIndex(start);
    // cout << "after branch_index " << branch_index << endl;
    DisconnectBranch(parent_of_start, branch_index);
    FreeNode(start);
    // cout << "after disconnect" << endl;
    deleteEmptyNodesUp(parent_of_start);
}




RTREE_TEMPLATE
bool RTREE_QUAL::is_leaf_in_OL(Node *a_leaf, vector<overlapping_leaf_tbas> *overlapping_leaves){
    // cout << "in is-leaf-in-tba " << endl;
    for(int folan = 0; folan < overlapping_leaves->size(); folan++){
        // cout << "is-leaf-in-tba loop, " << folan << endl;
        // cout << "a leaf id " << a_leaf->m_id << endl;
        // cout << "folan id " << (overlapping_leaves->at(folan)).this_leaf->m_id << endl;
        if(a_leaf->m_id == (overlapping_leaves->at(folan)).this_leaf->m_id ) return true;
    }
    return false;
}


RTREE_TEMPLATE
void RTREE_QUAL::shuffle_right(vector<overlapping_leaf_tbas> *overlapping_leaves, Node *this_leaf, int this_tba_start, int k){

    if(k <= 0) {
        // cout << "k was less than 0" << endl;
        return;
    }

    // cout << "shuffle right on " << this_leaf->m_id << " L, h, R: " << this_leaf->m_L << " " << this_leaf->m_holes << " " << this_leaf->m_R << endl;
    // if(this_leaf->m_id == 24867){cout << "shuffle right called on 24867!!!!!!!!!!" << endl;}
    // if(this_leaf->m_id == 1192){cout << "shuffle right called on 1192!!!!!!!!!!" << endl;}

    int folan, filan;

    if(k <= this_leaf->m_holes){
        // cout << "case 23" << endl;
        // \textbf{arr}[this-leaf.L: this-leaf.L + k] = this-tbi\;
        // for(int ii = 0; ii < k; ii++){
        //     cout << "tba " << this_tba_ids[ii] << " : ";
        //     for(int jj = 0; jj < NUMDIMS; jj++) cout << this_tba_mins[ii][jj] << " " << this_tba_maxes[ii][jj] << " , ";
        //     cout << endl;
        // }
        // for(int ii = 0; ii < k; ii++){
        //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
        //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
        //     cout << endl;
        // }
        memcpy(m_data_arr_mins[this_leaf->m_L], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
        memcpy(m_data_arr_maxes[this_leaf->m_L], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
        memcpy(m_data_arr_ids + this_leaf->m_L, tba_ids + this_tba_start, k*sizeof(int));

        // for(folan =0; folan< k; folan++){
        //     memcpy(m_data_arr_mins[this_leaf->m_L + folan], this_tba_mins + folan*NUMDIMS, NUMDIMS*sizeof(float));
        //     memcpy(m_data_arr_maxes[this_leaf->m_L + folan], this_tba_maxes + folan*NUMDIMS, NUMDIMS*sizeof(float));
        //     // memcpy(m_data_arr_ids + this_leaf->m_L, this_tba_ids, k*sizeof(int));
        // }

        // for(int ii = 0; ii < k; ii++){
        //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
        //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
        //     cout << endl;
        // }
        // for(folan = 0; folan < k; folan++){
        //     for(filan = 0; filan < NUMDIMS; filan++){
        //         m_data_arr_mins[this_leaf->m_L + folan][filan] = this_tba_mins[folan][filan];
        //         m_data_arr_maxes[this_leaf->m_L + folan][filan] = this_tba_maxes[folan][filan];
        //     }
        //     m_data_arr_ids[this_leaf->m_L + folan] =  this_tba_ids[folan];
        // }
        // this-leaf.L += k
        // this-leaf.left-sibling.R+= k\;
        // this-leaf.holes -= k\;
        // return\;
        this_leaf->m_L += k;
        (this_leaf->m_older_brother)->m_R += k;
        this_leaf->m_holes -= k;
        this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
        return;
    }

    else if(this_leaf->m_id == m_youngest_leaf->m_id){
        // cout << "case 24" << endl;


        if(k > (this_leaf->m_R - this_leaf->m_L)){
            // \textbf{arr}[this-leaf.L + k: this-leaf.R + k - this-leaf.holes] = \textbf{arr}[this-leaf.L + this-leaf.holes: this-leaf.R]\;
            memcpy(m_data_arr_mins[this_leaf->m_L + k], m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_L + k], m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_data_arr_ids + this_leaf->m_L + k, m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes, (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(int));


            // for(folan = 0; folan < (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes); folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_data_arr_mins[this_leaf->m_L + k + folan][filan] = m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes + folan][filan];
            //         m_data_arr_maxes[this_leaf->m_L + k + folan][filan] = m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes + folan][filan];
            //     }
            //     m_data_arr_ids[this_leaf->m_L + k + folan] =  m_data_arr_ids[this_leaf->m_L + this_leaf->m_holes + folan];
            // }
            
        }
        else{
            // % holes < k < size
            // \textbf{arr}[this-leaf.R: this-leaf.R + k - this-leaf.holes] = \textbf{arr}[this-leaf.L + this-leaf.holes: this-leaf.L + k]\;
            memcpy(m_data_arr_mins[this_leaf->m_R], m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(k - this_leaf->m_holes)*sizeof(float));            
            memcpy(m_data_arr_maxes[this_leaf->m_R], m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(k - this_leaf->m_holes)*sizeof(float));            
            memcpy(m_data_arr_ids + this_leaf->m_R, m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes, (k - this_leaf->m_holes)*sizeof(int));


            // for(folan = 0; folan < (k - this_leaf->m_holes); folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_data_arr_mins[this_leaf->m_R + folan][filan] = m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes + folan][filan];
            //         m_data_arr_maxes[this_leaf->m_R + folan][filan] = m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes + folan][filan];
            //     }
            //     m_data_arr_ids[this_leaf->m_R + folan] =  m_data_arr_ids[this_leaf->m_L + this_leaf->m_holes + folan];
            // }
            // 
            

        }
        
        // \textbf{arr}[this-leaf.L: this-leaf.L + k] = this-tbi\;
        // for(int ii = 0; ii < k; ii++){
        //     cout << "tba " << this_tba_ids[ii] << " : ";
        //     for(int jj = 0; jj < NUMDIMS; jj++) cout << this_tba_mins[ii][jj] << " " << this_tba_maxes[ii][jj] << " , ";
        //     cout << endl;
        // }
        // for(int ii = 0; ii < k; ii++){
        //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
        //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
        //     cout << endl;
        // }
        memcpy(m_data_arr_mins[this_leaf->m_L], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
        memcpy(m_data_arr_maxes[this_leaf->m_L], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
        memcpy(m_data_arr_ids + this_leaf->m_L, tba_ids + this_tba_start, k*sizeof(int));

        // for(folan = 0; folan < k; folan++){
        //     memcpy(m_data_arr_mins[this_leaf->m_L +  folan], this_tba_mins + folan*NUMDIMS, NUMDIMS*sizeof(float));
        //     memcpy(m_data_arr_maxes[this_leaf->m_L + folan], this_tba_maxes + folan*NUMDIMS, NUMDIMS*sizeof(float));
        // }

        // for(int ii = 0; ii < k; ii++){
        //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
        //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
        //     cout << endl;
        // }
        // for(folan = 0; folan < k; folan++){
        //     for(filan = 0; filan < NUMDIMS; filan++){
        //         m_data_arr_mins[this_leaf->m_L + folan][filan] = this_tba_mins[folan][filan];
        //         m_data_arr_maxes[this_leaf->m_L + folan][filan] = this_tba_maxes[folan][filan];
        //     }
        //     m_data_arr_ids[this_leaf->m_L + folan] =  this_tba_ids[folan];
        // }
        // this-leaf.L += k\;
        this_leaf->m_L += k;
        // this-leaf.R += k - this-leaf.holes\;
        this_leaf->m_R += (k - this_leaf->m_holes);
        // this-leaf.left-sibling.R += k\;
        (this_leaf->m_older_brother)->m_R += k;
        // this-leaf.holes = 0\;
        this_leaf->m_holes = 0;
        this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
        // return\;
        return;
    }
    // if this leaf is not in overlapping leaves
    // if(find(this_tba->begin(), this_tba->end(), this_leaf) == this_tba->end()){
    else if(!is_leaf_in_OL(this_leaf, overlapping_leaves)){

        // cout << "case 14" << endl;
        if(k >= (this_leaf->m_R - this_leaf->m_L)){
            // cout << "case 15" << endl;
            // pending-insertions += \textbf{arr}[this-leaf.L + this-leaf.holes: this-leaf.R]\;
            memcpy(m_pending_insertions_mins[m_pending_insertions_count], m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_pending_insertions_maxes[m_pending_insertions_count], m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_pending_insertions_ids + m_pending_insertions_count, m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes, (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(int));


            // for(folan = 0; folan < (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes); folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_pending_insertions_mins[m_pending_insertions_count + folan][filan] = m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes + folan][filan];
            //         m_pending_insertions_maxes[m_pending_insertions_count + folan][filan] = m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes + folan][filan];
            //     }
            //     m_pending_insertions_ids[m_pending_insertions_count + folan] = m_data_arr_ids[this_leaf->m_L + this_leaf->m_holes + folan];
            //     // if(m_pending_insertions_ids[m_pending_insertions_count + folan] == 100352){
            //     //     cout << "PUTTING 100351 INTO PENDING!!!!!!!!!!!!!" << endl;
            //     // }
            //     // cout << "putting into pending " << m_pending_insertions_ids[m_pending_insertions_count + folan]  << endl;
                
            // }
            m_pending_insertions_count += (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes);

            //  \textbf{arr}[this-leaf.L: this-leaf.R] = this-tbi[:(this-leaf.R - this-leaf.L)]
            // for(int ii = 0; ii < k; ii++){
            //     cout << "tba " << this_tba_ids[ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << this_tba_mins[ii][jj] << " " << this_tba_maxes[ii][jj] << " , ";
            //     cout << endl;
            // }
            // for(int ii = 0; ii < k; ii++){
            //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
            //     cout << endl;
            // }
            memcpy(m_data_arr_mins[this_leaf->m_L], tba_mins[this_tba_start], (this_leaf->m_R - this_leaf->m_L)*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_L], tba_maxes[this_tba_start], (this_leaf->m_R - this_leaf->m_L)*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_ids + this_leaf->m_L, tba_ids + this_tba_start, (this_leaf->m_R - this_leaf->m_L)*sizeof(int));

            // for(folan = 0; folan < (this_leaf->m_R - this_leaf->m_L); folan++){
            //     memcpy(m_data_arr_mins[this_leaf->m_L + folan], this_tba_mins + folan*NUMDIMS, NUMDIMS*sizeof(float));
            //     memcpy(m_data_arr_maxes[this_leaf->m_L + folan], this_tba_maxes+ folan*NUMDIMS, NUMDIMS*sizeof(float));
            // }

            // for(int ii = 0; ii < k; ii++){
            //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
            //     cout << endl;
            // }
            // for(folan = 0; folan < (this_leaf->m_R - this_leaf->m_L); folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_data_arr_mins[this_leaf->m_L + folan][filan] = this_tba_mins[folan][filan];
            //         m_data_arr_maxes[this_leaf->m_L + folan][filan] = this_tba_maxes[folan][filan];
            //     }
            //     m_data_arr_ids[this_leaf->m_L + folan] =  this_tba_ids[folan];
            // }

            // this-leaf.left-sibling.right-sibling = this-leaf.right-sibling\;
            // this-leaf.right-sibling.left-sibling = this-leaf.leaf-sibling\;

            // (this_leaf->m_older_brother)->m_younger_brother = this_leaf->m_younger_brother;
            // (this_leaf->m_younger_brother)->m_older_brother = this_leaf->m_older_brother;

            // this-leaf.left-sibling.R = this-leaf.R\;
            (this_leaf->m_older_brother)->m_R = this_leaf->m_R;

            Node *next_leaf = this_leaf->m_younger_brother; 
            int this_leaf_size = this_leaf->m_R - this_leaf->m_L;

            // remove this leaf from parent
            int branch_index = getNodeBranchIndex(this_leaf);
            Node* parent = this_leaf->m_parent;
            DisconnectBranch(parent, branch_index);
            FreeNode(this_leaf);
            deleteEmptyNodesUp(parent);

            // TODO: think about this, I do not remember exactly what disconnct branch does
            // shuffle_right(overlapping_leaves, next_leaf, this_tba[std::slice(this_leaf_size, k - this_leaf_size, 1)], k - this_leaf_size);
            // vector<Branch> stupid(this_tba.begin() + this_leaf_size, this_tba.end());
            // shuffle_right(overlapping_leaves, next_leaf, this_tba_mins + NUMDIMS * this_leaf_size, this_tba_maxes + NUMDIMS * this_leaf_size, this_tba_ids + NUMDIMS * this_leaf_size, k - this_leaf_size);
            shuffle_right(overlapping_leaves, next_leaf, this_tba_start + this_leaf_size, k - this_leaf_size);
        }
        else{
            // cout << "case 16" << endl;
            // % holes < k < size
            // pending-insertions += \textbf{arr}[this-leaf.L + this-leaf.holes: this-leaf.L + k]\;
            // cout << "m_p_count: " << m_pending_insertions_count << endl;
            // for(int ii = 0; ii < (k - this_leaf->m_holes); ii++){
            //     cout << m_data_arr_ids[this_leaf->m_L + this_leaf->m_holes + ii] << " " << m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes + ii][0] << " " << m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes + ii][0];
            //     cout << ", " << m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes + ii][1] << " " << m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes + ii][1] << endl;
            // }
            memcpy(m_pending_insertions_mins[m_pending_insertions_count], m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(k - this_leaf->m_holes)*sizeof(float));
            memcpy(m_pending_insertions_maxes[m_pending_insertions_count], m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(k - this_leaf->m_holes)*sizeof(float));
            memcpy(m_pending_insertions_ids + m_pending_insertions_count, m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes, (k - this_leaf->m_holes)*sizeof(int));
            // cout << "after" << endl;
            // for(int ii = 0; ii < (k - this_leaf->m_holes); ii++){
            //     cout << m_pending_insertions_ids[m_pending_insertions_count + ii] << " " << m_pending_insertions_mins[m_pending_insertions_count + ii][0] << " " << m_pending_insertions_maxes[m_pending_insertions_count+ ii][0];
            //     cout << ", " << m_pending_insertions_mins[m_pending_insertions_count+ ii][1] << " " << m_pending_insertions_maxes[m_pending_insertions_count + ii][1] << endl;
            // }


            // for(folan = 0; folan < (k - this_leaf->m_holes); folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_pending_insertions_mins[m_pending_insertions_count + folan][filan] = m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes + folan][filan];
            //         m_pending_insertions_maxes[m_pending_insertions_count + folan][filan] = m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes + folan][filan];
            //     }
            //     m_pending_insertions_ids[m_pending_insertions_count + folan] = m_data_arr_ids[this_leaf->m_L + this_leaf->m_holes + folan];
            //     // if(m_pending_insertions_ids[m_pending_insertions_count + folan] == 100352){
            //     //     cout << "PUTTING 100351 INTO PENDING!!!!!!!!!!!!!" << endl;
            //     // }
            //     // cout << "putting into pending " << m_pending_insertions_ids[m_pending_insertions_count + folan]  << endl;
                
            // }
            m_pending_insertions_count += (k - this_leaf->m_holes);
            // \textbf{arr}[this-leaf.L: this-leaf.L + k] = this-tbi\;
            // cout << "after replacement ____________________" << endl;
            // for(folan = 0; folan < k; folan++){
            //     cout << "tba has item " << this_tba_ids[folan] << "  with cover ";
            //     for(filan = 0; filan < NUMDIMS; filan++) cout << this_tba_mins[folan][filan] << " " << this_tba_maxes[folan][filan] << ", ";
            //     cout << endl; 
            // }
            memcpy(m_data_arr_mins[this_leaf->m_L], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_L], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_ids + this_leaf->m_L, tba_ids + this_tba_start, k*sizeof(int));
            // for(folan = 0; folan < k; folan++){
            //     memcpy(m_data_arr_mins[this_leaf->m_L + folan], this_tba_mins+folan*NUMDIMS, NUMDIMS*sizeof(float));
            //     memcpy(m_data_arr_maxes[this_leaf->m_L + folan], this_tba_maxes+folan*NUMDIMS, NUMDIMS*sizeof(float));
            // }
            // for(folan = 0; folan < k; folan++){
            //     cout << "just copied item " << m_data_arr_ids[this_leaf->m_L + folan] << "  with cover ";
            //     for(filan = 0; filan < NUMDIMS; filan++) cout << m_data_arr_mins[this_leaf->m_L + folan][filan] << " " << m_data_arr_maxes[this_leaf->m_L + folan][filan] << ", ";
            //     cout << endl; 
            //     //  cout << "float " << sizeof(float) << " k = " << k << " m_data_arr_mins[0] " << sizeof(m_data_arr_mins[0]) << " tba_mins[0] " << sizeof(this_tba_mins[0]) << " moving " << k*sizeof(m_data_arr_mins[0]) << endl;
            // }
            // for(folan = 0; folan < k; folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_data_arr_mins[this_leaf->m_L + folan][filan] = this_tba_mins[folan][filan];
            //         m_data_arr_maxes[this_leaf->m_L + folan][filan] = this_tba_maxes[folan][filan];
            //     }
            //     m_data_arr_ids[this_leaf->m_L + folan] =  this_tba_ids[folan];
            //     // cout << m_data_arr_ids[this_leaf->m_L + folan] << endl;
            // }
            // cout << "____________________" << endl;
            // this-leaf.L += k\;
            this_leaf->m_L += k;
            // cout << "changed leaf id: " << this_leaf->m_id << " L,h, R: " << this_leaf->m_L << " " << this_leaf->m_holes << " " << this_leaf->m_R << endl;
            // this-leaf.left-sibling.R += k\;
            (this_leaf->m_older_brother)->m_R += k;
            // cout << "changed leaf id: " << (this_leaf->m_older_brother)->m_id << " L,h, R: " << (this_leaf->m_older_brother)->m_L << " " << (this_leaf->m_older_brother)->m_holes << " " << (this_leaf->m_older_brother)->m_R << endl;
            // cout << "the id of this-leaf.olderbrother.youngerbrother " << ((this_leaf->m_older_brother)->m_younger_brother)->m_id << " this.older.older " << ((this_leaf->m_older_brother)->m_older_brother)->m_id << endl;

            // this-leaf.holes = 0\;
            this_leaf->m_holes = 0;
            this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
            // return\;
            // if leaf is regular, recalculate the MBR
            return;
        }
    
    }
    else{
        // cout << "case 25" << endl;

        if(k >= (this_leaf->m_R - this_leaf->m_L)){
            // cout << "case 17" << endl;
            // next-leaf = this-leaf.right-sibling\;
            Node * next_leaf = this_leaf->m_younger_brother;
            // this-leaf.left-sibling.R += (this-leaf.R - this-leaf.L)\;
            (this_leaf->m_older_brother)->m_R += (this_leaf->m_R - this_leaf->m_L);
            // Move this-leaf to end of \textbf{arr} and fix all the pointers\;
            int this_old_L = this_leaf->m_L, this_old_R = this_leaf->m_R;
            int last_filled_spot_in_arr = m_youngest_leaf->m_R;

            memcpy(m_data_arr_mins[last_filled_spot_in_arr], m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_data_arr_maxes[last_filled_spot_in_arr], m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_data_arr_ids + last_filled_spot_in_arr, m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes, (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(int));


            // for(folan = 0; folan < (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes); folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_data_arr_mins[last_filled_spot_in_arr + folan][filan] = m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes + folan][filan];
            //         m_data_arr_maxes[last_filled_spot_in_arr + folan][filan] = m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes + folan][filan];

            //     }
            //     m_data_arr_ids[last_filled_spot_in_arr + folan] = m_data_arr_ids[this_leaf->m_L + this_leaf->m_holes + folan];
            // }

            folan = (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes);
            this_leaf->m_L = last_filled_spot_in_arr;
            this_leaf->m_R = last_filled_spot_in_arr + folan;

            
            (this_leaf->m_older_brother)->m_younger_brother = this_leaf->m_younger_brother;
            (this_leaf->m_younger_brother)->m_older_brother = this_leaf->m_older_brother;

            this_leaf->m_older_brother = m_youngest_leaf;
            m_youngest_leaf->m_younger_brother = this_leaf;
            this_leaf->m_younger_brother = NULL;
            this_leaf->m_holes = 0;
            this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
            m_youngest_leaf = this_leaf;

            // \textbf{arr}[this-leaf.L: this-leaf.R] = this-tbi[:(this-leaf.R - this-leaf.L)]\;
            // \CommentSty{// then the rest of tbi still belongs to the next-leaf.left-sib, cause we removed the one in between(moved to the end}\;

            // for(int ii = 0; ii < this_old_R - this_old_L; ii++){
            //     cout << "tba " << this_tba_ids[ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << this_tba_mins[ii][jj] << " " << this_tba_maxes[ii][jj] << " , ";
            //     cout << endl;
            // }
            // for(int ii = 0; ii < k; ii++){
            //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
            //     cout << endl;
            // }
            memcpy(m_data_arr_mins[this_old_L], tba_mins[this_tba_start], (this_old_R - this_old_L)*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_old_L], tba_maxes[this_tba_start], (this_old_R - this_old_L)*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_ids + this_old_L, tba_ids + this_tba_start, (this_old_R - this_old_L)*sizeof(int));
            // for(folan = 0; folan < (this_old_R - this_old_L); folan++){
            //     memcpy(m_data_arr_mins[this_old_L + folan], this_tba_mins + folan*NUMDIMS, NUMDIMS*sizeof(float));
            //     memcpy(m_data_arr_maxes[this_old_L + folan], this_tba_maxes + folan*NUMDIMS, NUMDIMS*sizeof(float));
            // }
            // for(int ii = 0; ii < this_old_R - this_old_L; ii++){
            //     cout << "copied " << m_data_arr_ids[this_old_L + ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_old_L+ ii][jj] << " " << m_data_arr_maxes[this_old_L + ii][jj] << " , ";
            //     cout << endl;
            // }
            // for(folan = 0; folan < (this_old_R - this_old_L); folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_data_arr_mins[this_old_L + folan][filan] = this_tba_mins[folan][filan];
            //         m_data_arr_maxes[this_old_L + folan][filan] = this_tba_maxes[folan][filan];
            //     }
            //     m_data_arr_ids[this_old_L + folan] =  this_tba_ids[folan];
            // }
            // \textbf{\name right shuffle}(overlapping-leaves, next-leaf, this-tbi[(this-leaf.R - this-leaf.L):], (k - (this-leaf.R - this-leaf.L)) )\;
            // shuffle_right(overlapping_leaves, next_leaf, this_tba[std::slice((this_old_R - this_old_L), k - (this_old_R - this_old_L), 1)], k - (this_old_R - this_old_L));
            // vector<Branch> stupid(this_tba.begin() + (this_old_R - this_old_L), this_tba.end());
            shuffle_right(overlapping_leaves, next_leaf, this_tba_start + this_old_R - this_old_L, k - (this_old_R - this_old_L));
        }
        else{
        // cout << "case 162" << endl;
        // % holes < k < size
        // tmp = \textbf{arr}[this-leaf.L + this-leaf.holes: this-leaf.L + k]\;
        // vector<Branch> temp(k - this_leaf->m_holes);
        // float temp_mins[k - this_leaf->m_holes][NUMDIMS];
        // float temp_maxes[k - this_leaf->m_holes][NUMDIMS];
        // int temp_ids[k - this_leaf->m_holes];

        float* temp_mins;
        temp_mins = (float*)malloc((k - this_leaf->m_holes) * NUMDIMS * sizeof(float));

        float* temp_maxes;
        temp_maxes = (float*)malloc((k - this_leaf->m_holes)  * NUMDIMS* sizeof(float));

        // for(folan = 0; folan < k - this_leaf->m_holes; folan++){
        //     temp_mins[folan] = (float*)malloc(NUMDIMS*sizeof(float));
        //     temp_maxes[folan] = (float*)malloc(NUMDIMS*sizeof(float));
        // }

        int* temp_ids = (int*)malloc((k - this_leaf->m_holes) * sizeof(int));




        // cout << " after malloc " << endl;
        memcpy(temp_mins,  m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes], (k - this_leaf->m_holes)*NUMDIMS*sizeof(float));
        memcpy(temp_maxes,  m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes], (k - this_leaf->m_holes)*NUMDIMS*sizeof(float));
        memcpy(temp_ids, m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes, (k - this_leaf->m_holes)*sizeof(int));
        // for(folan = 0; folan < (k - this_leaf->m_holes); folan++){
        //     memcpy(temp_mins+folan * NUMDIMS,  m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes +  folan], NUMDIMS*sizeof(float));
        //     memcpy(temp_maxes + folan * NUMDIMS,  m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes +  folan], NUMDIMS*sizeof(float));
        // }
        // for(folan = 0; folan < (k - this_leaf->m_holes); folan++){
        //     //  temp[folan] = Branch(NULL, Rect(m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes + folan], m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes + folan]), m_data_arr_ids[this_leaf->m_L + this_leaf->m_holes + folan]);

        //      for(int filan = 0; filan < NUMDIMS; filan++){
        //         temp_mins[folan][filan] = m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes + folan][filan];
        //         temp_maxes[folan][filan] = m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes + folan][filan];

        //      }
        //      temp_ids[folan] = m_data_arr_ids[this_leaf->m_L + this_leaf->m_holes + folan];
        //      // TODO: I have no cle if this will actually work
        // }
        // \textbf{arr}[this-leaf.L: this-leaf.L + k] = this-tbi\;
        // cout << "---------- k = " << k << endl;
        // for(int ii = 0; ii < k; ii++){
        //     cout << "tba " << this_tba_ids[ii] << " : ";
        //     for(int jj = 0; jj < NUMDIMS; jj++) cout << this_tba_mins[ii][jj] << " " << this_tba_maxes[ii][jj] << " , ";
        //     cout << endl;
        // }
        // for(int ii = 0; ii < k; ii++){
        //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
        //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
        //     cout << endl;
        // }
        memcpy(m_data_arr_mins[this_leaf->m_L], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
        memcpy(m_data_arr_maxes[this_leaf->m_L], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
        memcpy(m_data_arr_ids + this_leaf->m_L, tba_ids + this_tba_start, k*sizeof(int));

        // 7 sep
        // put temp back in tba, so the globalisation can still work
        memcpy(tba_mins[this_tba_start], temp_mins, (k - this_leaf->m_holes)*NUMDIMS*sizeof(float));
        memcpy(tba_maxes[this_tba_start], temp_maxes, (k - this_leaf->m_holes)*NUMDIMS*sizeof(float));
        memcpy(tba_ids + this_tba_start, temp_ids, (k - this_leaf->m_holes)*sizeof(int));
        // 7 sep
        // for(folan = 0; folan < k; folan++){
        //     memcpy(m_data_arr_mins[this_leaf->m_L +  folan], this_tba_mins+folan*NUMDIMS, NUMDIMS*sizeof(float));
        //     memcpy(m_data_arr_maxes[this_leaf->m_L + folan], this_tba_maxes+folan*NUMDIMS, NUMDIMS*sizeof(float));
            
        // }
        // for(int ii = 0; ii < k; ii++){
        //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
        //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
        //     cout << endl;
        // }
        // for(folan = 0; folan < k; folan++){
        //     for(filan = 0; filan < NUMDIMS; filan++){
        //         m_data_arr_mins[this_leaf->m_L + folan][filan] = this_tba_mins[folan][filan];
        //         m_data_arr_maxes[this_leaf->m_L + folan][filan] = this_tba_maxes[folan][filan];
        //     }
        //     m_data_arr_ids[this_leaf->m_L + folan] =  this_tba_ids[folan];
        // }
        // cout << "+++++++++++++" << endl;
        // this-leaf.L += k\;
        this_leaf->m_L += k;
        // this-leaf.left-sibling.R += k\;
        (this_leaf->m_older_brother)->m_R += k;
        // k -= this-leaf.holes\;
        k -= this_leaf->m_holes;
        // this-leaf.holes = 0\;
        this_leaf->m_holes = 0;
        this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
        // \textbf{\name right shuffle}(overlapping-leaves, this-leaf.right-sibling, tmp, k)\;
        shuffle_right(overlapping_leaves, this_leaf->m_younger_brother, this_tba_start, k);

        free(temp_ids);
        // for(folan = 0; folan < k - this_leaf->m_holes; folan++){
        //     free(temp_mins[folan]);
        //     free(temp_maxes[folan]);
        // }
        free(temp_mins);
        free(temp_maxes);
        }
    }

}

RTREE_TEMPLATE
void RTREE_QUAL::shuffle_left(vector<overlapping_leaf_tbas> *overlapping_leaves, Node *this_leaf, int this_tba_start, int k){
    if(k <= 0) return;

    int folan, filan;
    // cout << "shuffle left on " << this_leaf->m_id << " L, h, R: " << this_leaf->m_L << " " << this_leaf->m_holes << " " << this_leaf->m_R << endl;

    // if(this_leaf->m_id == 24867){cout << "shuffle left called on 24867!!!!!!!!!!" << endl;}
    // if(this_leaf->m_id == 1192){cout << "shuffle left called on 1192!!!!!!!!!!" << endl;}
    if(k <= this_leaf->m_holes){
        // cout << "case 133" << endl;
        // % move the k items at the end to the holes, then replace them with tbi
        // \textbf{arr}[this-leaf.L + this-leaf.holes - k: this-leaf.L + this-leaf.holes] = \textbf{arr}[this-leaf.R - k: this-leaf.R]\;
        int src = max(this_leaf->m_R - k, this_leaf->m_L + this_leaf->m_holes);
        int how_many_to_move = min(k, this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes);
        int dest = min(this_leaf->m_L + this_leaf->m_holes - how_many_to_move, this_leaf->m_R - k - how_many_to_move);

        memcpy(m_data_arr_mins[dest], m_data_arr_mins[src], NUMDIMS*how_many_to_move*sizeof(float));
        memcpy(m_data_arr_maxes[dest], m_data_arr_maxes[src], NUMDIMS*how_many_to_move*sizeof(float));
        memcpy(m_data_arr_ids + dest, m_data_arr_ids + src, how_many_to_move*sizeof(int));

        // for(folan = 0; folan < k; folan++){
        //     for(filan =0; filan < NUMDIMS; filan++){
        //         m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes - k + folan][filan] = m_data_arr_mins[this_leaf->m_R - k + folan][filan];
        //         m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes - k + folan][filan] = m_data_arr_maxes[this_leaf->m_R - k + folan][filan];
        //     }
        //     m_data_arr_ids[this_leaf->m_L + this_leaf->m_holes - k + folan] = m_data_arr_ids[this_leaf->m_R - k + folan];
        // }
        // \textbf{arr}[this-leaf.R - k: this-leaf.R] = this-tbi\;
        // for(int ii = 0; ii < k; ii++){
        //     cout << "tba " << this_tba_ids[ii] << " : ";
        //     for(int jj = 0; jj < NUMDIMS; jj++) cout << this_tba_mins[ii][jj] << " " << this_tba_maxes[ii][jj] << " , ";
        //     cout << endl;
        // }
        // for(int ii = 0; ii < k; ii++){
        //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
        //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
        //     cout << endl;
        // }
        // cout << "K = " << k << endl;
        memcpy(m_data_arr_mins[this_leaf->m_R - k], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
        memcpy(m_data_arr_maxes[this_leaf->m_R - k], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
        memcpy(m_data_arr_ids + this_leaf->m_R - k, tba_ids + this_tba_start, k*sizeof(int));
        // for(folan = 0; folan < k; folan++){
        //     memcpy(m_data_arr_mins[this_leaf->m_R - k + folan], this_tba_mins+folan*NUMDIMS, NUMDIMS*sizeof(float));
        //     memcpy(m_data_arr_maxes[this_leaf->m_R - k +  folan], this_tba_maxes+folan*NUMDIMS, NUMDIMS*sizeof(float));
        // }
        // for(int ii = 0; ii < k; ii++){
        //     cout << "copied " << m_data_arr_ids[this_leaf->m_R - k + ii] << " : ";
        //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_R - k + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_R - k + ii][jj] << " , ";
        //     cout << endl;
        // }
        // for(folan = 0; folan < k; folan++){
        //     for(filan = 0; filan < NUMDIMS; filan++){
        //         m_data_arr_mins[this_leaf->m_R - k + folan][filan] = this_tba_mins[folan][filan];
        //         m_data_arr_maxes[this_leaf->m_R - k + folan][filan] = this_tba_maxes[folan][filan];
        //     }
        //     m_data_arr_ids[this_leaf->m_R - k + folan] = this_tba_ids[folan];
        // }
        // this-leaf.R -= k\;
        // cout << "CaseD 1\n";
        // cout << "where we put the tba: " << this_leaf->m_R - k << endl;
        // cout << "k, src, dest, hmtm: " << k << " " << src << " " << dest << " " << how_many_to_move << endl;
        // cout << "holes were: " << this_leaf->m_holes << endl;
        this_leaf->m_holes = dest - this_leaf->m_L;
        // cout << "holes are: " << this_leaf->m_holes << endl;
        // cout << "Right was " << this_leaf->m_R << endl;
        this_leaf->m_R -= k;
        // cout << "Changed to " << this_leaf->m_R << endl;
        // this-leaf.right-sibling.L -= k\;
        // cout << "Left was " << (this_leaf->m_younger_brother)->m_L << endl;
        (this_leaf->m_younger_brother)->m_L -= k;
        // cout << "Changed to " << (this_leaf->m_younger_brother)->m_L << endl;
        // this-leaf.holes -= k\;
        // this_leaf->m_holes -= k;
        
        this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
        // set the whole for the leaf on the right
        (this_leaf->m_younger_brother)->m_holes = 0;
        // return\;
        return;
    }
    
    else if(this_leaf->m_id == m_oldest_leaf->m_id){
        // cout << "was the oldest " << endl;
        if(k > (this_leaf->m_R - this_leaf->m_L)){
            // not sue what
            // TODO: I have an idea
            // let's move this_leaf->left to the end with its tbi
            // then add the space to the holes of this_leaf->left->left
            // then the oldest will stay as it is 
            // and we would be done
            // move this-leaf->left to the end
            //  cout << "CaseD 2\n";
            int right_leaf_size = ((this_leaf->m_younger_brother)->m_R - (this_leaf->m_younger_brother)->m_L);
            Node* right_leaf = this_leaf->m_younger_brother;
            int last_filled_spot_in_arr = m_youngest_leaf->m_R;

            memcpy(m_data_arr_mins[last_filled_spot_in_arr], m_data_arr_mins[right_leaf->m_L + right_leaf->m_holes], NUMDIMS*(right_leaf->m_R - right_leaf->m_L - right_leaf->m_holes)*sizeof(float));
            memcpy(m_data_arr_maxes[last_filled_spot_in_arr], m_data_arr_maxes[right_leaf->m_L + right_leaf->m_holes], NUMDIMS*(right_leaf->m_R - right_leaf->m_L - right_leaf->m_holes)*sizeof(float));
            memcpy(m_data_arr_ids + last_filled_spot_in_arr, m_data_arr_ids + right_leaf->m_L + right_leaf->m_holes, (right_leaf->m_R - right_leaf->m_L - right_leaf->m_holes)*sizeof(int));

 
            // for(folan = 0; folan < (right_leaf->m_R - right_leaf->m_L - right_leaf->m_holes); folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_data_arr_mins[last_filled_spot_in_arr + folan][filan] = m_data_arr_mins[right_leaf->m_L + right_leaf->m_holes + folan][filan];
            //         m_data_arr_maxes[last_filled_spot_in_arr + folan][filan] = m_data_arr_maxes[right_leaf->m_L + right_leaf->m_holes + folan][filan];

            //     }
            //     m_data_arr_ids[last_filled_spot_in_arr + folan] = m_data_arr_ids[right_leaf->m_L + right_leaf->m_holes + folan];
            // }

            folan = (right_leaf->m_R - right_leaf->m_L - right_leaf->m_holes);
            // cout << "Right_leaf LEFT was " << right_leaf->m_L << endl;
            right_leaf->m_L = last_filled_spot_in_arr;
            // cout << "Changed to " << right_leaf->m_L << endl;
            // cout << "Right_leaf RIGHT was " << right_leaf->m_R << endl;
            right_leaf->m_R = last_filled_spot_in_arr + folan;
            // cout << "Changed to " << right_leaf->m_R  << endl;

            // the tbi too
            // for(int ii = 0; ii < k; ii++){
            //     cout << "tba " << this_tba_ids[ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << this_tba_mins[ii][jj] << " " << this_tba_maxes[ii][jj] << " , ";
            //     cout << endl;
            // }
            // for(int ii = 0; ii < k; ii++){
            //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
            //     cout << endl;
            // }
            memcpy(m_data_arr_mins[right_leaf->m_R], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[right_leaf->m_R], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_ids + right_leaf->m_R, tba_ids + this_tba_start, k*sizeof(int));
            // for(folan = 0; folan < k; folan++){
            //     memcpy(m_data_arr_mins[right_leaf->m_R +  folan], this_tba_mins+folan*NUMDIMS, NUMDIMS*sizeof(float));
            //     memcpy(m_data_arr_maxes[right_leaf->m_R +  folan], this_tba_maxes+folan*NUMDIMS, NUMDIMS*sizeof(float));
            
            // }
            // for(int ii = 0; ii < k; ii++){
            //     cout << "copied " << m_data_arr_ids[this_leaf->m_R + ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_R + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_R + ii][jj] << " , ";
            //     cout << endl;
            // }
            // for(folan = 0; folan < k; folan++){
            //     for(filan =0; filan < NUMDIMS; filan++){
            //         m_data_arr_mins[right_leaf->m_R + folan][filan] = this_tba_mins[folan][filan];
            //         m_data_arr_maxes[right_leaf->m_R + folan][filan] = this_tba_maxes[folan][filan];
            //     }
            //     m_data_arr_ids[right_leaf->m_R + folan] = this_tba_ids[folan];
            // }
            right_leaf->m_R += k;
            // cout << "Right_leaf RIGHT CHANGED AGAIN " << right_leaf->m_R << endl;
            
            (right_leaf->m_older_brother)->m_younger_brother = right_leaf->m_younger_brother;
            (right_leaf->m_younger_brother)->m_older_brother = right_leaf->m_older_brother;

            right_leaf->m_older_brother = m_youngest_leaf;
            m_youngest_leaf->m_younger_brother = right_leaf;
            right_leaf->m_younger_brother = NULL;
            right_leaf->m_holes = 0;
            right_leaf->m_count = right_leaf->m_R - right_leaf->m_L - right_leaf->m_holes;
            m_youngest_leaf = right_leaf;
            
            // add the holes to this-leaf->left  this is the old thisleaf->left->left, so it's correct
           
            // cout << "this_leaf Left was " << (this_leaf->m_younger_brother)->m_L << endl;
            (this_leaf->m_younger_brother)->m_L -= right_leaf_size;
            // cout << "Changed to " << (this_leaf->m_younger_brother)->m_L << endl;
            (this_leaf->m_younger_brother)->m_holes += (right_leaf_size);
            // ASSERT((this_leaf->m_younger_brother)->m_count == (this_leaf->m_younger_brother)->m_R - (this_leaf->m_younger_brother)->m_L - (this_leaf->m_younger_brother)->m_holes);
            return;
        }
        else{
            // cout << " k smaller than size " << endl;
            int this_L = this_leaf->m_L, this_R = this_leaf->m_R;
            ASSERT(this_L == 0);
            Node* right_leaf = this_leaf->m_younger_brother;
            (this_leaf->m_younger_brother)->m_older_brother = NULL;
            m_oldest_leaf = right_leaf;
            // move this leaf to end and fix all pointers

            // cout << "id , L, R, h:" << this_leaf->m_id << " " << this_L << " " << this_R << " " << this_leaf->m_holes << endl;

            int last_filled_spot_in_arr = m_youngest_leaf->m_R;
            // cout << "src, dext: " << this_leaf->m_L + this_leaf->m_holes << " " << last_filled_spot_in_arr << endl;

            folan = (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes);


            memcpy(m_data_arr_mins[last_filled_spot_in_arr], m_data_arr_mins[this_L + this_leaf->m_holes], NUMDIMS*(folan)*sizeof(float));
            memcpy(m_data_arr_maxes[last_filled_spot_in_arr], m_data_arr_maxes[this_L + this_leaf->m_holes], NUMDIMS*(folan)*sizeof(float));
            memcpy(m_data_arr_ids + last_filled_spot_in_arr, m_data_arr_ids + this_L + this_leaf->m_holes, (folan)*sizeof(int));

            // for(folan = 0; folan < (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes); folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_data_arr_mins[last_filled_spot_in_arr + folan][filan] = m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes + folan][filan];
            //         m_data_arr_maxes[last_filled_spot_in_arr + folan][filan] = m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes + folan][filan];

            //     }
            //     m_data_arr_ids[last_filled_spot_in_arr + folan] = m_data_arr_ids[this_leaf->m_L + right_leaf->m_holes + folan];
            // }
            // cout << "CaseD 3" << endl;
            // cout << "This_leaf LEFT was " << this_leaf->m_L << endl;
            this_leaf->m_L = last_filled_spot_in_arr;
            // cout << "Changed to " << this_leaf->m_L << endl;
            // cout << "This_leaf RIGHT was " << this_leaf->m_R << endl;
            this_leaf->m_R = last_filled_spot_in_arr + folan;
            // cout << "Changed to " << this_leaf->m_R << endl;

            this_leaf->m_younger_brother = NULL;
            this_leaf->m_older_brother = m_youngest_leaf;
            m_youngest_leaf->m_younger_brother = this_leaf;
            this_leaf->m_holes = 0;
            // ASSERT(this_leaf->m_count == (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes));
            m_youngest_leaf = this_leaf;
            // cout << "k, copying tba to " << k << " " << this_R - k << endl;
            // \textbf{arr}[this-R - k: this-R] = this-tbi\;
            // for(int ii = 0; ii < k; ii++){
            //     cout << "tba " << this_tba_ids[ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << this_tba_mins[ii][jj] << " " << this_tba_maxes[ii][jj] << " , ";
            //     cout << endl;
            // }
            // for(int ii = 0; ii < k; ii++){
            //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
            //     cout << endl;
            // }
            memcpy(m_data_arr_mins[this_R - k], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_R - k], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_ids + this_R - k, tba_ids + this_tba_start, k*sizeof(int));
            // for(folan = 0; folan < k; folan++){
            //     memcpy(m_data_arr_mins[this_R - k + folan], this_tba_mins+folan*NUMDIMS, NUMDIMS*sizeof(float));
            //     memcpy(m_data_arr_maxes[this_R - k +  folan], this_tba_maxes+folan*NUMDIMS, NUMDIMS*sizeof(float));
            // }
            // for(int ii = 0; ii < k; ii++){
            //     cout << "copied " << m_data_arr_ids[this_R - k + ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_R - k + ii][jj] << " " << m_data_arr_maxes[this_R - k + ii][jj] << " , ";
            //     cout << endl;
            // }
            // for(folan =0; folan < k; folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_data_arr_mins[this_R - k + folan][filan] = this_tba_mins[folan][filan];
            //         m_data_arr_maxes[this_R - k + folan][filan] = this_tba_maxes[folan][filan];
            //     }
            //     m_data_arr_ids[this_R - k + folan] = this_tba_ids[folan];
            // }

            right_leaf->m_L =0;
            right_leaf->m_holes = this_R - k;
            // cout << "right leaf, oldest leaf, L, R, h: " << right_leaf->m_L << " " << right_leaf->m_holes << " " << right_leaf->m_R << endl;
            ASSERT(right_leaf->m_id == m_oldest_leaf->m_id);
            return;
            

            
        }
    }

    // if(find(this_tba->begin(), this_tba->end(), this_leaf) == this_tba->end()){
    else if(!is_leaf_in_OL(this_leaf, overlapping_leaves)){
        // cout << "was not hot" << endl;
        if(k > (this_leaf->m_R - this_leaf->m_L)){

            // cout << "CaseD 4" << endl;
            // cout << "id, src, dest, k, holes" << this_leaf->m_id << " " << this_leaf->m_L + this_leaf->m_holes << " " << m_pending_insertions_count << " " << k << " " << this_leaf->m_holes << endl;
            // % let's remove this leaf completely into the pending-insertions and continue
            // pending-insertions += \textbf{arr}[this-leaf.L + this-leaf.holes: this-leaf.R]\;
            memcpy(m_pending_insertions_mins[m_pending_insertions_count], m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_pending_insertions_maxes[m_pending_insertions_count], m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_pending_insertions_ids + m_pending_insertions_count, m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes, (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(int));

            // for(folan = 0; folan < (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes); folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_pending_insertions_mins[m_pending_insertions_count + folan][filan] = m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes + folan][filan];
            //         m_pending_insertions_maxes[m_pending_insertions_count + folan][filan] = m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes + folan][filan];
            //     }
            //     m_pending_insertions_ids[m_pending_insertions_count + folan] = m_data_arr_ids[this_leaf->m_L + this_leaf->m_holes + folan];
            //     // if(m_pending_insertions_ids[m_pending_insertions_count + folan] == 100352){
            //     //     cout << "PUTTING 100351 INTO PENDING!!!!!!!!!!!!!" << endl;
            //     // }
            //     // cout << "putting into pending " << m_pending_insertions_ids[m_pending_insertions_count + folan]  << endl;
                
            // }
            m_pending_insertions_count+=(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes);

            //  \textbf{arr}[this-leaf.L: this-leaf.R] = this-tbi[:(this-leaf.R - this-leaf.L)]
            // for(int ii = 0; ii < this_leaf->m_R - this_leaf->m_L; ii++){
            //     cout << "tba " << this_tba_ids[ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << this_tba_mins[ii][jj] << " " << this_tba_maxes[ii][jj] << " , ";
            //     cout << endl;
            // }
            // for(int ii = 0; ii < k; ii++){
            //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
            //     cout << endl;
            // }
            memcpy(m_data_arr_mins[this_leaf->m_L], tba_mins[this_tba_start], (this_leaf->m_R - this_leaf->m_L)*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_L], tba_maxes[this_tba_start], (this_leaf->m_R - this_leaf->m_L)*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_ids + this_leaf->m_L, tba_ids + this_tba_start, (this_leaf->m_R - this_leaf->m_L)*sizeof(int));
            // for( folan = 0; folan < (this_leaf->m_R - this_leaf->m_L); folan++){
            //     memcpy(m_data_arr_mins[this_leaf->m_L + folan], this_tba_mins+folan*NUMDIMS, NUMDIMS*sizeof(float));
            //     memcpy(m_data_arr_maxes[this_leaf->m_L +  folan], this_tba_maxes+folan*NUMDIMS, NUMDIMS*sizeof(float));
                    
            // }
            // for(int ii = 0; ii < this_leaf->m_R - this_leaf->m_L; ii++){
            //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
            //     cout << endl;
            // }
            // for(folan = 0; folan < (this_leaf->m_R - this_leaf->m_L); folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_data_arr_mins[this_leaf->m_L + folan][filan] = this_tba_mins[folan][filan];
            //         m_data_arr_maxes[this_leaf->m_L + folan][filan] = this_tba_maxes[folan][filan];
            //     }
            //     m_data_arr_ids[this_leaf->m_L + folan] =  this_tba_ids[folan];
            // }

            // this-leaf.left-sibling.right-sibling = this-leaf.right-sibling\;
            // this-leaf.right-sibling.left-sibling = this-leaf.leaf-sibling\;

            // (this_leaf->m_older_brother)->m_younger_brother = this_leaf->m_younger_brother;
            // (this_leaf->m_younger_brother)->m_older_brother = this_leaf->m_older_brother;

            // this-leaf.left-sibling.R = this-leaf.R\;
            // cout << "This_leaf->younger was " << (this_leaf->m_younger_brother)->m_L << endl;
            (this_leaf->m_younger_brother)->m_L = this_leaf->m_L;
            // cout << "Changed to " << (this_leaf->m_younger_brother)->m_L << endl;

            Node *next_leaf = this_leaf->m_older_brother; 
            int this_leaf_size = this_leaf->m_R - this_leaf->m_L;

            // remove this leaf from parent
            int branch_index = getNodeBranchIndex(this_leaf);
            Node* parent = this_leaf->m_parent;
            DisconnectBranch(parent, branch_index);
            // FreeNode(this_leaf);
            deleteEmptyNodesUp(parent);
            
            // TODO: think about this, I do not remember exactly what disconnct branch does
            // shuffle_left(overlapping_leaves, next_leaf, this_tba[std::slice(this_leaf_size, k - this_leaf_size, 1)], k - this_leaf_size);
            // vector<Branch> stupid(this_tba.begin() + this_leaf_size, this_tba.end());
            shuffle_left(overlapping_leaves, next_leaf, this_tba_start + this_leaf_size, k - this_leaf_size);
        }
        else{
            // cout << " not hot, k larger than holes, smaller than size" << endl;
            // % holes < k < size
            // % first let's try to fill the holes
            // \textbf{arr}[this-leaf.L: this-leaf.L + this-leaf.holes] = \textbf{arr}[this-leaf.R - this-leaf.holes: this-leaf.R]\;
            // cout << "CaseD 5" << endl;
            int src = max(this_leaf->m_L + this_leaf->m_holes, this_leaf->m_R - this_leaf->m_holes);
            int how_many_to_move = min(this_leaf->m_holes,this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes);
            // cout << "id, src, dest, k, holes" << this_leaf->m_id << " " << src << " " << this_leaf->m_L << " " << k << " " << this_leaf->m_holes << " hmtm " << how_many_to_move << endl;

            memcpy(m_data_arr_mins[this_leaf->m_L], m_data_arr_mins[src], NUMDIMS*(how_many_to_move)*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_L], m_data_arr_maxes[src], NUMDIMS*(how_many_to_move)*sizeof(float));
            memcpy(m_data_arr_ids + this_leaf->m_L, m_data_arr_ids + src, how_many_to_move*sizeof(int));


            // for(folan = 0; folan < (this_leaf->m_holes); folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_data_arr_mins[this_leaf->m_L + folan][filan] = m_data_arr_mins[this_leaf->m_R - this_leaf->m_holes + folan][filan];
            //         m_data_arr_maxes[this_leaf->m_L + folan][filan] = m_data_arr_maxes[this_leaf->m_R - this_leaf->m_holes + folan][filan];
            //     }
            //     m_data_arr_ids[this_leaf->m_L + folan] = m_data_arr_ids[this_leaf->m_R - this_leaf->m_holes + folan];
            // }
            // pending-insertions += \textbf{arr}[this-leaf.R - k: this-leaf.R - this-leaf.holes]\;
            memcpy(m_pending_insertions_mins[m_pending_insertions_count], m_data_arr_mins[this_leaf->m_R - k], NUMDIMS*(k - this_leaf->m_holes)*sizeof(float));
            memcpy(m_pending_insertions_maxes[m_pending_insertions_count], m_data_arr_maxes[this_leaf->m_R - k], NUMDIMS*(k - this_leaf->m_holes)*sizeof(float));
            memcpy(m_pending_insertions_ids + m_pending_insertions_count, m_data_arr_ids + this_leaf->m_R - k, (k - this_leaf->m_holes)*sizeof(int));

            // cout << "moving from " << this_leaf->m_R - k << " to pending, this many " << (k - this_leaf->m_holes) << endl;


            // for(folan = 0; folan < (k - this_leaf->m_holes); folan++){
            //     for(filan = 0; filan < NUMDIMS; filan ++){
            //         m_pending_insertions_mins[m_pending_insertions_count + folan][filan] = m_data_arr_mins[this_leaf->m_R - k + folan][filan];
            //         m_pending_insertions_maxes[m_pending_insertions_count + folan][filan] = m_data_arr_maxes[this_leaf->m_R - k + folan][filan];
            //     }
            //     m_pending_insertions_ids[m_pending_insertions_count + folan] = m_data_arr_ids[this_leaf->m_R - k + folan];
            //     // if(m_pending_insertions_ids[m_pending_insertions_count + folan] == 100352){
            //     //     cout << "PUTTING 100351 INTO PENDING!!!!!!!!!!!!!" << endl;
            //     // }
            //     // cout << "putting into pending " << m_pending_insertions_ids[m_pending_insertions_count + folan]  << endl;
            //     // m_pending_insertions_count++;
            // }
            m_pending_insertions_count+=(k - this_leaf->m_holes);
            // \textbf{arr}[this-leaf.R - k: this-leaf.R] = this-tbi\;
            // for(int ii = 0; ii < k; ii++){
            //     cout << "tba " << this_tba_ids[ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << this_tba_mins[ii][jj] << " " << this_tba_maxes[ii][jj] << " , ";
            //     cout << endl;
            // }
            // for(int ii = 0; ii < k; ii++){
            //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
            //     cout << endl;
            // }
            // cout << "this_leaf->m_R = " << this_leaf->m_R << " k = " << k << endl;
            memcpy(m_data_arr_mins[this_leaf->m_R - k], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_R - k], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_ids + this_leaf->m_R - k, tba_ids + this_tba_start, k*sizeof(int));
            // for(folan = 0; folan < k; folan++){
            //     memcpy(m_data_arr_mins[this_leaf->m_R - k + folan], this_tba_mins+folan*NUMDIMS, NUMDIMS*sizeof(float));
            //     memcpy(m_data_arr_maxes[this_leaf->m_R - k +  folan], this_tba_maxes+folan*NUMDIMS, NUMDIMS*sizeof(float));
                    
            // }            
            // for(int ii = 0; ii < k; ii++){
            //     cout << "copied " << m_data_arr_ids[this_leaf->m_R - k + ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_R - k + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_R - k + ii][jj] << " , ";
            //     cout << endl;
            // }
            // for(folan = 0; folan < k; folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_data_arr_mins[this_leaf->m_R - k + folan][filan] = this_tba_mins[folan][filan];
            //         m_data_arr_maxes[this_leaf->m_R - k + folan][filan] = this_tba_maxes[folan][filan];
            //     }
            //     m_data_arr_ids[this_leaf->m_R - k + folan] = this_tba_ids[folan];
            // }
            // this-leaf.R -= k\;
            // cout << "CaseD 5\n";
            // cout << "This leaf RIGHT was " << this_leaf->m_R << endl;
            this_leaf->m_R -= k;
            // cout << "Changed to " << this_leaf->m_R << endl;
            // this-leaf.right-sibling.L -= k\;
            // cout << "This_leaf->older LEFT was " << (this_leaf->m_older_brother)->m_L << endl;
            (this_leaf->m_younger_brother)->m_L -= k;
            // cout << "Changed to " << (this_leaf->m_older_brother)->m_L << endl;
            // this-leaf.holes = 0\;
            this_leaf->m_holes = 0;
            this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;

            // set the holes for the leaf on the right
            (this_leaf->m_younger_brother)->m_holes = 0;
            // return\;
            return;
        }
    }
    else{
        if(k >= (this_leaf->m_R - this_leaf->m_L)){
            // cout << "k > R - L" << endl;
            Node *next_leaf = this_leaf->m_older_brother;
            (this_leaf->m_younger_brother)->m_L = this_leaf->m_L;
            (this_leaf->m_younger_brother)->m_holes = 0;
            // this-leaf.left-sibling.right-sibling = this-leaf.right-sibling\;
            // this-leaf.right-sibling.left-sibling = this-leaf.leaf-sibling\;
            
            // (this_leaf->m_older_brother)->m_younger_brother = this_leaf->m_younger_brother;
            // (this_leaf->m_younger_brother)->m_older_brother = this_leaf->m_older_brother;

            // Move this-leaf to end of \textbf{arr} and fix all the pointers\;
            
            int this_old_L = this_leaf->m_L, this_old_R = this_leaf->m_R;
            int last_filled_spot_in_arr = m_youngest_leaf->m_R;
            // cout << "this_old_L before: " << this_old_L << " " << this_old_R << endl;

            memcpy(m_data_arr_mins[last_filled_spot_in_arr], m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_data_arr_maxes[last_filled_spot_in_arr], m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_data_arr_ids + last_filled_spot_in_arr, m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes, (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(int));


            // for(folan = 0; folan < (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes); folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_data_arr_mins[last_filled_spot_in_arr + folan][filan] = m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes + folan][filan];
            //         m_data_arr_maxes[last_filled_spot_in_arr + folan][filan] = m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes + folan][filan];

            //     }
            //     m_data_arr_ids[last_filled_spot_in_arr + folan] = m_data_arr_ids[this_leaf->m_L + this_leaf->m_holes + folan];
            // }
            folan = (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes);
            // cout << "CaseD 6\n";
            // cout << "This_leaf LEFT was " << this_leaf->m_L << endl;
            this_leaf->m_L = last_filled_spot_in_arr;
            // cout << "Changed to " << this_leaf->m_L << endl;
            // cout << "This_leaf RIGHT was " << this_leaf->m_R << endl;
            this_leaf->m_R = last_filled_spot_in_arr + folan;
            // cout << "Changed to " << this_leaf->m_R << endl;

            
            (this_leaf->m_older_brother)->m_younger_brother = this_leaf->m_younger_brother;
            (this_leaf->m_younger_brother)->m_older_brother = this_leaf->m_older_brother;

            
            // cout << (this_leaf->m_younger_brother)->m_L << " is the L of right sibling" << endl;
            // cout << "younger bro id, L, R, left_bro_id " << (this_leaf->m_younger_brother)->m_id << " " << (this_leaf->m_younger_brother)->m_L << " " << (this_leaf->m_younger_brother)->m_R << " " << ((this_leaf->m_younger_brother)->m_older_brother)->m_id << endl;
            // cout << "next leaf id, L, R, right_bro_id " << next_leaf->m_id << " " << next_leaf->m_L << " " << next_leaf->m_R << " " << (next_leaf->m_younger_brother)->m_id << endl; 

            this_leaf->m_older_brother = m_youngest_leaf;
            m_youngest_leaf->m_younger_brother = this_leaf;
            this_leaf->m_younger_brother = NULL;
            this_leaf->m_holes = 0;
            this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
            m_youngest_leaf = this_leaf;

            

            // cout << "this_old_L after: " << this_old_L << " " << this_old_R << endl;

            // cout << "moving to " << this_old_L << endl; 
            // \textbf{arr}[this-leaf.L: this-leaf.R] = this-tbi[:(this-leaf.R - this-leaf.L)]\;
            // for(int ii = 0; ii < (this_old_R - this_old_L); ii++){
            //     cout << "tba " << this_tba_ids[ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << this_tba_mins[ii][jj] << " " << this_tba_maxes[ii][jj] << " , ";
            //     cout << endl;
            // }
            // for(int ii = 0; ii < k; ii++){
            //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
            //     cout << endl;
            // }
            memcpy(m_data_arr_mins[this_old_L], tba_mins[this_tba_start], (this_old_R - this_old_L)*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_old_L], tba_maxes[this_tba_start], (this_old_R - this_old_L)*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_ids + this_old_L, tba_ids + this_tba_start, (this_old_R - this_old_L)*sizeof(int));
            // for(folan = 0; folan < (this_old_R - this_old_L); folan++){
            //     memcpy(m_data_arr_mins[this_old_L + folan], this_tba_mins+folan*NUMDIMS, NUMDIMS*sizeof(float));
            //     memcpy(m_data_arr_maxes[this_old_L +  folan], this_tba_maxes+folan*NUMDIMS, NUMDIMS*sizeof(float));
                
            // }
            // for(int ii = 0; ii < (this_old_R - this_old_L); ii++){
            //     cout << "copied " << m_data_arr_ids[this_old_L + ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_old_L + ii][jj] << " " << m_data_arr_maxes[this_old_L + ii][jj] << " , ";
            //     cout << endl;
            // }
            // for(folan = 0; folan < (this_old_R - this_old_L); folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_data_arr_mins[this_old_L + folan][filan] = this_tba_mins[folan][filan];
            //         m_data_arr_maxes[this_old_L + folan][filan] = this_tba_maxes[folan][filan];
            //     }
            //     m_data_arr_ids[this_old_L + folan] =  this_tba_ids[folan];
            //     // cout << "moved " << m_data_arr_ids[this_old_L + folan] << " into the arr" << endl;
            // }
            // \CommentSty{// then the rest of tbi still belongs to the next-leaf.left-sib, cause we removed the one in between(moved to the end}\;
            // shuffle_right(overlapping_leaves, next_leaf, this_tba[std::slice((this_old_R - this_old_L), k - (this_old_R - this_old_L), 1)], k - (this_old_R - this_old_L));
            // vector<Branch> stupid(this_tba.begin() + (this_old_R - this_old_L), this_tba.end());
            shuffle_left(overlapping_leaves, next_leaf, this_tba_start + this_old_R - this_old_L, k - (this_old_R - this_old_L));

        }
        else{
            // cout << "normal shuffle to the next case" << endl;
            // fill the holes
            // \textbf{arr}[this-leaf.L:this-leaf.L + this-leaf.holes] = \textbf{arr}[this-leaf.R - this-leaf.holes: this-leaf.R]
            int how_many_to_move = min(this_leaf->m_holes, this_leaf->m_R - this_leaf->m_holes-this_leaf->m_L);
            memcpy(m_data_arr_mins[this_leaf->m_L], m_data_arr_mins[this_leaf->m_R - how_many_to_move], NUMDIMS*how_many_to_move*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_L], m_data_arr_maxes[this_leaf->m_R - how_many_to_move], NUMDIMS*how_many_to_move*sizeof(float));
            memcpy(m_data_arr_ids + this_leaf->m_L, m_data_arr_ids + this_leaf->m_R - how_many_to_move, how_many_to_move*sizeof(int));


            // for(folan = 0; folan < this_leaf->m_holes; folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_data_arr_mins[this_leaf->m_L + folan][filan] = m_data_arr_mins[this_leaf->m_R - this_leaf->m_holes + folan][filan];
            //         m_data_arr_maxes[this_leaf->m_L + folan][filan] = m_data_arr_maxes[this_leaf->m_R - this_leaf->m_holes + folan][filan];
            //     }
            //     m_data_arr_ids[this_leaf->m_L + folan] = m_data_arr_ids[this_leaf->m_R - this_leaf->m_holes + folan];
            // }
            // tmp = \textbf{arr}[this-leaf.R - k: this-leaf.R - this-leaf.holes]
            // vector<Branch> temp(k - this_leaf->m_holes);
            // float temp_mins[k - this_leaf->m_holes][NUMDIMS];
            // float temp_maxes[k - this_leaf->m_holes][NUMDIMS];
            // int temp_ids[k - this_leaf->m_holes];

            float* temp_mins;
            temp_mins = (float*)malloc((k - this_leaf->m_holes) * NUMDIMS * sizeof(float));

            float* temp_maxes;
            temp_maxes = (float*)malloc((k - this_leaf->m_holes) * NUMDIMS * sizeof(float));

            // for(folan = 0; folan < k - this_leaf->m_holes; folan++){
            //     temp_mins[folan] = (float*)malloc(NUMDIMS*sizeof(float));
            //     temp_maxes[folan] = (float*)malloc(NUMDIMS*sizeof(float));
            // }

            int* temp_ids = (int*)malloc((k - this_leaf->m_holes) * sizeof(int));

            memcpy(temp_mins, m_data_arr_mins[this_leaf->m_R - k], (k - this_leaf->m_holes)*NUMDIMS*sizeof(float));
            memcpy(temp_maxes, m_data_arr_maxes[this_leaf->m_R - k], (k - this_leaf->m_holes)*NUMDIMS*sizeof(float));
            memcpy(temp_ids, m_data_arr_ids + this_leaf->m_R - k, (k - this_leaf->m_holes)*sizeof(int));

            // for(folan = 0; folan < (k - this_leaf->m_holes); folan++){
            //     memcpy(temp_mins + folan * NUMDIMS, m_data_arr_mins[this_leaf->m_R - k + folan], NUMDIMS*sizeof(float));
            //     memcpy(temp_maxes + folan * NUMDIMS, m_data_arr_maxes[this_leaf->m_R - k + folan], NUMDIMS*sizeof(float));
                    
            // }


            // for(folan = 0; folan < (k - this_leaf->m_holes); folan++){
            //     // temp[folan] = Branch(NULL, Rect(m_data_arr_mins[this_leaf->m_R - k + folan], m_data_arr_maxes[this_leaf->m_R - k + folan]), m_data_arr_ids[this_leaf->m_R - k + folan]);

            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         temp_mins[folan][filan] = m_data_arr_mins[this_leaf->m_R - k + folan][filan];
            //         temp_maxes[folan][filan] = m_data_arr_maxes[this_leaf->m_R - k + folan][filan];
            //     }
            //     temp_ids[folan] = m_data_arr_ids[this_leaf->m_R - k + folan];
            //     // TODO: I have no cle if this will actually work
            // }

            // for(folan = 0; folan < k; folan++){
            //     cout << "tba -> " << this_tba_ids[folan] << " with cover ";
            //     for(filan = 0; filan < NUMDIMS; filan++) cout <<  this_tba_mins[folan][filan] << " " <<  this_tba_maxes[folan][filan] << " , ";
            //     cout << "\n";
            // }

            // \textbf{arr}[this-leaf.R - k: this-leaf.R] = this-tbi\;
            
            memcpy(m_data_arr_mins[this_leaf->m_R - k], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_R - k], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_ids + this_leaf->m_R - k, tba_ids + this_tba_start, k*sizeof(int));

            // 7 Sep
            // put temp back in tba

            memcpy(tba_mins[this_tba_start], temp_mins, (k - this_leaf->m_holes)*NUMDIMS*sizeof(float));
            memcpy(tba_maxes[this_tba_start], temp_maxes, (k - this_leaf->m_holes)*NUMDIMS*sizeof(float));
            memcpy(tba_ids + this_tba_start, temp_ids, (k - this_leaf->m_holes)*sizeof(int));

            // 7 Sep

            // for(folan = 0; folan < k; folan++){
                
            //     memcpy(m_data_arr_mins[this_leaf->m_R - k + folan], this_tba_mins+folan*NUMDIMS, NUMDIMS*sizeof(float));
            //     memcpy(m_data_arr_maxes[this_leaf->m_R - k + folan], this_tba_maxes+folan*NUMDIMS, NUMDIMS*sizeof(float));
            // }
            

            // for(folan = 0; folan < k; folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_data_arr_mins[this_leaf->m_R - k + folan][filan] = this_tba_mins[folan][filan];
            //         m_data_arr_maxes[this_leaf->m_R - k + folan][filan] = this_tba_maxes[folan][filan];
            //     }
            //     m_data_arr_ids[this_leaf->m_R - k + folan] =  this_tba_ids[folan];
            // }
            // for(folan = 0; folan < k; folan++){
            //     cout << "Copied from " << m_data_arr_ids[this_leaf->m_R - k + folan] << " with cover ";
            //     for(filan = 0; filan < NUMDIMS; filan++) cout <<  m_data_arr_mins[this_leaf->m_R - k + folan][filan] << " " <<  m_data_arr_maxes[this_leaf->m_R - k + folan][filan] << " , ";
            //     cout << "\n";
            // }

            // cout << "CaseD 7\n";
            // cout << "This_leaf RIGHT was " << this_leaf->m_R << endl;
            this_leaf->m_R -= k;
            // cout << "Changed to " << this_leaf->m_R << endl;
            // cout << "This_leaf->younger was " << (this_leaf->m_younger_brother)->m_L << endl;
            (this_leaf->m_younger_brother)->m_L -= k;
            // cout << "Changed to " << (this_leaf->m_younger_brother)->m_L << endl;
            k -= this_leaf->m_holes;
            this_leaf->m_holes = 0;
            this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
            // set the holes for the leaf on the right
            (this_leaf->m_younger_brother)->m_holes = 0;
            shuffle_left(overlapping_leaves, this_leaf->m_older_brother, this_tba_start, k);

            free(temp_ids);
            // for(folan = 0; folan < k - this_leaf->m_holes; folan++){
            //     free(temp_mins[folan]);
            //     free(temp_maxes[folan]);
            // }
            free(temp_mins);
            free(temp_maxes);
        }
    }
}


RTREE_TEMPLATE
void RTREE_QUAL::ripple_v2(vector<overlapping_leaf_tbas> *overlapping_leaves){

    Node* this_leaf;
    // vector<Branch> this_tba;
    // float* this_tba_mins; float* this_tba_maxes; int* this_tba_ids; 
    int this_tba_count;
    int this_tba_start;
    bool dir_is_right;
    int k;
    int folan, filan;

    // cout << "in ripple  v2, ol.size() " << overlapping_leaves->size() << endl;
    for(int index =0; index < overlapping_leaves->size(); index++){
        this_leaf = overlapping_leaves->at(index).this_leaf;
        // this_tba = overlapping_leaves->at(index).this_tba;
        // this_tba_mins = overlapping_leaves->at(index).this_tba_mins;
        // this_tba_maxes = overlapping_leaves->at(index).this_tba_maxes;
        // this_tba_ids = overlapping_leaves->at(index).this_tba_ids;
        this_tba_start = overlapping_leaves->at(index).this_tba_start;
        this_tba_count = overlapping_leaves->at(index).this_tba_count;



        // cout << "PRINT tba" << endl;        
        // if(this_tba_count > 0) {
        //     cout << "in ripple, leaf_id: " << this_leaf->m_id << " tba.size(): " << this_tba_count << "  start: " << this_tba_start << " L, h, R: " << this_leaf->m_L << " " << this_leaf->m_holes << " " << this_leaf->m_R << endl;
        //     cout << "inseting these: !!!!!!!!!!!!!!!!!" << endl;
        //     for(int ii = this_tba_start; ii < this_tba_start + this_tba_count; ii++){
        //         cout << tba_ids[ii] << " " << tba_mins[ii][0] << " " << tba_maxes[ii][0] << endl;
        //     }
        //     cout << "!!!!!!!!!!!!!!!!!!!!!!" << endl;
        // }

        // cout << "in loop, this-tba.size() " << this_tba.size() << endl;
        // k = this_tba.size();
        k = this_tba_count;
        // if(this_leaf->m_id == 14769){
        //     cout << "ripple called on 14769!!!!!!!!!! m_id[99019] " << m_data_arr_ids[99019] << endl;

        // }

        // cout << "K  = " << k << endl;

        // if(this_tba.size() == 0) continue;
        if(k == 0) continue;
        // if(this_leaf->m_id == 1191){cout << "ripple called on 1191!!!!!!!!!!" << endl;}
        // if(this_leaf->m_id == 1192){cout << "ripple called on 1192!!!!!!!!!!" << endl;}
        


        if(this_leaf->m_id == m_oldest_leaf->m_id){
            // cout << "case 11" << endl;
            
            if(k <= this_leaf->m_holes){
                // cout << "2 2" << endl;
                // // arr[leaf.L + leaf.holes - k: leaf.L + leaf.holes] = leaf-tbi[:leaf.holes]
                // for(int ii = 0; ii < k; ii++){
                //     cout << "tba " << this_tba_ids[ii] << " : ";
                //     for(int jj = 0; jj < NUMDIMS; jj++) cout << this_tba_mins[ii][jj] << " " << this_tba_maxes[ii][jj] << " , ";
                //     cout << endl;
                // }
                // for(int ii = 0; ii < k; ii++){
                //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
                //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
                //     cout << endl;
                // }
                // memcpy(m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes - k], this_tba_mins[0], k*NUMDIMS*sizeof(float));
                // memcpy(m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes - k], this_tba_maxes[0], k*NUMDIMS*sizeof(float));
                memcpy(m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes - k, tba_ids + this_tba_start, k*sizeof(int));
                // for(folan = 0; folan < k; folan++){
                //     memcpy(m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes - k + folan], this_tba_mins+folan*NUMDIMS, NUMDIMS*sizeof(float));
                //     memcpy(m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes - k + folan], this_tba_maxes+folan*NUMDIMS, NUMDIMS*sizeof(float));
                        
                // }
                memcpy(m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes - k], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
                memcpy(m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes - k], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));

                // for(int ii = 0; ii < k; ii++){
                //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + this_leaf->m_holes - k + ii] << " : ";
                //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L  + this_leaf->m_holes - k+ ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes - k + ii][jj] << " , ";
                //     cout << endl;
                // }
                // for(folan = 0; folan < k; folan++){
                //     for(filan = 0; filan < NUMDIMS; filan++){
                //         // m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes - k + folan][filan] = (this_tba.at(folan).m_rect).m_min[filan];
                //         // m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes - k + folan][filan] = (this_tba.at(folan).m_rect).m_max[filan];
                //         m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes - k + folan][filan] = this_tba_mins[folan][filan];
                //         m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes - k + folan][filan] = this_tba_maxes[folan][filan];
                //     }
                //     // m_data_arr_ids[this_leaf->m_L + this_leaf->m_holes - k + folan] =  this_tba.at(folan).m_data;
                //     m_data_arr_ids[this_leaf->m_L + this_leaf->m_holes - k + folan] =  this_tba_ids[folan];

                // }
                this_leaf->m_holes -= k;
                this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
                continue;
            }
            // otherwose move it down
            // move leaf to end with tbi
            // continue
            Node* right_leaf = this_leaf->m_younger_brother;
            int old_R = this_leaf->m_R;
            int last_filled_spot = m_youngest_leaf->m_R;

            memcpy(m_data_arr_mins[last_filled_spot], m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes] , NUMDIMS* (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_data_arr_maxes[last_filled_spot], m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes] , NUMDIMS* (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
            memcpy(m_data_arr_ids + last_filled_spot, m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes , (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(int));



            // for(folan = 0; folan < (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes); folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_data_arr_mins[last_filled_spot + folan][filan] = m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes + folan][filan];
            //         m_data_arr_maxes[last_filled_spot + folan][filan] = m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes + folan][filan];
            //     }
            //     m_data_arr_ids[last_filled_spot + folan] = m_data_arr_ids[this_leaf->m_L + this_leaf->m_holes + folan];
            // }
            folan = (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes);
            this_leaf->m_L = last_filled_spot;
            this_leaf->m_R = last_filled_spot + folan;

            // and the tbi
            // for(int ii = 0; ii < k; ii++){
            //     cout << "tba " << this_tba_ids[ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << this_tba_mins[ii][jj] << " " << this_tba_maxes[ii][jj] << " , ";
            //     cout << endl;
            // }
            // for(int ii = 0; ii < k; ii++){
            //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
            //     cout << endl;
            // }
            // memcpy(m_data_arr_mins[this_leaf->m_R], this_tba_mins[0], k*NUMDIMS*sizeof(float));
            // memcpy(m_data_arr_maxes[this_leaf->m_R], this_tba_maxes[0], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_ids + this_leaf->m_R, tba_ids + this_tba_start, k*sizeof(int));
            // for(folan =0; folan < k; folan++){
            //     memcpy(m_data_arr_mins[this_leaf->m_R + folan], this_tba_mins+folan*NUMDIMS, NUMDIMS*sizeof(float));
            //     memcpy(m_data_arr_maxes[this_leaf->m_R + folan], this_tba_maxes+folan*NUMDIMS, NUMDIMS*sizeof(float));
                    
            // }
            memcpy(m_data_arr_mins[this_leaf->m_R], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_R], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
            // for(int ii = 0; ii < k; ii++){
            //     cout << "copied " << m_data_arr_ids[this_leaf->m_R + ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_R + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_R + ii][jj] << " , ";
            //     cout << endl;
            // }
            // for(folan =0; folan < k; folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_data_arr_mins[this_leaf->m_R + folan][filan] = this_tba_mins[folan][filan];
            //         m_data_arr_maxes[this_leaf->m_R + folan][filan] =  this_tba_maxes[folan][filan];
            //     }
            //     m_data_arr_ids[this_leaf->m_R + folan] = this_tba_ids[folan];
            // }

            


            this_leaf->m_R += k;

            this_leaf->m_older_brother = m_youngest_leaf;
            this_leaf->m_younger_brother = NULL;
            m_youngest_leaf->m_younger_brother = this_leaf;
            m_youngest_leaf = this_leaf;
            this_leaf->m_holes = 0;

            right_leaf->m_older_brother = NULL;
            m_oldest_leaf = right_leaf;
            right_leaf->m_L = 0;
            right_leaf->m_holes += old_R;

            continue;
            

        }
        else if(this_leaf->m_id == m_youngest_leaf->m_id){
            // cout << "case 12" << endl;
            // append tbi to end of arr
            // this-leaf.R += k
            // continue
            
            // for(int ii = 0; ii < k; ii++){
            //     cout << "tba " << this_tba_ids[ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << this_tba_mins[ii][jj] << " " << this_tba_maxes[ii][jj] << " , ";
            //     cout << endl;
            // }
            // for(int ii = 0; ii < k; ii++){
            //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
            //     cout << endl;
            // }
            // memcpy(m_data_arr_mins[this_leaf->m_R], this_tba_mins[0], k*NUMDIMS*sizeof(float));
            // memcpy(m_data_arr_maxes[this_leaf->m_R], this_tba_maxes[0], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_ids + this_leaf->m_R, tba_ids + this_tba_start, k*sizeof(int));

            // for(folan = 0; folan < k ;folan++){
            //     memcpy(m_data_arr_mins[this_leaf->m_R + folan], this_tba_mins+folan*NUMDIMS, NUMDIMS*sizeof(float));
            //     memcpy(m_data_arr_maxes[this_leaf->m_R + folan], this_tba_maxes+folan*NUMDIMS, NUMDIMS*sizeof(float));
                    
            // }
            memcpy(m_data_arr_mins[this_leaf->m_R], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_R], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
            

            // for(int ii = 0; ii < k; ii++){
            //     cout << "copied " << m_data_arr_ids[this_leaf->m_R + ii] << " : ";
            //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_R + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_R + ii][jj] << " , ";
            //     cout << endl;
            // }
            // for(folan = 0; folan < k; folan++){
            //     for(filan = 0; filan < NUMDIMS; filan++){
            //         m_data_arr_mins[this_leaf->m_R + folan][filan] = this_tba_mins[folan][filan];
            //         m_data_arr_maxes[this_leaf->m_R + folan][filan] = this_tba_maxes[folan][filan];
            //     }
            //     m_data_arr_ids[this_leaf->m_R + folan] = this_tba_ids[folan];
            // }
            this_leaf->m_R += k;
            this_leaf->m_count += k;
            continue;
        }
        else{
            dir_is_right = true;

            // if(this_leaf->m_id == m_oldest_leaf->m_id || find(this_tba->begin(), this_tba->end(), this_leaf->m_younger_brother) != this_tba->end()){
            if(this_leaf->m_id == m_oldest_leaf->m_id || is_leaf_in_OL(this_leaf->m_younger_brother, overlapping_leaves)){
                // right sibling is in OL
                dir_is_right = false;
            }
            // cout << "chose dir, dir-is-right is " << dir_is_right << endl;
            
            // cout << "leaf.size " << (this_leaf->m_R - this_leaf->m_L) << endl;
            if(k < (this_leaf->m_R - this_leaf->m_L)){
                // cout << "1" << endl;
                if(k <= this_leaf->m_holes){
                    // cout << "2" << endl;
                    // arr[leaf.L + leaf.holes - k: leaf.L + leaf.holes] = leaf-tbi[:leaf.holes]
                    // for(int ii = 0; ii < k; ii++){
                    //     cout << "tba " << this_tba_ids[ii] << " : ";
                    //     for(int jj = 0; jj < NUMDIMS; jj++) cout << this_tba_mins[ii][jj] << " " << this_tba_maxes[ii][jj] << " , ";
                    //     cout << endl;
                    // }
                    // for(int ii = 0; ii < k; ii++){
                    //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
                    //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
                    //     cout << endl;
                    // }
                    // memcpy(m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes - k], this_tba_mins[0], k*NUMDIMS*sizeof(float));
                    // memcpy(m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes - k], this_tba_maxes[0], k*NUMDIMS*sizeof(float));
                    memcpy(m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes - k, tba_ids + this_tba_start, k*sizeof(int));

                    // for(folan = 0; folan < k; folan++){
                    //     memcpy(m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes - k + folan], this_tba_mins+folan*NUMDIMS, NUMDIMS*sizeof(float));
                    //     memcpy(m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes - k + folan], this_tba_maxes+folan*NUMDIMS, NUMDIMS*sizeof(float));
                            
                    // }
                    memcpy(m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes - k], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
                    memcpy(m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes - k], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));

                    // for(int ii = 0; ii < k; ii++){
                    //     cout << "copied " << m_data_arr_ids[this_leaf->m_L  + this_leaf->m_holes - k+ ii] << " : ";
                    //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes - k + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes - k + ii][jj] << " , ";
                    //     cout << endl;
                    // }
                    // for(folan = 0; folan < k; folan++){
                    //     for(filan = 0; filan < NUMDIMS; filan++){
                    //         m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes - k + folan][filan] = this_tba_mins[folan][filan];
                    //         m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes - k + folan][filan] = this_tba_maxes[folan][filan];
                    //     }
                    //     m_data_arr_ids[this_leaf->m_L + this_leaf->m_holes - k + folan] =  this_tba_ids[folan];
                    // }
                    this_leaf->m_holes -= k;
                    this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
                    continue;
                }
                // otherwise, fill the holes completely first
                // fill the holes

                // for(int ii = 0; ii < k; ii++){
                //     cout << "tba " << this_tba_ids[ii] << " : ";
                //     for(int jj = 0; jj < NUMDIMS; jj++) cout << this_tba_mins[ii][jj] << " " << this_tba_maxes[ii][jj] << " , ";
                //     cout << endl;
                // }
                // for(int ii = 0; ii < k; ii++){
                //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
                //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
                //     cout << endl;
                // }
                memcpy(m_data_arr_mins[this_leaf->m_L], tba_mins[this_tba_start], this_leaf->m_holes*NUMDIMS*sizeof(float));
                memcpy(m_data_arr_maxes[this_leaf->m_L], tba_maxes[this_tba_start], this_leaf->m_holes*NUMDIMS*sizeof(float));
                memcpy(m_data_arr_ids + this_leaf->m_L, tba_ids + this_tba_start, this_leaf->m_holes*sizeof(int));

                // for(folan = 0; folan < this_leaf->m_holes; folan++){
                //     memcpy(m_data_arr_mins[this_leaf->m_L + folan], this_tba_mins+folan*NUMDIMS, NUMDIMS*sizeof(float));
                //     memcpy(m_data_arr_maxes[this_leaf->m_L + folan], this_tba_maxes+folan*NUMDIMS, NUMDIMS*sizeof(float));
                        
                // }

                // for(int ii = 0; ii < k; ii++){
                //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
                //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
                //     cout << endl;
                // }
                // for(folan = 0; folan < this_leaf->m_holes; folan++){
                //     for(filan = 0; filan < NUMDIMS; filan++){
                //         m_data_arr_mins[this_leaf->m_L + folan][filan] = this_tba_mins[folan][filan];
                //         m_data_arr_maxes[this_leaf->m_L + folan][filan] = this_tba_maxes[folan][filan];
                //     }
                //     m_data_arr_ids[this_leaf->m_L + folan] = this_tba_ids[folan];
                // }
                #ifdef timebreak
                    auto start_time = clock();
                #endif
                if(dir_is_right){
                    // cout << "3" << endl;
                    // shuffle_right(overlapping_leaves, this_leaf->m_younger_brother, (*this_tba)[std::slice(this_leaf->m_holes, k - this_leaf->m_holes, 1)], k - this_leaf->m_holes);
                    // vector<Branch> stupid(this_tba.begin() + this_leaf->m_holes, this_tba.end());
                    // shuffle_right(overlapping_leaves, this_leaf->m_younger_brother, stupid, k - this_leaf->m_holes);
                    // shuffle_right(overlapping_leaves, this_leaf->m_younger_brother, this_tba_mins[this_leaf->m_holes], this_tba_maxes[this_leaf->m_holes], this_tba_ids[this_leaf->m_holes], k - this_leaf->m_holes);
                    // cout << "before shuffle right passing " << k - this_leaf->m_holes << endl;

                    shuffle_right(overlapping_leaves, this_leaf->m_younger_brother, this_tba_start + this_leaf->m_holes, k - this_leaf->m_holes);
                    // cout << "after shuffle right" << endl;
                    this_leaf->m_holes = 0;
                    // TODO: im not actually sure these slice things work...
                    // TODO: i also am not sure if the k - holes should have a + 1
                }
                else{
                    // cout << "4" << endl;
                    // shuffle_left(overlapping_leaves, this_leaf->m_older_brother, this_tba[std::slice(this_leaf->m_holes, k - this_leaf->m_holes, 1)], k - this_leaf->m_holes);
                    // vector<Branch> stupid(this_tba.begin() + this_leaf->m_holes, this_tba.end());
                    // shuffle_left(overlapping_leaves, this_leaf->m_older_brother, stupid, k - this_leaf->m_holes);
                    // cout << "before shuffle left" << endl;
                    shuffle_left(overlapping_leaves, this_leaf->m_older_brother, this_tba_start + this_leaf->m_holes, k - this_leaf->m_holes);
                    // cout << "after shuffle left" << endl;
                    // TODO: same as above

                }
                #ifdef timebreak
                    ripple_time += (clock() - start_time);
                #endif

                // this_leaf->m_holes = 0;
                this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;


            }
            else{
                // cout << "5" << endl;
                // cout << "@@@@@@@ " << this_leaf->m_L << " " << this_leaf->m_holes << " " << this_leaf->m_R << " " << k << endl;
                // Set leaf.size holes at the start of leaf.right-sibling\;
                // Move leaf completely to the end of \textbf{arr}, and put leaf-tbi there too\;

                (this_leaf->m_younger_brother)->m_L = this_leaf->m_L;
                (this_leaf->m_younger_brother)->m_holes += (this_leaf->m_R - this_leaf->m_L);

                int last_filled_spot_in_arr = m_youngest_leaf->m_R;

                memcpy(m_data_arr_mins[last_filled_spot_in_arr], m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
                memcpy(m_data_arr_maxes[last_filled_spot_in_arr], m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
                memcpy(m_data_arr_ids + last_filled_spot_in_arr, m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes, (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(int));


                // for(folan = 0; folan < (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes); folan++){
                //     for(filan = 0; filan < NUMDIMS; filan++){
                //         m_data_arr_mins[last_filled_spot_in_arr + folan][filan] = m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes + folan][filan];
                //         m_data_arr_maxes[last_filled_spot_in_arr + folan][filan] = m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes + folan][filan];

                //     }
                //     m_data_arr_ids[last_filled_spot_in_arr + folan] = m_data_arr_ids[this_leaf->m_L + this_leaf->m_holes + folan];
                // }

                folan = (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes);
                this_leaf->m_L = last_filled_spot_in_arr;
                last_filled_spot_in_arr += folan;

                // for(int ii = 0; ii < k; ii++){
                //     cout << "tba " << this_tba_ids[ii] << " : ";
                //     for(int jj = 0; jj < NUMDIMS; jj++) cout << this_tba_mins[ii][jj] << " " << this_tba_maxes[ii][jj] << " , ";
                //     cout << endl;
                // }
                // for(int ii = 0; ii < k; ii++){
                //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
                //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
                //     cout << endl;
                // }

                memcpy(m_data_arr_mins[last_filled_spot_in_arr], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
                memcpy(m_data_arr_maxes[last_filled_spot_in_arr], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
                memcpy(m_data_arr_ids + last_filled_spot_in_arr, tba_ids + this_tba_start, k*sizeof(int));

                // for(folan = 0; folan < k; folan++){
                //     memcpy(m_data_arr_mins[last_filled_spot_in_arr + folan], this_tba_mins+folan*NUMDIMS, NUMDIMS*sizeof(float));
                //     memcpy(m_data_arr_maxes[last_filled_spot_in_arr + folan], this_tba_maxes+folan*NUMDIMS, NUMDIMS*sizeof(float));
                    
                // }

                //  for(int ii = 0; ii < k; ii++){
                //     cout << "copied " << m_data_arr_ids[last_filled_spot_in_arr + ii] << " : ";
                //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[last_filled_spot_in_arr + ii][jj] << " " << m_data_arr_maxes[last_filled_spot_in_arr + ii][jj] << " , ";
                //     cout << endl;
                // }

                // for(folan = 0; folan < k; folan++){
                //     for(filan = 0; filan < NUMDIMS; filan++){
                //         m_data_arr_mins[last_filled_spot_in_arr + folan][filan] = this_tba_mins[folan][filan];
                //         m_data_arr_maxes[last_filled_spot_in_arr + folan][filan] = this_tba_maxes[folan][filan];
                //     }
                //     m_data_arr_ids[last_filled_spot_in_arr + folan] = this_tba_ids[folan];
                //     // cout << "ADDED " << m_data_arr_ids[last_filled_spot_in_arr + folan] << endl;
                // }

                this_leaf->m_R = last_filled_spot_in_arr + k;

                (this_leaf->m_older_brother)->m_younger_brother = this_leaf->m_younger_brother;
                (this_leaf->m_younger_brother)->m_older_brother = this_leaf->m_older_brother;

                this_leaf->m_older_brother = m_youngest_leaf;
                m_youngest_leaf->m_younger_brother = this_leaf;
                this_leaf->m_younger_brother = NULL;
                this_leaf->m_holes = 0;
                this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
                m_youngest_leaf = this_leaf;


            }


        }
    }
}





RTREE_TEMPLATE
void RTREE_QUAL::ripple_v3(vector<overlapping_leaf_tbas> *overlapping_leaves){

    Node* this_leaf;
    int this_tba_count;
    int this_tba_start;
    int k;
    int folan, filan;

    // cout << "in ripple  v2, ol.size() " << overlapping_leaves->size() << endl;
    for(int index =0; index < overlapping_leaves->size(); index++){
        this_leaf = overlapping_leaves->at(index).this_leaf;
        this_tba_start = overlapping_leaves->at(index).this_tba_start;
        this_tba_count = overlapping_leaves->at(index).this_tba_count;



        // cout << "PRINT tba" << endl;        
        // if(this_tba_count > 0) {
        //     cout << "in ripple, leaf_id: " << this_leaf->m_id << " tba.size(): " << this_tba_count << "  start: " << this_tba_start << " L, h, R: " << this_leaf->m_L << " " << this_leaf->m_holes << " " << this_leaf->m_R << endl;
        //     cout << "inseting these: !!!!!!!!!!!!!!!!!" << endl;
        //     for(int ii = this_tba_start; ii < this_tba_start + this_tba_count; ii++){
        //         cout << tba_ids[ii] << " " << tba_mins[ii][0] << " " << tba_maxes[ii][0] << endl;
        //     }
        //     cout << "!!!!!!!!!!!!!!!!!!!!!!" << endl;
        // }

        // cout << "in loop, this-tba.size() " << this_tba.size() << endl;
        // k = this_tba.size();
        k = this_tba_count;
        
        if(k == 0) continue;

        #ifdef stats
            if(this_leaf->isRegular()) count_regular_leaves--;
            else count_irregular_leaves--;
        #endif

        if(k < this_leaf->m_holes){
            // fill the holes
            memcpy(m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes - k, tba_ids + this_tba_start, k*sizeof(int));
            memcpy(m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes - k], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes - k], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
            this_leaf->m_holes -= k;
            this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
            #ifdef stats
                if(this_leaf->isRegular()) count_regular_leaves++;
                else count_irregular_leaves++;
            #endif
            continue;
        }

        if(this_leaf->m_id == m_youngest_leaf->m_id){
            // append to the end of this leaf
            memcpy(m_data_arr_ids + this_leaf->m_R, tba_ids + this_tba_start, k*sizeof(int));
            memcpy(m_data_arr_mins[this_leaf->m_R], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
            memcpy(m_data_arr_maxes[this_leaf->m_R], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
            this_leaf->m_R += k;
            this_leaf->m_count += k;
            #ifdef stats
                if(this_leaf->isRegular()) count_regular_leaves++;
                else count_irregular_leaves++;
            #endif
            continue;
        }

        // otherwise
        // then just move this leaf to the end with its tba
        (this_leaf->m_younger_brother)->m_L = this_leaf->m_L;
        (this_leaf->m_younger_brother)->m_holes += (this_leaf->m_R - this_leaf->m_L);

        int last_filled_spot_in_arr = m_youngest_leaf->m_R;

        memcpy(m_data_arr_mins[last_filled_spot_in_arr + DEFAULT_HOLES], m_data_arr_mins[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
        memcpy(m_data_arr_maxes[last_filled_spot_in_arr + DEFAULT_HOLES], m_data_arr_maxes[this_leaf->m_L + this_leaf->m_holes], NUMDIMS*(this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(float));
        memcpy(m_data_arr_ids + last_filled_spot_in_arr + DEFAULT_HOLES, m_data_arr_ids + this_leaf->m_L + this_leaf->m_holes, (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes)*sizeof(int));

        folan = (this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes) + DEFAULT_HOLES;
        this_leaf->m_L = last_filled_spot_in_arr;
        last_filled_spot_in_arr += folan;

        // for(int ii = 0; ii < k; ii++){
        //     cout << "tba " << this_tba_ids[ii] << " : ";
        //     for(int jj = 0; jj < NUMDIMS; jj++) cout << this_tba_mins[ii][jj] << " " << this_tba_maxes[ii][jj] << " , ";
        //     cout << endl;
        // }
        // for(int ii = 0; ii < k; ii++){
        //     cout << "copied " << m_data_arr_ids[this_leaf->m_L + ii] << " : ";
        //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[this_leaf->m_L + ii][jj] << " " << m_data_arr_maxes[this_leaf->m_L + ii][jj] << " , ";
        //     cout << endl;
        // }

        memcpy(m_data_arr_mins[last_filled_spot_in_arr], tba_mins[this_tba_start], k*NUMDIMS*sizeof(float));
        memcpy(m_data_arr_maxes[last_filled_spot_in_arr], tba_maxes[this_tba_start], k*NUMDIMS*sizeof(float));
        memcpy(m_data_arr_ids + last_filled_spot_in_arr, tba_ids + this_tba_start, k*sizeof(int));


        //  for(int ii = 0; ii < k; ii++){
        //     cout << "copied " << m_data_arr_ids[last_filled_spot_in_arr + ii] << " : ";
        //     for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[last_filled_spot_in_arr + ii][jj] << " " << m_data_arr_maxes[last_filled_spot_in_arr + ii][jj] << " , ";
        //     cout << endl;
        // }

        this_leaf->m_R = last_filled_spot_in_arr + k;

        if(this_leaf->m_id == m_oldest_leaf->m_id){
            (this_leaf->m_younger_brother)->m_older_brother = NULL;
            m_oldest_leaf = this_leaf->m_younger_brother;
            ASSERT((this_leaf->m_younger_brother)->m_L == 0);

        }
        else{
            (this_leaf->m_older_brother)->m_younger_brother = this_leaf->m_younger_brother;
            (this_leaf->m_younger_brother)->m_older_brother = this_leaf->m_older_brother;
        }

        this_leaf->m_older_brother = m_youngest_leaf;
        m_youngest_leaf->m_younger_brother = this_leaf;
        this_leaf->m_younger_brother = NULL;
        this_leaf->m_holes = DEFAULT_HOLES;
        this_leaf->m_count = this_leaf->m_R - this_leaf->m_L - this_leaf->m_holes;
        m_youngest_leaf = this_leaf;

        #ifdef stats
            if(this_leaf->isRegular()) count_regular_leaves++;
            else count_irregular_leaves++;
        #endif


    }





}



RTREE_TEMPLATE
bool RTREE_QUAL::CountDataPoints(Node *a_node, int &a_foundCount){
    if (a_node->m_parent == NULL) a_foundCount += m_pending_insertions_count;
    if (a_node->IsInternalNode()) {
        for (int index = 0; index < a_node->m_count; ++index) {
            CountDataPoints(a_node->m_branch[index].m_child, a_foundCount);
        }
    } else {
        // ASSERT(a_node->m_count == (a_node->m_R - a_node->m_L - a_node->m_holes));
        // if(a_node->m_count != (a_node->m_R - a_node->m_L - a_node->m_holes)){
        //     cout << "COUNT WAS WRONG!!!!!!!!!!!" << endl;
        //     cout << "id: " << a_node->m_id <<  " count: " << a_node->m_count << " real count(r - l - holes): " << (a_node->m_R - a_node->m_L - a_node->m_holes) << endl;
        //     ASSERT(a_node->m_count == (a_node->m_R - a_node->m_L - a_node->m_holes));
        // }
        // a_foundCount += a_node->m_count;
        a_foundCount += (a_node->m_R - a_node->m_L - a_node->m_holes);
    }
    return true; // Continue searching
}



RTREE_TEMPLATE
bool RTREE_QUAL::CountNodes(Node *a_node, int &a_foundCount){
    a_foundCount += 1;
    if (a_node->IsInternalNode()) {
        for (int index = 0; index < a_node->m_count; ++index) {
            CountNodes(a_node->m_branch[index].m_child, a_foundCount);
        }
    }
    return true; // Continue searching
}


RTREE_TEMPLATE
bool RTREE_QUAL::SumDataPointIDs(Node *a_node, int &a_foundCount){
    if (a_node->IsInternalNode()) {
        for (int index = 0; index < a_node->m_count; ++index) {
            SumDataPointIDs(a_node->m_branch[index].m_child, a_foundCount);
        }
    } else {
        for(int i= 0; i < a_node->m_count; i++){
            a_foundCount += a_node->m_branch[i].m_data;
        }
        //        a_foundCount += a_node->m_count;
    }
    return true; // Continue searching
}




RTREE_TEMPLATE
int RTREE_QUAL::ChooseAxis(Rect node_cover, Rect query_cover) {

    int ms[NUMDIMS];
    int current_min_m = 21474836; int current_min_index = 0;
    for(int i =0; i < NUMDIMS; i++){
        ms[i] = std::abs(2 * query_cover.m_max[i] - query_cover.m_min[i] - node_cover.m_max[i]) + std::abs(2 * query_cover.m_min[i] - query_cover.m_max[i] - node_cover.m_min[i]);
        if(ms[i] < current_min_m){
            current_min_m = ms[i];
            current_min_index = i;
        }
    }
    // int m_x = std::abs(2 * query_cover.m_max[0] - query_cover.m_min[0] - node_cover.m_max[0]) + std::abs(2 * query_cover.m_min[0] - query_cover.m_max[0] - node_cover.m_min[0]);
    // int m_y = std::abs(2 * query_cover.m_max[1] - query_cover.m_min[1] - node_cover.m_max[1]) + std::abs(2 * query_cover.m_min[1] - query_cover.m_max[1] - node_cover.m_min[1]);

    // if(m_x < m_y){return 0;}
    // return 1;
    if(current_min_index == -1){cout << "IN CHOOSEAXIS, STILL -1, VERY WRONG!!" << endl;}
    return current_min_index;

}

// return the axis that has largest dif between min and max
RTREE_TEMPLATE
int RTREE_QUAL::ChooseAxisAgain(Rect node_cover, Rect query_cover, int ignore_this_dim) {

    float this_ms;
    float current_max_m = std::abs(node_cover.m_max[0] - node_cover.m_min[0]);
    int current_max_index = 0;
    for(int i =1; i < NUMDIMS; i++){
        this_ms = std::abs(node_cover.m_max[i] - node_cover.m_min[i]);
        if(this_ms > current_max_m && i != ignore_this_dim){
            current_max_m = this_ms;
            current_max_index = i;
        }
    }
    return current_max_index;

}



RTREE_TEMPLATE
int RTREE_QUAL::ChooseAxis_weird(Rect node_cover, Rect query_cover) {

    //    int ms[NUMDIMS];
    //    int current_min_m = 21474836; int current_min_index = 0;
    //    for(int i =0; i < NUMDIMS; i++){
    //        ms[i] = std::abs(2 * query_cover.m_max[i] - query_cover.m_min[i] - node_cover.m_max[i]) + std::abs(2 * query_cover.m_min[i] - query_cover.m_max[i] - node_cover.m_min[i]);
    //        if(ms[i] < current_min_m){
    //            current_min_m = ms[i];
    //            current_min_index = i;
    //        }
    //    }


    float this_ms;
    float current_max_m = std::abs(query_cover.m_max[0] - query_cover.m_min[0]);
    //    cout << "in choose, first " << current_max_m << endl;
    int current_max_index = 0;
    for(int i =1; i < NUMDIMS; i++){
        this_ms = std::abs(query_cover.m_max[i] - query_cover.m_min[i]);
        //        cout << "in choose, next " << this_ms << endl;
        if(this_ms > current_max_m){
            current_max_m = this_ms;
            current_max_index = i;
        }
    }
    // int m_x = std::abs(2 * query_cover.m_max[0] - query_cover.m_min[0] - node_cover.m_max[0]) + std::abs(2 * query_cover.m_min[0] - query_cover.m_max[0] - node_cover.m_min[0]);
    // int m_y = std::abs(2 * query_cover.m_max[1] - query_cover.m_min[1] - node_cover.m_max[1]) + std::abs(2 * query_cover.m_min[1] - query_cover.m_max[1] - node_cover.m_min[1]);

    // if(m_x < m_y){return 0;}
    // return 1;
    //    if(current_min_index == -1){cout << "IN CHOOSEAXIS, STILL -1, VERY WRONG!!" << endl;}
    float mid_of_node_cover_axis = (float) ((node_cover.m_max[current_max_index] + node_cover.m_min[current_max_index]) / 2);

    if(std::abs(query_cover.m_min[current_max_index] - mid_of_node_cover_axis) <
    std::abs(query_cover.m_max[current_max_index] - mid_of_node_cover_axis)){
        //        cout << "chose min\n";
        return (2 * current_max_index);
    } else{
        //        cout << "chose max\n";
        return (2 * current_max_index + 1);
    }

    //    return current_min_index;

}



RTREE_TEMPLATE
int RTREE_QUAL::ChooseAxisLargestSideMid(Rect node_cover, Rect query_cover) {

    float this_ms;
    float current_max_m = std::abs(node_cover.m_max[0] - node_cover.m_min[0]);
    //    cout << "in choose, first " << current_max_m << endl;
    int current_max_index = 0;
    for(int i =1; i < NUMDIMS; i++){
        this_ms = std::abs(node_cover.m_max[i] - node_cover.m_min[i]);
        //        cout << "in choose, next " << this_ms << endl;
        if(this_ms > current_max_m){
            current_max_m = this_ms;
            current_max_index = i;
        }
    }
    //    cout << "in choose, last " << current_max_m << endl;

    // now we have to figure out if q_min[current_max_index] ids closer to mid of node_cover
    float mid_of_node_cover_axis = (float) ((node_cover.m_max[current_max_index] + node_cover.m_min[current_max_index]) / 2);
    //    cout << "mid of node cover axis " << mid_of_node_cover_axis << endl;
    //    Rect query_cover_in_node;
    //    for(int i = 0; i < NUMDIMS; i++){
    //        query_cover_in_node.m_min[i] = std::max(query_cover.m_min[i], node_cover.m_min[i]);
    //        query_cover_in_node.m_max[i] = std::min(query_cover.m_max[i], node_cover.m_max[i]);
    ////        cout << "query cover in node, axis " << i << " min " << query_cover_in_node.m_min[i] << " max " << query_cover_in_node.m_max[i] << endl;
    //    }

    //    if(std::abs(query_cover_in_node.m_min[current_max_index] - mid_of_node_cover_axis) <
    //       std::abs(query_cover_in_node.m_max[current_max_index] - mid_of_node_cover_axis)){
    ////        cout << "chose min\n";
    //        return (2 * current_max_index);
    //    } else{
    ////        cout << "chose max\n";
    //        return (2 * current_max_index + 1);
    //    }
    if(std::abs(query_cover.m_min[current_max_index] - mid_of_node_cover_axis) <
    std::abs(query_cover.m_max[current_max_index] - mid_of_node_cover_axis)){
//                cout << "chose min and axis " << current_max_index << endl;
        return (2 * current_max_index);
    } else{
//                cout << "chose max and axis " << current_max_index << endl;
        return (2 * current_max_index + 1);
    }
}



RTREE_TEMPLATE
int RTREE_QUAL::CrackOnAxisCompMin(int L, int R, const float& crack_value, int axis) {
    //    return std::partition(arr+L, arr+R, [crack_value, axis](const Branch& em){ return em.m_rect.m_min[axis] > crack_value; }) - m_data_arr;
    return MyPartition(L, R, crack_value, axis, true);
}


RTREE_TEMPLATE
int RTREE_QUAL::CrackOnAxisCompMax(int L, int R, const float& crack_value, int axis) {
    //    return std::partition(arr+L, arr+R, [crack_value, axis](const Branch& em){ return em.m_rect.m_max[axis] < crack_value; }) - m_data_arr;
    return MyPartition(L, R, crack_value, axis, false);
}

RTREE_TEMPLATE
int RTREE_QUAL::CrackOnAxisCompMin_v2(int L, int R, const float& crack_value, int axis, float *cover_mins_left, float *cover_maxes_left, float *cover_mins_right, float *cover_maxes_right) {
    return MyPartition_v2(L, R, crack_value, axis, true, cover_mins_left, cover_maxes_left, cover_mins_right, cover_maxes_right);
}


RTREE_TEMPLATE
int RTREE_QUAL::CrackOnAxisCompMax_v2( int L, int R, const float& crack_value, int axis, float *cover_mins_left, float *cover_maxes_left, float *cover_mins_right, float *cover_maxes_right) {
    return MyPartition_v2(L, R, crack_value, axis, false, cover_mins_left, cover_maxes_left, cover_mins_right, cover_maxes_right);
}

RTREE_TEMPLATE
int RTREE_QUAL::CrackOnAxisCompMin_v3( int L, int R, const float& crack_value, int axis, Rect *left_rect, Rect *right_rect) {
    return MyPartition_v3(L, R, crack_value, axis, true, left_rect, right_rect);
}


RTREE_TEMPLATE
int RTREE_QUAL::CrackOnAxisCompMax_v3( int L, int R, const float& crack_value, int axis, Rect *left_rect, Rect *right_rect) {
    return MyPartition_v3(L, R, crack_value, axis, false, left_rect, right_rect);
}

RTREE_TEMPLATE
int RTREE_QUAL::CrackOnAxisCompMin_v4( int L, int R, const float& crack_value, int axis, Rect *left_rect, Rect *right_rect) {
    return MyPartition_v4(L, R, crack_value, axis, true, left_rect, right_rect);
}


RTREE_TEMPLATE
int RTREE_QUAL::CrackOnAxisCompMax_v4( int L, int R, const float& crack_value, int axis, Rect *left_rect, Rect *right_rect) {
    return MyPartition_v4(L, R, crack_value, axis, false, left_rect, right_rect);
}


RTREE_TEMPLATE
int RTREE_QUAL::CrackOnAxisCompMin_v45( int L, int R, const float& crack_value, int axis, Rect *left_rect, Rect *right_rect) {
    return MyPartition_v45(L, R, crack_value, axis, true, left_rect, right_rect);
}


RTREE_TEMPLATE
int RTREE_QUAL::CrackOnAxisCompMax_v45( int L, int R, const float& crack_value, int axis, Rect *left_rect, Rect *right_rect) {
    return MyPartition_v45(L, R, crack_value, axis, false, left_rect, right_rect);
}



RTREE_TEMPLATE
void RTREE_QUAL::swap(float &a, float &b) {
    float t = a;
    a = b;
    b = t;
}

//RTREE_TEMPLATE
//void RTREE_QUAL::swap(Branch *a, Branch *b) {
//    Branch t = *a;
//    *a = *b;
//    *b = t;
//}

RTREE_TEMPLATE
void RTREE_QUAL::swap_index(int i, int j) {
    //    if(i != j) {
    //        cout << "BEFORE SWAP" << endl;
    //        cout << i << " " << j << " " << m_data_arr[i].m_data << " " << m_data_arr[i].m_rect.m_min[0] << " "
    //             << m_data_arr[j]->m_data << " " << m_data_arr[j]->m_rect.m_min[0] << endl;
    //    }

    // older than 28sep
//    Branch t;
//    CopyBranch(t, m_data_arr[i]);
//    CopyBranch(m_data_arr[i], m_data_arr[j]);
//    CopyBranch(m_data_arr[j], t);
    // older than 28sep
    //    if(i != j) {
    //        cout << "AFTER SWAP" << endl;
    //        cout << i << " " << j << " " << m_data_arr[i].m_data << " " << m_data_arr[i].m_rect.m_min[0] << " "
    //             << m_data_arr[j]->m_data << " " << m_data_arr[j]->m_rect.m_min[0] << endl;
    //    }
    //    m_data_arr[j] = t;

    // 28sep
//    Branch* t;
//    t = m_data_arr[i];
//    m_data_arr[i] = m_data_arr[j];
//    m_data_arr[j] = t;
    // 28sep

    // 5oct
//    int temp = m_data_arr_datainfo[i];
//    m_data_arr_datainfo[i] = m_data_arr_datainfo[j];
//    m_data_arr_datainfo[j] = temp;

    //6oct
    //id
    std::swap(m_data_arr_ids[i], m_data_arr_ids[j]);
    //mins
    std::swap(m_data_arr_mins[i], m_data_arr_mins[j]);
    //maxes
    std::swap(m_data_arr_maxes[i], m_data_arr_maxes[j]);




}




RTREE_TEMPLATE
int RTREE_QUAL::
MyPartition_v45(int low, int high, float pivot_value, int axis, bool min_or_max, Rect *left_rect, Rect *right_rect){
    int x1 = low; int x2 = high - 1;
//    int k = 0;
    Rect x1_rect, x2_rect;

//    for(int ashghal=0; ashghal < NUMDIMS; ashghal++){
//        if(left_rect->m_min[ashghal] > left_rect->m_max[ashghal]){
//            cout << "IN PARTITION LEFT COVER MIN " << left_rect->m_min[ashghal] << " max " << left_rect->m_max[ashghal] << endl;
//        }
//
//        if(right_rect->m_min[ashghal] > right_rect->m_max[ashghal]){
//            cout << "IN PARTITION RIGHT COVER MIN " << right_rect->m_min[ashghal] << " max " << right_rect->m_max[ashghal] << endl;
//        }
//    }

    if(min_or_max){
        while(x1 <= x2 && x2 > 0){
            if (m_data_arr_mins[x1][axis] <= pivot_value) {
                // update left cover

//                if(m_data_arr_mins[x1][axis] < left_rect->m_min[axis]){left_rect->m_min[axis] = m_data_arr_mins[x1][axis];}
                if(m_data_arr_maxes[x1][axis] > left_rect->m_max[axis]){left_rect->m_max[axis] = m_data_arr_maxes[x1][axis];}

                x1++;
            }
            else{
                while(x2 > 0 && x2 >= x1 && m_data_arr_mins[x2][axis] > pivot_value){
                    // update right cover

                    if(m_data_arr_mins[x2][axis] < right_rect->m_min[axis]){right_rect->m_min[axis] = m_data_arr_mins[x2][axis];}
                    if(m_data_arr_maxes[x2][axis] > right_rect->m_max[axis]){right_rect->m_max[axis] = m_data_arr_maxes[x2][axis];}

                    x2--;
                }

                if(x1 < x2){
                    swap_index(x1, x2);
                    // update left cover

//                    if(m_data_arr_mins[x1][axis] < left_rect->m_min[axis]){left_rect->m_min[axis] = m_data_arr_mins[x1][axis];}
                    if(m_data_arr_maxes[x1][axis] > left_rect->m_max[axis]){left_rect->m_max[axis] = m_data_arr_maxes[x1][axis];}

                    // update right cover

                    if(m_data_arr_mins[x2][axis] < right_rect->m_min[axis]){right_rect->m_min[axis] = m_data_arr_mins[x2][axis];}
                    if(m_data_arr_maxes[x2][axis] > right_rect->m_max[axis]){right_rect->m_max[axis] = m_data_arr_maxes[x2][axis];}

                    x2--; x1++;
                }
            }
        }
    }
    else{
        while(x1 <= x2 && x2 > 0){
            if (m_data_arr_maxes[x1][axis] < pivot_value) {
                // update left cover

                if(m_data_arr_mins[x1][axis] < left_rect->m_min[axis]){
                    left_rect->m_min[axis] = m_data_arr_mins[x1][axis];
//                    if(left_rect->m_min[axis] > left_rect->m_max[axis] && left_rect->m_min[axis] != 21474836 && left_rect->m_max[axis] != -21474836){
//                        cout << "IN PARTITION LEFT COVER MIN after update " << left_rect->m_min[axis] << " max " << left_rect->m_max[axis] << " axis  " << axis << endl;
//
//                    }
                }
                if(m_data_arr_maxes[x1][axis] > left_rect->m_max[axis]){
                    left_rect->m_max[axis] = m_data_arr_maxes[x1][axis];
//                    if(left_rect->m_min[axis] > left_rect->m_max[axis] && left_rect->m_min[axis] != 21474836 && left_rect->m_max[axis] != -21474836){
//                        cout << "IN PARTITION LEFT COVER MIN after update " << left_rect->m_min[axis] << " max " << left_rect->m_max[axis] << " axis  " << axis << endl;
//
//                    }
                }

                x1++;

            }
            else{
                while(x2 > 0 && x2 >= x1 && m_data_arr_maxes[x2][axis] >= pivot_value){
                    // update right cover
                    if(m_data_arr_mins[x2][axis] < right_rect->m_min[axis]){
                        right_rect->m_min[axis] = m_data_arr_mins[x2][axis];
//                        if(right_rect->m_min[axis] > right_rect->m_max[axis] && right_rect->m_min[axis] != 21474836 && right_rect->m_max[axis] != -21474836){
//                            cout << "IN PARTITION RIGHT COVER MIN after update " << right_rect->m_min[axis] << " max " << right_rect->m_max[axis] << " axis  " << axis << endl;
//                        }
                    }
//                    if(m_data_arr_maxes[x2][axis] > right_rect->m_max[axis]){
//                        right_rect->m_max[axis] = m_data_arr_maxes[x2][axis];
////                        if(right_rect->m_min[axis] > right_rect->m_max[axis] && right_rect->m_min[axis] != 21474836 && right_rect->m_max[axis] != -21474836){
////                            cout << "IN PARTITION RIGHT COVER MIN after update " << right_rect->m_min[axis] << " max " << right_rect->m_max[axis] << " axis  " << axis << endl;
////                        }
//                    }



                    x2--;
                }
                if(x1 < x2){
                    swap_index(x1, x2);
                    // update left cover

                    if(m_data_arr_mins[x1][axis] < left_rect->m_min[axis]){
                        left_rect->m_min[axis] = m_data_arr_mins[x1][axis];
//                        if(left_rect->m_min[axis] > left_rect->m_max[axis] && left_rect->m_min[axis] != 21474836 && left_rect->m_max[axis] != -21474836){
//                            cout << "IN PARTITION LEFT COVER MIN after update " << left_rect->m_min[axis] << " max " << left_rect->m_max[axis] << " axis  " << axis << endl;
//
//                        }
                    }
                    if(m_data_arr_maxes[x1][axis] > left_rect->m_max[axis]){
                        left_rect->m_max[axis] = m_data_arr_maxes[x1][axis];
//                        if(left_rect->m_min[axis] > left_rect->m_max[axis] && left_rect->m_min[axis] != 21474836 && left_rect->m_max[axis] != -21474836){
//                            cout << "IN PARTITION LEFT COVER MIN after update " << left_rect->m_min[axis] << " max " << left_rect->m_max[axis] << " axis  " << axis << endl;
//
//                        }
                    }

                    // update right cover
                    if(m_data_arr_mins[x2][axis] < right_rect->m_min[axis]){
                        right_rect->m_min[axis] = m_data_arr_mins[x2][axis];
//                        if(right_rect->m_min[axis] > right_rect->m_max[axis] && right_rect->m_min[axis] != 21474836 && right_rect->m_max[axis] != -21474836){
//                            cout << "IN PARTITION RIGHT COVER MIN after update " << right_rect->m_min[axis] << " max " << right_rect->m_max[axis] << " axis  " << axis << endl;
//                        }
                    }
//                    if(m_data_arr_maxes[x2][axis] > right_rect->m_max[axis]){
//                        right_rect->m_max[axis] = m_data_arr_maxes[x2][axis];
////                        if(right_rect->m_min[axis] > right_rect->m_max[axis] && right_rect->m_min[axis] != 21474836 && right_rect->m_max[axis] != -21474836){
////                            cout << "IN PARTITION RIGHT COVER MIN after update " << right_rect->m_min[axis] << " max " << right_rect->m_max[axis] << " axis  " << axis << endl;
////                        }
//                    }


                    x2--; x1++;


                }
            }
        }
    }

    return x1;
}



RTREE_TEMPLATE
bool RTREE_QUAL::SearchNN(Rect query_rect, Node **result) {
    priority_queue<pair<Node *, Rect>, vector<pair<Node *, Rect>>, CompareNode> Q;
    // maybe we should have done it with nodes, not branches
    Q.push(make_pair(m_root, query_rect));
    while (!Q.empty()) {
        pair<Node *, Rect> element = Q.top();
        Q.pop();
        //        if(element.first->IsLeaf())
        if (element.first->m_level == 1) {
            *result = element.first;
            //            Rect test = NodeCover(result);
            //            printf("IN SEARCHNN: %d, %d, %d, %d\n", test.m_min[0], test.m_min[1], test.m_max[0], test.m_max[1]);
            return true;
        } else {
            for(int i =0; i  < element.first->m_count; i++){
                //            for (typename std::vector<Branch>::iterator it = element.first->m_branch.begin();
                //                 it != element.first->m_branch.end(); it++) {
                // *it ~ Branch*
                //                Q.push(make_pair(it->m_child, query_rect));
                Q.push(make_pair(element.first->m_child, query_rect));
            }
        }
    }
    return false;

}

RTREE_TEMPLATE
bool RTREE_QUAL::SearchNNtester(const float a_min[NUMDIMS], const float a_max[NUMDIMS]) {
    Rect rect;

    for (int axis = 0; axis < NUMDIMS; ++axis) {
        rect.m_min[axis] = a_min[axis];
        rect.m_max[axis] = a_max[axis];
    }

    return SearchNN(rect, NULL);
}

RTREE_TEMPLATE
void RTREE_QUAL::CopyBranch(Branch &current, const Branch &other){
    current.m_rect = other.m_rect;
    current.m_child = other.m_child;
    current.m_data = other.m_data;
}


RTREE_TEMPLATE
int RTREE_QUAL::getNodeBranchIndex(Node *a_node) {
    //    ASSERT(a_node);
    //    cout << a_node->m_id << " has " << a_node->m_count << " children\n";
    //    cout << a_node->m_id << " parent children count in getnodebranchindex " << (a_node->m_parent)->m_branch.size() << ", " << (a_node->m_parent)->m_count <<"\n";
    //    cout << "checking levels: " << a_node->m_level << " " << (a_node->m_parent)->m_level << "\n";
    // JUST FOR CHECKING
    //    for(int i = 0; i < (a_node->m_parent)->m_branch.size(); i++) {
    //        cout << "in getnodebranchindex for checkkkkk " << ((a_node->m_parent)->m_branch[i].m_child)->m_id << "\n";
    //    }
    for(int i = 0; i < (a_node->m_parent)->m_count; i++) {
        //        cout << "in getnodebranchindex for " << ((a_node->m_parent)->m_branch[i].m_child)->m_id << "\n";
        if (a_node->m_id == ((a_node->m_parent)->m_branch[i].m_child)->m_id) {
            return i;
        }
    }
    return -1;
}

RTREE_TEMPLATE
void RTREE_QUAL::Insert_anylevel(const Branch &a_branch, Node *start, int a_level){
    //    ASSERT(start);
    //    ASSERT(start->m_level >= a_level);
    ASSERT(a_branch.m_child->m_count > 0);
    //    cout << "in insert_any \n";
    if(start->m_level > a_level){
        // printf("calling pb in insert_anylevel\n");
        int index = PickBranch(&a_branch.m_rect, start);

        Insert_anylevel(a_branch, start->m_branch[index].m_child, a_level);
    }
    else{
        // start->m_level == a_level
        Node * current_node = start;
        Branch current_branch;
        CopyBranch(current_branch, a_branch);
        bool is_split = false; // so the loop starts
        bool just_go_till_top = false; // to indicate that no more adds are required, but we should recurse up to update the rects
        while(true){
            //            cout << "in while\n";
            ASSERT(current_branch.m_child->m_count > 0);
            // last current_node was not root
            if(!just_go_till_top) {
                //                cout << "in while not just till top\n";
                Node *newNode = AllocNode();
                Rect rect1; Rect rect2;
                // 2023
                if(current_branch.m_child->m_count == 0){
                    cout << "THIS IS WHERE YOU FUCKED UP early." << endl;
                    cout << "parent, node id " << (current_node->m_parent)->m_id << " " << current_branch.m_child->m_id << " parent.count " << (current_node->m_parent)->m_count << " currentBranch.count " << current_branch.m_child->m_count << endl;
                }
                // 2023
                is_split = AddBranch(&current_branch, current_node, &newNode, rect1, rect2);

                if (is_split) {
                    //                    cout << "in while was split\n";
                    // was split
                    if (current_node->m_parent != NULL) {
                        // 2023
                        if(newNode->m_count == 0){
                            cout << "THIS IS WHERE YOU FUCKED UP." << endl;
                            cout << "parent, node id " << (current_node->m_parent)->m_id << " " << newNode->m_id << " parent.count " << (current_node->m_parent)->m_count << " currentBranch.count " << current_branch.m_child->m_count << endl;
                        }
                        ASSERT(newNode->m_count > 0);
                        
                        // 2023
                        //                        cout << "in while is split, not root\n";
                        // update nodecover rectangles
                        //                        cout << "CHECK POINT 1\n";
                        int index = getNodeBranchIndex(current_node);
                        //                        (current_node->m_parent)->m_branch[index].m_rect = NodeCover(current_node);
                        (current_node->m_parent)->m_branch[index].m_rect = rect1;
                        // set parameters to move up the tree
                        current_branch.m_child = newNode;
                        
                        //                        current_branch.m_rect = NodeCover(newNode);
                        current_branch.m_rect = rect2;
                        current_node = current_node->m_parent;
                    } else {
                        //                        cout << "in while is split and root\n";
                        // current_node is root and it was split
                        // now we have to make a new root

                        Node *newRoot = AllocNode();
                        newRoot->m_level = current_node->m_level + 1;
                        newRoot->m_parent = NULL;
                        //                        newRoot->name = "new root in my insert";

                        current_branch.m_child = newNode;
                        //                        current_branch.m_rect = NodeCover(newNode);
                        current_branch.m_rect = rect2;

                        AddBranch(&current_branch, newRoot, NULL);

                        current_branch.m_child = current_node;
                        //                        current_branch.m_rect = NodeCover(current_node);
                        current_branch.m_rect = rect1;

                        AddBranch(&current_branch, newRoot, NULL);

                        m_root = newRoot;
                        break;
                    }
                } else {
                    //                    cout << "in while not split\n";
                    // was not split
                    if (current_node->m_parent != NULL) {
                        //                        cout << "in while not split and not root\n";
                        // we are not at the root
                        //                        cout << "CHECK POINT 2 " << current_node->m_id << " with parent:: " << (current_node->m_parent)->m_id << "\n" ;
                        int index = getNodeBranchIndex(current_node);
                        (current_node->m_parent)->m_branch[index].m_rect = CombineRect(&(current_branch.m_rect),
                                                                                       &((current_node->m_parent)->m_branch[index].m_rect));
                        current_node = current_node->m_parent;
                        just_go_till_top = true;
                    } else {

                        //                        cout << "THIS SHOULD NOT HAPPEN, or maybe it could\n";
                        break;
                    }
                }
            } else{
                //                cout << "in while just till top\n";
                if (current_node->m_parent != NULL) {
                    //                    cout << "in while just till top not root\n";
                    //                    cout << "CHECK POINT 3\n";
                    int index = getNodeBranchIndex(current_node);
                    //                    cout << "after get index\n";
                    (current_node->m_parent)->m_branch[index].m_rect = CombineRect(&(current_branch.m_rect),
                                                                                   &((current_node->m_parent)->m_branch[index].m_rect));
                    //                    (current_node->m_parent)->m_branch[index].m_rect = NodeCover(current_node);

                    current_node = current_node->m_parent;
                } else{
                    //                    cout << "in while just till top root\n";
                    break;
                }
            }

        }

    }
}

RTREE_TEMPLATE
void RTREE_QUAL::printTree(string file_name) {
    //    string file_name= "tree.txt";
    ofstream myfile;
    myfile.open(file_name.c_str());
    if(myfile.is_open()){
        //        myfile <<
        printTreeRec(m_root, myfile);
    }
}

RTREE_TEMPLATE
void RTREE_QUAL::printTreeRec(Node* a_node, ofstream &myfile){
    if (a_node->IsInternalNode()) {

        for (int index = 0; index < a_node->m_count; ++index) {
            myfile << a_node->m_branch[index].m_child->m_level << " " << a_node->m_branch[index].m_rect.m_min[0] << " " << a_node->m_branch[index].m_rect.m_min[1] << " " << a_node->m_branch[index].m_rect.m_max[0] << " " << a_node->m_branch[index].m_rect.m_max[1] << "\n";
            printTreeRec(a_node->m_branch[index].m_child, myfile);
        }
    }

}

RTREE_TEMPLATE
void RTREE_QUAL::printTreeDebugRec(Node* a_node, ofstream &myfile){
    if (a_node->IsInternalNode()) {
        myfile << a_node->m_id << " has " << a_node->m_count << "children" << endl;
        for (int index = 0; index < a_node->m_count; ++index) {
            printTreeDebugRec(a_node->m_branch[index].m_child, myfile);
        }
    } else{
        if(a_node->m_parent != NULL) {
            myfile << a_node->m_id << " son of " << (a_node->m_parent)->m_id << " : "
            << a_node->m_L << " " << a_node->m_R << "\n";
        } else{
            myfile << a_node->m_id << " is root : "
            << a_node->m_L << " " << a_node->m_R << "\n";
        }
    }

}


RTREE_TEMPLATE
void RTREE_QUAL::PrintLeafSizesRec(Node* a_node, ofstream &myfile){
    if (a_node->IsInternalNode()) {
        //        myfile << a_node->m_id << " has " << a_node->m_count << "children" << endl;
        for (int index = 0; index < a_node->m_count; ++index) {
            PrintLeafSizesRec(a_node->m_branch[index].m_child, myfile);
        }
    } else{
        myfile << a_node->m_id << " " << a_node->m_count << "\n";
    }

}

RTREE_TEMPLATE
void RTREE_QUAL::printTreeDebug(string file_name) {
    //    string file_name= "tree.txt";
    ofstream myfile;
    myfile.open(file_name.c_str());
    if(myfile.is_open()){
        //        myfile <<
        printTreeDebugRec(m_root, myfile);
    }
}

RTREE_TEMPLATE
void RTREE_QUAL::printLeafSizes(string file_name) {
    //    string file_name= "tree.txt";
    ofstream myfile;
    myfile.open(file_name.c_str());
    if(myfile.is_open()){
        //        myfile <<
        PrintLeafSizesRec(m_root, myfile);
    }
}


RTREE_TEMPLATE
void RTREE_QUAL::printLeafLinkedList(){
    cout << "LEAF LINKED LIST ---------------" << endl;
    Node *current_node = m_oldest_leaf;
    while(current_node != NULL){
        cout << "node " << current_node->m_id << " L: " << current_node->m_L << " R: " << current_node->m_R << " holes: " << current_node->m_holes << " count: " << current_node->m_count << endl;
        current_node = current_node->m_younger_brother;
    }
    cout << "END of LEAF LINKED LIST ---------------" << endl;
}

RTREE_TEMPLATE
void RTREE_QUAL::printIDs(){
    cout << "pending insertion ids: ... " << endl;
    for(int i = 0; i < m_pending_insertions_count; i++){
        cout << m_pending_insertions_ids[i] << " " ;
    }
    cout << endl << "----------" << endl;
    cout << "data arr ids: ... " << endl;
    for(int i = 0; i < 2 * DATA_COUNT; i++){
        cout << m_data_arr_ids[i] << " ";
    }
    cout << endl << "----------" << endl;

}


RTREE_TEMPLATE
void RTREE_QUAL::printLeafBoundsOfData(int data_index){
    printLeafBoundsOfDataRec(m_root, data_index);
}

RTREE_TEMPLATE
bool RTREE_QUAL::printLeafBoundsOfDataRec(Node *a_node, int data_index){
    if(a_node->IsInternalNode()){
        for(int i = 0; i < a_node->m_count; i++){
            bool print_or_not = printLeafBoundsOfDataRec(a_node->m_branch[i].m_child, data_index);
            if(print_or_not){
                cout << "level: " << a_node->m_level << " : " ;
                cout << "child's level: " << a_node->m_branch[i].m_child->m_level << " : ";
                for(int j = 0; j < NUMDIMS; j++){
                    cout << a_node->m_branch[i].m_rect.m_min[j] << " " << a_node->m_branch[i].m_rect.m_max[j] << " , ";
                }
                cout << endl;
                // // right brother rect
                // Node *right_leaf = a_node->m_branch[i].m_child->m_younger_brother;
                // if(right_leaf != NULL){
                //     int right_leaf_index = getNodeBranchIndex(right_leaf);
                //     Rect right_leaf_rect = (right_leaf->m_parent)->m_branch[right_leaf_index].m_rect;
                //     cout << "child's right brother bounds: ";
                //     for(int j = 0; j < NUMDIMS; j++){
                //         cout << right_leaf_rect.m_min[j] << " " << right_leaf_rect.m_max[j] << " , ";
                //     }
                //     cout << endl;
                // }

                // // left brother rect
                // Node *left_leaf = a_node->m_branch[i].m_child->m_older_brother;
                // if(left_leaf!= NULL){
                //     int left_leaf_index = getNodeBranchIndex(left_leaf);
                //     Rect left_leaf_rect = (left_leaf->m_parent)->m_branch[left_leaf_index].m_rect;
                //     cout << "child's left brother bounds: ";
                //     for(int j = 0; j < NUMDIMS; j++){
                //         cout << left_leaf_rect.m_min[j] << " " << left_leaf_rect.m_max[j] << " , ";
                //     }
                //     cout << endl;
                // }

                // cout << endl;
                return true;
            }
        }
    }
    else{
        for(int j = a_node->m_L + a_node->m_holes; j < a_node->m_R; j++){
            if(m_data_arr_ids[j] == data_index){
                cout << "data was: " << endl;
                for(int jj = 0; jj < NUMDIMS; jj++) cout << m_data_arr_mins[j][jj] << " " << m_data_arr_maxes[j][jj] << " ,";
                cout << endl; 
                return true;
            }
        }
        return false;
    }
    return false;
}



RTREE_TEMPLATE
void RTREE_QUAL::findAllLeavesWithDataInThem(int data_index){
    cout << "calling falwdit" << endl;
    Node* this_leaf = m_oldest_leaf;
    while(this_leaf != NULL){
        // if(this_leaf->m_id == 15372){
        //     cout << "cheacking leaf 15372  with bounds: L, h, R " << this_leaf->m_L << " " << this_leaf->m_holes << " " << this_leaf->m_R << endl;
        //     cout << "items in here: .........." << endl;
        //     for(int o = this_leaf->m_L + this_leaf->m_holes; o < this_leaf->m_R; o++){
        //         cout << m_data_arr_ids[o] << endl;
        //     }
        //     cout << "..........." << endl;
        // }
        // do some stuff
        for(int i = this_leaf->m_L + this_leaf->m_holes; i < this_leaf->m_R; i++){
            if(m_data_arr_ids[i] == data_index){
                cout << "data in leaf " << this_leaf->m_id << " is item " << i << " and L, holes, R: " << this_leaf->m_L << " " << this_leaf->m_holes << " " << this_leaf->m_R << endl;
                cout << "data bounds: " << m_data_arr_mins[i][0] << " " << m_data_arr_maxes[i][0] << ", " << m_data_arr_mins[i][1] << " " << m_data_arr_maxes[i][1] << endl;
                if(this_leaf->m_parent != NULL){
                    int branch_index = getNodeBranchIndex(this_leaf);
                    Rect cover = (this_leaf->m_parent)->m_branch[branch_index].m_rect;
                    cout << "leaf bound: " << cover.m_min[0] << " " << cover.m_max[0] << " " << cover.m_min[1]  << " " << cover.m_max[1] << endl;
                }
            }
        }
        // move on to the next
        this_leaf = this_leaf->m_younger_brother;
    }

    // also look in pending
    for(int i = 0; i < m_pending_insertions_count; i++){
        if(m_pending_insertions_ids[i] == data_index){
            cout << "data is object " << i << " in pending" << endl;
            cout << "data bounds: " << m_pending_insertions_mins[i][0] << " " << m_pending_insertions_maxes[i][0] << ", " << m_pending_insertions_mins[i][1] << " " << m_pending_insertions_maxes[i][1] << endl;

        }
    }
    cout << "FINISHED\n";
}

RTREE_TEMPLATE
void RTREE_QUAL::lookForEmptyNodesRec(Node *a_node){

    if(a_node->IsInternalNode()){
        if(a_node->m_count > 0){
            for(int i = 0; i < a_node->m_count; i++){
                    lookForEmptyNodesRec((a_node->m_branch[i]).m_child);
            }
        }
        else{
            cout << "node " << a_node->m_id << " was empty. very very bad. not recommended. BOOHOO" << endl;
            cout << "info: count = " << a_node->m_count << " level: " << a_node->m_level << " MAX-lvl " << m_root->m_level << endl;
            // return;
            exit(3);
        }
    }
    return;
    


}

RTREE_TEMPLATE
void RTREE_QUAL::lookForEmptyNodes(){
    lookForEmptyNodesRec(m_root);
}

RTREE_TEMPLATE
void RTREE_QUAL::findNodeCount(int node_id, Node* this_node){
    if(this_node->m_id == node_id){
        cout << "found node, has " << this_node->m_count << " count, lvl: " << this_node->m_level << endl;
        return;
    }
    if(this_node->IsInternalNode()){
        if(this_node->m_count > 0){
            for(int i = 0; i < this_node->m_count; i++){
                    findNodeCount(node_id, (this_node->m_branch[i]).m_child);
            }
        }
    }

    // cout << "finished searching for node based on id" << endl;
}


RTREE_TEMPLATE
std::pair<int, int> RTREE_QUAL::CountLeaves(){
    int reg_counter = 0; int irreg_counter = 0;
    Node* this_leaf = m_oldest_leaf;
    while(this_leaf != NULL){
        if(this_leaf->isRegular())  reg_counter++;
        else irreg_counter++;
        // go to next leaf
        this_leaf = this_leaf->m_younger_brother;
    }
    return std::make_pair(reg_counter, irreg_counter);
}



RTREE_TEMPLATE
void RTREE_QUAL::getLeafArea(ofstream &file, Node* a_node)
{
  if(a_node->m_level > 1)
  {
    for(int i = 0; i < a_node->m_count;i++)
    {
      getLeafArea(file, a_node->m_branch[i].m_child);
    }
  }
  else if(a_node->m_level == 1)
  {
    for(int i = 0; i < a_node->m_count;i++)
    {
      double area = (a_node->m_branch[i].m_rect.m_max[0] - a_node->m_branch[i].m_rect.m_min[0])*(a_node->m_branch[i].m_rect.m_max[1] - a_node->m_branch[i].m_rect.m_min[1]);
      file << area << endl;
    }
  }
  return;
}

RTREE_TEMPLATE
void RTREE_QUAL::getLeafArea(string filename)
{
  ofstream f;
  f.open(filename.c_str()); 
  getLeafArea(f, m_root);
  f.close();
}


#undef RTREE_TEMPLATE
#undef RTREE_QUAL

#endif //RTREE_723OCT21_H

