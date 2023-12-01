#ifndef RTREE_723OCT21_H
#define RTREE_723OCT21_H

// copied from 722oct on 19oct
// no longer checks the far side of the query piece in the partition function
// also has the HybridCover

// this one changes the choice cover after failed crack


// NOTE This file compiles under MSVC 6 SP5 and MSVC .Net 2003 it may not work on other compilers without modification.

// NOTE These next few lines may be win32 specific, you may need to modify them to compile on other platform
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

#define ASSERT assert // RTree uses ASSERT( condition )
#ifndef Min
#define Min std::min
#endif //Min
#ifndef Max
#define Max std::max
#endif //Max

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
//#define DATA_COUNT 70380191
// #define DATA_COUNT 64000000
//#define DATA_COUNT 153206662
//#define DATA_COUNT 2000000
#define CRACK_SWITCH_SIZE_THRESHOLD 8000


#ifdef stats
    int count_internal_nodes = 0;
    int count_regular_leaves = 0;
    int count_irregular_leaves = 0;
#endif

// #define NUMDIMS 3

#ifdef data_size6M
    #define DATA_COUNT 6000000
#elif data_size3_8M
    #define DATA_COUNT 3800000
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
#elif data_size9_5M
    #define DATA_COUNT 9500000
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
#elif data_size19M
    #define DATA_COUNT 19000000
#elif data_size12M
    #define DATA_COUNT 12000000
#elif data_size24M
    #define DATA_COUNT 24000000
#elif data_size48M
    #define DATA_COUNT 48000000
#elif data_size52_5M
    #define DATA_COUNT 52500000
#elif data_size66_5M
    #define DATA_COUNT 66500000
#elif data_size6_45M
    #define DATA_COUNT 6450000
#elif data_size8_17M
    #define DATA_COUNT 8170000
#elif data_size86_25M
    #define DATA_COUNT 86250000
#elif data_size109_25M
    #define DATA_COUNT 109250000
#elif data_size114_75M
    #define DATA_COUNT 114750000
#endif  

#ifdef dim3
    #define NUMDIMS 3
#else
    #define NUMDIMS 2
#endif

static float m_data_arr_mins[DATA_COUNT][NUMDIMS];
static float m_data_arr_maxes[DATA_COUNT][NUMDIMS];
static int m_data_arr_ids[DATA_COUNT];
//static int m_data_arr_datainfo[DATA_COUNT];

//float **m_data_arr_mins = new float*[DATA_COUNT];
//float **m_data_arr_maxes = new float*[DATA_COUNT];
//int m_data_arr_ids[DATA_COUNT];


//float** m_data_arr_mins = new float*[DATA_COUNT];
//float** m_data_arr_maxes = new float*[DATA_COUNT];
//int* m_data_arr_ids = new int[DATA_COUNT];

// FATEMEH OUT



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

                    /// Remove all entries from tree
                    void RemoveAll();

                    // print tree structure
                    void findObject(int data_index);

                    // print for debug
                    void findObjectRec(Node* a_node, int data_index);

                    /// Count the data elements in this container.  This is slow as no internal counter is maintained.
                    int Count();

                    void findAllLeavesWithDataInThem(int data_index);



                    /// Create 2d-crack defined in adaptive tree creation, temporarily public
                    //    void TwoDCrack(Node* node_parent_of_start, int parent_of_start_branch_index, Node* start, const float a_min[NUMDIMS], const float a_max[NUMDIMS]);

                    // tester
                    bool SearchNNtester(const float a_min[NUMDIMS], const float a_max[NUMDIMS]);

                    // print out the tree structure
                    void printTree(string file_name);

                    // print for debug
                    void printTreeDebug(string file_name);

                    void printLeafSizes(string file_name);

                    float getLeafAreaAverage(); 

                    // testing, should not be public
                    //    bool Insert_anylevel(const Branch &a_branch, Node *start, int a_level);

                    /// Iterator is not remove safe.
                    class Iterator {
                    private:

                        enum {
                            MAX_STACK = 32
                        }; //  Max stack size. Allows almost n^32 where n is number of branches in node

                        struct StackElement {
                            Node *m_node;
                            int m_branchIndex;
                        };

                    public:

                        Iterator() { Init(); }

                        ~Iterator() {}

                        /// Is iterator invalid
                        bool IsNull() { return (m_tos <= 0); }

                        /// Is iterator pointing to valid data
                        bool IsNotNull() { return (m_tos > 0); }

                        /// Access the current data element. Caller must be sure iterator is not NULL first.
                        DATATYPE &operator*() {
                            //                    ASSERT(IsNotNull());
                            StackElement &curTos = m_stack[m_tos - 1];
                            return curTos.m_node->m_branch[curTos.m_branchIndex].m_data;
                        }

                        /// Access the current data element. Caller must be sure iterator is not NULL first.
                        const DATATYPE &operator*() const {
                            //                    ASSERT(IsNotNull());
                            StackElement &curTos = m_stack[m_tos - 1];
                            return curTos.m_node->m_branch[curTos.m_branchIndex].m_data;
                        }

                        /// Find the next data element
                        bool operator++() { return FindNextData(); }

                        /// Get the bounds for this node
                        void GetBounds(float a_min[NUMDIMS], float a_max[NUMDIMS]) {
                            //                    ASSERT(IsNotNull());
                            StackElement &curTos = m_stack[m_tos - 1];
                            Branch &curBranch = curTos.m_node->m_branch[curTos.m_branchIndex];

                            for (int index = 0; index < NUMDIMS; ++index) {
                                a_min[index] = curBranch.m_rect.m_min[index];
                                a_max[index] = curBranch.m_rect.m_max[index];
                            }
                        }

                    private:

                        /// Reset iterator
                        void Init() { m_tos = 0; }

                        /// Find the next data element in the tree (For internal use only)
                        bool FindNextData() {
                            for (;;) {
                                if (m_tos <= 0) {
                                    return false;
                                }
                                StackElement curTos = Pop(); // Copy stack top cause it may change as we use it

                                if (curTos.m_node->IsLeaf()) {
                                    // Keep walking through data while we can
                                    if (curTos.m_branchIndex + 1 < curTos.m_node->m_count) {
                                        // There is more data, just point to the next one
                                        Push(curTos.m_node, curTos.m_branchIndex + 1);
                                        return true;
                                    }
                                    // No more data, so it will fall back to previous level
                                } else {
                                    if (curTos.m_branchIndex + 1 < curTos.m_node->m_count) {
                                        // Push sibling on for future tree walk
                                        // This is the 'fall back' node when we finish with the current level
                                        Push(curTos.m_node, curTos.m_branchIndex + 1);
                                    }
                                    // Since cur node is not a leaf, push first of next level to get deeper into the tree
                                    Node *nextLevelnode = curTos.m_node->m_branch[curTos.m_branchIndex].m_child;
                                    Push(nextLevelnode, 0);

                                    // If we pushed on a new leaf, exit as the data is ready at TOS
                                    if (nextLevelnode->IsLeaf()) {
                                        return true;
                                    }
                                }
                            }
                        }

                        /// Push node and branch onto iteration stack (For internal use only)
                        void Push(Node *a_node, int a_branchIndex) {
                            m_stack[m_tos].m_node = a_node;
                            m_stack[m_tos].m_branchIndex = a_branchIndex;
                            ++m_tos;
                            //                    ASSERT(m_tos <= MAX_STACK);
                        }

                        /// Pop element off iteration stack (For internal use only)
                        StackElement &Pop() {
                            //                    ASSERT(m_tos > 0);
                            --m_tos;
                            return m_stack[m_tos];
                        }

                        StackElement m_stack[MAX_STACK];              ///< Stack as we are doing iteration instead of recursion
                        int m_tos;                                    ///< Top Of Stack index

                        friend class RTree; // Allow hiding of non-public functions while allowing manipulation by logical owner
                    };

                    /// Get 'first' for iteration
                    void GetFirst(Iterator &a_it) {
                        a_it.Init();
                        Node *first = m_root;
                        //        printf("m_root children count = %d\n", first->m_count);
                        //        printf("m_root children count = %d\n", first->m_branch.size());
                        for (typename std::vector<Branch>::iterator it = first->m_branch.begin(); it != first->m_branch.end(); it++) {
                            //            printf("bounds: %d, %d, %d, %d\n", it->m_rect.m_min[0],  it->m_rect.m_min[1], it->m_rect.m_max[0],it->m_rect.m_max[1]);
                        }
                        while (first) {
                            //            printf("in getFirst, while\n");
                            if (first->IsInternalNode() && first->m_count > 1) {
                                //                printf("in getFirst, first was internal\n");
                                a_it.Push(first, 1); // Descend sibling branch later
                            } else if (first->IsLeaf()) {
                                //                printf("in getFirst, first was leaf\n");
                                if (first->m_count) {
                                    //                    printf("in getFirst, first was leaf and had child\n");
                                    a_it.Push(first, 0);
                                }
                                break;
                            }
                            first = first->m_branch[0].m_child;
                        }
                    }

                    /// Get Next for iteration
                    void GetNext(Iterator &a_it) {
                        //        printf("in getNext, tos = %d\n", a_it.m_tos);
                        ++a_it;
                    }

                    /// Is iterator NULL, or at end?
                    bool IsNull(Iterator &a_it) { return a_it.IsNull(); }

                    /// Get object at iterator position
                    DATATYPE &GetAt(Iterator &a_it) { return *a_it; }

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
                    };

                    /// Node for each branch level
                    struct Node {
                        bool IsInternalNode() { return (m_level > 0); } // Not a leaf, but a internal node
                        bool IsLeaf() { return (m_level == 0); } // A leaf, contains data
                        bool isRegular() { return m_regular; } // leaf is regular is m_regular is true, irregular ow

                        // for debug
                        std::string name;
                        int m_count;                                  ///< Count
                        int m_level;                                  ///< Leaf is zero, others positive
                        bool m_regular;                               ///< Some leaves can be irregular. Leaf is irregular if this is false, regular ow
                        //        std::vector<Branch> m_branch;                 ///< Branch
                        Branch m_branch[MAXNODES];                    ///< Branch
//                        Branch m_data_points[MAXDATAPOINTS];            ///< data points, couldnt keep them in m_branch since the limits differ
                        int m_L;                                        ///< left side of data_points interval in m_data_arr
                        int m_R;                                        ///< right side of data_points interval in m_data_arr
                        bool m_transferred;                             ///< whether data has been copied to m_data_points, if true has been copied, if false data is in m_data_arr
                        Node* m_parent;
                        int m_id;

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

                    struct LTA_v2{
                        Node* this_leaf;
                        int Ls[2*NUMDIMS + 1 + 1];
                        int Rs[2*NUMDIMS + 1 + 1];
                        float crack_covers_mins[2*NUMDIMS + 1 + 1][NUMDIMS];
                        float crack_covers_maxes[2*NUMDIMS + 1 + 1][NUMDIMS];
                        int how_many_created = 0;

                    };

                    typedef vector<LTA_v2> LeavesToAdd_v2;


                    // END ELEGANT

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

                    void Add_ltas(LeavesToAdd *all_lta);
                    // to work with lta_v2
                    void Add_ltas_v2(LeavesToAdd_v2 *all_lta);




                    // FATEMEH OUT


                    // FATEMEH
                    bool CountDataPoints(Node *a_node, int &a_foundCount);
                    bool SumDataPointIDs(Node *a_node, int &a_foundCount);
                    bool CountNodes(Node *a_node, int &a_foundCount);
                    void PrintLeafSizesRec(Node *a_node, ofstream &myfile);
                    void getLeafAreaAverage(Node *a_node, float &sum, int &count);
                    // FATEMEH OUT

                    bool SearchNN(Rect query_rect, Node **result);

                    void RemoveAllRec(Node *a_node);

                    void Reset();

                    void CountRec(Node *a_node, int &a_count);


                    void CopyRec(Node *current, Node *other);

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



                    BranchList m_pendings;
                    Rect root_cover;

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
}

RTREE_TEMPLATE
RTREE_QUAL::RTree(std::string data_file_name){
    std::ifstream data_file(data_file_name.c_str());
    if(data_file.is_open()){cout << data_file_name << endl;}
    // int data_count = std::count(std::istreambuf_iterator<char>(data_file),
    //                             std::istreambuf_iterator<char>(), '\n');
    // cout << "FOR FATEMEH " << data_count << endl;
    // if(data_count < DATA_COUNT){
    //     cout << "NOT ENOUGH DATA!!\n";
    //     exit(3);
    // }


    data_file.clear();
    data_file.seekg(0, ios::beg);

    // float min_x, min_y, max_x, max_y;
    float min[NUMDIMS]; float max[NUMDIMS];

    for(int i = 0; i < DATA_COUNT; i++){
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
        m_data_arr_ids[i] = i;
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
    m_root->m_R = DATA_COUNT;

    root_cover = LeafCover(0, DATA_COUNT);

    #ifdef stats
        count_irregular_leaves = 1;
    #endif

    //    cout << "testing " << m_data_arr[0].m_rect.m_min[1] << endl;

}

RTREE_TEMPLATE
RTREE_QUAL::RTree(const RTree &other) : RTree() {
    CopyRec(m_root, other.m_root);
}


RTREE_TEMPLATE
RTREE_QUAL::~RTree() {
    Reset(); // Free, or reset node memory
}


RTREE_TEMPLATE
void RTREE_QUAL::Insert(const float a_min[NUMDIMS], const float a_max[NUMDIMS], const DATATYPE &a_dataId) {


    Branch branch;
    branch.m_data = a_dataId;
    branch.m_child = NULL;

    for (int axis = 0; axis < NUMDIMS; ++axis) {
        branch.m_rect.m_min[axis] = a_min[axis];
        branch.m_rect.m_max[axis] = a_max[axis];
    }

    InsertRect(branch, &m_root, 0);
}


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

    Search_2026(m_root, &rect, foundCount, &all_lta);
    Add_ltas_v2(&all_lta);

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
void RTREE_QUAL::CopyRec(Node *current, Node *other) {
    //    cout << "in copyrec \n";
    current->m_level = other->m_level;
    current->m_count = other->m_count;

    if (current->IsInternalNode())  // not a leaf node
        {
        //        cout << "in copyrec not leaf m_count =  " << current->m_count<<"\n";
        for (int index = 0; index < current->m_count; ++index) {
            //            Branch *currentBranch;
            //            current->m_branch.push_back(*currentBranch);
            Branch *currentBranch = &current->m_branch[index];
            Branch *otherBranch = &other->m_branch[index];
            //            Branch *otherBranch;
            // TODO check the pointers later
            //            other->m_branch.push_back(*otherBranch);
            //            Branch *otherBranch = &other->m_branch[index];

            std::copy(otherBranch->m_rect.m_min,
                      otherBranch->m_rect.m_min + NUMDIMS,
                      currentBranch->m_rect.m_min);

            std::copy(otherBranch->m_rect.m_max,
                      otherBranch->m_rect.m_max + NUMDIMS,
                      currentBranch->m_rect.m_max);

            currentBranch->m_child = AllocNode();
            CopyRec(currentBranch->m_child, otherBranch->m_child);
        }
        } else // A leaf node
        {
        for (int index = 0; index < current->m_count; ++index) {


            //            Branch *currentBranch;
            //            current->m_branch.push_back(*currentBranch);
            Branch *currentBranch = &current->m_branch[index];
            Branch *otherBranch = &other->m_branch[index];
            //            Branch *otherBranch;
            // TODO check the pointers later
            //            other->m_branch.push_back(*otherBranch);
            //            Branch *otherBranch = &other->m_branch[index];


            //            for (int j = 0; j < NUMDIMS; j++) {
            //                currentBranch->m_rect.m_min[j] = otherBranch->m_rect.m_min[j];
            //            }
            std::copy(otherBranch->m_rect.m_min,
                      otherBranch->m_rect.m_min + NUMDIMS,
                      currentBranch->m_rect.m_min);


            std::copy(otherBranch->m_rect.m_max,
                      otherBranch->m_rect.m_max + NUMDIMS,
                      currentBranch->m_rect.m_max);


            currentBranch->m_data = otherBranch->m_data;
            //            current->m_branch.push_back(*currentBranch);
        }
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
        Rect for_me = NodeCover(a_node);
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
            a_node->m_branch[index].m_rect = Rect(a_coverRect.m_min, a_coverRect.m_max);
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
    if( !a_node->IsLeaf()) {
        //        cout << "NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO\n";
        rect = a_node->m_branch[0].m_rect;
        for (int index = 1; index < a_node->m_count; ++index) {
            rect = CombineRect(&rect, &(a_node->m_branch[index].m_rect));
        }
    } else{
//        rect = m_data_arr[a_node->m_L]->m_rect;
//         rect = Rect(m_data_arr_mins[a_node->m_L], m_data_arr_maxes[a_node->m_L]);
//         for (int index = a_node->m_L + 1; index < a_node->m_R; ++index) {
// //            rect = CombineRect(&rect, &(m_data_arr[index]->m_rect));
//             alaki_rect = Rect(m_data_arr_mins[index], m_data_arr_maxes[index]);
//             rect = CombineRect(&rect, &(alaki_rect));
//         }
        cout << "NODE COVER CALLED ON LEAF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        exit(4);
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
//         cout << "REMEMBER YOU HAVE NOT FIXED THIS!! ADDBRACNH FOR IRREG LEAF\n";
//         cout << "NEVER MIND SHOULD NEVER HAVE HAPPENED!\n";
//         return false;
//     }
// }





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


RTREE_TEMPLATE
bool RTREE_QUAL::AddBranch(const Branch *a_branch, Node *a_node, Node **a_newNode) {
    //            ASSERT(a_branch);
    //            ASSERT(a_node);
    #ifdef stats
        if(a_node->m_level > 1) count_internal_nodes++;
        else{
            if(a_branch->m_child->isRegular()) count_regular_leaves++;
            else count_irregular_leaves++;
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


// 2023
RTREE_TEMPLATE
bool RTREE_QUAL::AddBranch(const Branch *a_branch, Node *a_node, Node **a_newNode, Rect &a_coverRect, Rect &a_newCoverRect) {
    //            ASSERT(a_branch);
    //            ASSERT(a_node);
    #ifdef stats
        if(a_node->m_level > 1) count_internal_nodes++;
        else{
            if(a_branch->m_child->isRegular()) count_regular_leaves++;
            else count_irregular_leaves++;
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
        else{
            if(a_node->m_branch[a_index].m_child->isRegular()) count_regular_leaves--;
            else count_irregular_leaves--;
        }
    #endif

    // Remove element by swapping with the last element to prevent gaps in array
    a_node->m_branch[a_index] = a_node->m_branch[a_node->m_count - 1];
    //    if(a_node->m_branch[a_index].m_child != NULL){
    //        std::vector<Branch>().swap((a_node->m_branch[a_index].m_child)->m_branch);
    //    }
    //    a_node->m_branch[a_index] = a_node->m_branch.back();
    //    a_node->m_branch.pop_back();
    a_node->m_count -= 1;
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

    for (int index = 0; index < a_node->m_count; ++index) {
        Rect *curRect = &a_node->m_branch[index].m_rect;
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
                if(a_parVars->m_branchBuf[index].m_child->isRegular()) count_regular_leaves--;
                else count_irregular_leaves--;
            }
        #endif
    }
    #ifdef stats
        if(a_nodeA->m_level > 1){
            count_internal_nodes -= a_parVars->m_total;
        }
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
                // if(m_data_arr_ids[index] == 45152){
                //     cout << "in regular, checking overlap for 45152, result:  " << Overlap(a_rect, m_data_arr_mins[index], m_data_arr_maxes[index]) << endl;
                // }
                if (Overlap(a_rect, m_data_arr_mins[index], m_data_arr_maxes[index])) {
                    // cout << m_data_arr_ids[index] << "\n" ;
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

                // // 2023 DEBUGGING
                // cout << "before crack call, other piece rect: " << endl;
                // for(int ii = 0; ii < NUMDIMS; ii++){
                //     cout << "( " << other_piece.m_min[ii] << ", " << other_piece.m_max[ii] << ") ";
                // }
                // cout << endl;
                // cout << "before crack call, query piece rect: " << endl;
                // for(int ii = 0; ii < NUMDIMS; ii++){
                //     cout << "( " << potential_query_piece.m_min[ii] << ", " << potential_query_piece.m_max[ii] << ") ";
                // }
                // cout << endl;

                // // 2023 DEBUGGING
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

                // // 2023 DEBUGGING
                // cout << "after crack call, other piece rect: " << endl;
                // for(int ii = 0; ii < NUMDIMS; ii++){
                //     cout << "( " << other_piece.m_min[ii] << ", " << other_piece.m_max[ii] << ") ";
                // }
                // cout << endl;
                // cout << "after crack call, query piece rect: " << endl;
                // for(int ii = 0; ii < NUMDIMS; ii++){
                //     cout << "( " << potential_query_piece.m_min[ii] << ", " << potential_query_piece.m_max[ii] << ") ";
                // }
                // cout << endl;

                // // 2023 DEBUGGING

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
                    // 2023
                    // potential_query_piece = Rect(this_piece_rect.m_min, this_piece_rect.m_max);
                    // 2023
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
                    // cout << m_data_arr_ids[folan] << endl;
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

                    // this_lta.crack_covers_mins[this_lta.how_many_created][ashghal] = this_piece_rect.m_min[ashghal];
                    // this_lta.crack_covers_maxes[this_lta.how_many_created][ashghal] = this_piece_rect.m_max[ashghal];
                    this_lta.crack_covers_mins[this_lta.how_many_created][ashghal] = potential_query_piece.m_min[ashghal];
                    this_lta.crack_covers_maxes[this_lta.how_many_created][ashghal] = potential_query_piece.m_max[ashghal];
                }
                this_lta.how_many_created++;
            }
            // aval bayad bozorgtarin tikke ro peyda konim
            ///// START OF STOCH PART
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
        ////// END OF STOCH PART

            all_lta->push_back(this_lta);
        }
    }
    return true; // Continue searching
}




RTREE_TEMPLATE
void RTREE_QUAL::Add_ltas_v2(LeavesToAdd_v2 *all_lta){
    // LOGGER 1 sep 2021
//    int count_pieces_created = 0;
//    for(int i = 0; i < all_lta->size(); i++){
//        count_pieces_created += ((all_lta->at(i)).how_many_created);
//    }
//    cout << "pieces_created " << count_pieces_created << endl;
    // END LOGGER
    // for(int i =0; i < all_lta->size(); i++){
    //     cout << "all_lta[" << i << "] for node with id: " << all_lta->at(i).this_leaf->m_id << "cracked into " << all_lta->at(i).how_many_created <<endl;
    //     for(int j = 0; j < all_lta->at(i).how_many_created; j++){
    //         cout << "from " << all_lta->at(i).Ls[j] << " til " << all_lta->at(i).Rs[j] << endl;
    //     }
    // }
    for(int i = 0; i < all_lta->size(); i++){
        if((all_lta->at(i).this_leaf)->m_parent != NULL){
            Node* parent_of_start = (all_lta->at(i).this_leaf)->m_parent;
            int branch_index = getNodeBranchIndex((all_lta->at(i).this_leaf));
            DisconnectBranch(parent_of_start, branch_index);
            for(int j = 0 ; j < ((all_lta->at(i)).how_many_created); j++){
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

                // 2023 DEBUGGING
                // Rect checker = LeafCover(this_branch.m_child->m_L, this_branch.m_child->m_R);
                // for(int ii = 0; ii < NUMDIMS; ii++){
                //     if(this_branch.m_rect.m_min[ii] > checker.m_min[ii] || this_branch.m_rect.m_max[ii] < checker.m_max[ii]){
                //         cout << "COVER WAS WRONG THIS IS WHRE THE SHIT GOES DOWN..." << endl;
                //         cout << this_branch.m_rect.m_min[ii] << " " << checker.m_min[ii] << " " << this_branch.m_rect.m_max[ii] << " " << checker.m_max[ii] << endl;
                //         exit(5);
                //     }
                // }
//                for(int ashghal = 0; ashghal < NUMDIMS; ashghal++){
//                    if(this_branch.m_rect.m_min[ashghal] > this_branch.m_rect.m_max[ashghal]){
//                        cout << "TOOYE ADDLTA, MIN>MAX " << ashghal << " min " << this_branch.m_rect.m_min[ashghal] << " max " << this_branch.m_rect.m_max[ashghal] << endl;
//                    }
//                }



                Insert_anylevel(this_branch, parent_of_start, 1);
            }
        }
        else{
            m_root->m_count = 0;
            m_root->m_level++;
            m_root->m_regular = true;
            m_root->m_parent = NULL;

            #ifdef stats
                count_irregular_leaves--;
                count_internal_nodes++;
            #endif

            for(int j = 0 ; j < ((all_lta->at(i)).how_many_created); j++){
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
                // 2023 DEBUGGING
                // this_branch.m_rect = LeafCover(this_branch.m_child->m_L, this_branch.m_child->m_R);
                // Rect checker = LeafCover(this_branch.m_child->m_L, this_branch.m_child->m_R);
                // for(int ii = 0; ii < NUMDIMS; ii++){
                //     if(this_branch.m_rect.m_min[ii] > checker.m_min[ii] || this_branch.m_rect.m_max[ii] < checker.m_max[ii]){
                //         cout << "COVER WAS WRONG THIS IS WHRE THE SHIT GOES DOWN..." << endl;
                //         cout << this_branch.m_rect.m_min[ii] << " " << checker.m_min[ii] << " " << this_branch.m_rect.m_max[ii] << " " << checker.m_max[ii] << endl;
                //         exit(5);
                //     }
                // }
//                for(int ashghal = 0; ashghal < NUMDIMS; ashghal++){
//                    if(this_branch.m_rect.m_min[ashghal] > this_branch.m_rect.m_max[ashghal]){
//                        cout << "TOOYE ADDLTA, MIN>MAX " << ashghal << " min " << this_branch.m_rect.m_min[ashghal] << " max " << this_branch.m_rect.m_max[ashghal] << endl;
//                    }
//                }

                InsertRect(this_branch, &m_root, 1);
            }
        }
    }
}


RTREE_TEMPLATE
bool RTREE_QUAL::CountDataPoints(Node *a_node, int &a_foundCount){
    if (a_node->m_parent == NULL) a_foundCount += m_pendings.branch_list.size();
    if (a_node->IsInternalNode()) {
        for (int index = 0; index < a_node->m_count; ++index) {
            CountDataPoints(a_node->m_branch[index].m_child, a_foundCount);
        }
    } else {
        a_foundCount += a_node->m_count;
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
    //    cout << "in get index"<< (a_node->m_parent)->m_branch.size() << "\n";
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
    //    cout << "in insert_any \n";
    if(start->m_level > a_level){
        printf("calling pb in insert_anylevel\n");
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
            // last current_node was not root
            if(!just_go_till_top) {
                //                cout << "in while not just till top\n";
                Node *newNode = AllocNode();
                Rect rect1; Rect rect2;
                is_split = AddBranch(&current_branch, current_node, &newNode, rect1, rect2);

                if (is_split) {
                    //                    cout << "in while was split\n";
                    // was split
                    if (current_node->m_parent != NULL) {
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
void RTREE_QUAL::findObject(int data_index) {
    findObjectRec(m_root, data_index);
}

RTREE_TEMPLATE
void RTREE_QUAL::findObjectRec(Node* a_node, int data_index){
     if (a_node->IsInternalNode()) {
        for (int index = 0; index < a_node->m_count; ++index) {
            // myfile << a_node->m_branch[index].m_child->m_level << " " << a_node->m_branch[index].m_rect.m_min[0] << " " << a_node->m_branch[index].m_rect.m_min[1] << " " << a_node->m_branch[index].m_rect.m_max[0] << " " << a_node->m_branch[index].m_rect.m_max[1] << "\n";
            findObjectRec(a_node->m_branch[index].m_child, data_index);
        }
    }
    else
    {
        for(int i = a_node->m_L; i < a_node->m_R; i++){
            if(m_data_arr_ids[i] == data_index){
                cout << "data in leaf " << a_node->m_id << " is item " << i << " and L, R: " << a_node->m_L << " " << a_node->m_R << endl;
                cout << "data bounds: " << m_data_arr_mins[i][0] << " " << m_data_arr_maxes[i][0] << ", " << m_data_arr_mins[i][1] << " " << m_data_arr_maxes[i][1] << endl;
                 if(a_node->m_parent != NULL){
                    int branch_index = getNodeBranchIndex(a_node);
                    Rect cover = (a_node->m_parent)->m_branch[branch_index].m_rect;
                    cout << "leaf bound: " << cover.m_min[0] << " " << cover.m_max[0] << " " << cover.m_min[1]  << " " << cover.m_max[1] << endl;
                }
            }
        }
    }
  return;

}


RTREE_TEMPLATE
void RTREE_QUAL::getLeafAreaAverage(Node* a_node, float &sum, int &count)
{
  if(a_node->m_level > 1)
  {
    for(int i = 0; i < a_node->m_count;i++)
    {
        // file << a_node->m_level - 1 << " " << a_node->m_branch[i].m_rect.m_min[0] << " " << a_node->m_branch[i].m_rect.m_max[0] << " " << a_node->m_branch[i].m_rect.m_min[1] << " " << a_node->m_branch[i].m_rect.m_max[1] << endl;
        getLeafAreaAverage(a_node->m_branch[i].m_child, sum, count);
    }
  }
  else if(a_node->m_level == 1)
  {
    for(int i = 0; i < a_node->m_count;i++)
    {
        double area = (a_node->m_branch[i].m_rect.m_max[0] - a_node->m_branch[i].m_rect.m_min[0])*(a_node->m_branch[i].m_rect.m_max[1] - a_node->m_branch[i].m_rect.m_min[1]);
        // file << a_node->m_level - 1 << " " << a_node->m_branch[i].m_rect.m_min[0] << " " << a_node->m_branch[i].m_rect.m_max[0] << " " << a_node->m_branch[i].m_rect.m_min[1] << " " << a_node->m_branch[i].m_rect.m_max[1] << endl;
        sum += area;
        count++;
        // file << area << endl;
    }
  }
  return;
}

RTREE_TEMPLATE
float RTREE_QUAL::getLeafAreaAverage()
{
  float sum = 0.0; int count =0;
  getLeafAreaAverage(m_root, sum, count);
  return sum / count;
}



#undef RTREE_TEMPLATE
#undef RTREE_QUAL

#endif //RTREE_723OCT21_H

