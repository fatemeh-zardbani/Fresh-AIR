#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <stdio.h>
#include <sys/time.h>
// #include "RTree_723oct21s_L2.h"
// #include "RTree_723oct21s_stochastic_15June23.h"
// #include "RTree_723oct21s_stochastic_14dec22.h"
// #include "RTree_723oct21s_stochastic_24July23.h"
// #include "RTree_723oct21s_stochastic_25July23.h"
// #include "RTree_723oct21s_stochastic_2Aug23.h"

#ifdef version1
   #include "RTree_723oct21s_stochastic_24July23.h"
#elif version2
    #include "RTree_723oct21s_stochastic_25July23.h"
#elif version3
    #include "RTree_723oct21s_stochastic_2Aug23.h"
#elif version32
    #include "RTree_723oct21s_stochastic_2Aug23_v2.h"
#elif version4
    #include "RTree_723oct21s_stochastic_25July23_v2.h"
#elif version5
    #include "RTree_723oct21s_stochastic_25July23_v3.h"
#elif version6
    #include "RTree_723oct21s_stochastic_14Aug23.h"
#elif version7
    #include "RTree_723oct21s_stochastic_15Aug23.h"
#elif version72
    #include "RTree_723oct21s_stochastic_15Aug23_v2.h"
#elif version8
    #include "RTree_723oct21s_stochastic_7Sep23.h"
#elif version9
    #include "RTree_723oct21s_stochastic_122Sep23.h"
#elif version10
    #include "RTree_723oct21s_stochastic_12Sep23.h"
#elif version11
    #include "RTree_723oct21s_stochastic_13Sep23.h"
#elif version12
    #include "RTree_723oct21s_stochastic_132Sep23.h"
#elif version122
    #include "RTree_723oct21s_stochastic_14Sep23.h"
#elif version122R
    #include "RTree_723oct21s_stochastic_14Sep23R.h"
#elif versionCS
    #include "RTree_723oct21s_stochastic_21Sep23_CS.h"
#elif versionCSR
    #include "RTree_723oct21s_stochastic_21Sep23_CSR.h"
#elif versionCSv11
    #include "RTree_723oct21s_stochastic_21Sep23_CSv11.h"
#elif versionCS2
    #include "RTree_723oct21s_stochastic_21Sep23_CS2.h"
#elif versionCS2R
    #include "RTree_723oct21s_stochastic_21Sep23_CS2_ripple.h"
#elif versionCS2RC
    #include "RTree_723oct21s_stochastic_21Sep23_CS2RC.h"
#elif versionCS2RP
    #include "RTree_723oct21s_stochastic_21Sep23_CS2RP.h"
#elif versionCS2RM
    #include "RTree_723oct21s_stochastic_21Sep23_CS2RM.h"
#elif versionGQ
    #include "RTree_723oct21s_stochastic_21Sep23_GQ.h"
#elif versionGQ1
    #include "RTree_723oct21s_stochastic_21Sep23_GQ1.h"
#elif versionGQ2
    #include "RTree_723oct21s_stochastic_21Sep23_GQ2.h"
#elif versionGQ3
    #include "RTree_723oct21s_stochastic_21Sep23_GQ3.h"
#elif versionGQ3M
    #include "RTree_723oct21s_stochastic_21Sep23_GQ3M.h"
#elif versionGQ3MR
    #include "RTree_723oct21s_stochastic_21Sep23_GQ3MR.h"
#elif versionGS
    #include "RTree_723oct21s_stochastic_21Sep23_GS.h"
#elif versionGS2
    #include "RTree_723oct21s_stochastic_21Sep23_GS2.h"
#elif versionGS3
    #include "RTree_723oct21s_stochastic_21Sep23_GS3.h"
#elif versionGS3M
    #include "RTree_723oct21s_stochastic_21Sep23_GS3M.h"
#elif versionGS4
    #include "RTree_723oct21s_stochastic_21Sep23_GS4.h"
#elif versionGS4M
    #include "RTree_723oct21s_stochastic_21Sep23_GS4M.h"
#elif versionGS4Mv11
    #include "RTree_723oct21s_stochastic_21Sep23_GS4Mv11.h"
#elif versionGS4Mv11R
    #include "RTree_723oct21s_stochastic_21Sep23_GS4Mv11R.h"
#elif versionGS4Mv11RP
    #include "RTree_723oct21s_stochastic_21Sep23_GS4Mv11RP.h"
#elif versionGS4Mv11RS
    #include "RTree_723oct21s_stochastic_21Sep23_GS4Mv11RS.h"
#elif versionGS4Mv11RSR
    #include "RTree_723oct21s_stochastic_21Sep23_GS4Mv11RSR.h"
#elif versionGS4Mv11RSRP
    #include "RTree_723oct21s_stochastic_21Sep23_GS4Mv11RSRP.h"
#elif versionGS4Mv11RSRC
    #include "RTree_723oct21s_stochastic_21Sep23_GS4Mv11RSRC.h"
#elif versionGS3M
    #include "RTree_723oct21s_stochastic_21Sep23_GS3M.h"
#elif versionCQGS
    #include "RTree_723oct21s_stochastic_21Sep23_CQ_GS.h"
#elif versionCQGS2
    #include "RTree_723oct21s_stochastic_21Sep23_CQ_GS2.h"
#elif versionCQGS3
    #include "RTree_723oct21s_stochastic_21Sep23_CQ_GS3.h"
#elif versionGSTerminal
    #include "RTree_723oct21s_stochastic_21Sep23_GS_terminal.h"
#elif versionGSL
    #include "RTree_723oct21s_stochastic_21Sep23_GSL.h"
#elif versionGSTSpecial
    #include "RTree_723oct21s_stochastic_21Sep23_GST_Special.h"
#endif



#ifdef fo8
    const int fo = 8;
#elif fo16
    const int fo = 16;
#elif fo32
    const int fo = 32;
#endif


#ifdef maxt64
    const int maxt = 64;
#elif maxt256
    const int maxt = 256;
#elif maxt512
    const int maxt = 512;
#elif maxt1024
    const int maxt = 1024;
#elif maxt2048
    const int maxt = 2048;
#endif


#ifdef mint32
    const int mint = 32;
#elif mint48
    const int mint = 48;
#elif mint128
    const int mint = 128;
#elif mint192
    const int mint = 192;
#elif mint256
    const int mint = 256;
#elif mint384
    const int mint = 384;
#elif mint512
    const int mint = 512;
#elif mint768
    const int mint = 768;
#elif mint1024
    const int mint = 1024;
#elif mint1536
    const int mint = 1536;
#endif







using namespace std;


typedef int ValueType;
typedef float ELEMTYPE;
typedef double ELEMTYPEREAL;

struct Rect
{
    Rect()  {}

    //Rect(ELEMTYPE a_minX, ELEMTYPE a_minY, ELEMTYPE a_maxX, ELEMTYPE a_maxY)
    //{
    //    m_min[0] = a_minX;
    //    m_min[1] = a_minY;
    //    m_max[0] = a_maxX;
    //    m_max[1] = a_maxY;
    //}

    Rect(ELEMTYPE a_min[NUMDIMS], ELEMTYPE a_max[NUMDIMS]){
        for(int i = 0; i < NUMDIMS; i++) m_min[i] = a_min[i];
        for(int i = 0; i < NUMDIMS; i++) m_max[i] = a_max[i];
    }


    ELEMTYPE m_min[NUMDIMS];
    ELEMTYPE m_max[NUMDIMS];
};


bool Overlap(Rect *a_rectA, Rect *a_rectB);

bool Overlap(Rect *a_rectA, Rect *a_rectB) {
            ASSERT(a_rectA && a_rectB);

    for (int index = 0; index < 2; ++index) {
        if (a_rectA->m_min[index] > a_rectB->m_max[index] ||
            a_rectB->m_min[index] > a_rectA->m_max[index]) {
            return false;
        }
    }
    return true;
}



double timing(){
    static struct timeval t1, t2;
    gettimeofday(&t2,NULL);
    double ret = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) * 1e-6;
    t1 = t2;
    return ret;
}


int main(int argc, char **argv){


     auto seed = time(NULL); // 1695893019
    // auto seed = 0;
//    auto seed = 1697617670;
//    cout << seed << endl;
    // auto seed = 1692015686; // crashes at 81199
    srand(seed);

    string data_file_name = argv[1];
    int query_size = stoi(argv[2]);
    string query_file_name = argv[3];
    string time_file_name = argv[4];
    int ratio = stoi(argv[5]);
    int insert_count = stoi(argv[6]);
    
    std::ifstream query_file(query_file_name.c_str());
    //if(query_file.is_open()){cout << "it is open:)\n";}

    #ifdef stats
        ofstream stats_file;
        string stats_file_name = argv[7];
        stats_file.open(stats_file_name.c_str());   
    #endif

    ofstream times_file;
    times_file.open(time_file_name.c_str());

    
    typedef RTree<ValueType, fo, fo/2, maxt, mint> MyTree;
    MyTree tree(data_file_name);

    ELEMTYPE min[NUMDIMS]; ELEMTYPE max[NUMDIMS];
    int this_id;

  //  cout << "INSERTION COMPLETE\n";

    query_file.clear();
    query_file.seekg(0, ios::beg);
   // cout << "cleared q file\n";
    //cout << "query size " << query_size << endl;
    Rect queries[query_size];
   // cout << "created rect array\n";

    vector<Rect> data(DATA_COUNT);
    std::ifstream data_file(data_file_name.c_str());

    for(int i = 0; i < DATA_COUNT; i++){
        data[i] = Rect();
        for(int j = 0; j < NUMDIMS; j++) data_file >> data[i].m_min[j];
        for(int j = 0; j < NUMDIMS; j++) data_file >> data[i].m_max[j];
    }
    

    for(int i = 0; i < query_size; i++){
        queries[i] = Rect();
        for(int j = 0; j < NUMDIMS; j++) query_file >> queries[i].m_min[j];
        for(int j = 0; j < NUMDIMS; j++) query_file >> queries[i].m_max[j];
       
    }
    //cout << "LOADED QUERIES\n";
    int found_count; int sum_ids; int total_count;
    int lin_count;
    clock_t q_time = 0;
    int rightestright = 0;
    int last_id = DATA_COUNT -1;
    // timing();
    clock_t insert_time, start_insert_time;
    clock_t start_query_time;
    // cout << "DATA COUNT: " << DATA_COUNT << " data.size() " << data.size() << endl;
    // exit(4);
    for(int i = 0; i < query_size; i++){
    // for(int i = 0; i < 2; i++){
        cerr << "QUERY " << i << "...\n";
        insert_time = 0.0;
        if(i % ratio == 0 && i != 0){
            for(int o =0; o < insert_count; o++){
                // insert one
                for(int j = 0; j < NUMDIMS; j++) data_file >> min[j];
                for(int j = 0; j < NUMDIMS; j++) data_file >> max[j];
                this_id = ++last_id;
                // timing();
                start_insert_time = clock();
                // #ifdef versionGS 
                //     tree.Insert_gradual(min, max, this_id);
                // #elif versionGS2
                //     tree.Insert_gradual(min, max, this_id);
                // #elif versionGS3
                //     tree.Insert_gradual(min, max, this_id);
                // #elif versionGS3M
                //     tree.Insert_gradual(min, max, this_id);
                // #elif versionGSTerminal
                //     tree.Insert_gradual(min, max, this_id);
                // #elif versionGSTSpecial
                //     tree.Insert_gradual(min, max, this_id);
                // #elif versionCQGS
                //     tree.Insert_gradual(min, max, this_id);
                // #elif versionCQGS2
                //     tree.Insert_gradual(min, max, this_id);
                // #elif versionCQGS3
                //     tree.Insert_gradual(min, max, this_id);
               #if defined(versionGS) || defined(versionGS2) || defined(versionGS3) || defined(versionGS4) || defined(versionGS4M) || defined(versionGS4Mv11) || defined(versionGS4Mv11R) || defined(versionGS4Mv11RP) || defined(versionGS4Mv11RS) || defined(versionGS4Mv11RSR) || defined(versionGS4Mv11RSRC) || defined(versionGS4Mv11RSRP) || defined(versionGS3M)  || defined(versionGSTerminal) || defined(versionGSL) || defined(versionGSTSpecial) || defined(CQGS) || defined(CQGS2) || defined(CQGS3)
                    tree.Insert_gradual(min, max, this_id);
                    // cout << "after insert gradual " << this_id << endl;
                    // #ifdef stats
                    //     int rightestright = tree.getRightestRight();
                    //     // int holes = tree.getCountOfHoles();
                    //     int height = tree.TreeHeight();
                    //     // int pending_count = tree.getPendingCount();
                    //     // auto leaf_count = tree.CountLeaves();
                    //     // int node_count = tree.CountNodes();
                    //     // stats_file << node_count << " " << leaf_count.first << " " << leaf_count.second << " " << pending_count << " " << height << " " << rightestright << endl;
                    //     // stats_file << rightestright << " " << holes << " " << push_back_count << " " << pushed_back_items_count << " " << fill_holes_count << " " << split_leaf_count << endl;
                    //     // stats_file << "i " << rightestright << " " << holes << " " << count_total_holes  << endl;
                    //     stats_file << "i " << rightestright << " " << count_total_holes << " " << height << " " << count_internal_nodes << " " << count_regular_leaves << " " << count_irregular_leaves << endl;
                    //     // if(count_irregular_leaves != leaf_count.second || count_regular_leaves != leaf_count.first)
                    //     // {
                    //     //     cout << "leaf count is wrong! " << count_regular_leaves << " " << leaf_count.first << " " << count_irregular_leaves << " " << leaf_count.second  << endl;
                    //     //     exit(3);
                    //     // }
                    //     // if(count_internal_nodes != node_count - (leaf_count.first + leaf_count.second ))
                    //     // {
                    //     //     cout << "node count is wrong " << count_internal_nodes << " " <<  node_count -(leaf_count.first + leaf_count.second) << endl;
                    //     //     exit(3);
                    //     // }
                    //     // if(holes != count_total_holes) {
                    //     //     stats_file << "FREAK OUT!!!!!!!!!!!!!" << endl;
                    //     //     exit(3);
                    //     // }
                    //     push_back_count = 0;
                    //     fill_holes_count = 0;
                    //     split_leaf_count = 0;
                    //     pushed_back_items_count = 0;
                        
                    // #endif
                #elif defined(versionCS) || defined(versionCSR) || defined(versionCS2) || defined(versionCSv11) || defined(versionCS2R) || defined(versionCS2RC) || defined(versionCS2RM) || defined(versionCS2RP)
                    // cout << "inserting item " << this_id << endl;
                    tree.Insert_full(min, max, this_id);
                    // cout << "afer insert" << endl;
                    // rightestright = tree.getRightestRight();
                    // int holes = tree.getCountOfHoles();
                    // cout << "11 rr  = " << rightestright << " holes = " << holes << endl; 
                    // #ifdef stats
                    //     int rightestright = tree.getRightestRight();
                    //     // int holes = tree.getCountOfHoles();
                    //     int height = tree.TreeHeight();
                    //     // int pending_count = tree.getPendingCount();
                    //     // auto leaf_count = tree.CountLeaves();
                    //     // int node_count = tree.CountNodes();
                    //     // if(count_irregular_leaves != leaf_count.second || count_regular_leaves != leaf_count.first)
                    //     // {
                    //     //     cout << "leaf count is wrong! " << count_regular_leaves << " " << leaf_count.first << " " << count_irregular_leaves << " " << leaf_count.second  << endl;
                    //     //     exit(3);
                    //     // }
                    //     // if(count_internal_nodes != node_count - (leaf_count.first + leaf_count.second ))
                    //     // {
                    //     //     cout << "node count is wrong " << count_internal_nodes << " " <<  node_count -(leaf_count.first + leaf_count.second) << endl;
                    //     //     exit(3);
                    //     // }
                    //     // stats_file << node_count << " " << leaf_count.first << " " << leaf_count.second << " " << pending_count << " " << height << " " << rightestright << endl;
                    //     // stats_file << rightestright << " " << holes << " " << push_back_count << " " << pushed_back_items_count << " " << fill_holes_count << " " << split_leaf_count << endl;
                    //     // stats_file << "i " << rightestright << " " << holes << " " << count_total_holes  << endl;
                    //     stats_file << "i " << rightestright << " " << count_total_holes << " " << height << " " << count_internal_nodes << " " << count_regular_leaves << " " << count_irregular_leaves << endl;                     
                    //     // if(holes != count_total_holes) {
                    //         // stats_file << "FREAK OUT!!!!!!!!!!!!!" << endl;
                    //         // exit(3);
                    //     // }
                    //     push_back_count = 0;
                    //     fill_holes_count = 0;
                    //     split_leaf_count = 0;
                    //     pushed_back_items_count = 0;
                        
                    // #endif
                #else
                    tree.Insert(min, max, this_id);
                #endif
                // insert_time += timing();
                insert_time += clock() - start_insert_time;
                // cout << "inserted " << this_id << " " << min[0] << " " << max[0] << " , " << min[1] << " " << max[1] << endl;
                data.push_back(Rect(min, max));
            }
        }

        // cout << "after all insertions before query" << endl;


        // cout << "before insertion and query." << endl;
        // if(i == 82960 || i == 82961 || i == 82962) tree.printLeafLinkedList();

        // for(int j = 0; j < NUMDIMS; j++) data_file >> min[j];
        // for(int j = 0; j < NUMDIMS; j++) data_file >> max[j];
        // this_id = ++last_id;
        // tree.Insert(min, max, this_id);
        // // cout << "inserted " << this_id << endl;
        // data.push_back(Rect(min, max));

        // timing();
        start_query_time = clock();
        // #ifdef versionGS 
        //     found_count = tree.QueryAdaptive_v15(queries[i].m_min, queries[i].m_max);
        // #elif versionGS2
        //     found_count = tree.QueryAdaptive_v15(queries[i].m_min, queries[i].m_max);
        // #elif versionGS3
        //     found_count = tree.QueryAdaptive_v15(queries[i].m_min, queries[i].m_max);
        // #elif versionGS3M
        //     found_count = tree.QueryAdaptive_v15(queries[i].m_min, queries[i].m_max);
        // #elif versionGSTerminal
        //     found_count = tree.QueryAdaptive_v15(queries[i].m_min, queries[i].m_max);
        // #elif versionGSTSpecial
        //     // cout << "calling QA15..." << endl;
        //     found_count = tree.QueryAdaptive_v15(queries[i].m_min, queries[i].m_max);
        // #elif versionCS
        //     found_count = tree.QueryAdaptive_v15(queries[i].m_min, queries[i].m_max);
        // #elif versionCS2
        //     found_count = tree.QueryAdaptive_v15(queries[i].m_min, queries[i].m_max);
         #if defined(versionGS) || defined(versionGS2) || defined(versionGS3) || defined(versionGS4) || defined(versionGS4M) || defined(versionGS4Mv11) || defined(versionGS4Mv11R)|| defined(versionGS4Mv11RP) || defined(versionGS4Mv11RS) || defined(versionGS4Mv11RSRC) || defined(versionGS4Mv11RSR) || defined(versionGS4Mv11RSRP) || defined(versionGS3M) || defined(versionGSTerminal) || defined(versionGSL) || defined(versionGSTSpecial) || defined(versionCS) || defined(versionCSR) || defined(versionCS2) || defined(versionCSv11) || defined(versionCS2R) || defined(versionCS2RC) || defined(versionCS2RM) || defined(versionCS2RP)
            found_count = tree.QueryAdaptive_v15(queries[i].m_min, queries[i].m_max);
            // cout << "after qa " << i << endl; 
            #ifdef stats
                int rightestright = tree.getRightestRight();
                // int holes = tree.getCountOfHoles();
                int height = tree.TreeHeight();
                // int pending_count = tree.getPendingCount();
                // auto leaf_count = tree.CountLeaves();
                // int node_count = tree.CountNodes();
                // stats_file << node_count << " " << leaf_count.first << " " << leaf_count.second << " " << pending_count << " " << height << " " << rightestright << endl;
                // stats_file << rightestright << " " << holes << " " << push_back_count << " " << pushed_back_items_count << " " << fill_holes_count << " " << split_leaf_count << endl;
                // stats_file << "q " << rightestright << " " << holes  << endl;
                // stats_file << "q " << rightestright << " " << count_total_holes  << endl;
                stats_file << "i " << rightestright << " " << count_total_holes << " " << height << " " << count_internal_nodes << " " << count_regular_leaves << " " << count_irregular_leaves << endl;
                // if(count_irregular_leaves != leaf_count.second || count_regular_leaves != leaf_count.first)
                // {
                //     cout << "leaf count is wrong! " << count_regular_leaves << " " << leaf_count.first << " " << count_irregular_leaves << " " << leaf_count.second  << endl;
                //     exit(3);
                // }
                // if(count_internal_nodes != node_count - (leaf_count.first + leaf_count.second ))
                // {
                //     cout << "node count is wrong " << count_internal_nodes << " " <<  node_count -(leaf_count.first + leaf_count.second) << endl;
                //     exit(3);
                // }
//                push_back_count = 0;
 //               fill_holes_count = 0;
  //              split_leaf_count = 0;
    //            pushed_back_items_count = 0;
                
            #endif
        #else
            found_count = tree.QueryAdaptive_v14(queries[i].m_min, queries[i].m_max);
            #ifdef stats
                int rightestright = tree.getRightestRight();
                // int holes = tree.getCountOfHoles();
                int height = tree.TreeHeight();
                // int pending_count = tree.getPendingCount();
                // auto leaf_count = tree.CountLeaves();
                // int node_count = tree.CountNodes();
                // stats_file << node_count << " " << leaf_count.first << " " << leaf_count.second << " " << pending_count << " " << height << " " << rightestright << endl;
                // stats_file << rightestright << " " << holes << " " << push_back_count << " " << pushed_back_items_count << " " << fill_holes_count << " " << split_leaf_count << endl;
                // stats_file << "q " << rightestright << " " << holes  << endl;
                // stats_file << "q " << rightestright << " " << count_total_holes  << endl;
                stats_file << "i " << rightestright << " " << count_total_holes << " " << height << " " << count_internal_nodes << " " << count_regular_leaves << " " << count_irregular_leaves << endl;
                // if(count_irregular_leaves != leaf_count.second || count_regular_leaves != leaf_count.first)
                // {
                //     cout << "leaf count is wrong! " << count_regular_leaves << " " << leaf_count.first << " " << count_irregular_leaves << " " << leaf_count.second  << endl;
                //     exit(3);
                // }
                // if(count_internal_nodes != node_count - (leaf_count.first + leaf_count.second ))
                // {
                //     cout << "node count is wrong " << count_internal_nodes << " " <<  node_count -(leaf_count.first + leaf_count.second) << endl;
                //     exit(3);
                // }
                
            #endif
        #endif
        q_time = clock() - start_query_time;
        
        // if(i == 0)
        // {
        //     tree.printLeafLinkedList();
        // }
        // q_time = timing();
        // times_file << found_count << " " << q_time <<  "\n";
        times_file << found_count << " " << (double)insert_time/CLOCKS_PER_SEC << " " << (double)q_time/CLOCKS_PER_SEC <<  "\n";


        // cout << "after insertion and query." << endl;
        // if(i == 82960 || i == 82961 || i == 82962) tree.printLeafLinkedList();

        #ifdef linearcheck
        // linear scan
            lin_count = 0;
        
        // // // timing();
        // // // // cout << "data.size() " << data.size() << endl;  
            for(int ii =0; ii < data.size(); ii++){
                // if( ii == 32000694)
                // {
                //     cout << "Query LINEAR \n";
                //     for(int jj = 0; jj < NUMDIMS; jj++)
                //     {
                //         cout << queries[i].m_min[jj] << " " << queries[i].m_max[jj] << ",";
                //     }
                //     cout << "\n";
                //     cout << "Data object 32000694 LINEAR \n";
                //     for(int jj = 0; jj < NUMDIMS; jj++)
                //     {
                //         cout << data[32000694].m_min[jj] << " " << data[32000694].m_max[jj] << ",";
                //     }
                //     cout << "\n";

                // }
                if(Overlap(&(queries[i]), &(data[ii]))) {
                    #ifdef debug
                        cerr << ii << endl;
                    #endif
                    // cerr << ii << " " << data[ii].m_min[0] << " " << data[ii].m_min[1] << " " << data[ii].m_max[0] << " " << data[ii].m_max[1] << endl;
                    lin_count++;
                }
            }
        #endif
        // q_time = timing();
        // times_file << lin_count << " " << q_time <<  "\n";
        // cerr << "QUERY " << i << " " << found_count << endl;
        //cout << "QUERY " << i << " found count: " << found_count << " truth: " << lin_count << endl;
	cout << found_count << endl;
        // cerr << tree.getLeafAreaAverage() << endl;

        // if(i > 61630)
        //     tree.printLeafLinkedList();


        // total_count = tree.CountDataPoints();
        // cout << "COUNT " << total_count << " pending insertions: " << tree.getPendingCount() << "\n";

        // rightestright = tree.getRightestRight();
        // cout << "########## rightest right: " << rightestright << endl;
        // tree.findTrashInAllPlaces();
        // tree.lookForEmptyNodes();
        // cout << "after empty check" << endl;
        // if(rightestright > 3*DATA_COUNT){return 3;}
        // tree.findAllLeavesWithDataInThem(152010);
        // tree.findDataInAllPlaces(14125118);
        // tree.printLeafBoundsOfData(32000694);

        // cout << "Query \n";
        // for(int ii = 0; ii < NUMDIMS; ii++)
        // {
        //     cout << queries[i].m_min[ii] << " " << queries[i].m_max[ii] << ",";
        // }
        // cout << "\n";
        // if(i > 1100){
        // cout << "Data object 32000694 \n";
        //     for(int ii = 0; ii < NUMDIMS; ii++)
        //     {
        //         cout << data[32000694].m_min[ii] << " " << data[32000694].m_max[ii] << ",";
        //     }
        //     cout << "\n";

        // }
        #ifdef linearcheck
            if(lin_count != found_count){
                cout << "WRONG ANSWER, TRAGEDY, DISASTER, ALL THAT..." << endl;
                cout << "Query \n";
                for(int ii = 0; ii < NUMDIMS; ii++)
                {
                    cout << queries[i].m_min[ii] << " " << queries[i].m_max[ii] << ",";
                }
                cout << "\n";

                cout << "Data object 502037 \n";
                for(int ii = 0; ii < NUMDIMS; ii++)
                {
                    cout << data[502037].m_min[ii] << " " << data[502037].m_max[ii] << ",";
                }
                cout << "\n";

            // // //     tree.printLeafBoundsOfData(648515);
            // //     // tree.findAllLeavesWithDataInThem(519955);

            // // //     // tree.printLeafLinkedList();
                return 3;
            }
        #endif


        // #ifdef stats
        //     int rightestright = tree.getRightestRight();
        //     int holes = tree.getCountOfHoles();
        //     // int height = tree.TreeHeight();
        //     // int pending_count = tree.getPendingCount();
        //     // auto leaf_count = tree.CountLeaves();
        //     // int node_count = tree.CountNodes();
        //     // stats_file << node_count << " " << leaf_count.first << " " << leaf_count.second << " " << pending_count << " " << height << " " << rightestright << endl;
        //     stats_file << rightestright << " " << holes << " " << push_back_count << " " << pushed_back_items_count << " " << fill_holes_count << " " << split_leaf_count << endl;
        //     push_back_count = 0;
        //     fill_holes_count = 0;
        //     split_leaf_count = 0;
        //     pushed_back_items_count = 0;
            
        // #endif

        


    }
    // tree.getLeafArea("leafArea_freshAir_CSmain_2m_25k_2r_40ic.txt");


    query_file.close();
    // times_file.close();
    cout << "DONE!!\n";
    return 0;
}

