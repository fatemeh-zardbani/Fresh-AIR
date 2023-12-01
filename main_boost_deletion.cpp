//  g++ -I ./boost_1_61_0/ main.cpp -std=c++14 -o folan && ./folan
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>

#include <boost/geometry/index/rtree.hpp>

// to store queries results
#include <vector>

// just for output
#include <iostream>
#include <boost/foreach.hpp>

// FATEMEHi
#include <stdlib.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;
//#define NUMDIMS 8


double timing(){
    static struct timeval t1, t2;
    gettimeofday(&t2,NULL);
    double ret = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) * 1e-6;
    t1 = t2;
    return ret;
}
// FATEMEH OUT

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

int main(int argc, char **argv)
{


    //  auto seed = time(NULL); // 1695893019
    auto seed = 0;
//    auto seed = 1697617670;
//    cout << seed << endl;
    // auto seed = 1692015686; // crashes at 81199
    srand(seed);


    
    typedef bg::model::point<float, 2, bg::cs::cartesian> point;
    typedef bg::model::box<point> box;
    typedef std::pair<box, unsigned> value;

    // FATEMEH 
    //int data_size = 19349214;
    //int data_size = 7550261;
    // int data_size = 70380191;
    // int data_size = 500000;

    //int data_size = 64000000;
    // in ghalate  bara 3d, az file e main_3d.cpp estefade kon
    // string data_file_name = "../rtree_data/edges_mbrs_float.txt";
    //string data_file_name = "../rtree_data/server_uniform_shape_data_exp_64m_float.txt";
    
    // int query_size = 100000;
    //int query_size = 10000;
    //uniform_3d_data_points_double.txt 500000 ../rtree_data/uniform_3d_query_500k.txt
    // string query_file_name = "../rtree_data/usa_e-2%_uniform_1m_other_float.txt";
    // server_uniform_shape_data_exp_64m_float_double_points.txt 500000 ../rtree_data/server_uniform_query_float4.txt
    //string query_file_name = "../rtree_data/server_uniform_query_float4_partialmatch_10k.txt";


    string data_file_name = argv[1];
    int query_size = stoi(argv[2]);
    string query_file_name = argv[3];
    string time_file_name = argv[4];
    int ratio = stoi(argv[5]);
    int insert_count = stoi(argv[6]);
    int delete_count = stoi(argv[7]);
    int data_size = stoi(argv[8]);
    string delete_file_name = argv[9];
    
    std::ifstream data_file(data_file_name.c_str());
    if(data_file.is_open()){cout << "it is open:)"<< data_file_name <<"\n";}
    std::ifstream query_file(query_file_name.c_str());
    if(query_file.is_open()){cout << "it is open:)\n";}


    int query_count = std::count(std::istreambuf_iterator<char>(query_file),
                                 std::istreambuf_iterator<char>(), '\n');
    int data_count = std::count(std::istreambuf_iterator<char>(data_file),
                                 std::istreambuf_iterator<char>(), '\n');
    cout << data_count << " " << query_count << "\n";
    if(query_count < query_size) {cout<<"NOT ENOUGH QUERIES IN THE FILE\n"; return 1;}
    if(data_count < data_size) {cout<<"NOT ENOUGH DATA IN THE FILE\n"; return 1;}

    data_file.clear();
    data_file.seekg(0, ios::beg);
    query_file.clear();
    query_file.seekg(0, ios::beg);

    // read all data and put into vector
    std::vector<value> all_data;
    //float min_values[NUMDIMS]; float max_values[NUMDIMS];
    float min_x, min_y, max_x, max_y;
    for(int i = 0; i < data_size; i++){
        data_file >> min_x >> min_y >> max_x >> max_y;
        box b(point(min_x, min_y), point(max_x, max_y));
        all_data.push_back(std::make_pair(b, i));
        if(i == 3229830){
            cout << "3229830 box: " << min_x << " " << min_y << " " << max_x << " " << max_y << endl;
        }
    }


    ofstream times_file;
    times_file.open(time_file_name.c_str());


    timing();
    bgi::rtree< value, bgi::rstar<16> > rtree(all_data);
    double indexing_time = timing();

    std::vector<box> queries;

    for(int i = 0; i < query_size; i++){
        query_file >> min_x >> min_y >> max_x >> max_y;
        box b(point(min_x, min_y), point(max_x, max_y));
	    queries.push_back(b);
    }
    // FATEMEH OUT



    int how_many_deleted_till_now = 0;
    std::ifstream delete_file(delete_file_name.c_str());
    std::vector<value> tbd;
    int trash_id;
    for(int i = 0; i < delete_count*(query_size/ratio - 1); i++){
        delete_file >> min_x >> min_y >> max_x >> max_y >> trash_id;
        box b(point(min_x, min_y), point(max_x, max_y));
        tbd.push_back(std::make_pair(b, trash_id));
    }



    // FATEMEH
    // run all the queries and time them

    int total_data_count = data_size;
    int last_id = data_size -1;
    int this_id;
    float min[2]; float max[2];
    clock_t q_time = 0;
    std::vector<value> result_s;
    int item_index_to_be_deleted_choice;
    box b;
    clock_t insert_time, start_insert_time, delete_time, start_delete_time;
    clock_t start_query_time;
    // timing();
    for(int i = 0; i < query_size; i++){
       insert_time = 0.0;
       delete_time = 0.0;
        if(i % ratio == 0 && i != 0){
            for(int o =0; o < insert_count; o++){
                // insert one
                for(int j = 0; j < 2; j++) data_file >> min[j];
                for(int j = 0; j < 2; j++) data_file >> max[j];
                this_id = ++last_id;
                b = box(point(min[0], min[1]), point(max[0], max[1]));
                // timing();
                start_insert_time = clock();
                rtree.insert(make_pair(b, this_id));
                insert_time += clock() - start_insert_time;
                all_data.push_back(std::make_pair(b, this_id));
                // cout << "inserted " << this_id << " " << min[0] << " " << max[0] << " , " << min[1] << " " << max[1] << endl;
                // data.push_back(Rect(min, max));
            }

            for(int o =0; o < delete_count; o++){
                // item_index_to_be_deleted_choice = (rand() % all_data.size());
                item_index_to_be_deleted_choice = tbd[how_many_deleted_till_now].second;
                start_delete_time = clock();
                if(rtree.remove(all_data[item_index_to_be_deleted_choice]) == 0){
                    cout << "item " << item_index_to_be_deleted_choice << " to be deleted not found in tree" << endl;
                    exit(5);
                }
                // rtree.remove(tbd[how_many_deleted_till_now]);
                // if(rtree.remove(tbd[how_many_deleted_till_now]) == 0){
                //     cout << "item " << tbd[how_many_deleted_till_now].second << " with box "  << tbd[how_many_deleted_till_now].first.min_corner().get<0>() <<  " " << tbd[how_many_deleted_till_now].first.min_corner().get<1>() << " to be deleted not found in tree, last id: " << last_id << " " << this_id << endl;
                //     exit(5);
                // }
                how_many_deleted_till_now++;
                delete_time += clock() - start_delete_time;
                // all_data.erase(all_data.begin() + item_index_to_be_deleted_choice);
            }
        }
        //cout << "QUERY " << i << "\n";
        result_s.clear();
        // timing();
        start_query_time = clock();
        rtree.query(bgi::intersects(queries[i]), std::back_inserter(result_s));
        q_time = clock() - start_query_time;
        times_file << (double)insert_time/CLOCKS_PER_SEC << " " << (double)delete_time/CLOCKS_PER_SEC << " " << (double)q_time/CLOCKS_PER_SEC << " " << indexing_time << "\n";
	    cout << "QUERY " << i << " " << result_s.size() << endl;
    }

    query_file.close();
    data_file.close();
    times_file.close();

    cout << "DONE!!\n";

    // FATEMEH OUT


    
    return 0;
}
