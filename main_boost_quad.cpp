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
    int data_size = stoi(argv[7]);
    
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
    //float min_z, max_z;
    for(int i = 0; i < data_size; i++){
	//point p1, p2;
	//for(int j=0; j < NUMDIMS; j++){ 
	//	data_file >> min_values[j];
	//	p1.set<j>(min_values[j]);
	//}
	//for(int j=0; j < NUMDIMS; j++){
	//       	data_file >> max_values[j];
	//	p2.set<j>(max_values[j]);
	//}
        data_file >> min_x >> min_y >> max_x >> max_y;
        
	//data_file >> min_x >> min_y>> min_z >> max_x >> max_y >> max_z;

	box b(point(min_x, min_y), point(max_x, max_y));
        //box b(point(min_x, min_y, min_z), point(max_x, max_y, max_z));
	//box b(p1, p2);
	all_data.push_back(std::make_pair(b, i));
    }

    // FATEMEH OUT

    // FATEMEH
    //string time_file_name= "../rtree_results/march2021/uniform_shape_exp/srtree16_uniform_shape_times_rerun.txt";
    // string time_file_name = argv[1];
    ofstream times_file;
    times_file.open(time_file_name.c_str());

    // FATEMEH OUT

    // FATEMEH
    // create the rtree using default constructor
    timing();
    bgi::rtree< value, bgi::quadratic<16> > rtree(all_data);
    double indexing_time = timing();
    // times_file << indexing_time << "\n";
    // cout << "SIZE " << sizeof rtree << endl;
    // FATEMEH OUT



    // FATEMEH
    // read all queries and put them in a vector
    // i do this because when I do it at the same time, it freaks out
    // not important enough to fix:))
    std::vector<box> queries;

    for(int i = 0; i < query_size; i++){
        query_file >> min_x >> min_y >> max_x >> max_y;
//        printf("%d, %d, %d, %d \n", min_x, min_y, max_x, max_y);
        box b(point(min_x, min_y), point(max_x, max_y));
        //point p1, p2;
	//for(int j=0; j < NUMDIMS; j++){
	//      data_file >> min_values[j];
	//      p1.set<j>(min_values[j]);
	//}
	//for(int j=0; j < NUMDIMS; j++){
	//       data_file >> max_values[j];
	//       p2.set<j>(max_values[j]);
	//}
	//box b(p1, p2);
	queries.push_back(b);
    }
    // FATEMEH OUT



    // FATEMEH
    // run all the queries and time them
    int last_id = data_size -1;
    int this_id;
    float min[2]; float max[2];
    clock_t q_time = 0;
    std::vector<value> result_s;
    box b;
    clock_t insert_time, start_insert_time;
    clock_t start_query_time;
    // timing();
    for(int i = 0; i < query_size; i++){
       insert_time = 0.0;
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
                // cout << "inserted " << this_id << " " << min[0] << " " << max[0] << " , " << min[1] << " " << max[1] << endl;
                // data.push_back(Rect(min, max));
            }
        }
        //cout << "QUERY " << i << "\n";
        result_s.clear();
        // timing();
        start_query_time = clock();
        rtree.query(bgi::intersects(queries[i]), std::back_inserter(result_s));
        q_time = clock() - start_query_time;
        times_file << (double)insert_time/CLOCKS_PER_SEC << " " << (double)q_time/CLOCKS_PER_SEC << " " << indexing_time << "\n";
	    cout << "QUERY " << i << " " << result_s.size() << endl;
    }

    query_file.close();
    data_file.close();
    times_file.close();

    cout << "DONE!!\n";

    // FATEMEH OUT


    
    return 0;
}
