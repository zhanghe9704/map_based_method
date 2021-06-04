#include <iostream>
#include "../include/map.hpp"
#include "../include/tracking.hpp"

/******
* In this example, we read the truncated map from an output of COSY Infinity 9.x
* and track the particle with a given initial coordinates using the map.
* Then we symplificate the truncated map and track the same particle using the symplificated map.
* Use the python script plot.py to plot the tracking result in phase space and save them in x-xp.pdf and y-yp.pdf.
******/
int main() {
    int dim = 2;
    int da_order = 4;  //This should be 1 above the order of the map for gf calculation.
    int da_dim = 2*dim;
    int da_scratch = 1000;
    da_init(da_order, da_dim, da_scratch);  //Initialize the da environment.
    da_change_order(da_order-1);            //Reduce the order by one.

    //Read the truncated map from file.
    vector<Map> trunc_map(da_dim, Map(dim, da_full_length()));
    char filename[100] = "4D_2ndOrder_COSY.txt";
    readCOSYMap(dim, filename, trunc_map);

    //Print the map to screen.
    std::cout<<"Load map: "<<std::endl;
    for(auto &v: trunc_map)
        std::cout<<v<<std::endl;

    char output[100] = "track_truncated_map_cosy.txt";

    double xi[4] = {3.3000000000000000E-02,   0.000000000000000,  3.3000000000000000E-02,   0.000000000000000};
    double xf[4] = {0,0,0,0};

    int n_turns = 1000000;
    bool output_to_file = true;
    int dn = 1000;
    mapTrack(dim, xi, trunc_map, n_turns, xf, output_to_file, output, dn);


    //Convert the Map type into DA
    vector<DAVector> da_map(da_dim);
    for(int i=0; i<da_dim; ++i) {
            map2da(trunc_map.at(i), da_map.at(i));
            da_map.at(i) = da_map.at(i) - da_map.at(i).con();   //Remove the constant part!
    }


    int gf_type = 1; //Choose the type of the generating function.
    map_symp(da_map, dim, gf_type, da_map);

    //Print the map to screen again as DA vectors.
    std::cout<<"Symp DA map: "<<std::endl;
    for(auto&v: da_map) v.print();

    //Convert DA type into Map.
    vector<Map> symp_map(da_dim, Map(dim, da_full_length()));
    for(int i=0; i<da_dim; ++i) {
        da2map(da_map.at(i), symp_map.at(i));
    }

    //Tracking use the symplectic map.
    char output_symp[100] = "track_symplectic_map_type_1_cosy.txt";
    mapTrack(dim, xi, symp_map, n_turns, xf, output_to_file, output_symp, dn);
    return 0;
}

