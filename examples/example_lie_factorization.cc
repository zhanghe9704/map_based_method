#include "../tpsa_lib/include/da.h"
#include "../include/da_functions.h"
#include <iostream>
#include "../include/map.hpp"
#include "../include/gf.hpp"

using std::vector;
using std::cout;
using std::endl;

/******
* In this example, we read the truncated map from an output of COSY Infinity 9.x
* and perform the Lie factorization for it.
******/

int main() {
    int dim = 2;
    int da_order = 10;
    int da_dim = 2*dim;
    int da_scratch = 5000;
    //Initialize the da environment. Reserve one more order for the following Lie factorization calculation.
    da_init(da_order+1, da_dim, da_scratch);

    //Read the truncated map from file.
    vector<Map> trunc_map(da_dim, Map(dim, da_full_length()));
    char filename[100] = "test_Liemap_10th_order_4_dim.txt"; //Truncated Taylor map
    readCOSYMap(dim, filename, trunc_map);

    //Print the map to screen.
    for(auto &v: trunc_map)
        cout<<v<<endl;

    //Convert the map into a group of da vectors. (For demonstration. No necessary to perform tracking.)
    vector<DAVector> da_map(da_dim);
    for(int i=0; i<da_dim; ++i) {
            map2da(trunc_map.at(i), da_map.at(i));
    }

    for(auto&v: da_map) v.print();

    vector<double> c;
    vector<DAVector> l, fo;
    lie_factorization(da_map, da_dim, da_order, c, l , fo);

    cout<<"Lie factorization: "<<endl;
    cout<<"Constants: "<<endl;
    for(int i=0; i<da_dim; ++i) {
    cout<<c.at(i)<<endl;
    }
    cout<<"Linear map: "<<endl;
    for(auto&v : l) v.print();
    cout<<"f_i: "<<endl;
    for(auto&v : fo) v.print();

    lie_factorization_inverse_order(da_map, da_dim, da_order, c, l , fo);
    cout<<"Lie factorization in verse order: "<<endl;
    cout<<"Constants: "<<endl;
    for(int i=0; i<da_dim; ++i) {
    cout<<c.at(i)<<endl;
    }
    cout<<"Linear map: "<<endl;
    for(auto&v : l) v.print();
    cout<<"f_i: "<<endl;
    for(auto&v : fo) v.print();

    da_set_eps(1e-15);
    lie_factorization(da_map, da_dim, da_order, c, l , fo ,true);

    cout<<"Lie factorization: "<<endl;
    cout<<"Constants: "<<endl;
    for(int i=0; i<da_dim; ++i) {
    cout<<c.at(i)<<endl;
    }
    cout<<"Linear map: "<<endl;
    for(auto&v : l) v.print();
    cout<<"f_i: "<<endl;
    for(auto&v : fo) v.print();

    lie_factorization_inverse_order(da_map, da_dim, da_order, c, l , fo, true);
    cout<<"Lie factorization in verse order: "<<endl;
    cout<<"Constants: "<<endl;
    for(int i=0; i<da_dim; ++i) {
    cout<<c.at(i)<<endl;
    }
    cout<<"Linear map: "<<endl;
    for(auto&v : l) v.print();
    cout<<"f_i: "<<endl;
    for(auto&v : fo) v.print();

    return 0;
}
