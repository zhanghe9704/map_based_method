#ifndef GF_HPP_INCLUDED
#define GF_HPP_INCLUDED

#include "map.hpp"
#include "../tpsa_lib/include/da.h"

void map2da(Map& m, DAVector& d);
void da2map(DAVector& d, Map& m);
void generating_function(std::vector<DAVector>& ivecs, int dim, const int type, DAVector& gf) ;
void gf_eqns(std::vector<DAVector>& map, int dim, const int type, std::vector<DAVector>& eqns);
void gf2map(DAVector& gf, int dim, const int type, std::vector<DAVector>& m);
void map_symp(std::vector<DAVector>& truc_map, const int dim, const int gf_type, std::vector<DAVector>& sym_map);
#endif // GF_HPP_INCLUDED
