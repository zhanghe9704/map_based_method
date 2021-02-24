#ifndef DA_FUNCTIONS_H_INCLUDED
#define DA_FUNCTIONS_H_INCLUDED

#include "..\tpsa_lib\include\da.h"

DAVector potential(std::vector<DAVector> &field, unsigned int nv);
DAVector da_gmd(DAVector& f, std::vector<DAVector>& g, unsigned int nv);
DAVector da_flow(std::vector<DAVector>& f, DAVector& g, unsigned int nv);
DAVector lie_exp(DAVector& f, DAVector& g, unsigned int nv);
void lie_factorization(std::vector<DAVector>& m, unsigned int nv, int no, std::vector<double>& c,
                       std::vector<DAVector>& l, std::vector<DAVector>& f, bool combine = false);
void lie_factorization_inverse_order(std::vector<DAVector>& m, unsigned int nv, int no, std::vector<double>& c,
                       std::vector<DAVector>& l, std::vector<DAVector>& f, bool combine = false);

#endif // GF_HPP_INCLUDED
