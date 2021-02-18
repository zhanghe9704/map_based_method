#include "../include/da_functions.h"
#include <assert.h>


DAVector potential(std::vector<DAVector>& field, unsigned int nv) {
    DAVector ptl;
    for(unsigned int i=0; i<nv; ++i) {
        DAVector ff = field.at(i);
        for(unsigned int j=i+1; j<nv; ++j) {
            DAVector s;
            da_substitute_const(ff, j, 0, s);
            ff = s;
        }
        ptl = ptl + da_int(ff, i);
    }
    return ptl;
}

DAVector possion_bracket(DAVector& x, DAVector& y, unsigned int nv) {
    DAVector res = 0;
    for(unsigned int i=0; i<nv; i+=2) {
        res += da_der(x, i)*da_der(y,i+1) - da_der(x, i+1)*da_der(y, i);
    }
    return res;
}

DAVector da_gmd(DAVector& f, std::vector<DAVector>& g, unsigned int nv) {
    DAVector res = 0;
    for(unsigned int i=0; i<nv; ++i) {
        res += da_der(f,i)*g.at(i);
    }
    return res;
}

DAVector da_flow(std::vector<DAVector>& f, DAVector& g, unsigned int nv, double eps) {
    DAVector dg = da_gmd(g, f, nv);
    DAVector res = g + dg;
    double cnt = 1.0;

    while(dg.norm()>eps) {
        cnt += 1;
        dg /= cnt;
        dg = da_gmd(dg, f, nv);
        res += dg;
        if(cnt>100) {
            assert(false&&"Error in da_flow, no convergence after 100 iterations.");
        }
    }
    return res;
}

DAVector lie_exp(DAVector& f, DAVector& g, unsigned int nv, double eps) {
    std::vector<DAVector> x(nv);
    for(unsigned int i=0; i<nv; i+=2) {
        x.at(i) = -da_der(f, i+1);
        x.at(i+1) = da_der(f, i);
    }
    return da_flow(x, g, nv, eps);
}

/** \brief Calculate the Lie factorization from a Taylor map.
 * Calculate the Lie factorization M(x) =n (L exp(:f3:)exp(:f4:)...exp(:f(n+1):)) x + C from a Taylor map,
 * where L is the linear map, C is the constants.
 * Because the function f will be calculated up to (n+1)-th order for an n-th order Taylor map, one needs to initiallzie
 * the DA environment with at least the (n+1)-th order to get the correct result.
 * \param m The truncated Taylor map.
 * \param nv Number of DA variables.
 * \param no Order of the truncated map.
 * \param c Constants of the map.
 * \param l Linear map.
 * \param f The fi of the Lie factorization.
 * \param eps If the absolute value of a coefficient in a DA vector is less than eps, it will be set to zero.
 * \return none.
 *
 */
void lie_factorization(std::vector<DAVector>& m, unsigned int nv, int no, std::vector<double>& c,
                       std::vector<DAVector>& l, std::vector<DAVector>& f, double eps) {
    c.clear();
    l.clear();
    f.clear();
    std::vector<DAVector> n;
    for(unsigned int i=0; i<nv; ++i) {
        c.push_back(m.at(i).con());     //The constants of the map m.
        n.push_back(m.at(i)-c.at(i));
    }
    da_change_order(1);
    for(unsigned int i=0; i<nv; ++i) l.push_back(n.at(i));  //Linear map.

    std::vector<DAVector> ln(nv);
    inv_map(l, nv, ln);  //Inverse linear map: ln.

    da_restore_order();

    std::vector<DAVector> nl(nv);
    da_composition(n, ln, nl);

    for(int i=1; i<no; ++i) {
        da_change_order(i+1);
        for(unsigned int j=0; j<nv; j+=2) {
            ln.at(j) = nl.at(j+1) - da[j+1];
            ln.at(j+1) = -nl.at(j) + da[j];
        }

        da_change_order(i+2);
        DAVector fi = potential(ln, nv);
        fi.clean(eps);
        f.push_back(fi);

        fi *= -1;
        da_restore_order();
        for(unsigned int j=0; j<nv; ++j) {
            nl.at(j) = lie_exp(fi, nl.at(j), nv, eps);
        }
    }

    da_restore_order();
}
