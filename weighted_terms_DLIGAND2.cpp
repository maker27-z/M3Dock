//
// Created by 91686 on 2023/11/17.
//

#include "weighted_terms_DLIGAND2.h"
#include "molecule.h"

fl weighted_terms_DLIGAND2::eval(sz t1, sz t2, fl r2) const {
    fl total_e = 0;
    r2 = r2*r2;
    if(t1 < 0 || t1 >= matype_pro) return 0;
    int k1 = iskind[t1];
    if(k1 < 0) return 0;

    if(t2 < 0 || t2 < matype_pro) return 0;
    int k2 = iskind[t2];
//    k2 += matype_pro;
    if(k2 < 0) return 0;


    if(r2 > Rcut*Rcut) return 0;
    int b = r2bin(sqrt(r2));

    total_e = fl(edfire[k1][k2][b]);
    return total_e;

}