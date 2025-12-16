//
// Created by 91686 on 2023/11/16.
//

#include "weighted_terms_xscore.h"



fl weighted_terms_xscore::VDW_term(sz t1, sz t2, fl r) const {
    float d0, d, cutoff, tmp, tmp1, tmp2, asum, sum;
    cutoff = DIST_CUTOFF;

    d0 = XSCORE_radius(t1) + XSCORE_radius(t2);
    d = r;
    if(d > cutoff) return 0.0;

    tmp1=d0/d;
    tmp1=tmp1*tmp1*tmp1*tmp1; tmp2=tmp1*tmp1;
    tmp=tmp2-2.00*tmp1;

    asum = tmp;
    asum*=(-1.00);
    if(asum<0.00) return 0.0;

    sum = asum;

    return sum;

}

fl weighted_terms_xscore::HB_term(sz t1, sz t2, fl r) const {

}
fl weighted_terms_xscore::HM_term(sz t1, sz t2, fl r) const {

}
fl weighted_terms_xscore::RT_term(sz t1, sz t2, fl r) const {

}

fl weighted_terms_xscore::eval(sz t1, sz t2, fl r) const {
    fl total_e = 0;
    total_e = VDW_term(t1,t2,r) + HB_term(t1,t2,r) + HM_term(t1,t2,r) + RT_term(t1,t2,r);
    return total_e;

}