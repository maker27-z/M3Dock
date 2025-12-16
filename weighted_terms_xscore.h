//
// Created by 91686 on 2023/11/16.
//

#include "common.h"
#include "atom.h"
#include "xtools.h"
#ifndef LSHADE_ADAM_FINAL_WEIGHTED_TERMS_XSCORE_H
#define LSHADE_ADAM_FINAL_WEIGHTED_TERMS_XSCORE_H


struct weighted_terms_xscore  {
    weighted_terms_xscore(const flv& weights); // does not own t
    atom_type::t atom_typing_used() const { return atom_typing_used_; }
    fl cutoff() const { return cutoff_; }
    fl eval(sz t1, sz t2, fl r) const; // intentionally not checking for cutoff
    fl eval(const atom& a,const  atom& b, fl r) const;
    fl VDW_term(sz t1, sz t2, fl r) const;
    fl HB_term(sz t1, sz t2, fl r) const;
    fl HM_term(sz t1, sz t2, fl r) const;
    fl RT_term(sz t1, sz t2, fl r) const;


private:
    weighted_terms_xscore() {}

    flv weights;
    atom_type::t atom_typing_used_;
    fl cutoff_;
    szv enabled_usable_terms;
};



#endif //LSHADE_ADAM_FINAL_WEIGHTED_TERMS_XSCORE_H
