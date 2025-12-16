//
// Created by 91686 on 2023/11/17.
//

#ifndef LSHADE_ADAM_FINAL_WEIGHTED_TERMS_DLIGAND2_H
#define LSHADE_ADAM_FINAL_WEIGHTED_TERMS_DLIGAND2_H
#include "common.h"
#include "atom.h"
template<typename T>
class two_dim_array{
    std::vector<std::vector<T>> m_data;
    sz m_dim;

public:
    two_dim_array() {}
    two_dim_array(sz rows, sz columns, const T& filler_val) :  m_dim(rows+columns) {
        m_data.resize(rows);
        for(auto& row : m_data) {
            row.assign(columns,filler_val); // 设置每一行的列数
        }
    }

    T& operator()(sz i, sz j){
        return m_data[i][j];
    }
    const T& operator()(sz i, sz j) const{
        return m_data[i][j];
    }

};

struct weighted_terms_DLIGAND2  {
    weighted_terms_DLIGAND2() {
        cutoff_ = 15;
    }

    weighted_terms_DLIGAND2(const flv& weights); // does not own t
    atom_type::t atom_typing_used() const { return atom_typing_used_; }
    fl cutoff() const { return cutoff_; }
    fl eval(sz t1, sz t2, fl r) const; // intentionally not checking for cutoff
    fl eval(const atom& a,const  atom& b, fl r) const;


    int get_DLIGAND_lig_atomtype_nums() const {return DLIGAND_lig_atomtype_nums;}
    int get_DLIGAND_grid_atomtype_nums() const{return DLIGAND_grid_atomtype_nums;}


private:

    flv weights;
    atom_type::t atom_typing_used_;
    fl cutoff_;
    szv enabled_usable_terms;
    int DLIGAND_grid_atomtype_nums = 167;
    int DLIGAND_lig_atomtype_nums = 11;
};


#endif //LSHADE_ADAM_FINAL_WEIGHTED_TERMS_DLIGAND2_H
