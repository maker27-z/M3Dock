#pragma once
//
// Created by 86137 on 2023/4/12.
//
#include "model.h"
#include "conf.h"
#include "incrementable.h"
#include <algorithm>
#include "random.h"
#include "energy.h"
#include "unordered_map"
#ifndef LSHADE_ADAM_FINAL_NSGA2_H
#define LSHADE_ADAM_FINAL_NSGA2_H



struct nsga2 {
    sz num_saved_mins;//最后out的容量
    int num_steps;//迭代次数
    fl min_rmsd;//评估最后的结果
    nsga2() : num_steps(1500), num_saved_mins(20), min_rmsd(1.0) { }
    void operator()(model& m, model& m_ad4,output_container_nsga2& out, const precalculate& p, const precalculate& prec_ad4, const precalculate& prec_vdw, const precalculate& prec_forcefield, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator, const scoring_function& sf,const scoring_function& sf_ad4,const boost::optional<model>& ref) const;
    void init_individual1(model& m, energy_cal& energy, rng& generator, const vec& corner1, const vec& corner2, const scoring_function& sf, output_type_nsga2& tmp) const;
};

struct nsga2_aux {
public:
    int popsize;//种群规模
    int dim;//问题的维度，7+n
    int gen;//种群迭代次数
    std::vector<output_type_nsga2> P;//第一代种群赋给P
    std::vector<output_type_nsga2> Q;//生成第一代种群的子代Q
    std::vector<output_type_nsga2> R;//父子种群合并后的R
    output_type_nsga2 init_tmp;
    nsga2_aux() :popsize(100), dim(7) {}
    nsga2_aux(const model& m) {
        dim = m.ligand_degrees_of_freedom(0) + 6;//只考虑只有一个配体的情况
        popsize = 300;
        P.resize(popsize);
        //Q.resize(popsize);
        R.resize(popsize * 2);
    }
    std::vector<output_type_nsga2> make_new_pop(model& m, const precalculate& p, const igrid& ig, vec v, std::vector<output_type_nsga2> pop, energy_cal& energy, const vec& corner1, const vec& corner2, rng& generator, const scoring_function& sf,const boost::optional<model>& ref,change &g,const scoring_function& sf_ad4);//生成新种群
    std::unordered_map<int, std::vector<output_type_nsga2>> fast_nondominated_sort(std::vector<output_type_nsga2>& pop);//快速非支配排序
    void calculate_objectives(model& m, const precalculate& p, const igrid& ig, vec v, output_type_nsga2& indivisual, const scoring_function& sf, const scoring_function& sf_ad4, energy_cal& energy,const boost::optional<model>& ref,change &g);
    void calculate_crowding_distance(std::vector<output_type_nsga2>& rank_i_vector);//计算拥挤距离
    std::vector<output_type_nsga2> crossover_mutation(model& m, const precalculate& p, const igrid& ig, vec v, energy_cal& energy, rng& generator, const vec& corner1, const vec& corner2, output_type_nsga2 parent1, output_type_nsga2 parent2, const scoring_function& sf,const boost::optional<model>& ref,change &g,const scoring_function& sf_ad4);//交叉变异算子
    output_type_nsga2 binary_tournament(output_type_nsga2 indivisual1, output_type_nsga2 indivisual2);//二元竞标赛选择
    void init_individual(model& m, energy_cal& energy, rng& generator, const vec& corner1, const vec& corner2, const scoring_function& sf, output_type_nsga2& tmp,const boost::optional<model>& ref,const scoring_function& sf_ad4);



};





#endif //LSHADE_ADAM_FINAL_NSGA2_H
