//
// Created by 91686 on 2023/7/10.
//

#ifndef LSHADE_ADAM_FINAL_MPSOD_H
#define LSHADE_ADAM_FINAL_MPSOD_H



#include "model.h"
#include "conf.h"
#include "incrementable.h"
#include <algorithm>
#include "random.h"
#include "energy.h"
#include "unordered_map"

struct MPSOD {
    sz num_saved_mins;//最后out的容量
    int num_steps;//迭代次数
    fl min_rmsd;//评估最后的结果
    MPSOD() : num_steps(1000), num_saved_mins(20), min_rmsd(1.0) { }
    void operator()(model& m, output_container_MPSOD& out, const precalculate& p, const precalculate& prec_forcefield, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator, const scoring_function& sf,const boost::optional<model>& ref) const;
};


struct MPSOD_aux {
public:
    int popsize;//种群规模
    int dim;//问题的维度，7+n
    std::vector<output_type_MPSOD> pop;

    fl ** lamda;
    int T_size;
    fl delta;//delta: probability that parent solutions are selected from neighbourhood


    std::vector<output_type_MPSOD> indArray;

    int ** neighborhood; // neighborhood[i][T]接近第i个λ的有T_个权重向量
    int neighborhood_size;

    int nr;
    fl * z; // z相当于一个解，vector(objectives.size())
    std::vector<output_type_MPSOD> next;


    //MOPSO
    fl w;
    fl r1;
    fl r2;
    fl c1;
    fl c2;

    int archiveMaxSize;
    std::vector<output_type_MPSOD> gbest;
    std::vector<output_type_MPSOD> pbest;
    std::vector<output_type_MPSOD> velocity;


    MPSOD_aux() :popsize(100), archiveMaxSize(100), dim(7) {}
    MPSOD_aux(const model& m,const vec& corner1, const vec& corner2) {
        dim = m.ligand_degrees_of_freedom(0) + 6;//只考虑只有一个配体的情况
        popsize = dim * 8;
        pop.resize(popsize);
        next.resize(popsize);
        archiveMaxSize = dim * 8;
        w = 0.1;

        lamda = new fl*[popsize];
        neighborhood = new int*[popsize];
        z = new double[pop[0].objectives.size()];
        T_size = 0.1 * popsize;
        delta = 0.9;
//        nr = 2;


        //MOPSO
        gbest.resize(popsize);
        pbest.resize(popsize);
        velocity.resize(popsize);



    }
    void initUniformWeight();
    void initNeighborhood();
    void initPopulation(model& m,  const precalculate& p, const precalculate& prec_forcefield, const igrid& ig,
                        const precalculate& p_widened,  const vec& corner1, const vec& corner2,
                        rng& generator,const scoring_function& sf,const boost::optional<model>& ref);
    void initpopVector(model& m,  const precalculate& p, const precalculate& prec_forcefield, const igrid& ig,
                        const precalculate& p_widened,  const vec& corner1, const vec& corner2,
                        rng& generator,const scoring_function& sf,const boost::optional<model>& ref,std::vector<output_type_MPSOD>& popObj);
    void initIdealPoint();
    void initNext(model& m,  const precalculate& p, const precalculate& prec_forcefield, const igrid& ig,
                        const precalculate& p_widened,  const vec& corner1, const vec& corner2,
                        rng& generator,const scoring_function& sf,const boost::optional<model>& ref);
    std::vector<output_type_MPSOD> matingSelection(model& m,  const precalculate& p, const precalculate& prec_forcefield, const igrid& ig,
                                                   const precalculate& p_widened,  const vec& corner1, const vec& corner2,
                                                   rng& generator,const scoring_function& sf,const boost::optional<model>& ref,energy_cal& energy);
    void updateReference(const std::vector<output_type_MPSOD>& ind);
    void updateProblem(output_type_MPSOD  indiv, int id, int type,rng& generator);
    fl fitnessFunction(output_type_MPSOD  individual, fl * lamda_k);
    output_type_MPSOD DECrossover(model& m, energy_cal& energy, rng& generator, const vec& corner1, const vec& corner2,std::vector<output_type_MPSOD> parents,const scoring_function& sf,const boost::optional<model>& ref);
    void mutation(const vec& corner1, const vec& corner2,rng& generator,output_type_MPSOD& child);
    void calculate_objectives(model& m, const precalculate& p, const igrid& ig, vec v, output_type_MPSOD& indivisual, const scoring_function& sf, energy_cal& energy,const boost::optional<model>& ref,change &g);
    void calculate_crowding_distance(std::vector<output_type_MPSOD>& pop);
    std::unordered_map<int, std::vector<output_type_MPSOD>> fast_nondominated_sort(std::vector<output_type_MPSOD>& pop);//快速非支配排序
    void classification(model& m,  const precalculate& p, const precalculate& prec_forcefield, const igrid& ig,
                        const precalculate& p_widened,  const vec& corner1, const vec& corner2,
                        rng& generator,const scoring_function& sf,const boost::optional<model>& ref,energy_cal& energy);
    std::vector<std::vector<fl>>  divide();

    std::vector<output_type_MPSOD> TournamentSelection(int K,int N,rng& generator);
    //MOPSO
    //初始化每个粒子的速度
    void initVelocity(model& m,  const precalculate& p, const precalculate& prec_forcefield, const igrid& ig,
                      const precalculate& p_widened,  const vec& corner1, const vec& corner2,
                      rng& generator,const scoring_function& sf,const boost::optional<model>& ref);

    void initpbest(model& m,  const precalculate& p, const precalculate& prec_forcefield, const igrid& ig,
                   const precalculate& p_widened,  const vec& corner1, const vec& corner2,
                   rng& generator,const scoring_function& sf,const boost::optional<model>& ref);

    void computeSpeed(model& m,  const precalculate& p, const precalculate& prec_forcefield, const igrid& ig,
                      const precalculate& p_widened,  const vec& corner1, const vec& corner2,
                      rng& generator,const scoring_function& sf,const boost::optional<model>& ref);

    void computeNewPositions(const vec& corner1, const vec& corner2,std::vector<output_type_MPSOD>& parents);

    void DE_Correction(const vec& corner1, const vec& corner2,std::vector<output_type_MPSOD>& parents);

    void mutation(const vec& corner1, const vec& corner2,rng& generator,std::vector<output_type_MPSOD>& parents);
};

#endif //LSHADE_ADAM_FINAL_MPSOD_H
