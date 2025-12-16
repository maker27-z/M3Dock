// //
//// Created by 91686 on 2023/7/6.
////
//
//#ifndef LSHADE_ADAM_FINAL_MOEAD_H
//#define LSHADE_ADAM_FINAL_MOEAD_H
//
//
//#include "model.h"
//#include "conf.h"
//#include "incrementable.h"
//#include <algorithm>
//#include "random.h"
//#include "energy.h"
//#include "unordered_map"
//
//struct MOEAD {
//    sz num_saved_mins;//最后out的容量
//    int num_steps;//迭代次数
//    fl min_rmsd;//评估最后的结果
//    MOEAD() : num_steps(100000), num_saved_mins(20), min_rmsd(1.0) { }
//    void operator()(model& m, output_container_MOEAD& out, const precalculate& p, const precalculate& prec_forcefield, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator, const scoring_function& sf,const boost::optional<model>& ref) const;
//};
//
//
//struct MOEAD_aux {
//public:
//    int popsize;//种群规模
//    int dim;//问题的维度，7+n
//    std::vector<output_type_MOEAD> pop;
//
//    fl ** lamda;
//    int T_size;
//    fl delta;//delta: probability that parent solutions are selected from neighbourhood
//
//
//    std::vector<output_type_MOEAD> indArray;
//
//    int ** neighborhood; // neighborhood[i][T]接近第i个λ的有T_个权重向量
//    int neighborhood_size;
//
//    int nr;
//    fl * z; // z相当于一个解，vector(objectives.size())
//
//    MOEAD_aux() :popsize(100), dim(7) {}
//    MOEAD_aux(const model& m,const vec& corner1, const vec& corner2) {
//        dim = m.ligand_degrees_of_freedom(0) + 6;//只考虑只有一个配体的情况
//        popsize = dim * 8;
//        pop.resize(popsize);
//        lamda = new fl*[popsize];
//        neighborhood = new int*[popsize];
//        z = new double[pop[0].objectives.size()];
//        T_size = 0.1 * popsize;
//        delta = 0.9;
////        nr = 2;
//    }
//    void initUniformWeight();
//    void initNeighborhood();
//    void initPopulation(model& m,  const precalculate& p, const precalculate& prec_forcefield, const igrid& ig,
//                        const precalculate& p_widened,  const vec& corner1, const vec& corner2,
//                        rng& generator,const scoring_function& sf,const boost::optional<model>& ref);
//    void initIdealPoint();
//    void matingSelection(std::vector<int> &list, int cid, int size, int type,rng& generator);
//    void updateReference(const output_type_MOEAD& ind);
//    void updateProblem(output_type_MOEAD  indiv, int id, int type,rng& generator);
//    fl fitnessFunction(output_type_MOEAD  individual, fl * lamda_k);
//    output_type_MOEAD DECrossover(model& m, energy_cal& energy, rng& generator, const vec& corner1, const vec& corner2,std::vector<output_type_MOEAD> parents,const scoring_function& sf,const boost::optional<model>& ref);
//    void mutation(const vec& corner1, const vec& corner2,rng& generator,output_type_MOEAD& child);
//    void calculate_objectives(model& m, const precalculate& p, const igrid& ig, vec v, output_type_MOEAD& indivisual, const scoring_function& sf, energy_cal& energy,const boost::optional<model>& ref,change &g);
//    void deleteParams();
//};
//
//
//#endif //LSHADE_ADAM_FINAL_MOEAD_H
