//
// Created by 91686 on 2023/7/3.
//


#include "model.h"
#include "conf.h"
#include "incrementable.h"
#include <algorithm>
#include "random.h"
#include "energy.h"
#include "unordered_map"

#ifndef LSHADE_ADAM_FINAL_SMPSO_H
#define LSHADE_ADAM_FINAL_SMPSO_H

struct SMPSO {
    sz num_saved_mins;//最后out的容量
    int num_steps;//迭代次数
    fl min_rmsd;//评估最后的结果
    SMPSO() : num_steps(1000), num_saved_mins(20), min_rmsd(1.0) { }
    void operator()(model& m, model& m_ad4,output_container_SMPSO& out, const precalculate& p, const precalculate& prec_ad4, const precalculate& prec_vdw, const precalculate& prec_forcefield, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator, const scoring_function& sf,const scoring_function& sf_ad4,const boost::optional<model>& ref) const;
};


struct SMPSO_aux {
public:
    int popsize;//种群规模
    int dim;//问题的维度，7+n
    fl w;
    fl r1;
    fl r2;
    fl c1;
    fl c2;
    fl *deltaMax;
    fl *deltaMin;
    fl *upperLimit;
    fl *lowerLimit;
    int archiveMaxSize;
    output_type_SMPSO gbest;
    std::vector<output_type_SMPSO> pbest;
    std::vector<output_type_SMPSO> velocity;

    std::vector<output_type_SMPSO> pop;
    std::vector<output_type_SMPSO> archive;
    output_type_SMPSO init_tmp;

    SMPSO_aux() :popsize(100), archiveMaxSize(100), dim(7) {}
    SMPSO_aux(const model& m,const vec& corner1, const vec& corner2) {
        dim = m.ligand_degrees_of_freedom(0) + 6;//只考虑只有一个配体的情况
        popsize = 150;
        archiveMaxSize = 150;
        w = 0.1;
        pop.resize(popsize);
        archive.resize(archiveMaxSize);
        pbest.resize(popsize);
        velocity.resize(popsize);
        lowerLimit = new double[dim];
        upperLimit = new double[dim];
        deltaMax = new double[dim];
        deltaMin = new double[dim];
        for (int i = 0; i < dim; i++) {
            if(i < 3){
                lowerLimit[i] = corner1[i];
                upperLimit[i] = corner2[i];
                deltaMax[i] = (upperLimit[i] - lowerLimit[i])/2.0;
                deltaMin[i] = -deltaMax[i];
            }
            else if(i >= 3 && i < 6){
                if(i==3){
                    lowerLimit[i] = 0;
                    upperLimit[i] = 2 * pi;
                    deltaMax[i] = (upperLimit[i] - lowerLimit[i])/2.0;
                    deltaMin[i] = -deltaMax[i];
                }else if(i==4){
                    lowerLimit[i] = 0;
                    upperLimit[i] = pi;
                    deltaMax[i] = (upperLimit[i] - lowerLimit[i])/2.0;
                    deltaMin[i] = -deltaMax[i];
                }else{
                    lowerLimit[i] = 0;
                    upperLimit[i] = 2 * pi;
                    deltaMax[i] = (upperLimit[i] - lowerLimit[i])/2.0;
                    deltaMin[i] = -deltaMax[i];
                }
            }
            else{
                lowerLimit[i] = -pi;
                upperLimit[i] = pi;
                deltaMax[i] = (upperLimit[i] - lowerLimit[i])/2.0;
                deltaMin[i] = -deltaMax[i];
            }
        }
    }
    void initpop(model& m,  model& m_ad4, const precalculate& p, const precalculate& prec_ad4, const precalculate& prec_vdw, const precalculate& prec_forcefield, const igrid& ig,
                 const precalculate& p_widened,  const vec& corner1, const vec& corner2,
                 rng& generator,const scoring_function& sf,const scoring_function& sf_ad4,const boost::optional<model>& ref);

    void initVelocity(model& m,  model& m_ad4, const precalculate& p, const precalculate& prec_ad4, const precalculate& prec_vdw, const precalculate& prec_forcefield, const igrid& ig,
                      const precalculate& p_widened,  const vec& corner1, const vec& corner2,
                      rng& generator,const scoring_function& sf,const scoring_function& sf_ad4,const boost::optional<model>& ref);

    void initArchive(model& m,  model& m_ad4, const precalculate& p, const precalculate& prec_ad4, const precalculate& prec_vdw, const precalculate& prec_forcefield, const igrid& ig,
                     const precalculate& p_widened,  const vec& corner1, const vec& corner2,
                     rng& generator,const scoring_function& sf,const scoring_function& sf_ad4,const boost::optional<model>& ref);

    void initpbest(model& m,  model& m_ad4, const precalculate& p, const precalculate& prec_ad4, const precalculate& prec_vdw, const precalculate& prec_forcefield, const igrid& ig,
                   const precalculate& p_widened,  const vec& corner1, const vec& corner2,
                   rng& generator,const scoring_function& sf,const scoring_function& sf_ad4,const boost::optional<model>& ref);
    void calculate_objectives(model& m, const precalculate& p, const igrid& ig, vec v, output_type_SMPSO& indivisual, const scoring_function& sf, const scoring_function& sf_ad4, energy_cal& energy,const boost::optional<model>& ref,change &g);

    void computeSpeed(model& m,  const precalculate& p, const precalculate& prec_forcefield, const igrid& ig,
                      const precalculate& p_widened,  const vec& corner1, const vec& corner2,
                      rng& generator,const scoring_function& sf,const boost::optional<model>& ref);

    void computeNewPositions(const vec& corner1, const vec& corner2);

    void mutation(const vec& corner1, const vec& corner2,rng& generator);

    void updateArchive();

    void update_particle_memory(rng& generator);

    bool add(output_type_SMPSO x);

    void calculate_crowding_distance(std::vector<output_type_SMPSO>& pop);

    int indexWorst();
    fl constrictionCoefficient(double c1, double c2);

    fl velocityConstriction(fl v, int dim_index);
};
#endif //LSHADE_ADAM_FINAL_SMPSO_H
