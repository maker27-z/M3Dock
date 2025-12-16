/*
parallel_mc.h
   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#ifndef VINA_PARALLEL_MC_H
#define VINA_PARALLEL_MC_H

#include "monte_carlo.h"
#include "lshade.h"
#include "nsga2.h"
#include "MODEA.h"
#include "SMPSO.h"
#include "MOEAD.h"
#include "MPSOD.h"
struct parallel_mc {
    lshade ls;
    nsga2 ns;
    MODEA modea;
//    SMPSO smpso;
//    MOEAD moead;
//    MPSOD mpsod;
    sz num_tasks;
    sz num_threads;
    bool display_progress;
    parallel_mc() : num_tasks(8), num_threads(1), display_progress(true) {}
    void operator()(const model& m, const model& m_ad4, output_container_nsga2 & out, const precalculate& p, const precalculate& prec_ad4, const precalculate& prec_vdw, const precalculate& prec_forcefield,  const precalculate& prec_DLIGAND2, const igrid& ig, const igrid& DLIGAND2_grid, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, rng& generator,const scoring_function& sf, const scoring_function& sf_ad4, const boost::optional<model>& ref) const;
    void operator()(const model& m, const model& m_ad4, output_container_MODEA & out, const precalculate& p, const precalculate& prec_ad4, const precalculate& prec_vdw, const precalculate& prec_forcefield,  const precalculate& prec_DLIGAND2, const igrid& ig, const igrid& DLIGAND2_grid, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, rng& generator,const scoring_function& sf, const scoring_function& sf_ad4, const boost::optional<model>& ref) const;
    void operator()(const model& m, const model& m_ad4, output_container_SMPSO & out, const precalculate& p, const precalculate& prec_ad4, const precalculate& prec_vdw, const precalculate& prec_forcefield,  const precalculate& prec_DLIGAND2, const igrid& ig, const igrid& DLIGAND2_grid, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, rng& generator,const scoring_function& sf, const scoring_function& sf_ad4, const boost::optional<model>& ref) const;
    void operator()(const model& m, const model& m_ad4, output_container_MOEAD & out, const precalculate& p, const precalculate& prec_ad4, const precalculate& prec_vdw, const precalculate& prec_forcefield,  const precalculate& prec_DLIGAND2, const igrid& ig, const igrid& DLIGAND2_grid, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, rng& generator,const scoring_function& sf, const scoring_function& sf_ad4, const boost::optional<model>& ref) const;
    void operator()(const model& m, const model& m_ad4, output_container_MPSOD & out, const precalculate& p, const precalculate& prec_ad4, const precalculate& prec_vdw, const precalculate& prec_forcefield,  const precalculate& prec_DLIGAND2, const igrid& ig, const igrid& DLIGAND2_grid, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, rng& generator,const scoring_function& sf, const scoring_function& sf_ad4, const boost::optional<model>& ref) const;
};

#endif
