/*
coords.h
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

#ifndef VINA_COORDS_H
#define VINA_COORDS_H

#include "conf.h"
#include "atom.h" // for atomv

fl rmsd_upper_bound(const vecv& a, const vecv& b);
std::pair<sz, fl> find_closest(const vecv& a, const output_container& b);
std::pair<sz, fl> find_closest(const vecv& a, const output_container_nsga2 & b);
void add_to_output_container(output_container& out, const output_type& t, fl min_rmsd, sz max_size);
void add_to_output_container(output_container_nsga2 & out, const output_type_nsga2& t, fl min_rmsd, sz max_size);
void add_to_output_container(output_container_MODEA & out, const output_type_MODEA& t, fl min_rmsd, sz max_size);
void add_to_output_container_(output_container_MODEA & out, const output_type_MODEA& t, fl min_rmsd, sz max_size);
void add_to_output_container(output_container_SMPSO & out, const output_type_SMPSO& t, fl min_rmsd, sz max_size);
void add_to_output_container(output_container_MOEAD & out, const output_type_MOEAD& t, fl min_rmsd, sz max_size);
void add_to_output_container(output_container_MPSOD & out, const output_type_MPSOD& t, fl min_rmsd, sz max_size);

#endif
