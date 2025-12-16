/*
non_cache
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

#include "non_cache.h"
#include "curl.h"

non_cache::non_cache(const model& m, const grid_dims& gd_,  precalculate* p_, fl slope_) : sgrid(m, szv_grid_dims(gd_), p_->cutoff_sqr()), gd(gd_), p(p_), slope(slope_) {}
non_cache::non_cache(const grid_dims& gd_, precalculate* p_,Molecule_DLIGAND2 mol, fl slope_) : sgrid(szv_grid_dims(gd_), *p_, mol, p_->cutoff_sqr()), gd(gd_), p(p_), slope(slope_) {}

// fl non_cache::set_outTypeModea(output_type_MODEA out){
//    outputTypeModea = out;
//}

fl non_cache::eval      (const model& m, fl v) const { // clean up
	fl e = 0;
	const fl cutoff_sqr = p->cutoff_sqr();

	sz n = num_atom_types(p->atom_typing_used());
	VINA_FOR(i, m.num_movable_atoms()) {
		fl this_e = 0;
		fl out_of_bounds_penalty = 0;
		const atom& a = m.atoms[i];
		sz t1 = a.get(p->atom_typing_used());
		if(t1 >= n) continue;
		const vec& a_coords = m.coords[i];
		vec adjusted_a_coords; adjusted_a_coords = a_coords;
		VINA_FOR_IN(j, gd) {
			if(gd[j].n > 0) {
				if     (a_coords[j] < gd[j].begin) { adjusted_a_coords[j] = gd[j].begin; out_of_bounds_penalty += std::abs(a_coords[j] - gd[j].begin); }
				else if(a_coords[j] > gd[j].end  ) { adjusted_a_coords[j] = gd[j]  .end; out_of_bounds_penalty += std::abs(a_coords[j] - gd[j]  .end); }
			}
		}
		out_of_bounds_penalty *= slope;

		const szv& possibilities = sgrid.possibilities(adjusted_a_coords);

		VINA_FOR_IN(possibilities_j, possibilities) {
			const sz j = possibilities[possibilities_j];
			const atom& b = m.grid_atoms[j];
			sz t2 = b.get(p->atom_typing_used());
			if(t2 >= n) continue;
			vec r_ba; r_ba = adjusted_a_coords - b.coords; // FIXME why b-a and not a-b ?
			fl r2 = sqr(r_ba);
			if(r2 < cutoff_sqr) {
				sz type_pair_index = get_type_pair_index(p->atom_typing_used(), a, b);
				this_e +=  p->eval_fast(type_pair_index, r2);
			}
		}
		curl(this_e, v);
		e += this_e + out_of_bounds_penalty;
	}
	return e;
}

bool non_cache::within(const model& m, fl margin) const {
	VINA_FOR(i, m.num_movable_atoms()) {
		if(m.atoms[i].is_hydrogen()) continue;
		const vec& a_coords = m.coords[i];
		VINA_FOR_IN(j, gd)
			if(gd[j].n > 0)
				if(a_coords[j] < gd[j].begin - margin || a_coords[j] > gd[j].end + margin) 
					return false;
	}
	return true;
}

fl non_cache::eval_deriv(      model& m, fl v) const { // clean up
	fl e = 0;
    //首先确定截断距离，超过截断距离的原子不计算能量
	const fl cutoff_sqr = p->cutoff_sqr();
    //原子类型
	sz n = num_atom_types(p->atom_typing_used());
    int count = 0;
    //遍历所有配体原子
	VINA_FOR(i, m.num_movable_atoms()) {
        //用于累加能量
		fl this_e = 0;
        //梯度
		vec deriv(0, 0, 0);
        //梯度惩罚项
		vec out_of_bounds_deriv(0, 0, 0);
        //能量惩罚项
		fl out_of_bounds_penalty = 0;
        //得到当前配体原子
		const atom& a = m.atoms[i];
        //确定当前原子类型
		sz t1 = a.get(p->atom_typing_used());
        //如果当前原子类型大于规定原子类型，则将梯度置为0，且跳过当前原子
		if(t1 >= n) { m.minus_forces[i].assign(0); continue; }
        //确定当前原子坐标
		const vec& a_coords = m.coords[i];
        //临时原子坐标
		vec adjusted_a_coords; adjusted_a_coords = a_coords;
        //根据gd网格范围得到梯度惩罚项和能量惩罚项
		VINA_FOR_IN(j, gd) {
			if(gd[j].n > 0) {
				if     (a_coords[j] < gd[j].begin) { adjusted_a_coords[j] = gd[j].begin; out_of_bounds_deriv[j] = -1; out_of_bounds_penalty += std::abs(a_coords[j] - gd[j].begin); }
				else if(a_coords[j] > gd[j].end  ) { adjusted_a_coords[j] = gd[j]  .end; out_of_bounds_deriv[j] =  1; out_of_bounds_penalty += std::abs(a_coords[j] - gd[j]  .end); }
			}
		}
		out_of_bounds_penalty *= slope;
		out_of_bounds_deriv *= slope;

        //根据配体所有原子得到可能发生作用的蛋白原子索引
		const szv& possibilities = sgrid.possibilities(adjusted_a_coords);

        //根据蛋白索引遍历
		VINA_FOR_IN(possibilities_j, possibilities) {
			const sz j = possibilities[possibilities_j];
            //得到当前蛋白原子
			const atom& b = m.grid_atoms[j];
            //得到当前蛋白原子类型
			sz t2 = b.get(p->atom_typing_used());
//            printf("possibilities.size=%d\n",possibilities.size());
//            printf("possibilities_j=%d,j=%d,t2=%d\n",possibilities_j,j,t2);
            //若当前蛋白原子大于原子类型，则跳过
			if(t2 >= n) continue;
            //得到配体与蛋白原子之间的原子距离
			vec r_ba; r_ba = adjusted_a_coords - b.coords; // FIXME why b-a and not a-b ?
            //得到配体与蛋白原子之间的原子距离的平方
			fl r2 = sqr(r_ba);
            //若距离平方小于截断值的平方，则计算
			if(r2 < cutoff_sqr) {
                //得到蛋白与配体的原子对索引
				sz type_pair_index = get_type_pair_index(p->atom_typing_used(), a, b);
                //根据原子对索引，在预计算对象网格中用线性插值计算能量
				pr e_dor =  p->eval_deriv(type_pair_index, r2);
                //将计算得到的能量进行累加
				this_e += e_dor.first;
                //将计算得到的梯度进行累加
				deriv += e_dor.second * r_ba;
                count++;

            }
		}
		curl(this_e, deriv, v);
		m.minus_forces[i] = deriv + out_of_bounds_deriv;
		e += this_e + out_of_bounds_penalty;
	}
	return e;
}
fl non_cache::eval_deriv_ad4(model& m, fl v) const{
    return 0;
}
fl non_cache::eval_deriv_DLIGAND2(const output_type_MODEA& outputTypeModea, fl v) const { // clean up
//        struct timeval t11,t22;
//        double timeuse;
//        gettimeofday(&t11,NULL);
    fl e = 0;
    //得到截断距离的平方
    const fl cutoff_sqr = p->cutoff_sqr();
//    for(int i = 0; i <  outputTypeModea.lig_num; i++){
//        p->minus_forces[i].assign(0);
//    }
    //得到截断距离的平方
    p->minus_forces.resize(outputTypeModea.lig_num);
    sz n = matype_pro;
    int count = 0;
    //从配体原子中的初始索引值开始遍历
    for(int i = outputTypeModea.natm - outputTypeModea.lig_num; i < outputTypeModea.natm; i++) {
        fl this_e = 0;
        vec deriv(0, 0, 0);
        vec out_of_bounds_deriv(0, 0, 0);
        fl out_of_bounds_penalty = 0;
        const Atom_DLIGAND2& a = outputTypeModea.atoms[i];
        sz t1 = a.type;
        if(t1 < n) { p->minus_forces[i - outputTypeModea.natm + outputTypeModea.lig_num].assign(0); continue; }
        vec a_coords;
        VINA_FOR_IN(j, a_coords)
            a_coords[j] = outputTypeModea.getx(i)[j];
        vec adjusted_a_coords; adjusted_a_coords = a_coords;
        VINA_FOR_IN(j, gd) {
            if(gd[j].n > 0) {
                if     (a_coords[j] < gd[j].begin) { adjusted_a_coords[j] = gd[j].begin; out_of_bounds_deriv[j] = -1; out_of_bounds_penalty += std::abs(a_coords[j] - gd[j].begin); }
                else if(a_coords[j] > gd[j].end  ) { adjusted_a_coords[j] = gd[j]  .end; out_of_bounds_deriv[j] =  1; out_of_bounds_penalty += std::abs(a_coords[j] - gd[j]  .end); }
            }
        }
        out_of_bounds_penalty *= slope;
        out_of_bounds_deriv *= slope;

        const szv& possibilities = sgrid.possibilities(adjusted_a_coords);

        VINA_FOR_IN(possibilities_j, possibilities) {
            const sz j = possibilities[possibilities_j];
            const Atom_DLIGAND2& b = outputTypeModea.atoms[j];
            vec b_coords;
            VINA_FOR_IN(k, b_coords)
                b_coords[k] = outputTypeModea.getx(j)[k];
            sz t2 = b.type;
            if(t2 >=  n) continue;
            vec r_ba; r_ba = adjusted_a_coords - b_coords; // FIXME why b-a and not a-b ?
            fl r2 = sqr(r_ba);
            if(r2 < cutoff_sqr) {

//                sz type_pair_index = get_type_pair_index(p.atom_typing_used(), a, b);
                pr e_dor =  p->eval_deriv_DLIGAND2(t2, t1, r2);
                this_e += e_dor.first;
                deriv += e_dor.second * r_ba;
                count++;

            }
        }

        curl(this_e, deriv, v);
        p->minus_forces[i - outputTypeModea.natm + outputTypeModea.lig_num] = deriv + out_of_bounds_deriv;
        e += this_e + out_of_bounds_penalty;
    }
//    gettimeofday(&t22,NULL);
//    timeuse = (t22.tv_sec - t11.tv_sec) + (double)(t22.tv_usec - t11.tv_usec)/1000000.0;
//
//    cout<<"DLIGAND_time = "<<timeuse<<endl;  //���ʱ�䣨��λ����
    return e;
}


fl non_cache::eval_inter_ligand_receptor_deriv(model& m, fl v) const {
    fl e = 0;
    const fl cutoff_sqr = p->cutoff_sqr();

    sz n = num_atom_types(p->atom_typing_used());

    VINA_FOR(i, m.num_movable_atoms()) {
        if (!m.is_atom_in_ligand(i)) continue; // we only want flex-rigid interactions

        fl this_e = 0;
        vec deriv(0, 0, 0);
        vec out_of_bounds_deriv(0, 0, 0);
        fl out_of_bounds_penalty = 0;
        const atom &a = m.atoms[i];
        sz t1 = a.get(p->atom_typing_used());
        if (t1 >= n) {
            m.minus_forces[i].assign(0);
            continue;
        }

        const vec &a_coords = m.coords[i];
        vec adjusted_a_coords;
        adjusted_a_coords = a_coords;

        VINA_FOR_IN(j, gd) {
            if (gd[j].n > 0) {
                if (a_coords[j] < gd[j].begin) {
                    adjusted_a_coords[j] = gd[j].begin;
                    out_of_bounds_deriv[j] = -1;
                    out_of_bounds_penalty += std::abs(a_coords[j] - gd[j].begin);
                }
                else if (a_coords[j] > gd[j].end) {
                    adjusted_a_coords[j] = gd[j].end;
                    out_of_bounds_deriv[j] = 1;
                    out_of_bounds_penalty += std::abs(a_coords[j] - gd[j].end);
                }
            }
        }

        out_of_bounds_penalty *= slope;
        out_of_bounds_deriv *= slope;

        const szv &possibilities = sgrid.possibilities(adjusted_a_coords);

        VINA_FOR_IN(possibilities_j, possibilities) {
            const sz j = possibilities[possibilities_j];
            const atom &b = m.grid_atoms[j];
            sz t2 = b.get(p->atom_typing_used());
            if (t2 >= n) continue;
            vec r_ba;
            r_ba = adjusted_a_coords - b.coords; // FIXME why b-a and not a-b ?
            fl r2 = sqr(r_ba);
            if (r2 < cutoff_sqr) {
                sz type_pair_index = get_type_pair_index(p->atom_typing_used(), a, b);
                pr e_dor = p->eval_deriv(type_pair_index, r2);
                this_e += e_dor.first;
                deriv += e_dor.second * r_ba;
            }
        }
        curl(this_e, deriv, v);
        m.minus_forces[i] = deriv + out_of_bounds_deriv;
        e += this_e + out_of_bounds_penalty;
    }
    return e;
}

fl non_cache::eval_intra_deriv(model& m, fl v) const {
    fl e = 0;
    const fl cutoff_sqr = p->cutoff_sqr();

    sz n = num_atom_types(p->atom_typing_used());

    VINA_FOR(i, m.num_movable_atoms()) {
        if (m.is_atom_in_ligand(i)) continue; // we only want flex-rigid interactions

        fl this_e = 0;
        vec deriv(0, 0, 0);
        vec out_of_bounds_deriv(0, 0, 0);
        fl out_of_bounds_penalty = 0;
        const atom &a = m.atoms[i];
        sz t1 = a.get(p->atom_typing_used());
        if (t1 >= n) {
            m.minus_forces[i].assign(0);
            continue;
        }

        const vec &a_coords = m.coords[i];
        vec adjusted_a_coords;
        adjusted_a_coords = a_coords;

        VINA_FOR_IN(j, gd) {
            if (gd[j].n > 0) {
                if (a_coords[j] < gd[j].begin) {
                    adjusted_a_coords[j] = gd[j].begin;
                    out_of_bounds_deriv[j] = -1;
                    out_of_bounds_penalty += std::abs(a_coords[j] - gd[j].begin);
                }
                else if (a_coords[j] > gd[j].end) {
                    adjusted_a_coords[j] = gd[j].end;
                    out_of_bounds_deriv[j] = 1;
                    out_of_bounds_penalty += std::abs(a_coords[j] - gd[j].end);
                }
            }
        }

        out_of_bounds_penalty *= slope;
        out_of_bounds_deriv *= slope;

        const szv &possibilities = sgrid.possibilities(adjusted_a_coords);

        VINA_FOR_IN(possibilities_j, possibilities) {
            const sz j = possibilities[possibilities_j];
            const atom &b = m.grid_atoms[j];
            sz t2 = b.get(p->atom_typing_used());
            if (t2 >= n) continue;
            vec r_ba;
            r_ba = adjusted_a_coords - b.coords; // FIXME why b-a and not a-b ?
            fl r2 = sqr(r_ba);
            if (r2 < cutoff_sqr) {
                sz type_pair_index = get_type_pair_index(p->atom_typing_used(), a, b);
                pr e_dor = p->eval_deriv(type_pair_index, r2);
                this_e += e_dor.first;
                deriv += e_dor.second * r_ba;
            }
        }
        curl(this_e, deriv, v);
        m.minus_forces[i] = deriv + out_of_bounds_deriv;
        e += this_e + out_of_bounds_penalty;
    }
    return e;
}


