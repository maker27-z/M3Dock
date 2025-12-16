//
// Created by 91686 on 2023/10/16.
//

#include <sys/time.h>
#include "MODEA.h"
#include"adam_cal.h"
#include "matplotlibcpp.h"
#include "tee.h"
//#include <boost/python.hpp>
const fl beta1 = 0.5;

const fl beta2 = 0.999;
const fl sigma = 0.00000001;
const fl lr = 0.01;

void MyNormalize(vec& angles) {
    while (angles[0] > 2 * pi) {
        angles[0] -= 2 * pi;
    }
    while (angles[0] < 0) {
        angles[0] += 2 * pi;
    }
    while (angles[1] > 2 * pi) {
        angles[1] -= 2 * pi;
    }
    while (angles[1] < -2 * pi) {
        angles[1] += 2 * pi;
    }
    if (angles[1] > pi && angles[1] < 2 * pi) {
        angles[1] = 2 * pi - angles[1];
    }
    if (angles[1] > -2 * pi && angles[1] < -pi) {
        angles[1] = 2 * pi + angles[1];
    }
    if (angles[1] > -pi && angles[1] < 0) {
        angles[1] = -angles[1];
    }
    while (angles[2] > 2 * pi) {
        angles[2] -= 2 * pi;
    }
    while (angles[2] < 0) {
        angles[2] += 2 * pi;
    }

    if(!(angles[0] > 0 && angles[0] < 2 * pi)){
        std::cout<<"angles[0]: "<<angles[0]<<endl;
    }
    if(!(angles[1] > 0 && angles[1] < pi)){
        std::cout<< std::fixed << std::setprecision(16) << "angles[1]: "<<angles[1]<<endl;
    }
    if(!(angles[2] > 0 && angles[2] < 2 * pi)){
        std::cout<<"angles[2]: "<<angles[2]<<endl;
    }


}

void quaternion_to_3angles(const output_type_MODEA& x, vec& three_angle) {
    double theta = 2 * acos(x.c.ligands[0].rigid.orientation.R_component_1());
    if (theta < 0.00000001) {
        three_angle[0] = 0;
        three_angle[1] = 0;
        three_angle[2] = 0;
    }
    else {
        vec quaternion_v;//v_x,v_y,v_z
        quaternion_v[0] = x.c.ligands[0].rigid.orientation.R_component_2() / sin(theta / 2); //v_x
        quaternion_v[1] = x.c.ligands[0].rigid.orientation.R_component_3() / sin(theta / 2); //v_y
        quaternion_v[2] = x.c.ligands[0].rigid.orientation.R_component_4() / sin(theta / 2); //v_z
        fl alpha = acos(quaternion_v[2]);
        fl beta = 0;
        fl sin_beta = quaternion_v[1] / sin(alpha);
        fl cos_bata = quaternion_v[0] / sin(alpha);
        if (sin_beta > 0)
            beta = acos(cos_bata);
        else if (sin_beta < 0)
            beta = 2 * pi - acos(cos_bata);
        three_angle[0] = theta;
        three_angle[1] = alpha;
        three_angle[2] = beta;
    }
    assert(three_angle[0] >= 0 && three_angle[0] <= 2 * pi);
    assert(three_angle[1] >= 0 && three_angle[1] <= pi);
    assert(three_angle[2] >= 0 && three_angle[2] <= 2 * pi);
}

void threeangles_to_quaternion(const vec& three_angle, std::vector<double>& quater) {
    if (fabs(three_angle[0] - 0.0) < 1e-6) {
        quater[0] = 0;
        quater[1] = 0;
        quater[2] = 0;
        quater[3] = 0;
    }
    else {
        quater[0] = three_angle[0] / 2;//    ֮theta
        quater[1] = sin(three_angle[1]) * cos(three_angle[2]);
        quater[2] = sin(three_angle[1]) * sin(three_angle[2]);
        quater[3] = cos(three_angle[1]);
    }
}


void MODEA_aux::mutate_cross(model& m, energy_cal& energy, rng& generator, const vec& corner1, const vec& corner2) {
	conf_size s = m.get_size();
	change g(s);
	fl best_e = max_fl;
	best_index = -1;
	int po_size_cur = nowpopulation.size();
        std::vector<int> index;
    //为每个个体执行变异交叉操作
	for (int i = 0; i < po_size_cur; i++) {

//		int pbest = random_int(0, p * po_size_cur - 1, generator);//          Ⱥ         pbest
        //先生成三个随机个体
		int index1, index2, index3;//   ɲ ͬ  x_r1  x_r2     ±
		index1 = random_int(0, po_size_cur - 1, generator);

		index2 = random_int(0, po_size_cur - 1, generator);
		bool isA = false;
		while (index1 == index2) {//r2       Ⱥ   ⲿ 浵A Ĳ        ѡȡ
			index2 = random_int(0, po_size_cur - 1, generator);
		}
		if (index2 >= po_size_cur) {//      Ⱥ  ģ Ĳ  ־    ⲿ 浵A е Ԫ
			index2 -= po_size_cur;
			isA = true;
		}
        index3 = random_int(0, po_size_cur - 1, generator);
        while (index1 == index3 || index2 == index3) {//  ǰ   岻 ܺ͸ø         ͬ
            index3 = random_int(0, po_size_cur - 1, generator);
        }
        index.push_back(index1);
        index.push_back(index2);
        index.push_back(index3);
        std::vector<output_type_MODEA> tmp_random;
        tmp_random.push_back(nowpopulation[index1]);
        tmp_random.push_back(nowpopulation[index2]);
        tmp_random.push_back(nowpopulation[index3]);
        int flag = 0;
        float maxe = -std::numeric_limits<float>::infinity();
        fast_nondominated_sort(tmp_random);
        calculate_crowding_distance(tmp_random);
        for(int m = 0; m < tmp_random.size(); m++){
            if(tmp_random[m].rank==1){
                if(tmp_random[m].crowding_distance>maxe){
                    flag = m;
                    maxe = tmp_random[m].crowding_distance;
                }
            }
        }

//        std::sort(index.begin(), index.end());

        index1 = index[flag];
        index.erase(index.begin() + flag);
        index2 = index[0];
        index3 = index[1];
		std::vector<double> quater;
		quater.resize(4);
		vec angle_de;//û   õ
		bool is_mutate_qt = false;
		int dim_j = random_int(0, dim - 1, generator);

        //先遍历个体的每个维度，如果符合交叉条件，则进行交叉
		for (int j = 0; j < dim; j++) {
			double p_cr = random_fl(0, 1.0, generator);
			if (p_cr < CR || j == dim_j) { //          Ҫ    ,      һ  ά    ԭ

					if (j < 3) {

						fl t = 1;
						mutatepopulation[i].c.ligands[0].rigid.position[j] = nowpopulation[index1](j) +
							+ F * (nowpopulation[index2](j) - nowpopulation[index3](j));
						if (mutatepopulation[i].c.ligands[0].rigid.position[j] < corner1[j]) {
							mutatepopulation[i].c.ligands[0].rigid.position[j] = corner1[j];
						}
						else if (mutatepopulation[i].c.ligands[0].rigid.position[j] > corner2[j]) {
							mutatepopulation[i].c.ligands[0].rigid.position[j] = corner2[j];
						}
					}
					else if (j >= 3 && j < 6) {
						mutatepopulation[i].rotor_angle[j - 3] = nowpopulation[index1].rotor_angle[j - 3]
								+ F * (nowpopulation[index2].rotor_angle[j - 3] - nowpopulation[index3].rotor_angle[j - 3]);

						//MyNormalize(angle_de);
						//threeangles_to_quaternion(angle_de, quater);
					}
					else {
						mutatepopulation[i].c.ligands[0].torsions[j - 6] = nowpopulation[index1](j)
							+ F * (nowpopulation[index2](j) - nowpopulation[index3](j));
						normalize_angle(mutatepopulation[i].c.ligands[0].torsions[j - 6]);
					}

			}
			else {//     㽻    ʣ ֱ Ӽ̳ ԭ    Ļ
				if (j < 3) {
					mutatepopulation[i].c.ligands[0].rigid.position[j] = nowpopulation[i](j);
				}
				else if (j >= 3 && j < 6) {
					//is_mutate_qt = false;
					mutatepopulation[i].rotor_angle[j - 3] = nowpopulation[i].rotor_angle[j - 3];
				}
				else {
					mutatepopulation[i].c.ligands[0].torsions[j - 6] = nowpopulation[i](j);
					normalize_angle(mutatepopulation[i].c.ligands[0].torsions[j - 6]);
				}
			}
		}
        //将三个角度转换为四元数orientation
		MyNormalize(mutatepopulation[i].rotor_angle);
		threeangles_to_quaternion(mutatepopulation[i].rotor_angle, quater);
		fl cos_theta = cos(quater[0]);
		fl sin_theta = sin(quater[0]);
		qt q(cos_theta, quater[1] * sin_theta, quater[2] * sin_theta, quater[3] * sin_theta);
		mutatepopulation[i].c.ligands[0].rigid.orientation = q;
//		mutatepopulation[i].e = energy(mutatepopulation[i].c, g);
//		if (mutatepopulation[i].e < best_e) {
//			best_e = mutatepopulation[i].e;
//			best_index = i;
//		}
	}
}

void MODEA_aux::select(model& m, energy_cal& energy, rng& generator) {

//	std::vector<double> fitness_improvement;
//	int count_zero = 0;
//	int po_size_cur = nowpopulation.size();
//	for (int i = 0; i < po_size_cur; i++) {
//		if (nowpopulation[i].e <= mutatepopulation[i].e) {
//			continue;
//		}
//		else {//
//			A.push_back(nowpopulation[i]);// ⲿ 鵵A    ʧ ܸ
//			SCR.push_back(nowpopulation[i].CR);
//			SF.push_back(F);
//			if (SCR.back() == 0) {
//				count_zero++;
//			}
//			fitness_improvement.push_back(fabs(nowpopulation[i].e - mutatepopulation[i].e));
//			nowpopulation[i] = mutatepopulation[i];
//		}
//	}

}

void MODEA_aux::calculate_objectives(model& m, const precalculate& p, const igrid& ig, vec v, output_type_MODEA& indivisual, const scoring_function& sf, const scoring_function& sf_ad4, energy_cal& energy,const boost::optional<model>& ref,change &g, change &vdw_g, change &vina_g, change &DLIGAND2_g){
    conf_size s = m.get_size();
//    change DLIGAND2_g(s);
    fl intramolecular_energy;
    for(int k = 0; k<indivisual.objectives.size();k++){
        if(k==0){
//                tmp.objectives[k] = energy.energy_cal_force(tmp.c);
            //首先计算vina分子内能量
            intramolecular_energy = energy.energy_cal_intra(indivisual.c,g);

//            indivisual.objectives[k] = energy.energy_cal_vina(sf, indivisual.c, intramolecular_energy,g);

//            intramolecular_energy = energy.AD4_energy_cal_intra(indivisual.c);
//            indivisual.objectives[k] = energy.energy_cal_forcefield(indivisual.c);
//            indivisual.objectives[k] = intramolecular_energy;
//            indivisual.objectives[k] = energy.energy_cal_autodock4(sf_ad4,indivisual.c);
//            indivisual.objectives[k] = energy.energy_cal_intra(indivisual.c,g);
//            indivisual.objectives[k] = energy.energy_cal_vdw(sf_ad4, indivisual.c);
//            indivisual.objectives[k] = energy.energy_cal_vdw(sf_ad4, indivisual.c,vdw_g);
            //根据vina分子内能量计算vina评分函数
            indivisual.objectives[k] = energy.energy_cal_vina_deriv(sf, indivisual.c, intramolecular_energy,vina_g);
//            indivisual.binding_energy = energy.energy_cal_vina_binding(sf, indivisual.c, ttt,intramolecular_energy,vina_g);

//              printf("indivisual.binding_energy:%f,indivisual.vina_energy:%f\n",indivisual.binding_energy,ttt);
//            printf("obj1vina:%f\n",ttt);
//            indivisual.objectives[k] = energy.energy_cal_autodock4(sf_ad4, indivisual.c,vina_g);
//            int b = 0;
//            float bind_score = energy.energy_cal_xscore(&indivisual);
//            indivisual.objectives[k] = bind_score*(-1.364);
        }else if(k==1){
//                tmp.objectives[k] = energy.energy_cal_hydrophobic(tmp.c);
//            indivisual.objectives[k] = energy.AD4_energy_cal_inter(sf, indivisual.c, intramolecular_energy);
//            float bind_score = energy.energy_cal_xscore(&indivisual);
//            indivisual.objectives[k] = bind_score*(-1.364);
//            indivisual.objectives[k] = energy.energy_cal_DLIGAND2(indivisual);
            //计算DLIGAND2能量评分函数
            indivisual.objectives[k] = energy.energy_cal_DLIGAND2(indivisual,DLIGAND2_g);
//            vec d_position = DLIGAND2_g.ligands[0].rigid.position;
//            vec d_orientation = DLIGAND2_g.ligands[0].rigid.orientation;
//            flv d_torsions = DLIGAND2_g.ligands[0].torsions;
//            indivisual.objectives[k] = energy.energy_cal_vdw(sf_ad4, indivisual.c,vdw_g);
//            int a = 0;
        }else if(k==2){
//                        indivisual.objectives[k] = energy.energy_cal_autodock4(sf_ad4,indivisual.c);
//                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
//            indivisual.objectives[k] = energy.energy_cal_rmsd(ref, indivisual.c);
//            indivisual.objectives[k] = energy.energy_cal_DLIGAND2(indivisual,DLIGAND2_g);
//            indivisual.objectives[k] = energy.energy_cal_vdw(sf_ad4, indivisual.c);
            //计算vdw能量评分函数
            indivisual.objectives[k] = energy.energy_cal_vdw(sf_ad4, indivisual.c,vdw_g);
//            indivisual.objectives[k] = energy.energy_cal_autodock4(sf_ad4, indivisual.c,vdw_g);
//int a = 0;
//            float a = indivisual.objectives[k];
//            float b = a;
        } else if (k == 3) {
//                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
//                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
//                tmp.objectives[k] = energy.energy_cal_Smog2016(tmp);
            indivisual.objectives[k] = energy.energy_cal_vdw(sf_ad4, indivisual.c,vdw_g);

        }
    }

}

std::unordered_map<int, std::vector<output_type_MODEA>> MODEA_aux::fast_nondominated_sort(std::vector<output_type_MODEA>& pop) {
    std::unordered_map<int, std::vector<output_type_MODEA>> F;
    for (auto& p : pop) {
        //p.Slaves.clear();
        p.Slaves.clear();

        p.dominated_count = 0;
        for (auto& q : pop) {
            if (p < q) {
                //      ǰ    p֧  q,  ô    ֧  ĸ   q    Sp
                p.Slaves.push_back(&q);
            }
            else if (q < p) {
                //       q֧  p,  ô  np  ֵ  һ
                p.dominated_count += 1;
            }
        }
        //ѡ     ʼ  Ⱥ  rank1 ĸ  壬   ҽ rank1 еĸ      F[1]
        if (p.dominated_count == 0) {
            p.rank = 1;
            F[1].push_back(p);
        }

    }

    int i = 1;
    std::vector<output_type_MODEA> Q;
    while (!F[i].empty()) {
        Q.clear();
        //ѭ    ǰF[i]           Ž⼯    rankn еĽ⼯
        for (auto& p : F[i]) {
            //ѭ    ǰ         Ž⼯ е Sp   ϣ     ǰ     ֧  Ľ
            for (auto& q : p.Slaves) {
                //  Sp е ÿ       np  -1    ֧ 䵱ǰ    Sp  ĸ       -1
                q->dominated_count = q->dominated_count - 1;
                //  ж Sp     np Ƿ Ϊ0,   Ϊ0,   ʾ  ǰSp еĸø    岻   κ           ֧
                //     ǰSp  npΪ0 Ľ    rank+1 У       Q У   Ϊ  ֧
                if (q->dominated_count == 0) {
                    q->rank = i + 1;
                    Q.push_back(*q);
                }
            }
        }
        i = i + 1;
        F[i] = Q;
    }
    return F;
}

void MODEA_aux::calculate_crowding_distance(std::vector<output_type_MODEA>& rank_i_vector) {
    int length = rank_i_vector.size();
    for (int i = 0; i < length; i++) {
        rank_i_vector[i].crowding_distance = 0.0;
    }
    //    1          Ϊ βΣ         ʱҲ   Կ  ɽ    ݵ ʵ ΰ󶨸              ں      ڶ            һ в      п   Ӱ 쵽       ݵ ʵ Ρ        ϣ       ں         ֻ    ,   Ե    Ǽ         ϣ        ֻ   ľͱ     const
    //Ϊʲô  ֱ  ֵ     أ  ȷʵ     ǵ            ʱֵ   ݾ     һ     ⣬ Ǿ      ܿ  ܻ    Ӱ 졣    ֪  ֵ    ʵ ʾ           һ ݸ     ʹ ã   ô    һЩ   ӵ  ࣬       string     ÿһ ο          ĺܶ  ʱ 䣬  ôͨ     ô  ξͺ  б Ҫ ˡ
    // ܵ   ˵  Ϊ           󴫲 ʱ     ܣ     Ҫ     ã   Ϊ            ϣ    ֻ            const
    for (int i = 0; i < rank_i_vector[0].objectives.size(); i++) {
        std::sort(rank_i_vector.begin(), rank_i_vector.end(), [i](const output_type_MODEA& a, const output_type_MODEA& b) {
            return a.objectives[i] < b.objectives[i];
        });
        rank_i_vector[0].crowding_distance = rank_i_vector[length - 1].crowding_distance = std::numeric_limits<fl>::infinity();
        fl fmax = rank_i_vector[length - 1].objectives[i];
        fl fmin = rank_i_vector[0].objectives[i];
        if (fmax != fmin) {
            for (int j = 1; j < length - 1; j++) {
                rank_i_vector[j].crowding_distance += (rank_i_vector[j + 1].objectives[i] - rank_i_vector[j - 1].objectives[i]) / (fmax - fmin);
            }
        }

    }

}


bool sort_comparator(const output_type_MODEA& x,const output_type_MODEA& y){
    if(x.rank != y.rank){
        if(x.rank < y.rank){
            return true;
        }else{
            return false;
        }
    } else if(x.crowding_distance != y.crowding_distance){
        if(x.crowding_distance > y.crowding_distance){
            return true;
        }else{
            return false;
        }
    } else{
        return true;
    }
}

//ӵ   Ⱦ   Ƚ             ӵ       ͬһrank ĸ
bool crowding_distance_comparator(const output_type_MODEA& x, const output_type_MODEA& y) {
    return x.crowding_distance > y.crowding_distance;
}
void MODEA_aux::nd_cd_sort(std::vector<output_type_MODEA> &pop) {
    std::unordered_map<int, std::vector<output_type_MODEA>> F = fast_nondominated_sort(pop);

    std::vector<output_type_MODEA> P_t_1;
    int i = 1;
    //    µĸ     Ⱥ  С    N    ô    ѭ
    while (P_t_1.size() + F[i].size() < nowpopulation.size()) {
        calculate_crowding_distance(F[i]);
        P_t_1.insert(P_t_1.end(), F[i].begin(), F[i].end());
        i = i + 1;
    }
    calculate_crowding_distance(F[i]);
    std::sort(F[i].begin(), F[i].end(), crowding_distance_comparator);
    int left = nowpopulation.size() - P_t_1.size();
    //P_t_1.insert(P_t_1.end(), F[i].begin(), F[i].begin() + aux.P.size() - P_t_1.size());
    for (int a = 0; a < left; a++) {
        P_t_1.push_back(F[i][a]);

    }
    nowpopulation = P_t_1;
    for (auto &s: nowpopulation) {
//        s.e = m_ad4.eval(p,ig,v,s.c);
//        s.e = s.objectives[1]/(1+0.05846*(m_ad4.ligand_degrees_of_freedom(0)));
        s.e = s.objectives[0] + s.objectives[1] + s.objectives[2];
    }



}
//
//float get_objective_distance(output_type_MODEA outt){
//        float d,tmpx,tmpy,tmpz;
//        extern output_type_MODEA std_tmp;
//        tmpx=(outt.objectives[0]-std_tmp.objectives[0])*(outt.objectives[0]-std_tmp.objectives[0]);
//        tmpy=(outt.objectives[1]-std_tmp.objectives[1])*(outt.objectives[1]-std_tmp.objectives[1]);
//        tmpz=(outt.objectives[2]-std_tmp.objectives[2])*(outt.objectives[2]-std_tmp.objectives[2]);
//
//        d=sqrt(tmpx+tmpy+tmpz);
//
//        return d;
//}
bool areAlmostEqual(float a, float b, float tolerance) {
    return std::fabs(a - b) < tolerance;
}

fl dotProduct(const std::vector<fl>& gradient1, const std::vector<fl>& gradient2){
    fl dot = 0.0;
    VINA_FOR(i,gradient1.size()){
        dot += gradient1[i] * gradient2[i];
    }
    return dot;
}

std::vector<fl> subtract(const std::vector<fl>& gradient1, const std::vector<fl>& gradient2){
    std::vector<fl> result(gradient1.size());
    VINA_FOR(i, gradient1.size()){
        result[i] = gradient1[i] - gradient2[i];
    }
    return result;
}

fl calculate_alpha(const std::vector<fl>& gradient1, const std::vector<fl>& gradient2){
    std::vector<fl> diff = subtract(gradient1,gradient2);
    fl numerator = dotProduct(diff, gradient2);
    fl denominator = dotProduct(diff, diff);
    if(denominator == 0){
        std::cout << "Error: Denominator is zero." << std::endl;
        return 0;
    }
    fl alpha = -numerator / denominator;

    alpha = std::max(0.0,std::min(1.0,alpha));

    return alpha;
}

void increment_position(conf& c, change& g, fl beta1, fl beta2, vec& V_d_position, vec& S_d_position, fl lr, fl sigma, unsigned t) {
    vec d_position = g.ligands[0].rigid.position;

    V_d_position = beta1 * V_d_position + (1 - beta1) * d_position;
    S_d_position = beta2 * S_d_position + (1 - beta2) * elementwise_product(d_position, d_position);

    vec V_d_position_corr = V_d_position / (1 - pow(beta1, t));
    vec S_d_position_corr = S_d_position / (1 - pow(beta2, t));
    vec_sqrt(S_d_position_corr);
    S_d_position_corr += sigma;

    c.ligands[0].rigid.position = c.ligands[0].rigid.position - lr * vec_division(V_d_position_corr, S_d_position_corr);
}


void increment_orientation(conf& c, change& g, fl beta1, fl beta2, vec& V_d_orientation, vec& S_d_orientation, fl lr, fl sigma, unsigned t) {
    vec d_orientation = g.ligands[0].rigid.orientation;

    V_d_orientation = beta1 * V_d_orientation + (1 - beta1) * d_orientation;
    S_d_orientation = beta2 * S_d_orientation + (1 - beta2) * elementwise_product(d_orientation, d_orientation);

    vec V_d_orientation_corr = V_d_orientation / (1 - pow(beta1, t));
    vec S_d_orientation_corr = S_d_orientation / (1 - pow(beta2, t));
    vec_sqrt(S_d_orientation_corr);
    S_d_orientation_corr += sigma;

    vec rotation = -lr * vec_division(V_d_orientation_corr, S_d_orientation_corr);
    quaternion_increment(c.ligands[0].rigid.orientation, rotation);
}

void increment_torsions(conf& c, change& g, fl beta1, fl beta2, flv& V_d_torsions, flv& S_d_torsions, fl lr, fl sigma, unsigned t) {
    flv d_torsions = g.ligands[0].torsions;

    V_d_torsions = torsions_add_torsions(s_mul_torsions(beta1, V_d_torsions), s_mul_torsions(1 - beta1, d_torsions));
    S_d_torsions = torsions_add_torsions(s_mul_torsions(beta2, S_d_torsions), s_mul_torsions(1 - beta2, torsions_sqr(d_torsions)));

    flv V_d_torsions_corr = torsions_div_s(V_d_torsions, 1 - pow(beta1, t));
    flv S_d_torsions_corr = torsions_div_s(S_d_torsions, 1 - pow(beta2, t));
    torsions_sqrt(S_d_torsions_corr);
    torsions_add_s(S_d_torsions_corr, sigma);

    flv torsions_ = s_mul_torsions(lr, torsions_div_torsions(V_d_torsions_corr, S_d_torsions_corr));

    for (sz i = 0; i < torsions_.size(); i++) {
        c.ligands[0].torsions[i] -= normalized_angle(torsions_[i]);
        normalize_angle(c.ligands[0].torsions[i]);
    }
}
//std::vector<float> convertPyObjectToVector(boost::python::object py_obj) {
//    std::vector<float> cpp_result;
//    boost::python::ssize_t len = boost::python::len(py_obj);
//    for (boost::python::ssize_t i = 0; i < len; ++i) {
//        cpp_result.push_back(boost::python::extract<float>(py_obj[i]));
//    }
//    return cpp_result;
//}
//std::vector<float> convert_to_vector(boost::python::object obj) {
//    // 确保 obj 是一个可迭代对象
//    boost::python::stl_input_iterator<float> begin(obj), end;
//    return std::vector<float>(begin, end);
//}
//boost::python::list vectorToList(const std::vector<std::vector<float>>& vec) {
//    boost::python::list list;
//    for (const auto& subVec : vec) {
//        boost::python::list sublist;
//        for (float i : subVec) {
//            sublist.append(i);
//        }
//        list.append(sublist);
//    }
//    return list;
//}

fl calculateL2Norms(const std::vector<fl>& gradients) {
    fl sum = 0.0;
        for (fl value : gradients) {
            sum += value * value;
        }

    return std::sqrt(sum);
}

void MODEA::operator()(model &m, model &m_ad4, output_container_MODEA &out, const precalculate &p,
                       const precalculate &prec_ad4, const precalculate &prec_vdw, const precalculate &prec_forcefield, const precalculate& prec_DLIGAND2,
                       const igrid &ig, const igrid& DLIGAND2_grid, const precalculate &p_widened, const igrid &ig_widened, const vec &corner1,
                       const vec &corner2, incrementable *increment_me, rng &generator, const scoring_function &sf,
                       const scoring_function &sf_ad4, const boost::optional<model> &ref) const {

    vec v(10, 10, 10);
    //用于初始化梯度对象change
    conf_size s = m.get_size();
    //创建能量计算对象energy_cal
    energy_cal energy(&m, &m_ad4, &p, &prec_ad4, &prec_vdw, &prec_forcefield,&prec_DLIGAND2,&ig, &DLIGAND2_grid, v);
    //创建算法的辅助对象，用于进化操作
    MODEA_aux aux(m);
    //创建梯度对象
    change g(s);
    //确定评分函数个数
    std::vector<fl> obj(3, 0.0);
    extern Ligand * liganda;
//    extern OBMolecule obl;

    extern string lig_name;
    string pdb_name = "./"+lig_name+"/"+lig_name+"_protein.pdb";
    string mol2_name = "./"+lig_name+"/"+lig_name+"_ligand.mol2";


//    extern string lig_name;
//    string pdb_name = "../example/"+lig_name+"/"+lig_name+"_protein.pdb";
//    string mol2_name = "../example/"+ lig_name+"/"+lig_name+"_ligand.mol2";

    //读取DLIGAND2函数所需要用到的蛋白质和配体文件中相关信息
    Molecule_DLIGAND2 *mol = new Molecule_DLIGAND2(pdb_name);
    mol->rdmol2(mol2_name);
    //用于初始化种群的临时个体
    output_type_MODEA tmp(s, obj,*liganda,*mol);


//    atomv m_atoms = m.get_ligand_atoms();
//    int count = 0;
//    VINA_FOR(i, m_atoms.size()){
//        VINA_FOR(j, tmp.num_atom){
//
//            if(areAlmostEqual(m.get_ligand_coords()[i][0], tmp.atom[j].coor[0], 0.01) && areAlmostEqual(m.get_ligand_coords()[i][1], tmp.atom[j].coor[1], 0.01) && areAlmostEqual(m.get_ligand_coords()[i][2], tmp.atom[j].coor[2], 0.01)){
//                tmp.atom[j].coor[0] = m_ad4.coords[i][0];
//                tmp.atom[j].coor[1] = m_ad4.coords[i][1];
//                tmp.atom[j].coor[2] = m_ad4.coords[i][2];
//            }
//        }
//    }
//    float bind_score = energy.energy_cal_xscore(&tmp);
//    float a = bind_score * (-1.364);
//    float b = a;





//    extern output_type_MODEA std_tmp;
//    extern model std_model;
//
//    std_tmp = tmp;
////    float intra;
////    intra = energy.energy_cal_intra(std_tmp.c,g);
////    float b  = energy.energy_cal_autodock4_std(sf_ad4, std_tmp.c,std_model);
//
//
//    for (int k = 0; k < std_tmp.objectives.size(); k++) {
//        if (k == 0) {
//
//
////        int count = 0;
//
//            //  ҵ ƥ   
//
//
////                intra = energy.energy_cal_intra(tmp.c,g);
////                tmp.objectives[k] = energy.energy_cal_force(tmp.c);
////                intra = energy.AD4_energy_cal_intra(tmp.c);
//            //tmp.objectives[k] = intra;
////                tmp.objectives[k] = energy.energy_cal_forcefield(tmp.c);
////                intra = energy.energy_cal_intra(tmp.c,g);
////                tmp.objectives[k] = energy.energy_cal_vina(sf, tmp.c, intra,g);
//            std_tmp.objectives[k] = energy.energy_cal_autodock4_std(sf_ad4, tmp.c,std_model);
////                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4, tmp.c);
////                float bind_score = energy.energy_cal_xscore(&tmp);
////                tmp.objectives[k] = bind_score * (-1.364);
//        } else if (k == 1) {
////                tmp.objectives[k] = energy.energy_cal_hydrophobic(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_inter(sf, tmp.c, intra,g)/(1+0.05846*(m_ad4.ligand_degrees_of_freedom(0)));
////                tmp.objectives[k] = energy.AD4_energy_cal_inter(sf, tmp.c, intra);
//            float bind_score = energy.energy_cal_xscore(&tmp);
//            std_tmp.objectives[k] = bind_score * (-1.364);
////                tmp.objectives[k] = energy.energy_cal_DLIGAND2(tmp);
//
//
//        } else if (k == 2) {
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
//            std_tmp.objectives[k] = energy.energy_cal_DLIGAND2(tmp);
//
//        } else if (k == 3) {
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
////                tmp.objectives[k] = energy.energy_cal_Smog2016(tmp);
//            std_tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4, std_tmp.c);
//
//        }
//
//    }
//
//    int ss = 0;



//    output_type_MODEA ns(s,obj);//      Ϊ           Ա
//    init_individual1(m, energy, generator, corner1, corner2, sf, ns);
//    output_type_MODEA tmp(s, obj);
    //确定临时个体的ad_number
    atomv m_atoms = m_ad4.get_ligand_atoms();
    int count = 0;
    VINA_FOR(i, m_atoms.size()){
        VINA_FOR(j, tmp.num_atom){

            if(areAlmostEqual(m_ad4.get_ligand_coords()[i][0], tmp.atom[j].coor[0], 0.01) && areAlmostEqual(m_ad4.get_ligand_coords()[i][1], tmp.atom[j].coor[1], 0.01) && areAlmostEqual(m_ad4.get_ligand_coords()[i][2], tmp.atom[j].coor[2], 0.01)){
                tmp.atom[j].ad_number = m_atoms[i].ad_number;
                count++;
                break;
            }
        }
    }
    //确定临时个体的ad_num
    int num = 0;
    VINA_FOR(i, m_atoms.size()){
        if(num == tmp.lig_num){
            break;
        }
        for(int j =  tmp.atoms.size()-1; j > 0; j--){
            if(areAlmostEqual(m_ad4.get_ligand_coords()[i][0], tmp.atoms[j].x[0], 0.01) && areAlmostEqual(m_ad4.get_ligand_coords()[i][1], tmp.atoms[j].x[1], 0.01) && areAlmostEqual(m_ad4.get_ligand_coords()[i][2], tmp.atoms[j].x[2], 0.01)){
                tmp.atoms[j].ad_num = m_atoms[i].ad_number;
                num++;
                break;
            }
        }
    }
    //初始化临时个体
    aux.init_tmp = tmp;
    //初始化临时个体的position、orientation、rotor_angle
    for (int i = 0; i < aux.popsize; i++) {
        tmp.c.randomize(corner1, corner2, generator);

        for (int j = 0; j < 3; j++) {
            if (j == 1) {
                tmp.rotor_angle[j] = random_fl(0, pi, generator);
            } else {
                tmp.rotor_angle[j] = random_fl(0, 2 * pi, generator);
            }
        }
        //    Ҫ    ʼ   õ      Ƕ ת  Ϊ  Ԫ    Ȼ   ٽ    ת   õ      ʼ    Ԫ      tmp  orientation
        //初始化临时个体tmp的orientation
        std::vector<double> quater;
        quater.resize(4);
        threeangles_to_quaternion(tmp.rotor_angle, quater);
        fl cos_theta = cos(quater[0]);
        fl sin_theta = sin(quater[0]);
        qt q(cos_theta, quater[1] * sin_theta, quater[2] * sin_theta, quater[3] * sin_theta);
        quaternion_normalize(q);//  һ
        tmp.c.ligands[0].rigid.orientation = q;
        //初始化总能量
        tmp.e = 0;
        tmp.rm = 0.0;
        fl intra;
        //初始化三个评分函数的梯度对象
        change DLIGAND2_g(s);
        change vdw_g(s);
        change vina_g(s);
        //利用临时个体tmp初始化整个初始种群和变异种群
        for (int k = 0; k < tmp.objectives.size(); k++) {
            if (k == 0) {
                intra = energy.energy_cal_intra(tmp.c,g);
//                tmp.objectives[k] = energy.energy_cal_force(tmp.c);
//                intra = energy.AD4_energy_cal_intra(tmp.c);
                //tmp.objectives[k] = intra;
//                tmp.objectives[k] = energy.energy_cal_forcefield(tmp.c);
//                intra = energy.energy_cal_intra(tmp.c,g);
//                tmp.objectives[k] = energy.energy_cal_vina(sf, tmp.c, intra,g);
//                tmp.objectives[k] = energy.energy_cal_autodock4(sf_ad4,tmp.c);
//                tmp.objectives[k] = energy.energy_cal_intra(tmp.c,g);
//                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4, tmp.c);
//                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4, tmp.c,vdw_g);
                  tmp.objectives[k] = energy.energy_cal_vina_deriv(sf, tmp.c, intra,vina_g);
//                tmp.objectives[k] = energy.energy_cal_autodock4(sf_ad4, tmp.c,vina_g);
//                float bind_score = energy.energy_cal_xscore(&tmp);
//                tmp.objectives[k] = bind_score * (-1.364);

            } else if (k == 1) {
//                tmp.objectives[k] = energy.energy_cal_hydrophobic(tmp.c);
//                tmp.objectives[k] = energy.energy_cal_inter(sf, tmp.c, intra,g)/(1+0.05846*(m_ad4.ligand_degrees_of_freedom(0)));
//                tmp.objectives[k] = energy.AD4_energy_cal_inter(sf, tmp.c, intra);
//                float bind_score = energy.energy_cal_xscore(&tmp);
//                tmp.objectives[k] = bind_score * (-1.364);
//                tmp.objectives[k] = energy.energy_cal_DLIGAND2(tmp);
                tmp.objectives[k] = energy.energy_cal_DLIGAND2(tmp,DLIGAND2_g);
//                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4, tmp.c,vdw_g);

            } else if (k == 2) {
//                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
//                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
//                tmp.objectives[k] = energy.energy_cal_DLIGAND2(tmp);
//                tmp.objectives[k] = energy.energy_cal_DLIGAND2(tmp,DLIGAND2_g);
//                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4, tmp.c);
                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4, tmp.c,vdw_g);
//                tmp.objectives[k] = energy.energy_cal_autodock4(sf_ad4, tmp.c,vdw_g);

//                tmp.objectives[k] = energy.energy_cal_autodock4(sf_ad4,tmp.c);

            } else if (k == 3) {
//                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
//                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
//                tmp.objectives[k] = energy.energy_cal_Smog2016(tmp);
                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4, tmp.c,vdw_g);

            }

            aux.nowpopulation[i] = tmp;
            aux.mutatepopulation[i] = tmp;
        }
    }
int sss = 0;
    //初始化MODEA算法的~pop种群
    //~pop
    for (int i = 0; i < aux.popsize; i++) {

        for (int j = 0; j < aux.dim; j++) {

                if (j < 3) {
                    tmp.c.ligands[0].rigid.position[j] = corner1[j] + corner2[j] + aux.nowpopulation[i](j);
                }
                else if (j >= 3 && j < 6) {
                   if(j==3 || j==5){
                       tmp.rotor_angle[j - 3] = 2 * pi + aux.nowpopulation[i].rotor_angle[j - 3];
                   }else if(j==4){
                       tmp.rotor_angle[j - 3] = pi + aux.nowpopulation[i].rotor_angle[j - 3];
                   }
                }
                else {
                    tmp.c.ligands[0].torsions[j - 6] = aux.nowpopulation[i](j);

                    normalize_angle(tmp.c.ligands[0].torsions[j - 6]);
                }

            }
        std::vector<double> quater;
        quater.resize(4);
        threeangles_to_quaternion(tmp.rotor_angle, quater);
        fl cos_theta = cos(quater[0]);
        fl sin_theta = sin(quater[0]);
        qt q(cos_theta, quater[1] * sin_theta, quater[2] * sin_theta, quater[3] * sin_theta);
        quaternion_normalize(q);//  һ
        tmp.c.ligands[0].rigid.orientation = q;

        fl intra;
        change DLIGAND2_g(s);
        change vdw_g(s);
        change vina_g(s);

        for (int k = 0; k < tmp.objectives.size(); k++) {
            if (k == 0) {
                intra = energy.energy_cal_intra(tmp.c,g);
//                tmp.objectives[k] = energy.energy_cal_force(tmp.c);
//                intra = energy.AD4_energy_cal_intra(tmp.c);
                //tmp.objectives[k] = intra;
//                tmp.objectives[k] = energy.energy_cal_forcefield(tmp.c);
//                intra = energy.energy_cal_intra(tmp.c,g);
//                tmp.objectives[k] = energy.energy_cal_vina(sf, tmp.c, intra,g);
//                tmp.objectives[k] = energy.energy_cal_autodock4(sf_ad4,tmp.c);
//                tmp.objectives[k] = energy.energy_cal_intra(tmp.c,g);
//                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4, tmp.c);
//                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4, tmp.c,vdw_g);
                tmp.objectives[k] = energy.energy_cal_vina_deriv(sf, tmp.c, intra,vina_g);
//                tmp.objectives[k] = energy.energy_cal_autodock4(sf_ad4, tmp.c,vina_g);
//                float bind_score = energy.energy_cal_xscore(&tmp);
//                tmp.objectives[k] = bind_score * (-1.364);

            } else if (k == 1) {
//                tmp.objectives[k] = energy.energy_cal_hydrophobic(tmp.c);
//                tmp.objectives[k] = energy.energy_cal_inter(sf, tmp.c, intra,g)/(1+0.05846*(m_ad4.ligand_degrees_of_freedom(0)));
//                tmp.objectives[k] = energy.AD4_energy_cal_inter(sf, tmp.c, intra);
//                float bind_score = energy.energy_cal_xscore(&tmp);
//                tmp.objectives[k] = bind_score * (-1.364);
//                tmp.objectives[k] = energy.energy_cal_DLIGAND2(tmp);
                tmp.objectives[k] = energy.energy_cal_DLIGAND2(tmp,DLIGAND2_g);
//                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4, tmp.c,vdw_g);

            } else if (k == 2) {
//                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
//                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
//                tmp.objectives[k] = energy.energy_cal_DLIGAND2(tmp,DLIGAND2_g);
//                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4, tmp.c);
                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4, tmp.c,vdw_g);
//                tmp.objectives[k] = energy.energy_cal_autodock4(sf_ad4, tmp.c,vdw_g);
//                tmp.objectives[k] = energy.energy_cal_autodock4(sf_ad4,tmp.c);


            } else if (k == 3) {
//                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
//                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
//                tmp.objectives[k] = energy.energy_cal_Smog2016(tmp);
                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4, tmp.c,vdw_g);

            }

            aux.oppositepopulation[i] = tmp;

        }
    }

//    for (int i = 0; i < 20; i++) {
//
//        for (int j = 0; j < aux.dim; j++) {
//
//            if (j < 3) {
//                tmp.c.ligands[0].rigid.position[j] = corner1[j] + corner2[j] + aux.nowpopulation[i](j);
//            }
//            else if (j >= 3 && j < 6) {
//                if(j==3 || j==5){
//                    tmp.rotor_angle[j - 3] = 2 * pi + aux.nowpopulation[i].rotor_angle[j - 3];
//                }else if(j==4){
//                    tmp.rotor_angle[j - 3] = pi + aux.nowpopulation[i].rotor_angle[j - 3];
//                }
//            }
//            else {
//                tmp.c.ligands[0].torsions[j - 6] = aux.nowpopulation[i](j);
//
//                normalize_angle(tmp.c.ligands[0].torsions[j - 6]);
//            }
//
//        }
//        std::vector<double> quater;
//        quater.resize(4);
//        threeangles_to_quaternion(tmp.rotor_angle, quater);
//        fl cos_theta = cos(quater[0]);
//        fl sin_theta = sin(quater[0]);
//        qt q(cos_theta, quater[1] * sin_theta, quater[2] * sin_theta, quater[3] * sin_theta);
//        quaternion_normalize(q);//  һ
//        tmp.c.ligands[0].rigid.orientation = q;
//
//        fl intra;
//        change DLIGAND2_g(s);
//        change vdw_g(s);
//        change vina_g(s);
//
//        for (int k = 0; k < tmp.objectives.size(); k++) {
//            if (k == 0) {
//                intra = energy.energy_cal_intra(tmp.c,g);
////                tmp.objectives[k] = energy.energy_cal_force(tmp.c);
////                intra = energy.AD4_energy_cal_intra(tmp.c);
//                //tmp.objectives[k] = intra;
////                tmp.objectives[k] = energy.energy_cal_forcefield(tmp.c);
////                intra = energy.energy_cal_intra(tmp.c,g);
////                tmp.objectives[k] = energy.energy_cal_vina(sf, tmp.c, intra,g);
////                tmp.objectives[k] = energy.energy_cal_autodock4(sf_ad4,tmp.c);
////                tmp.objectives[k] = energy.energy_cal_intra(tmp.c,g);
////                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4, tmp.c);
////                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4, tmp.c,vdw_g);
//                tmp.objectives[k] = energy.energy_cal_vina_deriv(sf, tmp.c, intra,vina_g);
////                tmp.objectives[k] = energy.energy_cal_autodock4(sf_ad4, tmp.c,vina_g);
////                float bind_score = energy.energy_cal_xscore(&tmp);
////                tmp.objectives[k] = bind_score * (-1.364);
//
//            } else if (k == 1) {
////                tmp.objectives[k] = energy.energy_cal_hydrophobic(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_inter(sf, tmp.c, intra,g)/(1+0.05846*(m_ad4.ligand_degrees_of_freedom(0)));
////                tmp.objectives[k] = energy.AD4_energy_cal_inter(sf, tmp.c, intra);
////                float bind_score = energy.energy_cal_xscore(&tmp);
////                tmp.objectives[k] = bind_score * (-1.364);
////                tmp.objectives[k] = energy.energy_cal_DLIGAND2(tmp);
//                tmp.objectives[k] = energy.energy_cal_DLIGAND2(tmp,DLIGAND2_g);
//
//
//            } else if (k == 2) {
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
////                tmp.objectives[k] = energy.energy_cal_DLIGAND2(tmp,DLIGAND2_g);
////                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4, tmp.c);
//                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4, tmp.c,vdw_g);
////                tmp.objectives[k] = energy.energy_cal_autodock4(sf_ad4, tmp.c,vdw_g);
////                tmp.objectives[k] = energy.energy_cal_autodock4(sf_ad4,tmp.c);
//
//
//            } else if (k == 3) {
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
////                tmp.objectives[k] = energy.energy_cal_Smog2016(tmp);
//                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4, tmp.c,vdw_g);
//
//            }
//
//            aux.archive_gradient[i] = tmp;
//
//        }
//    }
    //将初始化的初始种群和opposite种群合并，再通过快速非支配排序，得到最终的初始种群
    aux.R.clear();
    //    һ   ĸ     Ⱥ ϲ
    aux.R.insert(aux.R.end(), aux.nowpopulation.begin(), aux.nowpopulation.end());
    aux.R.insert(aux.R.end(), aux.oppositepopulation.begin(), aux.oppositepopulation.end());
    aux.nd_cd_sort(aux.R);
    //可无视
    std::vector<float> iterations;
    std::vector<float> best_fitness;
    std::vector<float> best_fitness_vdw;
    std::vector<float> best_fitness_vina;
    std::vector<float> best_fitness_DLIGAND2;
    extern float min_vina;
    extern float min_DLIGAND2;
    extern float min_vdw;
    extern float min_e;
    extern double deviceTime;
    //开始种群迭代
    for (int g1 = 0; g1 < num_steps; g1++) {
        if (increment_me)
            ++(*increment_me);
        //gradient
        change DLIGAND2_g(s);
        change vdw_g(s);
        change vina_g(s);
        change total_g(s);


                struct timeval t1,t2;
        double timeuse;
        gettimeofday(&t1,NULL);
//        cout<<"============================================="<<endl;
//        cout<<"aux.nowpopulation[10].objectives[0]: "<<aux.nowpopulation[10].objectives[0]<<endl;
//        cout<<"aux.nowpopulation[10].objectives[1]: "<<aux.nowpopulation[10].objectives[1]<<endl;
//        cout<<"aux.mutatepopulation[10].objectives[0]: "<<aux.mutatepopulation[10].objectives[0]<<endl;
//        cout<<"aux.mutatepopulation[10].objectives[1]: "<<aux.mutatepopulation[10].objectives[1]<<endl;
//        cout<<"start......."<<endl;
//        cout<<" "<<endl;
        //变异交叉操作
        aux.mutate_cross(m, energy, generator, corner1, corner2);
        //可无视
        aux.select(m, energy, generator);//select     ж nowpopulation
        //重点，计算能量评分函数
        for (int i = 0; i < aux.popsize; i++) {
            aux.calculate_objectives(m, p, ig, v, aux.mutatepopulation[i], sf, sf_ad4, energy, ref, g,vdw_g, vina_g, DLIGAND2_g);
        }
        gettimeofday(&t2,NULL);
        timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
        deviceTime+=timeuse;
//        std::cout<<"run time is: "<<timeuse<<std::endl;  //输出时间（单位：ｓ）
//                gettimeofday(&t2,NULL);
//        timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
//        cout<<"============================================="<<endl;
//        cout<<"aux.nowpopulation[10].objectives[0]: "<<aux.nowpopulation[10].objectives[0]<<endl;
//        cout<<"aux.nowpopulation[10].objectives[1]: "<<aux.nowpopulation[10].objectives[1]<<endl;
//        cout<<"aux.mutatepopulation[10].objectives[0]: "<<aux.mutatepopulation[10].objectives[0]<<endl;
//        cout<<"aux.mutatepopulation[10].objectives[1]: "<<aux.mutatepopulation[10].objectives[1]<<endl;
//        cout<<"end......."<<endl;
//        cout<<" "<<endl;
//        cout<<"vina_time = "<<timeuse<<endl;  //输出时间（单位：ｓ）

        ////
//        aux.calculate_crowding_distance(aux.mutatepopulation);
//        aux.calculate_crowding_distance(aux.nowpopulation);
//
//        std::sort(aux.mutatepopulation.begin(),aux.mutatepopulation.end(),crowding_distance_comparator);
//        std::sort(aux.nowpopulation.begin(),aux.nowpopulation.end(),crowding_distance_comparator);
        //得到变异种群中最优对象
        int flag = 0;
        float maxe = -std::numeric_limits<float>::infinity();
        std::unordered_map<int, std::vector<output_type_MODEA>> F = aux.fast_nondominated_sort(aux.mutatepopulation);
        aux.calculate_crowding_distance(aux.mutatepopulation);
        for(int i = 0; i < aux.mutatepopulation.size(); i++){
            if(aux.mutatepopulation[i].rank==1){
                if(aux.mutatepopulation[i].crowding_distance>maxe){
                    flag = i;
                    maxe = aux.mutatepopulation[i].crowding_distance;
                }
            }
        }
        output_type_MODEA mutation_best_tmp = aux.mutatepopulation[flag];
//        tee log;
//        log.init("./2vkm/energy.log");
//        log << "mode |   vina       |   DLIGAND2 |     vdw    |     e    |\n";
//        log << "     | (kcal/mol)  | (kcal/mol) | (kcal/mol) |  (float) |\n";
//        log << "-----+-------------+------------+-------------+----------\n";
//        log<<"before MGDA energy: \n"<<aux.mutatepopulation[flag].objectives[0]
//                << "    " << std::setw(8) << std::setprecision(3) <<aux.mutatepopulation[flag].objectives[1]
//                << "    " << std::setw(8) << std::setprecision(3) <<aux.mutatepopulation[flag].objectives[2]
//                << "    " << std::setw(8) << std::setprecision(3) << aux.mutatepopulation[flag].objectives[0]+aux.mutatepopulation[flag].objectives[1]+aux.mutatepopulation[flag].objectives[2]<<"\n";

        //MGDA
        if(1==0) {
            //gradient
            //得到配体的旋转键数目
            sz total_torsions = m_ad4.ligand_degrees_of_freedom(0);
            std::vector<fl> vdw_gradient(6 + total_torsions);
            std::vector<fl> DLIGAND2_gradient(6 + total_torsions);
            std::vector<fl> vina_gradient(6 + total_torsions);

            vec V_d_position(zero_vec);
            vec S_d_position(zero_vec);
            vec V_d_orientation(zero_vec);
            vec S_d_orientation(zero_vec);
            flv V_d_torsions(total_torsions, 0);
            flv S_d_torsions(total_torsions, 0);
            int t = 1;
            //对最优个体执行20次梯度下降
            for (int k = 0; k < 20; k++) {
                aux.calculate_objectives(m, p, ig, v, aux.mutatepopulation[flag], sf, sf_ad4, energy, ref, g, vdw_g,
                                         vina_g, DLIGAND2_g);
//            aux.mutatepopulation[0].e = energy.energy_cal_DLIGAND2(aux.mutatepopulation[0], DLIGAND2_g);
                //得到DLIGAND2梯度
                for (int i = 0; i < DLIGAND2_gradient.size(); i++) {
                    if (i < 3) {
                        DLIGAND2_gradient[i] = DLIGAND2_g.ligands[0].rigid.position[i];
                    } else if (i >= 3 && i < 6) {
                        DLIGAND2_gradient[i] = DLIGAND2_g.ligands[0].rigid.orientation[i - 3];
                    } else {
                        DLIGAND2_gradient[i] = DLIGAND2_g.ligands[0].torsions[i - 6];
                    }
                }
                //得到vdw梯度
                for (int i = 0; i < vdw_gradient.size(); i++) {
                    if (i < 3) {
                        vdw_gradient[i] = vdw_g.ligands[0].rigid.position[i];
                    } else if (i >= 3 && i < 6) {
                        vdw_gradient[i] = vdw_g.ligands[0].rigid.orientation[i - 3];
                    } else {
                        vdw_gradient[i] = vdw_g.ligands[0].torsions[i - 6];
                    }
                }
                //得到vina梯度
                for (int i = 0; i < vina_gradient.size(); i++) {
                    if (i < 3) {
                        vina_gradient[i] = vina_g.ligands[0].rigid.position[i];
                    } else if (i >= 3 && i < 6) {
                        vina_gradient[i] = vina_g.ligands[0].rigid.orientation[i - 3];
                    } else {
                        vina_gradient[i] = vina_g.ligands[0].torsions[i - 6];
                    }
                }
                fl vina_norm = calculateL2Norms(vina_gradient);
                fl DLIGAND2_norm = calculateL2Norms(DLIGAND2_gradient);
                fl vdw_norm = calculateL2Norms(vdw_gradient);
                //归一化DLIGAND2_gradient
                for (int i = 0; i < DLIGAND2_gradient.size(); i++) {
                    if (DLIGAND2_norm != 0) {
                        DLIGAND2_gradient[i] = DLIGAND2_gradient[i] / DLIGAND2_norm;
                    }
                }
                //归一化vina_gradient
                for (int i = 0; i < vina_gradient.size(); i++) {
                    if (vina_norm != 0) {
                        vina_gradient[i] = vina_gradient[i] / vina_norm;
                    }
                }
                //归一化vdw_gradient
                for (int i = 0; i < vdw_gradient.size(); i++) {
                    if (vdw_norm != 0) {
                        vdw_gradient[i] = vdw_gradient[i] / vdw_norm;
                    }
                }
                std::vector<std::vector<float>> total_gradients(3, std::vector<float>(DLIGAND2_gradient.size()));
                for (int i = 0; i < total_gradients.size(); i++) {
                    for (int j = 0; j < DLIGAND2_gradient.size(); j++) {
                        if(i==0){
                            total_gradients[i][j] =  DLIGAND2_gradient[j];
                        }else if(i==1){
                            total_gradients[i][j] =  vina_gradient[j];
                        }else if(i==2){
                            total_gradients[i][j] = vdw_gradient[j] ;
                        }
                    }
                }
                std::vector<float> alpha_vector(DLIGAND2_gradient.size(), 1.0 / DLIGAND2_gradient.size());
                std::vector<float> c_result;
//                for (std::vector<float> value: total_gradients) {
//                    for(float v:value){
//                        std::cout << v << " ";
//
//                    }
//                    cout <<endl;
//                }
                //根据FrankWolf优化器得到α系数
                MinNormSolverNumpy msn;
                std::pair<std::vector<float>, float> result;
                result = msn.find_min_norm_element(total_gradients);
//                                        for (float value: result.first) {
//                            std::cout << value << " ";
//                        }
//                        std::cout << std::endl;
                c_result = result.first;

//                try {
//                boost::python::object ignored = boost::python::exec("print(\"Hello, World\")");
                    // 加载 Python 模块

                    // 调用函数

//                    extern boost::python::object solver_class;
//                    extern boost::python::object find_min_norm_element;
//                    extern boost::python::list py_args;
//                    extern boost::python::object result;
//                    Py_Initialize();
//                    find_min_norm_element = solver_class.attr("find_min_norm_element");
//
//                    py_args = vectorToList(total_gradients);
//                    result = find_min_norm_element(py_args);
//
//                    // 检查 result 是否是一个元组，并且包含两个元素
//                    if (PyTuple_Check(result.ptr()) && PyTuple_Size(result.ptr()) == 2) {
//                        // 提取元组中的第一个元素（解向量）
//                        boost::python::object sol_vec_obj = result[0];
//                        // 提取元组中的第二个元素（最小范数）
//                        boost::python::object nd_obj = result[1];
//                        c_result = convert_to_vector(sol_vec_obj);
//
////                        // 打印结果（可选）
//                        for (float value: c_result) {
//                            std::cout << value << " ";
//                        }
//                        std::cout << std::endl;
//                        // 进一步处理 sol_vec_obj 和 nd_obj
//                        // 例如，将 sol_vec_obj 转换为 std::vector
//                        // 注意：这需要额外的类型检查和转换代码
//                        Py_Finalize();
//                    }
                    // 将结果转换为 C++ 类型（假设为 double）
                    // 将 Python 结果转换为 C++ vector
//            std::vector<float> c_result = convertPyObjectToVector(result);


//                } catch (boost::python::error_already_set const &) {
//                    PyErr_Print();
//            std::string perror_str = boost::python::p

//            std::cout << "Error in Python: " << perror_str << std::endl;
//                }
//            fl alpha = calculate_alpha(vdw_gradient, DLIGAND2_gradient);
//            std::vector<fl> w(DLIGAND2_gradient.size());
//            VINA_FOR(j, w.size()) {
//                w[j] = alpha * vdw_gradient[j] + (1 - alpha) * DLIGAND2_gradient[j];
//                if (j < 3) {
//                    total_g.ligands[0].rigid.position[j] = w[j];
//                } else if (j >= 3 && j < 6) {
//                    total_g.ligands[0].rigid.orientation[j - 3] = w[j];
//                } else {
//                    total_g.ligands[0].torsions[j - 6] = w[j];
//                }
//            }
                //根据α系数得到公共梯度方向
            VINA_FOR(j, DLIGAND2_gradient.size()) {

                if (j < 3) {
                    total_g.ligands[0].rigid.position[j] =   c_result[0] * DLIGAND2_gradient[j] + c_result[1] * vina_gradient[j] + c_result[2] * vdw_gradient[j];
                } else if (j >= 3 && j < 6) {
                    total_g.ligands[0].rigid.orientation[j - 3] = c_result[0] * DLIGAND2_gradient[j] + c_result[1] * vina_gradient[j] + c_result[2] * vdw_gradient[j];
                } else {
                    total_g.ligands[0].torsions[j - 6] = c_result[0] * DLIGAND2_gradient[j] + c_result[1] * vina_gradient[j] + c_result[2] * vdw_gradient[j];
                }
            }

                //archive gradient
//                if(aux.archive_gradient.size()==0){
//                    aux.archive_gradient[0] =   aux.mutatepopulation[flag];
//                }
//                aux.mutatepopulation[flag].total_g = total_g;
//                int sig = aux.add(aux.mutatepopulation[flag]);
//                if(sig==-1) {
//                    int a1 = random_int(0, aux.archive_gradient.size() - 1, generator);
//                    int a2 = random_int(0, aux.archive_gradient.size() - 1, generator);
//                    output_type_MODEA archive1 = aux.archive_gradient[a1];
//                    output_type_MODEA archive2 = aux.archive_gradient[a2];
//                    if (crowding_distance_comparator(archive1, archive2))
//                        total_g = archive1.total_g;
//                    else total_g = archive2.total_g;
//                }
                //根据公共梯度方向执行adam梯度下降算法
                increment_position(aux.mutatepopulation[flag].c, total_g, beta1, beta2, V_d_position, S_d_position, lr,
                                   sigma,
                                   t);
                increment_orientation(aux.mutatepopulation[flag].c, total_g, beta1, beta2, V_d_orientation,
                                      S_d_orientation, lr,
                                      sigma, t);
                increment_torsions(aux.mutatepopulation[flag].c, total_g, beta1, beta2, V_d_torsions, S_d_torsions, lr,
                                   sigma,
                                   t);
                t++;
            }


//            Py_Finalize();
//        aux.mutatepopulation[0].e = energy(aux.mutatepopulation[0].c, g);
            aux.calculate_objectives(m, p, ig, v, aux.mutatepopulation[flag], sf, sf_ad4, energy, ref, g, vdw_g, vina_g,
                                     DLIGAND2_g);
//            aux.calculate_crowding_distance(aux.mutatepopulation);
            quaternion_to_3angles(aux.mutatepopulation[flag], aux.mutatepopulation[flag].rotor_angle);
            //std::cout << "energy:" << std::setw(15) << energy(mutatepopulation[best_index].c, g) << std::endl;
//            if (aux.nowpopulation[0].crowding_distance < aux.mutatepopulation[0].crowding_distance) {
//                aux.nowpopulation[0] = aux.mutatepopulation[0];
//            }

            if(!(aux.mutatepopulation[flag]<mutation_best_tmp)){
                aux.mutatepopulation[flag] = mutation_best_tmp;
            }
//            log<<"after MGDA energy: \n"<<aux.mutatepopulation[flag].objectives[0]
//               <<"    " << std::setw(8) << std::setprecision(3)<<aux.mutatepopulation[flag].objectives[1]
//               <<"    " << std::setw(8) << std::setprecision(3)<<aux.mutatepopulation[flag].objectives[2]
//               << "    " << std::setw(8) << std::setprecision(3) << aux.mutatepopulation[flag].objectives[0]+aux.mutatepopulation[flag].objectives[1]+aux.mutatepopulation[flag].objectives[2]<<"\n";
//            log.endl();
        }

        //PCGrad
        if(1==0){
            //gradient
            sz total_torsions = m_ad4.ligand_degrees_of_freedom(0);
            std::vector<fl> vdw_gradient(6 + total_torsions);
            std::vector<fl> DLIGAND2_gradient(6 + total_torsions);
            std::vector<fl> vina_gradient(6 + total_torsions);

            vec V_d_position(zero_vec);
            vec S_d_position(zero_vec);
            vec V_d_orientation(zero_vec);
            vec S_d_orientation(zero_vec);
            flv V_d_torsions(total_torsions, 0);
            flv S_d_torsions(total_torsions, 0);
            int t = 1;
            for (int k = 0; k < 20; k++) {
                aux.calculate_objectives(m, p, ig, v, aux.mutatepopulation[flag], sf, sf_ad4, energy, ref, g, vdw_g,
                                         vina_g, DLIGAND2_g);
//            aux.mutatepopulation[0].e = energy.energy_cal_DLIGAND2(aux.mutatepopulation[0], DLIGAND2_g);
                for (int i = 0; i < DLIGAND2_gradient.size(); i++) {
                    if (i < 3) {
                        DLIGAND2_gradient[i] = DLIGAND2_g.ligands[0].rigid.position[i];
                    } else if (i >= 3 && i < 6) {
                        DLIGAND2_gradient[i] = DLIGAND2_g.ligands[0].rigid.orientation[i - 3];
                    } else {
                        DLIGAND2_gradient[i] = DLIGAND2_g.ligands[0].torsions[i - 6];
                    }
                }
                for (int i = 0; i < vdw_gradient.size(); i++) {
                    if (i < 3) {
                        vdw_gradient[i] = vdw_g.ligands[0].rigid.position[i];
                    } else if (i >= 3 && i < 6) {
                        vdw_gradient[i] = vdw_g.ligands[0].rigid.orientation[i - 3];
                    } else {
                        vdw_gradient[i] = vdw_g.ligands[0].torsions[i - 6];
                    }
                }
                for (int i = 0; i < vina_gradient.size(); i++) {
                    if (i < 3) {
                        vina_gradient[i] = vina_g.ligands[0].rigid.position[i];
                    } else if (i >= 3 && i < 6) {
                        vina_gradient[i] = vina_g.ligands[0].rigid.orientation[i - 3];
                    } else {
                        vina_gradient[i] = vina_g.ligands[0].torsions[i - 6];
                    }
                }
                fl vina_norm = calculateL2Norms(vina_gradient);
                fl DLIGAND2_norm = calculateL2Norms(DLIGAND2_gradient);
                fl vdw_norm = calculateL2Norms(vdw_gradient);
                for (int i = 0; i < DLIGAND2_gradient.size(); i++) {
                    DLIGAND2_gradient[i] = DLIGAND2_gradient[i]/DLIGAND2_norm;
                }
                for (int i = 0; i < vina_gradient.size(); i++) {
                    vina_gradient[i] = vina_gradient[i]/vina_norm;
                }
                for (int i = 0; i < vdw_gradient.size(); i++) {
                    vdw_gradient[i] = vdw_gradient[i]/vdw_norm;
                }
                std::vector<std::vector<float>> total_gradients(3, std::vector<float>(DLIGAND2_gradient.size()));
                for (int i = 0; i < total_gradients.size(); i++) {
                    for (int j = 0; j < DLIGAND2_gradient.size(); j++) {
                        if(i==0){
                            total_gradients[i][j] =  DLIGAND2_gradient[j];
                        }else if(i==1){
                            total_gradients[i][j] =  vina_gradient[j];
                        }else if(i==2){
                            total_gradients[i][j] = vdw_gradient[j] ;
                        }
                    }
                }
             PCGrad pcg;
                std::vector<std::vector<float>> total_gradients_PC(3, std::vector<float>(DLIGAND2_gradient.size()));
                total_gradients_PC = total_gradients;
                pcg.project_gradients(total_gradients,total_gradients_PC);
                std::vector<float> PCG_total_gradient(DLIGAND2_gradient.size());
                VINA_FOR(cc,total_gradients.size()){
                    VINA_FOR(gg,PCG_total_gradient.size()){
                        PCG_total_gradient[gg] += total_gradients_PC[cc][gg];
                    }
                }

                VINA_FOR(j, DLIGAND2_gradient.size()) {

                    if (j < 3) {
                        total_g.ligands[0].rigid.position[j] =   PCG_total_gradient[j];
                    } else if (j >= 3 && j < 6) {
                        total_g.ligands[0].rigid.orientation[j - 3] = PCG_total_gradient[j];
                    } else {
                        total_g.ligands[0].torsions[j - 6] = PCG_total_gradient[j];
                    }
                }
                increment_position(aux.mutatepopulation[flag].c, total_g, beta1, beta2, V_d_position, S_d_position, lr,
                                   sigma,
                                   t);
                increment_orientation(aux.mutatepopulation[flag].c, total_g, beta1, beta2, V_d_orientation,
                                      S_d_orientation, lr,
                                      sigma, t);
                increment_torsions(aux.mutatepopulation[flag].c, total_g, beta1, beta2, V_d_torsions, S_d_torsions, lr,
                                   sigma,
                                   t);
                t++;
            }



            aux.calculate_objectives(m, p, ig, v, aux.mutatepopulation[flag], sf, sf_ad4, energy, ref, g, vdw_g, vina_g,
                                     DLIGAND2_g);
            quaternion_to_3angles(aux.mutatepopulation[flag], aux.mutatepopulation[flag].rotor_angle);

        }

        //MGDA+PCGrad
        if(1==0){
            //gradient
            sz total_torsions = m_ad4.ligand_degrees_of_freedom(0);
            std::vector<fl> vdw_gradient(6 + total_torsions);
            std::vector<fl> DLIGAND2_gradient(6 + total_torsions);
            std::vector<fl> vina_gradient(6 + total_torsions);

            vec V_d_position(zero_vec);
            vec S_d_position(zero_vec);
            vec V_d_orientation(zero_vec);
            vec S_d_orientation(zero_vec);
            flv V_d_torsions(total_torsions, 0);
            flv S_d_torsions(total_torsions, 0);
            int t = 1;
            for (int k = 0; k < 20; k++) {
                aux.calculate_objectives(m, p, ig, v, aux.mutatepopulation[flag], sf, sf_ad4, energy, ref, g, vdw_g,
                                         vina_g, DLIGAND2_g);
//            aux.mutatepopulation[0].e = energy.energy_cal_DLIGAND2(aux.mutatepopulation[0], DLIGAND2_g);
                for (int i = 0; i < DLIGAND2_gradient.size(); i++) {
                    if (i < 3) {
                        DLIGAND2_gradient[i] = DLIGAND2_g.ligands[0].rigid.position[i];
                    } else if (i >= 3 && i < 6) {
                        DLIGAND2_gradient[i] = DLIGAND2_g.ligands[0].rigid.orientation[i - 3];
                    } else {
                        DLIGAND2_gradient[i] = DLIGAND2_g.ligands[0].torsions[i - 6];
                    }
                }
                for (int i = 0; i < vdw_gradient.size(); i++) {
                    if (i < 3) {
                        vdw_gradient[i] = vdw_g.ligands[0].rigid.position[i];
                    } else if (i >= 3 && i < 6) {
                        vdw_gradient[i] = vdw_g.ligands[0].rigid.orientation[i - 3];
                    } else {
                        vdw_gradient[i] = vdw_g.ligands[0].torsions[i - 6];
                    }
                }
                for (int i = 0; i < vina_gradient.size(); i++) {
                    if (i < 3) {
                        vina_gradient[i] = vina_g.ligands[0].rigid.position[i];
                    } else if (i >= 3 && i < 6) {
                        vina_gradient[i] = vina_g.ligands[0].rigid.orientation[i - 3];
                    } else {
                        vina_gradient[i] = vina_g.ligands[0].torsions[i - 6];
                    }
                }
                fl vina_norm = calculateL2Norms(vina_gradient);
                fl DLIGAND2_norm = calculateL2Norms(DLIGAND2_gradient);
                fl vdw_norm = calculateL2Norms(vdw_gradient);
                for (int i = 0; i < DLIGAND2_gradient.size(); i++) {
                    DLIGAND2_gradient[i] = DLIGAND2_gradient[i]/DLIGAND2_norm;
                }
                for (int i = 0; i < vina_gradient.size(); i++) {
                    vina_gradient[i] = vina_gradient[i]/vina_norm;
                }
                for (int i = 0; i < vdw_gradient.size(); i++) {
                    vdw_gradient[i] = vdw_gradient[i]/vdw_norm;
                }
                std::vector<std::vector<float>> total_gradients(3, std::vector<float>(DLIGAND2_gradient.size()));
                for (int i = 0; i < total_gradients.size(); i++) {
                    for (int j = 0; j < DLIGAND2_gradient.size(); j++) {
                        if(i==0){
                            total_gradients[i][j] =  DLIGAND2_gradient[j];
                        }else if(i==1){
                            total_gradients[i][j] =  vina_gradient[j];
                        }else if(i==2){
                            total_gradients[i][j] = vdw_gradient[j] ;
                        }
                    }
                }

                std::vector<float> c_result;
//                for (std::vector<float> value: total_gradients) {
//                    for(float v:value){
//                        std::cout << v << " ";
//
//                    }
//                    cout <<endl;
//                }
                MinNormSolverNumpy msn;
                PCGrad pcg;
                std::pair<std::vector<float>, float> result;
                std::vector<std::vector<float>> total_gradients_PCG(3, std::vector<float>(DLIGAND2_gradient.size()));
                total_gradients_PCG = total_gradients;
                for(int i = 0; i < total_gradients_PCG.size();i++){
                    for(int j = 0; j < total_gradients.size();j++){
                        if(pcg.dot(total_gradients_PCG[i], total_gradients[j]) > pcg.dot(total_gradients_PCG[i], total_gradients_PCG[i])){
                            pcg.plus_vector(total_gradients_PCG[i],pcg.co_vector(
                                    pcg.dot(total_gradients_PCG[i],total_gradients[j])
                                    /
                                    pcg.dot(total_gradients[j],total_gradients[j])
                                    ,total_gradients[j]));
                        }
//                        if(pcg.dot(total_gradients_PCG[j], total_gradients[i]) > pcg.dot(total_gradients_PCG[j], total_gradients_PCG[j])){
//                            pcg.plus_vector(total_gradients_PCG[j],pcg.co_vector(
//                                    pcg.dot(total_gradients_PCG[j],total_gradients[i])
//                                    /
//                                    pcg.dot(total_gradients[i],total_gradients[i])
//                                    ,total_gradients[i]));
//                        }
                    }
                }

//                pcg.project_gradients(total_gradients,total_gradients_PCG);
//                std::vector<float> PCG_total_gradient(DLIGAND2_gradient.size());
//                VINA_FOR(cc,total_gradients.size()){
//                    VINA_FOR(gg,PCG_total_gradient.size()){
//                        PCG_total_gradient[gg] += total_gradients_PCG[cc][gg];
//                    }
//                }
                result = msn.find_min_norm_element(total_gradients_PCG);
//                                        for (float value: result.first) {
//                            std::cout << value << " ";
//                        }
//                        std::cout << std::endl;
                c_result = result.first;



                VINA_FOR(j, DLIGAND2_gradient.size()) {

                    if (j < 3) {
                        total_g.ligands[0].rigid.position[j] =   c_result[0] * DLIGAND2_gradient[j] + c_result[1] * vina_gradient[j] + c_result[2] * vdw_gradient[j];
                    } else if (j >= 3 && j < 6) {
                        total_g.ligands[0].rigid.orientation[j - 3] = c_result[0] * DLIGAND2_gradient[j] + c_result[1] * vina_gradient[j] + c_result[2] * vdw_gradient[j];
                    } else {
                        total_g.ligands[0].torsions[j - 6] = c_result[0] * DLIGAND2_gradient[j] + c_result[1] * vina_gradient[j] + c_result[2] * vdw_gradient[j];
                    }
                }
                increment_position(aux.mutatepopulation[flag].c, total_g, beta1, beta2, V_d_position, S_d_position, lr,
                                   sigma,
                                   t);
                increment_orientation(aux.mutatepopulation[flag].c, total_g, beta1, beta2, V_d_orientation,
                                      S_d_orientation, lr,
                                      sigma, t);
                increment_torsions(aux.mutatepopulation[flag].c, total_g, beta1, beta2, V_d_torsions, S_d_torsions, lr,
                                   sigma,
                                   t);
                t++;
            }


//            Py_Finalize();
//        aux.mutatepopulation[0].e = energy(aux.mutatepopulation[0].c, g);
            aux.calculate_objectives(m, p, ig, v, aux.mutatepopulation[flag], sf, sf_ad4, energy, ref, g, vdw_g, vina_g,
                                     DLIGAND2_g);
//            aux.calculate_crowding_distance(aux.mutatepopulation);
            quaternion_to_3angles(aux.mutatepopulation[flag], aux.mutatepopulation[flag].rotor_angle);
            //std::cout << "energy:" << std::setw(15) << energy(mutatepopulation[best_index].c, g) << std::endl;
//            if (aux.nowpopulation[0].crowding_distance < aux.mutatepopulation[0].crowding_distance) {
//                aux.nowpopulation[0] = aux.mutatepopulation[0];
//            }
        }

        aux.R.clear();
        //    һ   ĸ     Ⱥ ϲ
        aux.R.insert(aux.R.end(), aux.nowpopulation.begin(), aux.nowpopulation.end());
        aux.R.insert(aux.R.end(), aux.mutatepopulation.begin(), aux.mutatepopulation.end());
        aux.nd_cd_sort(aux.R);
//        count=0;
//        std::cout<<"popsize: "<<aux.nowpopulation.size()<<std::endl;
//        for (auto &p: aux.nowpopulation) {
//            if (p.rank == 1) {
//                count++;
//            }
//        }
//        std::cout<<"count: "<<count<<std::endl;
        extern output_container_MODEA out_300;
        extern output_container_MODEA out_500;
        extern output_container_MODEA out_800;
        extern output_container_MODEA out_1000;
//        extern output_container_MODEA out_process_pf;
//        int flag_pf = 0;
//        float maxe_pf = -std::numeric_limits<float>::infinity();
//        for(int i = 0; i < aux.nowpopulation.size(); i++){
//            if(aux.nowpopulation[i].rank==1){
//                if(aux.nowpopulation[i].crowding_distance>maxe_pf){
//                    flag_pf = i;
//                    maxe_pf = aux.nowpopulation[i].crowding_distance;
//                }
//            }
//        }
//        output_type_MODEA out_pf = aux.nowpopulation[flag_pf];
//        m_ad4.set(out_pf.c);
//        out_pf.coords = m_ad4.get_heavy_atom_movable_coords();
//        add_to_output_container_(out_process_pf, out_pf, min_rmsd, num_saved_mins);

        if(g1==299){
            for (auto &s: aux.nowpopulation) {

                s.e = s.objectives[0] + s.objectives[1] + s.objectives[2];
            }
            for (auto &p: aux.nowpopulation) {
                if (p.rank == 1) {
                    m_ad4.set(p.c);
                    p.coords = m_ad4.get_heavy_atom_movable_coords();
                    add_to_output_container_(out_300, p, min_rmsd, num_saved_mins);

                }

            }

        }else if(g1==499){
            for (auto &s: aux.nowpopulation) {

                s.e = s.objectives[0] + s.objectives[1] + s.objectives[2];
            }
            for (auto &p: aux.nowpopulation) {
                if (p.rank == 1) {
                    m_ad4.set(p.c);
                    p.coords = m_ad4.get_heavy_atom_movable_coords();
                    add_to_output_container_(out_500, p, min_rmsd, num_saved_mins);

                }

            }

        }else if(g1==799){
            for (auto &s: aux.nowpopulation) {

                s.e = s.objectives[0] + s.objectives[1] + s.objectives[2];
            }
            for (auto &p: aux.nowpopulation) {
                if (p.rank == 1) {
                    m_ad4.set(p.c);
                    p.coords = m_ad4.get_heavy_atom_movable_coords();
                    add_to_output_container_(out_800, p, min_rmsd, num_saved_mins);

                }

            }

        }else if(g1==999){
            for (auto &s: aux.nowpopulation) {

                s.e = s.objectives[0] + s.objectives[1] + s.objectives[2];
            }
            for (auto &p: aux.nowpopulation) {
                if (p.rank == 1) {
                    m_ad4.set(p.c);
                    p.coords = m_ad4.get_heavy_atom_movable_coords();
                    add_to_output_container_(out_1000, p, min_rmsd, num_saved_mins);

                }

            }

        }
        bool compare_e(const output_type_MODEA& a, const output_type_MODEA& b);
        std::sort(aux.nowpopulation.begin(),aux.nowpopulation.end(),compare_e);

        iterations.push_back(g1); //   ʵ ʵĵ        滻
//        if(aux.nowpopulation[0].e<min_e){
            best_fitness.push_back(aux.nowpopulation[0].e); //   ʵ ʵ       Ӧ ȷ    滻
//            min_e = aux.nowpopulation[0].e;
//        }else{
//            best_fitness.push_back(min_e); //   ʵ ʵ       Ӧ ȷ    滻

//        }
//        if(aux.nowpopulation[0].objectives[0]<min_vina){
            best_fitness_vina.push_back(aux.nowpopulation[0].objectives[0]); //   ʵ ʵ       Ӧ ȷ    滻
//            min_vina = aux.nowpopulation[0].objectives[0];
//        }else{
//            best_fitness_vina.push_back(min_vina); //   ʵ ʵ       Ӧ ȷ    滻
//
//        }
//        if(aux.nowpopulation[0].objectives[1]<min_DLIGAND2){
            best_fitness_DLIGAND2.push_back(aux.nowpopulation[0].objectives[1]); //   ʵ ʵ       Ӧ ȷ    滻
//            min_DLIGAND2 = aux.nowpopulation[0].objectives[1];
//        }else{
//            best_fitness_DLIGAND2.push_back(min_DLIGAND2); //   ʵ ʵ       Ӧ ȷ    滻
//
//        }
//        if(aux.nowpopulation[0].objectives[2]<min_vdw){
            best_fitness_vdw.push_back(aux.nowpopulation[0].objectives[2]); //   ʵ ʵ       Ӧ ȷ    滻
//            min_vdw = aux.nowpopulation[0].objectives[2];
//        }else{
//            best_fitness_vdw.push_back(min_vdw); //   ʵ ʵ       Ӧ ȷ    滻
//
//        }
//        best_fitness_vina.push_back(aux.nowpopulation[0].objectives[0]);
//        best_fitness_DLIGAND2.push_back(aux.nowpopulation[0].objectives[1]);
//        best_fitness_vdw.push_back(aux.nowpopulation[0].objectives[2]);
    }
//    std::cout<<"deviceTime: "<<deviceTime<<std::endl;





    //   ͼ

//    matplotlibcpp::plot(iterations, best_fitness, "r-"); // 'r'   7

//
////    matplotlibcpp::xlabel("iterations");
////    matplotlibcpp::ylabel("best_fitness");
////    matplotlibcpp::title("convergence");
////    matplotlibcpp::save("../csv/convergence_graph_"+lig_name+".png"); //     ͼ
//    matplotlibcpp::show();
//    //   һ   ļ
//    std::ofstream out_file("../csv/results_"+lig_name+"_total.csv");
//    std::ofstream out_filevdw("../csv/results_"+lig_name+"_vdw.csv");
//    std::ofstream out_filevina("../csv/results_"+lig_name+"_vina.csv");
//    std::ofstream out_fileDLIGAND2("../csv/results_"+lig_name+"_DLIGAND2.csv");
////
////
//////     ļ  Ƿ ɹ
//    if (out_file.is_open()||out_filevdw.is_open()||out_filevina.is_open()||out_fileDLIGAND2.is_open()) {
//        // д      У   ѡ
//        out_file << "Iteration,Fitness Value\n";
//        out_filevdw << "Iteration,Fitness Value\n";
//        out_filevina << "Iteration,Fitness Value\n";
//        out_fileDLIGAND2 << "Iteration,Fitness Value\n";
//
//        // д
//        for (size_t i = 0; i < iterations.size(); ++i) {
//            out_file << iterations[i] << ','
//                     << std::fixed << std::setprecision(2) << best_fitness[i] << '\n';
//            out_filevdw << iterations[i] << ','
//                     << std::fixed << std::setprecision(2) << best_fitness_vdw[i] << '\n';
//            out_filevina << iterations[i] << ','
//                     << std::fixed << std::setprecision(2) << best_fitness_vina[i] << '\n';
//            out_fileDLIGAND2 << iterations[i] << ','
//                     << std::fixed << std::setprecision(2) << best_fitness_DLIGAND2[i] << '\n';
//        }
//
//        //  ر  ļ
//        out_file.close();
//        out_filevdw.close();
//        out_filevina.close();
//        out_fileDLIGAND2.close();
//
//    } else {
//        std::cerr << "Unable to open the file." << std::endl;
//    }
    for (auto &s: aux.nowpopulation) {
//        s.e = m_ad4.eval(p,ig,v,s.c);
//        s.e = s.objectives[1]/(1+0.05846*(m_ad4.ligand_degrees_of_freedom(0)));
        s.e = s.objectives[0] + s.objectives[1] + s.objectives[2];
    }

//    std::sort(aux.P.begin(),aux.P.end(),sort_comparator);
//    for (int i = 0; i < aux.nowpopulation.size() - 1; i++) {
//        for (int j = 0; j < aux.nowpopulation.size() - 1 - i; j++) {
//            if (aux.nowpopulation[j + 1].e < aux.nowpopulation[j].e) {
//                output_type_MODEA tmp = aux.nowpopulation[j];
//                aux.nowpopulation[j] = aux.nowpopulation[j + 1];
//                aux.nowpopulation[j + 1] = tmp;
//            }
//        }
//    }

    for (auto &p: aux.nowpopulation) {
        if (p.rank == 1) {
            m_ad4.set(p.c);
            p.coords = m_ad4.get_heavy_atom_movable_coords();
            add_to_output_container_(out, p, min_rmsd, num_saved_mins);

        }

    }

    }






