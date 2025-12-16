////
//// Created by 86137 on 2023/4/12.
////
//
//#include "nsga2.h"
//#include <cmath>
//#include <fstream>
//#include "coords.h"
//#include "mutate.h"
//#include"tee.h"
//#include"adam_cal.h"
//#include <sstream>
//
//
//const int eta = 1;
//
//const fl beta1 = 0.5;
//
//const fl beta2 = 0.999;
//const fl sigma = 0.00000001;
//const fl lr = 0.01;
//void MyNormalize_nsga2(vec& angles) {
//    while (angles[0] > 2 * pi) {
//        angles[0] -= 2 * pi;
//    }
//    while (angles[0] < 0) {
//        angles[0] += 2 * pi;
//    }
//    while (angles[1] > 2 * pi) {
//        angles[1] -= 2 * pi;
//    }
//    while (angles[1] < -2 * pi) {
//        angles[1] += 2 * pi;
//    }
//    if (angles[1] > pi && angles[1] < 2 * pi) {
//        angles[1] = 2 * pi - angles[1];
//    }
//    if (angles[1] > -2 * pi && angles[1] < -pi) {
//        angles[1] = 2 * pi + angles[1];
//    }
//    if (angles[1] > -pi && angles[1] < 0) {
//        angles[1] = -angles[1];
//    }
//    while (angles[2] > 2 * pi) {
//        angles[2] -= 2 * pi;
//    }
//    while (angles[2] < 0) {
//        angles[2] += 2 * pi;
//    }
//    assert(angles[0] > 0 && angles[0] < 2 * pi);
//    assert(angles[1] > 0 && angles[1] < pi);
//    assert(angles[2] > 0 && angles[2] < 2 * pi);
//
//}
////拥挤度距离比较器，用来根据拥挤距离对同一rank的个体进行排序
//bool crowding_distance_comparator(const output_type_nsga2& x, const output_type_nsga2& y) {
//    return x.crowding_distance > y.crowding_distance;
//}
//bool sort_comparator(const output_type_nsga2& x,const output_type_nsga2& y){
//    if(x.rank != y.rank){
//        if(x.rank < y.rank){
//            return true;
//        }else{
//            return false;
//        }
//    } else if(x.crowding_distance != y.crowding_distance){
//        if(x.crowding_distance > y.crowding_distance){
//            return true;
//        }else{
//            return false;
//        }
//    } else{
//        return true;
//    }
//}
//
//bool sort_comparator_ptr(const output_type_nsga2* x, const output_type_nsga2* y) {
//    if (x->rank != y->rank) {
//        if (x->rank < y->rank) {
//            return true;
//        }
//        else {
//            return false;
//        }
//    }
//    else if (x->crowding_distance != y->crowding_distance) {
//        if (x->crowding_distance > y->crowding_distance) {
//            return true;
//        }
//        else {
//            return false;
//        }
//    }
//    else {
//        return true;
//    }
//}
//void threeangles_to_quaternion1(const vec& three_angle, std::vector<double>& quater) {
//    if (fabs(three_angle[0] - 0.0) < 1e-6) {
//        quater[0] = 0;
//        quater[1] = 0;
//        quater[2] = 0;
//        quater[3] = 0;
//    }
//    else {
//        quater[0] = three_angle[0] / 2;//二分之theta角
//        quater[1] = sin(three_angle[1]) * cos(three_angle[2]);
//        quater[2] = sin(three_angle[1]) * sin(three_angle[2]);
//        quater[3] = cos(three_angle[1]);
//    }
//}
//
//void quaternion_to_3angles1(const output_type_nsga2& x, vec& three_angle) {
//    double theta = 2 * acos(x.c.ligands[0].rigid.orientation.R_component_1());
//    if (theta < 0.00000001) {
//        three_angle[0] = 0;
//        three_angle[1] = 0;
//        three_angle[2] = 0;
//    }
//    else {
//        vec quaternion_v;//v_x,v_y,v_z
//        quaternion_v[0] = x.c.ligands[0].rigid.orientation.R_component_2() / sin(theta / 2); //v_x
//        quaternion_v[1] = x.c.ligands[0].rigid.orientation.R_component_3() / sin(theta / 2); //v_y
//        quaternion_v[2] = x.c.ligands[0].rigid.orientation.R_component_4() / sin(theta / 2); //v_z
//        fl alpha = acos(quaternion_v[2]);
//        fl beta = 0;
//        fl sin_beta = quaternion_v[1] / sin(alpha);
//        fl cos_bata = quaternion_v[0] / sin(alpha);
//        if (sin_beta > 0)
//            beta = acos(cos_bata);
//        else if (sin_beta < 0)
//            beta = 2 * pi - acos(cos_bata);
//        three_angle[0] = theta;
//        three_angle[1] = alpha;
//        three_angle[2] = beta;
//    }
//    assert(three_angle[0] >= 0 && three_angle[0] <= 2 * pi);
//    assert(three_angle[1] >= 0 && three_angle[1] <= pi);
//    assert(three_angle[2] >= 0 && three_angle[2] <= 2 * pi);
//}
//
//
//void increment_position(conf& c, change& g, fl beta1, fl beta2, vec& V_d_position, vec& S_d_position, fl lr, fl sigma, unsigned t) {
//    vec d_position = g.ligands[0].rigid.position;
//
//    V_d_position = beta1 * V_d_position + (1 - beta1) * d_position;
//    S_d_position = beta2 * S_d_position + (1 - beta2) * elementwise_product(d_position, d_position);
//
//    vec V_d_position_corr = V_d_position / (1 - pow(beta1, t));
//    vec S_d_position_corr = S_d_position / (1 - pow(beta2, t));
//    vec_sqrt(S_d_position_corr);
//    S_d_position_corr += sigma;
//
//    c.ligands[0].rigid.position = c.ligands[0].rigid.position - lr * vec_division(V_d_position_corr, S_d_position_corr);
//}
//
//
//void increment_orientation(conf& c, change& g, fl beta1, fl beta2, vec& V_d_orientation, vec& S_d_orientation, fl lr, fl sigma, unsigned t) {
//    vec d_orientation = g.ligands[0].rigid.orientation;
//
//    V_d_orientation = beta1 * V_d_orientation + (1 - beta1) * d_orientation;
//    S_d_orientation = beta2 * S_d_orientation + (1 - beta2) * elementwise_product(d_orientation, d_orientation);
//
//    vec V_d_orientation_corr = V_d_orientation / (1 - pow(beta1, t));
//    vec S_d_orientation_corr = S_d_orientation / (1 - pow(beta2, t));
//    vec_sqrt(S_d_orientation_corr);
//    S_d_orientation_corr += sigma;
//
//    vec rotation = -lr * vec_division(V_d_orientation_corr, S_d_orientation_corr);
//    quaternion_increment(c.ligands[0].rigid.orientation, rotation);
//}
//
//
//void increment_torsions(conf& c, change& g, fl beta1, fl beta2, flv& V_d_torsions, flv& S_d_torsions, fl lr, fl sigma, unsigned t) {
//    flv d_torsions = g.ligands[0].torsions;
//
//    V_d_torsions = torsions_add_torsions(s_mul_torsions(beta1, V_d_torsions), s_mul_torsions(1 - beta1, d_torsions));
//    S_d_torsions = torsions_add_torsions(s_mul_torsions(beta2, S_d_torsions), s_mul_torsions(1 - beta2, torsions_sqr(d_torsions)));
//
//    flv V_d_torsions_corr = torsions_div_s(V_d_torsions, 1 - pow(beta1, t));
//    flv S_d_torsions_corr = torsions_div_s(S_d_torsions, 1 - pow(beta2, t));
//    torsions_sqrt(S_d_torsions_corr);
//    torsions_add_s(S_d_torsions_corr, sigma);
//
//    flv torsions_ = s_mul_torsions(lr, torsions_div_torsions(V_d_torsions_corr, S_d_torsions_corr));
//
//    for (sz i = 0; i < torsions_.size(); i++) {
//        c.ligands[0].torsions[i] -= normalized_angle(torsions_[i]);
//        normalize_angle(c.ligands[0].torsions[i]);
//    }
//}
//
//
//std::unordered_map<int, std::vector<output_type_nsga2>> nsga2_aux::fast_nondominated_sort(std::vector<output_type_nsga2>& pop) {
//    std::unordered_map<int, std::vector<output_type_nsga2>> F;
//    for (auto& p : pop) {
//        //p.Slaves.clear();
//        p.Slaves.clear();
//
//        p.dominated_count = 0;
//        for (auto& q : pop) {
//            if (p < q) {
//                // 如果当前个体p支配q,那么将被支配的个体q放入Sp中
//                p.Slaves.push_back(&q);
//            }
//            else if (q < p) {
//                //否则如果q支配p,那么将np的值加一
//                p.dominated_count += 1;
//            }
//        }
//        //选择出初始种群的rank1的个体，并且将rank1中的个体放入F[1]中
//        if (p.dominated_count == 0) {
//            p.rank = 1;
//            F[1].push_back(p);
//        }
//
//    }
//
//    int i = 1;
//    std::vector<output_type_nsga2> Q;
//    while (!F[i].empty()) {
//        Q.clear();
//        //循环当前F[i]的帕累托最优解集，即rankn中的解集
//        for (auto& p : F[i]) {
//            //循环当前帕累托最优解集中的Sp集合，即当前解的所支配的解
//            for (auto& q : p.Slaves) {
//                //将Sp中的每个个体的np都-1，即支配当前个体Sp里的个体数量-1
//                q->dominated_count = q->dominated_count - 1;
//                // 判断Sp个体的np是否为0,如果为0,则表示当前Sp中的该个个体不被任何其他个体所支配
//                // 将当前Sp中np为0的解放入rank+1中，即放入Q中，成为非支配解
//                if (q->dominated_count == 0) {
//                    q->rank = i + 1;
//                    Q.push_back(*q);
//                }
//            }
//        }
//        i = i + 1;
//        F[i] = Q;
//    }
//    return F;
//}
//
//
//void nsga2_aux::calculate_objectives(model& m, const precalculate& p, const igrid& ig, vec v, output_type_nsga2& indivisual, const scoring_function& sf, const scoring_function& sf_ad4,energy_cal& energy,const boost::optional<model>& ref,change &g){
//    conf_size s = m.get_size();
////    change g(s);
//    fl intramolecular_energy;
//    for(int k = 0; k<indivisual.objectives.size();k++){
//        if(k==0){
////                tmp.objectives[k] = energy.energy_cal_force(tmp.c);
////            intramolecular_energy = energy.energy_cal_intra(indivisual.c,g);
////            intramolecular_energy = energy.AD4_energy_cal_intra(indivisual.c);
////            indivisual.objectives[k] = energy.energy_cal_forcefield(indivisual.c);
////            indivisual.objectives[k] = intramolecular_energy;
//            indivisual.objectives[k] = energy.energy_cal_autodock4(sf_ad4,indivisual.c);
//        }else if(k==1){
////                tmp.objectives[k] = energy.energy_cal_hydrophobic(tmp.c);
////            indivisual.objectives[k] = energy.energy_cal_inter(sf, indivisual.c, intramolecular_energy,g);
////            indivisual.objectives[k] = energy.AD4_energy_cal_inter(sf, indivisual.c, intramolecular_energy);
//            float bind_score = energy.energy_cal_xscore(&indivisual);
//            indivisual.objectives[k] = bind_score*(-1.364);
//        }else if(k==2){
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
////            indivisual.objectives[k] = energy.energy_cal_rmsd(ref, indivisual.c);
//            indivisual.objectives[k] = energy.energy_cal_Smog2016(indivisual);
//        }
//    }
//
//}
//void nsga2_aux::calculate_crowding_distance(std::vector<output_type_nsga2>& rank_i_vector) {
//    int length = rank_i_vector.size();
//    for (int i = 0; i < length; i++) {
//        rank_i_vector[i].crowding_distance = 0.0;
//    }
//    //关于1：当引用作为形参，函数调用时也可以看成将传递的实参绑定给它，这样我们在函数体内对这个引用做的一切操作都有可能影响到函数传递的实参。如果我们希望参数在函数体内是只读的,所以当我们加了引用有希望参数是只读的就必须加 const。
//    //为什么不直接值传递呢？ 确实，但是当参数是类对象时值传递就有了一个问题，那就是性能可能会大受影响。我们知道值传递实际就是向函数拷贝一份副本来使用，那么对于一些复杂的类，尤其是 string 这样每一次拷贝可能消耗很多的时间，那么通过引用传参就很有必要了。
//    //总的来说因为我想提高类对象传参时的性能，所以要用引用，因为用了引用我又希望它只读所以我用了const。
//    for (int i = 0; i < rank_i_vector[0].objectives.size(); i++) {
//        std::sort(rank_i_vector.begin(), rank_i_vector.end(), [i](const output_type_nsga2& a, const output_type_nsga2& b) {
//            return a.objectives[i] < b.objectives[i];
//            });
//        rank_i_vector[0].crowding_distance = rank_i_vector[length - 1].crowding_distance = std::numeric_limits<fl>::infinity();
//        fl fmax = rank_i_vector[length - 1].objectives[i];
//        fl fmin = rank_i_vector[0].objectives[i];
//        if (fmax != fmin) {
//            for (int j = 1; j < length - 1; j++) {
//                rank_i_vector[j].crowding_distance += (rank_i_vector[j + 1].objectives[i] - rank_i_vector[j - 1].objectives[i]) / (fmax - fmin);
//            }
//        }
//
//    }
//
//}
//void nsga2_aux::init_individual(model& m, energy_cal& energy, rng& generator, const vec& corner1, const vec& corner2, const scoring_function& sf, output_type_nsga2& tmp,const boost::optional<model>& ref,const scoring_function& sf_ad4) {
//    conf_size s = m.get_size();
//    change g(s);
//
//    tmp.c.randomize(corner1, corner2, generator);
//
//    for (int j = 0; j < 3; j++) {
//        if (j == 1) {
//            tmp.rotor_angle[j] = random_fl(0, pi, generator);
//        }
//        else {
//            tmp.rotor_angle[j] = random_fl(0, 2 * pi, generator);
//        }
//    }
//    //这里要将初始化好的三个角度转换为四元数，然后再将这个转换好的随机初始化四元数赋给tmp的orientation
//    std::vector<double> quater;
//    quater.resize(4);
//    threeangles_to_quaternion1(tmp.rotor_angle, quater);
//    fl cos_theta = cos(quater[0]);
//    fl sin_theta = sin(quater[0]);
//    qt q(cos_theta, quater[1] * sin_theta, quater[2] * sin_theta, quater[3] * sin_theta);
//    quaternion_normalize(q);//归一化
//    tmp.c.ligands[0].rigid.orientation = q;
//    //    tmp.e = energy(tmp.c, g);
////    tmp.rm = 0.0;
//    fl intra;
//    for(int k = 0; k<tmp.objectives.size();k++){
//        if(k==0){
////                intra = energy.energy_cal_intra(tmp.c,g);
////                tmp.objectives[k] = energy.energy_cal_force(tmp.c);
////                intra = energy.AD4_energy_cal_intra(tmp.c);
//            //tmp.objectives[k] = intra;
////                tmp.objectives[k] = energy.energy_cal_forcefield(tmp.c);
//                tmp.objectives[k] = energy.energy_cal_autodock4(sf_ad4,tmp.c);
////            tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4,tmp.c);
//
//        }else if(k==1){
////                tmp.objectives[k] = energy.energy_cal_hydrophobic(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_inter(sf, tmp.c, intra,g)/(1+0.05846*(m.ligand_degrees_of_freedom(0)));
////                tmp.objectives[k] = energy.AD4_energy_cal_inter(sf, tmp.c, intra);
//            float bind_score = energy.energy_cal_xscore(&tmp);
//            tmp.objectives[k] = bind_score*(-1.364);
//
//        }else if(k==2){
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
//            tmp.objectives[k] = energy.energy_cal_Smog2016(tmp);
//
//        }else if(k==3){
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
////                tmp.objectives[k] = energy.energy_cal_Smog2016(tmp);
//            tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4,tmp.c);
//
//        }
//    }
//}
//std::vector<output_type_nsga2> nsga2_aux::crossover_mutation(model& m, const precalculate& p, const igrid& ig, vec v, energy_cal& energy, rng& generator, const vec& corner1, const vec& corner2,
//    output_type_nsga2 parent1, output_type_nsga2 parent2, const scoring_function& sf,const boost::optional<model>& ref,change &g,const scoring_function& sf_ad4) {
//    conf_size s = m.get_size();
//    std::vector<fl> obj(3, 0.0);
//    output_type_nsga2 kid1(s, obj);
//    output_type_nsga2 kid2(s, obj);
//    kid1 = init_tmp;
//    kid2 = init_tmp;
//    init_individual(m, energy, generator, corner1, corner2, sf, kid1,ref,sf_ad4);
//    init_individual(m, energy, generator, corner1, corner2, sf, kid2,ref,sf_ad4);
//    fl beta;
//    std::vector<output_type_nsga2> temp;
//    for (int i = 0; i < dim; i++) {
//        fl p_cr = random_fl(0, 1.0, generator);
//        if (p_cr < 0.9) {
//            beta = pow((2 * p_cr), (1 / (1.0 + eta)));
//        }
//        else {
//            beta = pow((1 / (2 - p_cr * 2)), (1 / (1.0 + eta)));
//        }
//        if (i < 3) {
//            kid1.c.ligands[0].rigid.position[i] = 0.5 * ((1 + beta) * parent1.c.ligands[0].rigid.position[i] + (1 - beta) * parent2.c.ligands[0].rigid.position[i]);
//            kid2.c.ligands[0].rigid.position[i] = 0.5 * ((1 - beta) * parent1.c.ligands[0].rigid.position[i] + (1 + beta) * parent2.c.ligands[0].rigid.position[i]);
//            if (kid1.c.ligands[0].rigid.position[i] < corner1[i]) {
//                kid1.c.ligands[0].rigid.position[i] = corner1[i];
//            }
//            else if (kid1.c.ligands[0].rigid.position[i] > corner2[i]) {
//                kid1.c.ligands[0].rigid.position[i] = corner2[i];
//            }
//
//            if (kid2.c.ligands[0].rigid.position[i] < corner1[i]) {
//                kid2.c.ligands[0].rigid.position[i] = corner1[i];
//            }
//            else if (kid2.c.ligands[0].rigid.position[i] > corner2[i]) {
//                kid2.c.ligands[0].rigid.position[i] = corner2[i];
//            }
//        }
//        else if (i >= 3 && i < 6) {
//            kid1.rotor_angle[i - 3] = 0.5 * ((1 + beta) * parent1.rotor_angle[i - 3] + (1 - beta) * parent2.rotor_angle[i - 3]);
//            kid2.rotor_angle[i - 3] = 0.5 * ((1 - beta) * parent1.rotor_angle[i - 3] + (1 + beta) * parent2.rotor_angle[i - 3]);
//
//        }
//        else {
//            kid1.c.ligands[0].torsions[i - 6] = 0.5 * ((1 + beta) * parent1.c.ligands[0].torsions[i - 6] + (1 - beta) * parent2.c.ligands[0].torsions[i - 6]);
//            kid2.c.ligands[0].torsions[i - 6] = 0.5 * ((1 - beta) * parent1.c.ligands[0].torsions[i - 6] + (1 + beta) * parent2.c.ligands[0].torsions[i - 6]);
//            normalize_angle(kid1.c.ligands[0].torsions[i - 6]);
//            normalize_angle(kid2.c.ligands[0].torsions[i - 6]);
//
//        }
//
//    }
//
//    fl delta;
//    for (int j = 0; j < dim; j++) {
//        fl p_mut = random_fl(0, 1.0, generator);
//        if (p_mut < 1.0/(dim*5)) {
//            delta = pow((2 * p_mut), (1 / (eta + 1.0)));
//        }
//        else {
//            delta = 1 - pow((2 * (1 - p_mut)), (1 / (eta + 1.0)));
//        }
//        if (j < 3) {
//            kid1.c.ligands[0].rigid.position[j] = kid1.c.ligands[0].rigid.position[j] + delta*(corner2[j] - corner1[j]);
//
//            if (kid1.c.ligands[0].rigid.position[j] < corner1[j]) {
//                kid1.c.ligands[0].rigid.position[j] = corner1[j];
//            }
//            else if (kid1.c.ligands[0].rigid.position[j] > corner2[j]) {
//                kid1.c.ligands[0].rigid.position[j] = corner2[j];
//            }
//
//
//        }
//        else if (j >= 3 && j < 6) {
//
//            if(j==3){
//                kid1.rotor_angle[j - 3] = kid1.rotor_angle[j - 3] + delta*2*pi;
//            }else if(j==4){
//                kid1.rotor_angle[j - 3] = kid1.rotor_angle[j - 3] + delta*pi;
//            }else if(j==5){
//                kid1.rotor_angle[j - 3] = kid1.rotor_angle[j - 3] + delta*2*pi;
//            }
//
//        }
//        else {
//            kid1.c.ligands[0].torsions[j - 6] = kid1.c.ligands[0].torsions[j - 6] + delta*2*pi;
//
//            normalize_angle(kid1.c.ligands[0].torsions[j - 6]);
//
//        }
//
//    }
//    std::vector<double> quater1; // 初始化四元数变量
//    std::vector<double> quater2; // 初始化四元数变量
//
//    quater1.resize(4);
//    quater2.resize(4);
//
//    MyNormalize_nsga2(kid1.rotor_angle);
//    MyNormalize_nsga2(kid2.rotor_angle);
//    threeangles_to_quaternion1(kid1.rotor_angle, quater1);
//    threeangles_to_quaternion1(kid2.rotor_angle, quater2);
//
//    fl cos_theta1 = cos(quater1[0]);
//    fl sin_theta1 = sin(quater1[0]);
//    qt q1(cos_theta1, quater1[1] * sin_theta1, quater1[2] * sin_theta1, quater1[3] * sin_theta1);
////    quaternion_normalize(q1);//归一化
//
//    fl cos_theta2 = cos(quater2[0]);
//    fl sin_theta2 = sin(quater2[0]);
//    qt q2(cos_theta2, quater2[1] * sin_theta2, quater2[2] * sin_theta2, quater2[3] * sin_theta2);
////    quaternion_normalize(q2);//归一化
//
//    kid1.c.ligands[0].rigid.orientation = q1;
//    kid2.c.ligands[0].rigid.orientation = q2;
//
//    calculate_objectives(m, p, ig, v, kid1, sf,sf_ad4,energy,ref,g);
//    calculate_objectives(m, p, ig, v, kid2, sf,sf_ad4,energy,ref,g);
//
//    temp.push_back(kid1);
//    temp.push_back(kid2);
//    return temp;
//}
//
//output_type_nsga2 nsga2_aux::binary_tournament(output_type_nsga2 indivisual1, output_type_nsga2 indivisual2) {
//    if (indivisual1.rank != indivisual2.rank) {
//        if (indivisual1.rank < indivisual2.rank) {
//            return indivisual1;
//        }
//        else {
//            return indivisual2;
//        }
//    }
//    else if (indivisual1.crowding_distance != indivisual2.crowding_distance) {
//        if (indivisual1.crowding_distance > indivisual2.crowding_distance) {
//            return indivisual1;
//        }
//        else {
//            return indivisual2;
//        }
//    }
//    else {
//        return indivisual1;
//    }
//}
//
//std::vector<output_type_nsga2> nsga2_aux::make_new_pop(model& m, const precalculate& p, const igrid& ig, vec v, std::vector<output_type_nsga2> pop, energy_cal& energy, const vec& corner1, const vec& corner2, rng& generator, const scoring_function& sf,const boost::optional<model>& ref,change &g,const scoring_function& sf_ad4) {
////    std::vector<output_type_nsga2> empty;
////    Q.swap(empty);
//    std::vector<output_type_nsga2> Q;
//    for (int i = 0; i < popsize / 2; i++) {
//        fl r1 = random_int(0, popsize-1, generator);
//        fl r2 = random_int(0, popsize-1, generator);
//        output_type_nsga2 parent1 = binary_tournament(pop[r1], pop[r2]);
//
//        r1 = random_int(0, popsize-1, generator);
//        r2 = random_int(0, popsize-1, generator);
//        output_type_nsga2 parent2 = binary_tournament(pop[r1], pop[r2]);
//
//        while((parent1.c==parent2.c)){
//            r1 = random_int(0, popsize-1, generator);
//            r2 = random_int(0, popsize-1, generator);
//            parent2 = binary_tournament(pop[r1], pop[r2]);
//        }
//        std::vector<output_type_nsga2> kids = crossover_mutation(m, p, ig, v, energy, generator, corner1, corner2, parent1, parent2, sf,ref,g,sf_ad4);
//        Q.insert(Q.end(), kids.begin(), kids.end());
//
//
//    }
//    return Q;
//}
//
//void nsga2::init_individual1(model& m, energy_cal& energy, rng& generator, const vec& corner1, const vec& corner2, const scoring_function& sf, output_type_nsga2& tmp) const {
//    tmp.c.randomize(corner1, corner2, generator);
//
//    for (int j = 0; j < 3; j++) {
//        if (j == 1) {
//            tmp.rotor_angle[j] = random_fl(0, pi, generator);
//        }
//        else {
//            tmp.rotor_angle[j] = random_fl(0, 2 * pi, generator);
//        }
//    }
//    //这里要将初始化好的三个角度转换为四元数，然后再将这个转换好的随机初始化四元数赋给tmp的orientation
//    std::vector<double> quater;
//    quater.resize(4);
//    threeangles_to_quaternion1(tmp.rotor_angle, quater);
//    fl cos_theta = cos(quater[0]);
//    fl sin_theta = sin(quater[0]);
//    qt q(cos_theta, quater[1] * sin_theta, quater[2] * sin_theta, quater[3] * sin_theta);
//    quaternion_normalize(q);//归一化
//    tmp.c.ligands[0].rigid.orientation = q;
//    //    tmp.e = energy(tmp.c, g);
////    fl intra;
////    for(int k = 0; k<tmp.objectives.size();k++){
////        if(k==0){
////            intra = energy.energy_cal_intra(tmp.c);
//////                tmp.objectives[k] = energy.energy_cal_force(tmp.c);
//////                intra = energy.AD4_energy_cal_intra(tmp.c);
////            tmp.objectives[k] = intra;
////        }else if(k==1){
//////                tmp.objectives[k] = energy.energy_cal_hydrophobic(tmp.c);
////            tmp.objectives[k] = energy.energy_cal_inter(sf, tmp.c, intra);
//////                tmp.objectives[k] = energy.AD4_energy_cal_inter(sf, tmp.c, intra);
////        }else if(k==2){
//////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
////            tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
////        }
////    }
//}
//
//bool areAlmostEqual(float a, float b, float tolerance) {
//    return std::fabs(a - b) < tolerance;
//}
//void nsga2::operator()(model& m, model& m_ad4, output_container_nsga2& out, const precalculate& p, const precalculate& prec_ad4, const precalculate& prec_vdw, const precalculate& prec_forcefield, const igrid& ig,
//                       const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2,
//                       incrementable* increment_me, rng& generator, const scoring_function& sf,const scoring_function& sf_ad4,const boost::optional<model>& ref) const {
//    vec v(10, 10, 10);
//    conf_size s = m.get_size();
//    energy_cal energy(&m, &m_ad4, &p, &prec_ad4, &prec_vdw, &prec_forcefield,&ig, v);
//    nsga2_aux aux(m);
//    change g(s);
//    std::vector<fl> obj(3, 0.0);
//    extern Ligand * liganda;
//    extern OBMolecule obl;
//    output_type_nsga2 tmp(s, obj, obl.obatom ,*liganda);
//
////    output_type_nsga2 ns(s,obj);//用来作为最佳能量个体对比
////    init_individual1(m, energy, generator, corner1, corner2, sf, ns);
////    output_type_nsga2 tmp(s, obj);
//    atomv m_atoms = m.get_ligand_atoms();
//    int count = 0;
//    VINA_FOR(i, m_atoms.size()){
//        VINA_FOR(j, tmp.num_atom){
//
//            if(areAlmostEqual(m.get_ligand_coords()[i][0], tmp.atom[j].coor[0], 0.01) && areAlmostEqual(m.get_ligand_coords()[i][1], tmp.atom[j].coor[1], 0.01) && areAlmostEqual(m.get_ligand_coords()[i][2], tmp.atom[j].coor[2], 0.01)){
//                tmp.atom[j].ad_number = m_atoms[i].ad_number;
//                count++;
//                break;
//            }
//        }
//    }
//    aux.init_tmp = tmp;
//    for (int i = 0; i < aux.popsize; i++) {
//        tmp.c.randomize(corner1, corner2, generator);
//
//        for (int j = 0; j < 3; j++) {
//            if (j == 1) {
//                tmp.rotor_angle[j] = random_fl(0, pi, generator);
//            } else {
//                tmp.rotor_angle[j] = random_fl(0, 2 * pi, generator);
//            }
//        }
//        //这里要将初始化好的三个角度转换为四元数，然后再将这个转换好的随机初始化四元数赋给tmp的orientation
//        std::vector<double> quater;
//        quater.resize(4);
//        threeangles_to_quaternion1(tmp.rotor_angle, quater);
//        fl cos_theta = cos(quater[0]);
//        fl sin_theta = sin(quater[0]);
//        qt q(cos_theta, quater[1] * sin_theta, quater[2] * sin_theta, quater[3] * sin_theta);
//        quaternion_normalize(q);//归一化
//        tmp.c.ligands[0].rigid.orientation = q;
//        tmp.e = 0;
//        tmp.rm = 0.0;
//        fl intra;
//        for (int k = 0; k < tmp.objectives.size(); k++) {
//            if (k == 0) {
////                intra = energy.energy_cal_intra(tmp.c,g);
////                tmp.objectives[k] = energy.energy_cal_force(tmp.c);
////                intra = energy.AD4_energy_cal_intra(tmp.c);
//                //tmp.objectives[k] = intra;
////                tmp.objectives[k] = energy.energy_cal_forcefield(tmp.c);
//                tmp.objectives[k] = energy.energy_cal_autodock4(sf_ad4,tmp.c);
////                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4, tmp.c);
//
//            } else if (k == 1) {
////                tmp.objectives[k] = energy.energy_cal_hydrophobic(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_inter(sf, tmp.c, intra,g)/(1+0.05846*(m.ligand_degrees_of_freedom(0)));
////                tmp.objectives[k] = energy.AD4_energy_cal_inter(sf, tmp.c, intra);
//                float bind_score = energy.energy_cal_xscore(&tmp);
//                tmp.objectives[k] = bind_score * (-1.364);
//
//            } else if (k == 2) {
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
//                tmp.objectives[k] = energy.energy_cal_Smog2016(tmp);
//
//            } else if (k == 3) {
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
////                tmp.objectives[k] = energy.energy_cal_Smog2016(tmp);
//                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4, tmp.c);
//
//            }
//
//            aux.P[i] = tmp;
//
//        }
//    }
//
//
//    // 初始化的第一代种群需要先进行非支配排序，然后选择、交叉、变异
//        aux.fast_nondominated_sort(aux.P);
//
//        // 生成的初始化种群的子代Q_t
//
//        aux.Q = aux.make_new_pop(m, p, ig, v, aux.P, energy, corner1, corner2, generator, sf, ref, g, sf_ad4);
//
//        for (int g1 = 0; g1 < num_steps; g1++) {
//            if (increment_me)
//                ++(*increment_me);
//            aux.R.clear();
//            //将第一代的父子种群合并
//            aux.R.insert(aux.R.end(), aux.P.begin(), aux.P.end());
//            aux.R.insert(aux.R.end(), aux.Q.begin(), aux.Q.end());
//
//            //对合并后的父子种群进行快速非支配排序
//            std::unordered_map<int, std::vector<output_type_nsga2>> F = aux.fast_nondominated_sort(aux.R);
////        std::vector<output_type_nsga2> as = F[1];
//            //新的父代种群Pt+1
//            std::vector<output_type_nsga2> P_t_1;
//            int i = 1;
//            //如果新的父代种群大小大于N，那么跳出循环
//            while (P_t_1.size() + F[i].size() < aux.P.size()) {
//                aux.calculate_crowding_distance(F[i]);
//                P_t_1.insert(P_t_1.end(), F[i].begin(), F[i].end());
//                i = i + 1;
//            }
//            aux.calculate_crowding_distance(F[i]);
//            std::sort(F[i].begin(), F[i].end(), crowding_distance_comparator);
//            int left = aux.P.size() - P_t_1.size();
//            //P_t_1.insert(P_t_1.end(), F[i].begin(), F[i].begin() + aux.P.size() - P_t_1.size());
//            for (int a = 0; a < left; a++) {
//                P_t_1.push_back(F[i][a]);
//
//            }
//            aux.Q = aux.make_new_pop(m, p, ig, v, P_t_1, energy, corner1, corner2, generator, sf, ref, g, sf_ad4);
//
//            aux.P = P_t_1;
//
//
////            int t = 1;
////            for (auto &s: aux.P) {
//////        s.e = m.eval(p,ig,v,s.c);
//////        s.e = s.objectives[1]/(1+0.05846*(m.ligand_degrees_of_freedom(0)));
////                s.e = s.objectives[0] + s.objectives[1]+s.objectives[2];
////            }
//
////    std::sort(aux.P.begin(),aux.P.end(),sort_comparator);
////            for (int i = 0; i < aux.P.size() - 1; i++) {
////                for (int j = 0; j < aux.P.size() - 1 - i; j++) {
////                    if (aux.P[j + 1].e < aux.P[j].e) {
////                        output_type_nsga2 tmp = aux.P[j];
////                        aux.P[j] = aux.P[j + 1];
////                        aux.P[j + 1] = tmp;
////                    }
////                }
////            }
//            //std::cout << "energy:" << std::setw(15) << energy(mutatepopulation[best_index].c, g) << std::setw(15);
////            for (int k = 0; k < 20; k++) {
////                for (int a = 0; a < 1; a++) {
//////                aux.calculate_objectives(m, p, ig, v, aux.P[a], sf,energy,ref,g);
////                    fl intramolecular_energy;
////                    for (int b = 0; b < aux.P[a].objectives.size(); b++) {
////                        if (b == 0) {
////                            sz total_torsions = m.ligand_degrees_of_freedom(0);//旋转键的个数
////                            conf_size s = m.get_size();
////                            change g(s);
////
////                            //std::cout << "iter " << ++iters << std::endl;
////                            vec V_d_position(zero_vec);
////                            vec S_d_position(zero_vec);
////                            vec V_d_orientation(zero_vec);
////                            vec S_d_orientation(zero_vec);
////                            flv V_d_torsions(total_torsions, 0);
////                            flv S_d_torsions(total_torsions, 0);
////                            intramolecular_energy = energy.energy_cal_intra(aux.P[a].c, g);
////                            aux.P[a].objectives[b] = intramolecular_energy;
////                            increment_position(aux.P[a].c, g, beta1, beta2, V_d_position, S_d_position, lr, sigma, t);
////                            increment_orientation(aux.P[a].c, g, beta1, beta2, V_d_orientation, S_d_orientation, lr,
////                                                  sigma, t);
////                            increment_torsions(aux.P[a].c, g, beta1, beta2, V_d_torsions, S_d_torsions, lr, sigma, t);
////                            t++;
////                        } else if (b == 1) {
////                            sz total_torsions = m.ligand_degrees_of_freedom(0);//旋转键的个数
////                            conf_size s = m.get_size();
////                            change g(s);
////
////                            //std::cout << "iter " << ++iters << std::endl;
////                            vec V_d_position(zero_vec);
////                            vec S_d_position(zero_vec);
////                            vec V_d_orientation(zero_vec);
////                            vec S_d_orientation(zero_vec);
////                            flv V_d_torsions(total_torsions, 0);
////                            flv S_d_torsions(total_torsions, 0);
////                            aux.P[a].objectives[b] = energy.energy_cal_inter(sf, aux.P[a].c, intramolecular_energy, g);
////                            increment_position(aux.P[a].c, g, beta1, beta2, V_d_position, S_d_position, lr, sigma, t);
////                            increment_orientation(aux.P[a].c, g, beta1, beta2, V_d_orientation, S_d_orientation, lr,
////                                                  sigma, t);
////                            increment_torsions(aux.P[a].c, g, beta1, beta2, V_d_torsions, S_d_torsions, lr, sigma, t);
////                            t++;
////                        } else if (b == 2) {
////                            aux.P[a].objectives[b] = energy.energy_cal_rmsd(ref, aux.P[a].c);
////                        }
////                    }
//////            mutatepopulation[best_index].e = energy(mutatepopulation[best_index].c, g);
////
////                }
////
////            }
////        aux.calculate_objectives(m, p, ig, v, aux.P[0], sf,energy,ref,g);
////            quaternion_to_3angles1(aux.P[0], aux.P[0].rotor_angle);
////        if(sort_comparator(aux.P[0], ns) || out.size() < num_saved_mins){
////            m.set(aux.P[0].c);
////            aux.P[0].coords = m.get_heavy_atom_movable_coords();
////            add_to_output_container(out,aux.P[0],min_rmsd,num_saved_mins);
////            if(sort_comparator(aux.P[0],ns)){
//
////                ns = aux.P[0];
////            }
////        }
//        }
//        for (auto &s: aux.P) {
////        s.e = m.eval(p,ig,v,s.c);
////        s.e = s.objectives[1]/(1+0.05846*(m.ligand_degrees_of_freedom(0)));
//            s.e = s.objectives[0] + s.objectives[1]+s.objectives[2];
//        }
//
////    std::sort(aux.P.begin(),aux.P.end(),sort_comparator);
//        for (int i = 0; i < aux.P.size() - 1; i++) {
//            for (int j = 0; j < aux.P.size() - 1 - i; j++) {
//                if (aux.P[j + 1].e < aux.P[j].e) {
//                    output_type_nsga2 tmp = aux.P[j];
//                    aux.P[j] = aux.P[j + 1];
//                    aux.P[j + 1] = tmp;
//                }
//            }
//        }
//
//        for (auto &p: aux.P) {
//            if (p.rank == 1) {
//
//
//                m_ad4.set(p.c);
//                p.coords = m.get_heavy_atom_movable_coords();
//                add_to_output_container(out, p, min_rmsd, num_saved_mins);
//
//            }
//
//        }
//
//}