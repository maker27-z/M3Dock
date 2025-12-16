////
//// Created by 91686 on 2023/7/6.
////
//
//#include "MOEAD.h"
//#include "coords.h"
//void threeangles_to_quaternion(const vec& three_angle, std::vector<double>& quater) {
//    if (fabs(three_angle[0] - 0.0) < 1e-6) {
//        quater[0] = 0;
//        quater[1] = 0;
//        quater[2] = 0;
//        quater[3] = 0;
//    }
//    else {
//        quater[0] = three_angle[0] / 2;//?????theta??
//        quater[1] = sin(three_angle[1]) * cos(three_angle[2]);
//        quater[2] = sin(three_angle[1]) * sin(three_angle[2]);
//        quater[3] = cos(three_angle[1]);
//    }
//}
//
//
//void MyNormalize(vec& angles) {
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
//
//
//void MOEAD_aux::initUniformWeight() {
//    if((pop[0].objectives.size() == 2) && (popsize <= 300)){
//        for(int n = 0; n < popsize; n++){
//            lamda[n] = new fl[pop[0].objectives.size()];
//            fl a = 1.0 * n / (popsize - 1);
//            lamda[n][0] = a;
//            lamda[n][1] = 1 - a;
//        }
//    }
//    else{
//        std::ostringstream os;
//        os << "../" << "W" << pop[0].objectives.size()<< "D_"
//           << popsize << ".dat";
//        std::string dataFileName;
//        dataFileName = os.str();
//
//    // Open the file
//        std::ifstream in(dataFileName.c_str());
//        if( !in ) {
//            std::cout << "initUniformWeight: failed when reading from file: : " <<
//                 dataFileName << std::endl;
//            exit(-1);
//        } // if
//
//        //int numberOfObjectives = 0;
//        int i = 0;
//        int j = 0;
//        std::string aux;
//        while (getline(in, aux)) {
//            std::istringstream iss(aux);
//            j = 0;
//            // TODO: Check if number of tokens per line is equal to number of
//            //       objectives
//            lamda[i] = new double[pop[0].objectives.size()];
//            while (iss) {
//                std::string token;
//                iss >> token;
//                if (token.compare("")!=0) {
//                    double value = atof(token.c_str());
//                    lamda[i][j] = value;
//                    //cout << "lambda[" << i << "," << j << "] = " << value << endl;
//                    j++;
//                } // if
//            } // while
//            i++;
//        } // while
//        in.close();
//
//    }
//}
//
//void MOEAD_aux::initPopulation(model& m,  const precalculate& p, const precalculate& prec_forcefield, const igrid& ig,
//                               const precalculate& p_widened,  const vec& corner1, const vec& corner2,
//                               rng& generator,const scoring_function& sf,const boost::optional<model>& ref) {
//    vec v(10, 10, 10);
//    conf_size s = m.get_size();
//    energy_cal energy(&m, &p, &prec_forcefield,&ig, v);
//    change g(s);
//    std::vector<fl> obj(2, 0.0);
//    output_type_MOEAD tmp(s, obj);
//
//
//    for(int i = 0; i < popsize; i++){
//        tmp.c.randomize(corner1, corner2, generator);
//
//        for (int j = 0; j < 3; j++) {
//            if (j == 1) {
//                tmp.rotor_angle[j] = random_fl(0, pi, generator);
//            }
//            else {
//                tmp.rotor_angle[j] = random_fl(0, 2 * pi, generator);
//            }
//        }
//
//        //?????????????????????????????????????????????????????????????????tmp??orientation
//        std::vector<double> quater;
//        quater.resize(4);
//        threeangles_to_quaternion(tmp.rotor_angle, quater);
//        fl cos_theta = cos(quater[0]);
//        fl sin_theta = sin(quater[0]);
//        qt q(cos_theta, quater[1] * sin_theta, quater[2] * sin_theta, quater[3] * sin_theta);
//        quaternion_normalize(q);//?????
//        tmp.c.ligands[0].rigid.orientation = q;
//        tmp.e = 0;
//        tmp.rm = 0.0;
//
//        fl intra;
//        for(int k = 0; k<tmp.objectives.size();k++){
//            if(k==0){
//                intra = energy.energy_cal_intra(tmp.c,g);
////                tmp.objectives[k] = energy.energy_cal_force(tmp.c);
////                intra = energy.AD4_energy_cal_intra(tmp.c);
//                tmp.objectives[k] = intra;
////                tmp.objectives[k] = energy.energy_cal_forcefield(tmp.c);
//            }else if(k==1){
////                tmp.objectives[k] = energy.energy_cal_hydrophobic(tmp.c);
//                tmp.objectives[k] = energy.energy_cal_inter(sf, tmp.c, intra,g);
////                tmp.objectives[k] = energy.AD4_energy_cal_inter(sf, tmp.c, intra);
//            }else if(k==2){
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
//                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
//            }
//        }
//
//        pop[i] = tmp;
//    }
//}
//
//double distVector(double * vector1, double * vector2, int dim) {
//    //int dim = vector1.size();
//    double sum = 0;
//    for (int n = 0; n < dim; n++) {
//        sum += (vector1[n] - vector2[n]) * (vector1[n] - vector2[n]);
//    }
//    return sqrt(sum);
//} // distVector
//
//void minFastSort(double * x, int * idx, int n, int m) {
//
//    for (int i = 0; i < m; i++) {
//        for (int j = i + 1; j < n; j++) {
//            if (x[i] > x[j]) {
//                double temp = x[i];
//                x[i] = x[j];
//                x[j] = temp;
//                int id = idx[i];
//                idx[i] = idx[j];
//                idx[j] = id;
//            } // if
//        }
//    } // for
//
//} // minFastSort
//
//void randomPermutation(int * perm, int size,rng& generator) {
//    int * index = new int[size];
//    bool * flag = new bool[size];
//
//    for (int n = 0; n < size; n++) {
//        index[n] = n;
//        flag[n] = true;
//    }
//
//    int num = 0;
//    while (num < size) {
//        int start = random_int(0, size - 1, generator);
//        //int start = int(size*nd_uni(&rnd_uni_init));
//        while (true) {
//            if (flag[start]) {
//                perm[num] = index[start];
//                flag[start] = false;
//                num++;
//                break;
//            }
//            if (start == (size - 1)) {
//                start = 0;
//            } else {
//                start++;
//            }
//        }
//    } // while
//
//    delete[] index;
//    delete[] flag;
//
//}
//void MOEAD_aux::initNeighborhood() {
//    double * x = new double[popsize];
//    int * idx = new int[popsize];
//
//    for (int i = 0; i < popsize; i++) {
//        // calculate the distances based on weight vectors
//        for (int j = 0; j < popsize; j++) {
//            x[j] = distVector(lamda[i], lamda[j],
//                              pop[0].objectives.size());
//            //x[j] = dist_vector(population[i].namda,population[j].namda);
//            idx[j] = j;
//            // cout << "x[" << j << "]: " << x[j] << ". idx[" << j << "]: " <<
//            //    idx[j] << endl ;
//        } // for
//
//        // find 'niche' nearest neighboring subproblems
//        minFastSort(x, idx, popsize, T_size);
//        //minfastsort(x,idx,population.size(),niche);
//
//        neighborhood[i] = new int[T_size];
//        for (int k = 0; k < T_size; k++) {
//            neighborhood[i][k] = idx[k];
//            //cout << "neg[ << i << "," << k << "]: " << neighborhood_[i][k] << endl;
//        }
//    } // for
//
//    delete[] x;
//    delete[] idx;
//
//}
//
//void MOEAD_aux::updateReference(const output_type_MOEAD& ind) {
//    for (int n = 0; n < pop[0].objectives.size(); n++) {
//        if (ind.objectives[n] < z[n]) {
//            z[n] = ind.objectives[n];
//        }
//    }
//}
//void MOEAD_aux::initIdealPoint() {
//    for (int i = 0; i < pop[0].objectives.size(); i++) {
//        z[i] = 1.0e+30;
//    } // for
//
//    for (int i = 0; i < popsize; i++) {
//        updateReference(pop[i]);
//    } // for
//}
//
//void MOEAD_aux::matingSelection(std::vector<int> &list, int cid, int size, int type,rng& generator) {
//// list : the set of the indexes of selected mating parents
//    // cid  : the id of current subproblem
//    // size : the number of selected mating parents
//    // type : 1 - neighborhood; otherwise - whole population
//    int ss;
//    int r;
//    int p;
//
//    //ss = neighborhood_[cid].length;
//    ss = T_size;
//    while (list.size() < size) {
//        if (type == 1) {
//            r = random_int(0, ss - 1, generator);
//            p = neighborhood[cid][r];
//            //p = population[cid].table[r];
//        } else {
//            p = random_int(0, popsize - 1, generator);
//        }
//        bool flag = true;
//        for (int i = 0; i < list.size(); i++) {
//            if (list[i] == p) // p is in the list
//            {
//                flag = false;
//                break;
//            }
//        }
//
//        //if (flag) list.push_back(p);
//        if (flag) {
//            list.push_back(p);
//        }
//    }
//}
//
//fl MOEAD_aux::fitnessFunction(output_type_MOEAD individual, fl *lamda_k) {
//    fl fitness;
//    fitness = 0.0;
//
//        fl maxFun = -1.0e+30;
//
//        for (int n = 0; n < pop[0].objectives.size(); n++) {
//            fl diff = fabs(individual.objectives[n] - z[n]);
//
//            fl feval;
//            if (lamda_k[n] == 0) {
//                feval = 0.0001 * diff;
//            } else {
//                feval = diff * lamda_k[n];
//            }
//            if (feval > maxFun) {
//                maxFun = feval;
//            }
//        } // for
//
//        fitness = maxFun;
//    return fitness;
//}
//void MOEAD_aux::updateProblem(output_type_MOEAD indiv, int id, int type,rng& generator) {
//// indiv: child solution
//    // id: the id of current subproblem
//    // type: update solutions in neighborhood (1) or whole population (otherwise)
//    int size;
//    int time;
//
//    time = 0;
//
//    if (type == 1) {
//        //size = neighborhood_[id].length;
//        size = T_size;
//    } else {
//        //size = population_.size();
//        size = popsize;
//    }
//    int * perm = new int[size];
//
//    randomPermutation(perm, size,generator);
//
//    for (int i = 0; i < size; i++) {
//
//        int k;
//
//        if (type == 1) {
//            k = neighborhood[id][perm[i]];
//        } else {
//            k = perm[i];      // calculate the values of objective function regarding
//            // the current subproblem
//        }
//
//        double f1, f2;
//
//        f1 = fitnessFunction(pop[k], lamda[k]);
//        f2 = fitnessFunction(indiv, lamda[k]);
//
//        if (f2 < f1) {
//            pop.erase(pop.begin()+k);
//            pop.push_back(indiv);
//
//            time++;
//        }
//        // the maximal number of solutions updated is not allowed to exceed 'limit'
//        if (time >= nr) {
//            delete[] perm;
//            return;
//        }
//    }
//
//    delete[] perm;
//}
//
//void init_individual(model& m, energy_cal& energy, rng& generator, const vec& corner1, const vec& corner2, const scoring_function& sf, output_type_MOEAD& tmp,const boost::optional<model>& ref) {
//    conf_size s = m.get_size();
//    change g(s);
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
//    //?????????????????????????????????????????????????????????????????tmp??orientation
//    std::vector<double> quater;
//    quater.resize(4);
//    threeangles_to_quaternion(tmp.rotor_angle, quater);
//    fl cos_theta = cos(quater[0]);
//    fl sin_theta = sin(quater[0]);
//    qt q(cos_theta, quater[1] * sin_theta, quater[2] * sin_theta, quater[3] * sin_theta);
//    quaternion_normalize(q);//?????
//    tmp.c.ligands[0].rigid.orientation = q;
//    //    tmp.e = energy(tmp.c, g);
////    tmp.rm = 0.0;
//    fl intra;
//    for(int k = 0; k<tmp.objectives.size();k++){
//        if(k==0){
////                tmp.objectives[k] = energy.energy_cal_force(tmp.c);
//            intra = energy.energy_cal_intra(tmp.c,g);
////            intra = energy.AD4_energy_cal_intra(tmp.c);
//            tmp.objectives[k] = intra;
////            tmp.objectives[k] = energy.energy_cal_forcefield(tmp.c);
//        }else if(k==1){
////                tmp.objectives[k] = energy.energy_cal_hydrophobic(tmp.c);
//            tmp.objectives[k] = energy.energy_cal_inter(sf, tmp.c, intra,g);
////            tmp.objectives[k] = energy.AD4_energy_cal_inter(sf, tmp.c, intra);
//        }else if(k==2){
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
//            tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
//        }
//    }
//}
//output_type_MOEAD MOEAD_aux::DECrossover(model &m, energy_cal &energy, rng &generator, const vec &corner1, const vec &corner2,std::vector<output_type_MOEAD> parents,const scoring_function& sf,const boost::optional<model>& ref) {
//    int jrand;
//    conf_size s = m.get_size();
//    std::vector<fl> obj(2, 0.0);
//    output_type_MOEAD child(s, obj);
//    init_individual(m, energy, generator, corner1, corner2, sf, child,ref);
//    fl cr = 0.8;
//    fl F = 0.5;
//    jrand = random_int(0, dim - 1,generator);
//
//
//    for(int j = 0;j < dim; j++){
//        if(random_fl(0.0, 1.0, generator) < cr || j == jrand){
//            fl value;
//
//            if(j < 3){
//                value = parents[2](j) + F * (parents[0](j) - parents[1](j));
//                if (value < corner1[j]) {
//                    value = corner1[j];
//                }
//                else if (value > corner2[j]) {
//                    value = corner2[j];
//                }
//                child.c.ligands[0].rigid.position[j] = value;
//            }
//            else if(j >= 3 && j < 6){
//                value = parents[2].rotor_angle[j - 3] + F * (parents[0].rotor_angle[j - 3] - parents[1].rotor_angle[j - 3]);
//                child.rotor_angle[j - 3] = value;
//            }
//            else{
//                value = parents[2](j) + F * (parents[0](j) - parents[1](j));
//                child.c.ligands[0].torsions[j - 6] = value;
//                normalize_angle(value);
//            }
//
//        }
//        else{
//            child = parents[2];
//        }
//    }
//    std::vector<double> quater;
//    quater.resize(4);
//    MyNormalize(child.rotor_angle);
//    threeangles_to_quaternion(child.rotor_angle, quater);
//
//    fl cos_theta = cos(quater[0]);
//    fl sin_theta = sin(quater[0]);
//    qt q(cos_theta, quater[1] * sin_theta, quater[2] * sin_theta, quater[3] * sin_theta);
//
//    child.c.ligands[0].rigid.orientation = q;
//    return child;
//}
//
//void MOEAD_aux::mutation(const vec &corner1, const vec &corner2, rng &generator,output_type_MOEAD& child) {
//    double rnd, delta1, delta2, mut_pow, deltaq;
//    double y, yl, yu, val, xy, initalval,eta = 20.0;
//    fl p_mut = 1.0/(dim);
//
//        for (int j = 0; j < dim; j++) {
//            if (random_fl(0, 1.0, generator) <= p_mut ) {
//                if (j < 3) {
//                    y = child.c.ligands[0].rigid.position[j];
//                    initalval = y;
//                    yl = corner1[j];
//                    yu = corner2[j];
//                    delta1 = (y-yl)/(yu-yl);
//                    delta2 = (yu-y)/(yu-yl);
//                    mut_pow = 1.0/(eta+1.0);
//                    rnd = random_fl(0, 1.0, generator);
//                    if(rnd <= 0.5){
//                        xy     = 1.0-delta1;
//                        val    = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(eta+1.0)));
//                        deltaq = pow(val,mut_pow) - 1.0;
//                    }else {
//                        xy     = 1.0-delta2;
//                        val    = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(eta+1.0)));
//                        deltaq = 1.0 - (pow(val,mut_pow));
//                    }
//                    y = y + deltaq*(yu-yl);
//                    if (y<yl)
//                        y = yl;
//                    if (y>yu)
//                        y = yu;
//                    if(std::isnan(y)) // y can be nan result from the pow
//                        child.c.ligands[0].rigid.position[j] = initalval;
//                    else
//                        child.c.ligands[0].rigid.position[j] = y;
//                }
//                else if (j >= 3 && j < 6) {
//                    if(j == 3){
//                        yl = 0;
//                        yu = 2 * pi;
//                    }else if(j == 4){
//                        yl = 0;
//                        yu = pi;
//                    }else{
//                        yl = 0;
//                        yu = 2 * pi;
//                    }
//                    y = child.rotor_angle[j - 3];
//                    initalval = y;
//                    delta1 = (y-yl)/(yu-yl);
//                    delta2 = (yu-y)/(yu-yl);
//                    mut_pow = 1.0/(eta+1.0);
//                    rnd = random_fl(0, 1.0, generator);
//                    if(rnd <= 0.5){
//                        xy     = 1.0-delta1;
//                        val    = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(eta+1.0)));
//                        deltaq = pow(val,mut_pow) - 1.0;
//                    }else {
//                        xy     = 1.0-delta2;
//                        val    = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(eta+1.0)));
//                        deltaq = 1.0 - (pow(val,mut_pow));
//                    }
//                    y = y + deltaq*(yu-yl);
//                    if (y<yl)
//                        y = yl;
//                    if (y>yu)
//                        y = yu;
//                    if(std::isnan(y)) // y can be nan result from the pow
//                        child.rotor_angle[j - 3] = initalval;
//                    else
//                        child.rotor_angle[j - 3] = y;
//                }
//                else {
//                    y = child.c.ligands[0].torsions[j - 6];
//                    initalval = y;
//                    yl = -pi;
//                    yu = pi;
//                    delta1 = (y-yl)/(yu-yl);
//                    delta2 = (yu-y)/(yu-yl);
//                    mut_pow = 1.0/(eta+1.0);
//                    rnd = random_fl(0, 1.0, generator);
//                    if(rnd <= 0.5){
//                        xy     = 1.0-delta1;
//                        val    = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(eta+1.0)));
//                        deltaq = pow(val,mut_pow) - 1.0;
//                    }else {
//                        xy     = 1.0-delta2;
//                        val    = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(eta+1.0)));
//                        deltaq = 1.0 - (pow(val,mut_pow));
//                    }
//                    y = y + deltaq*(yu-yl);
//                    if (y<yl)
//                        y = yl;
//                    if (y>yu)
//                        y = yu;
//                    if(std::isnan(y)) // y can be nan result from the pow
//                        child.c.ligands[0].torsions[j - 6] = initalval;
//                    else
//                        child.c.ligands[0].torsions[j - 6] = y;
//                    normalize_angle(child.c.ligands[0].torsions[j - 6]);
//
//
//
//
//                }
//
//
//
//            }
//        } // for
//        std::vector<double> quater; // ??????????????
//        quater.resize(4);
//        MyNormalize(child.rotor_angle);
//        threeangles_to_quaternion(child.rotor_angle, quater);
//
//        fl cos_theta1 = cos(quater[0]);
//        fl sin_theta1 = sin(quater[0]);
//        qt q(cos_theta1, quater[1] * sin_theta1, quater[2] * sin_theta1, quater[3] * sin_theta1);
//
//        child.c.ligands[0].rigid.orientation = q;
//}
//
//void
//MOEAD_aux::calculate_objectives(model &m, const precalculate &p, const igrid &ig, vec v, output_type_MOEAD &indivisual,
//                                const scoring_function &sf, energy_cal &energy, const boost::optional<model> &ref,
//                                change &g) {
//    conf_size s = m.get_size();
//    fl intramolecular_energy;
//    for(int k = 0; k<indivisual.objectives.size();k++){
//        if(k==0){
////                tmp.objectives[k] = energy.energy_cal_force(tmp.c);
//            intramolecular_energy = energy.energy_cal_intra(indivisual.c,g);
////            intramolecular_energy = energy.AD4_energy_cal_intra(indivisual.c);
////            indivisual.objectives[k] = energy.energy_cal_forcefield(indivisual.c);
//            indivisual.objectives[k] = intramolecular_energy;
//        }else if(k==1){
////                tmp.objectives[k] = energy.energy_cal_hydrophobic(tmp.c);
//            indivisual.objectives[k] = energy.energy_cal_inter(sf, indivisual.c, intramolecular_energy,g);
////            indivisual.objectives[k] = energy.AD4_energy_cal_inter(sf, indivisual.c, intramolecular_energy);
//        }else if(k==2){
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
//            indivisual.objectives[k] = energy.energy_cal_rmsd(ref, indivisual.c);
//        }
//    }
//}
//void MOEAD::operator()(model& m, output_container_MOEAD& out, const precalculate& p, const precalculate& prec_forcefield, const igrid& ig,
//                       const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2,
//                       incrementable* increment_me, rng& generator, const scoring_function& sf,const boost::optional<model>& ref) const {
//
//    MOEAD_aux aux(m,corner1,corner2);
//    vec v(10, 10, 10);
//    conf_size s = m.get_size();
//    energy_cal energy(&m, &p, &prec_forcefield,&ig, v);
//    change g(s);
//
//    aux.initPopulation(m,p,prec_forcefield,ig,p_widened,corner1,corner2,generator,sf,ref);
//
//    aux.initUniformWeight();
//
//    aux.initNeighborhood();
//
//    aux.initIdealPoint();
//
//    int g1 = 0;
//    while(g1 < num_steps){
//
//        int * permutation = new int[aux.popsize];
//        randomPermutation(permutation, aux.popsize, generator);
//        for (int i = 0; i < aux.popsize; i++) {
//            if (increment_me)
//                ++(*increment_me);
//            int n = permutation[i]; // or int n = i;
//            //int n = i ; // or int n = i;
//            int type;
//            fl rnd = random_fl(0.0,1.0,generator);
//
//            // STEP 2.1. Mating selection based on probability
//            if (rnd < aux.delta) // if (rnd < realb)
//            {
//                type = 1;   // neighborhood
//            } else {
//                type = 2;   // whole population
//            }
//            std::vector<int> p1;
//            aux.matingSelection(p1, n, 2, type,generator);
//
//            // STEP 2.2. Reproduction
//            output_type_MOEAD child;
//            std::vector<output_type_MOEAD> parents;
//            parents.push_back(aux.pop[p1[0]]);
//            parents.push_back(aux.pop[p1[1]]);
//            parents.push_back(aux.pop[n]);
//
//
//
//
//            // Apply DE crossover
//
//
//            child = aux.DECrossover(m,energy,generator,corner1,corner2,parents,sf,ref);
//
//            // Apply mutation
//            aux.mutation(corner1,corner2,generator,child);
//            // Evaluation
//            aux.calculate_objectives(m, p, ig, v, child, sf,energy,ref,g);
//
//            g1++;
//
//            // STEP 2.3. Repair. Not necessary
//
//            // STEP 2.4. Update z_
//            aux.updateReference(child);
//
//            // STEP 2.5. Update of solutions
//            aux.updateProblem(child, n, type,generator);
//        } // for
//
//        delete[] permutation;
//    }
//
//    for(auto& s : aux.pop){
//        s.e = s.objectives[0]+s.objectives[1];
//    }
//
//    for (int i = 0; i < aux.pop.size()-1; i++) {
//        for (int j = 0; j < aux.pop.size() - 1 - i; j++) {
//            if (aux.pop[j+1].e < aux.pop[j].e) {
//                output_type_MOEAD tmp = aux.pop[j];
//                aux.pop[j] = aux.pop[j + 1];
//                aux.pop[j + 1] = tmp;
//            }
//        }
//    }
//
//    for(auto& p : aux.pop){
//        m.set(p.c);
//        p.coords = m.get_heavy_atom_movable_coords();
//        add_to_output_container(out,p,min_rmsd,num_saved_mins);
//
//    }
//}