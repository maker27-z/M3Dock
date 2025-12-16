////
//// Created by 91686 on 2023/7/10.
////
//
//#include "MPSOD.h"
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
//void MPSOD_aux::initUniformWeight() {
//    if((2 == 2) && (popsize <= 300)){
//        for(int n = 0; n < popsize; n++){
//            lamda[n] = new fl[2];
//            fl a = 1.0 * n / (popsize - 1);
//            lamda[n][0] = a;
//            lamda[n][1] = 1 - a;
//        }
//    }
//    else{
//        std::ostringstream os;
//        os << "../" << "W" << 2<< "D_"
//           << popsize << ".dat";
//        std::string dataFileName;
//        dataFileName = os.str();
//
//        // Open the file
//        std::ifstream in(dataFileName.c_str());
//        if( !in ) {
//            std::cout << "initUniformWeight: failed when reading from file: : " <<
//                      dataFileName << std::endl;
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
//            lamda[i] = new double[2];
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
//void MPSOD_aux::initPopulation(model& m,  const precalculate& p, const precalculate& prec_forcefield, const igrid& ig,
//                               const precalculate& p_widened,  const vec& corner1, const vec& corner2,
//                               rng& generator,const scoring_function& sf,const boost::optional<model>& ref) {
//    vec v(10, 10, 10);
//    conf_size s = m.get_size();
//    energy_cal energy(&m, &p, &prec_forcefield,&ig, v);
//    change g(s);
//    std::vector<fl> obj(2, 0.0);
//    output_type_MPSOD tmp(s, obj);
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
//
//void MPSOD_aux::initVelocity(model& m,  const precalculate& p, const precalculate& prec_forcefield, const igrid& ig,
//                             const precalculate& p_widened,  const vec& corner1, const vec& corner2,
//                             rng& generator,const scoring_function& sf,const boost::optional<model>& ref) {
//    vec v(10, 10, 10);
//    conf_size s = m.get_size();
//    energy_cal energy(&m, &p, &prec_forcefield,&ig, v);
////    MPSOD_aux aux(m);
//    change g(s);
//    std::vector<fl> obj(2, 0.0);
//    output_type_MPSOD tmp(s, obj);
//
//    for(int i = 0;i < popsize;i++){
//        tmp.c.randomize(corner1, corner2, generator);
//        for(int j = 0;j < dim;j++){
//            if (j < 3) {
//                tmp.c.ligands[0].rigid.position[j] = 0.0;
//            }
//            else if (j >= 3 && j < 6) {
//                tmp.rotor_angle[j - 3] = 0.0;
//            }
//            else{
//                tmp.c.ligands[0].torsions[j - 6] = 0.0;
//            }
//        }
//        //这里要将初始化好的三个角度转换为四元数，然后再将这个转换好的随机初始化四元数赋给tmp的orientation
//        std::vector<double> quater;
//        quater.resize(4);
//        threeangles_to_quaternion(tmp.rotor_angle, quater);
//        fl cos_theta = cos(quater[0]);
//        fl sin_theta = sin(quater[0]);
//        qt q(cos_theta, quater[1] * sin_theta, quater[2] * sin_theta, quater[3] * sin_theta);
//        quaternion_normalize(q);//归一化
//        tmp.c.ligands[0].rigid.orientation = q;
//        tmp.e = 0;
//        tmp.rm = 0.0;
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
//        velocity[i] = tmp;
//    }
//
//}
//
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
//void MPSOD_aux::initNeighborhood() {
//    double * x = new double[popsize];
//    int * idx = new int[popsize];
//
//    for (int i = 0; i < popsize; i++) {
//        // calculate the distances based on weight vectors
//        for (int j = 0; j < popsize; j++) {
//            x[j] = distVector(lamda[i], lamda[j],
//                              2);
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
//void MPSOD_aux::updateReference(const std::vector<output_type_MPSOD>& ind) {
//    for(int i = 0; i < ind.size(); i++){
//        for (int n = 0; n < ind[0].objectives.size(); n++) {
//            if (ind[i].objectives[n] < z[n]) {
//                z[n] = ind[i].objectives[n];
//            }
//        }
//    }
//
//}
//void MPSOD_aux::initIdealPoint() {
//    for (int i = 0; i < 2; i++) {
//        z[i] = 1.0e+30;
//    } // for
//
//    for (int i = 0; i < popsize; i++) {
//        updateReference(pop);
//    } // for
//}
//
//std::vector<output_type_MPSOD> MPSOD_aux::TournamentSelection(int K, int N,rng& generator) {
//    std::vector<output_type_MPSOD> parents;
//    for(int i = 0; i< N; i++){
//        fl r1 = random_int(0, N-1, generator);
//        fl r2 = random_int(0, N-1, generator);
//        if (pop[r1].rank != pop[r2].rank) {
//            if (pop[r1].rank < pop[r2].rank) {
//                parents.push_back(pop[r1]);
//            }
//            else {
//                parents.push_back(pop[r2]);
//            }
//        }
//        else if (pop[r1].crowding_distance != pop[r2].crowding_distance) {
//            if (pop[r1].crowding_distance > pop[r2].crowding_distance) {
//                parents.push_back(pop[r1]);
//            }
//            else {
//                parents.push_back(pop[r2]);
//            }
//        }
//        else {
//            parents.push_back(pop[r1]);
//        }
//    }
//    return parents;
//}
//std::vector<output_type_MPSOD> MPSOD_aux::matingSelection(model& m,  const precalculate& p, const precalculate& prec_forcefield, const igrid& ig,
//                                                          const precalculate& p_widened,  const vec& corner1, const vec& corner2,
//                                                          rng& generator,const scoring_function& sf,const boost::optional<model>& ref,energy_cal& energy) {
//// list : the set of the indexes of selected mating parents
//    // cid  : the id of current subproblem
//    // size : the number of selected mating parents
//    // type : 1 - neighborhood; otherwise - whole population
////    int ss;
////    int r;
////    int p;
////
////    //ss = neighborhood_[cid].length;
////    ss = T_size;
////    while (list.size() < size) {
////        if (type == 1) {
////            r = random_int(0, ss - 1, generator);
////            p = neighborhood[cid][r];
////            //p = population[cid].table[r];
////        } else {
////            p = random_int(0, popsize - 1, generator);
////        }
////        bool flag = true;
////        for (int i = 0; i < list.size(); i++) {
////            if (list[i] == p) // p is in the list
////            {
////                flag = false;
////                break;
////            }
////        }
////
////        //if (flag) list.push_back(p);
////        if (flag) {
////            list.push_back(p);
////        }
////    }
//
//        // PopObj = PopObj - repmat(Z,size(PopObj,1),1);
//    std::vector<output_type_MPSOD> popObj;
//    popObj.resize(pop.size());
//    initpopVector(m,p,prec_forcefield,ig,p_widened,corner1,corner2,generator,sf,ref,popObj);
//
//        for(int i = 0; i < popsize; i++){
//            for(int j = 0; j < 2; j++){
//                popObj[i].objectives[j] = pop[i].objectives[j] - z[j];
//            }
//        }
//
//        //Parent = TournamentSelection(2,N,-CrowdingDistance(PopObj));
//    std::vector<output_type_MPSOD> parents = TournamentSelection(2,popsize,generator);
//
//    std::vector<output_type_MPSOD> P;
////    initpopVector(m,p,prec_forcefield,ig,p_widened,corner1,corner2,generator,sf,ref,P);
//        for(int m = 0; m < popsize; m++){
//            if(random_fl(0.0, 1.0, generator) < 0.9){
//                P.resize(T_size);
//                for(int n = 0;n < T_size; n++){
//                    P[n] = pop[neighborhood[m][n]];
//                }
//            }else{
//                P = pop;
//            }
//            pbest[m] = P[random_int(0,P.size() - 1,generator)];
//
//            std::vector<fl> W;
//            for(int a = 0; a < popObj[0].objectives.size(); a++){
//                fl wa = 0.0;
//                for(int b = 0; b < P.size(); b++){
//                    wa += lamda[b][a];
//                }
//                W.push_back(1.0*wa/P.size());
//            }
//
//            std::vector<fl> sum_normal;
//            for(int c = 0; c < P.size(); c++){
//                fl normal = 0.0;
//                for(int d = 0; d < popObj[0].objectives.size(); d++){
//                    normal += P[c].objectives[d] * P[c].objectives[d];
//                }
//                sum_normal.push_back(pow(normal,0.6));
//            }
//
//            std::vector<fl> best_fl;
//            for(int e = 0; e < P.size(); e++){
//                fl temp = 0.0;
//                for(int f = 0; f < popObj[0].objectives.size(); f++){
//                    temp += P[e].objectives[f] * W[f];
//                }
//                best_fl.push_back(1.0 * temp /sum_normal[e]);
//
//            }
//            auto maxIterator = std::max_element(best_fl.begin(), best_fl.end());
//            int maxIndex = std::distance(best_fl.begin(), maxIterator);
//            gbest[m] = P[maxIndex];
//        }
//    return parents;
//}
//
//fl MPSOD_aux::fitnessFunction(output_type_MPSOD individual, fl *lamda_k) {
//    fl fitness;
//    fitness = 0.0;
//
//    fl maxFun = -1.0e+30;
//
//    for (int n = 0; n < 2; n++) {
//        fl diff = fabs(individual.objectives[n] - z[n]);
//
//        fl feval;
//        if (lamda_k[n] == 0) {
//            feval = 0.0001 * diff;
//        } else {
//            feval = diff * lamda_k[n];
//        }
//        if (feval > maxFun) {
//            maxFun = feval;
//        }
//    } // for
//
//    fitness = maxFun;
//    return fitness;
//}
//void MPSOD_aux::updateProblem(output_type_MPSOD indiv, int id, int type,rng& generator) {
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
//void init_individual(model& m, energy_cal& energy, rng& generator, const vec& corner1, const vec& corner2, const scoring_function& sf, output_type_MPSOD& tmp,const boost::optional<model>& ref) {
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
//output_type_MPSOD MPSOD_aux::DECrossover(model &m, energy_cal &energy, rng &generator, const vec &corner1, const vec &corner2,std::vector<output_type_MPSOD> parents,const scoring_function& sf,const boost::optional<model>& ref) {
//    int jrand;
//    conf_size s = m.get_size();
//    std::vector<fl> obj(2, 0.0);
//    output_type_MPSOD child(s, obj);
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
//void MPSOD_aux::mutation(const vec& corner1, const vec& corner2,rng& generator,std::vector<output_type_MPSOD>& parents) {
//    double rnd, delta1, delta2, mut_pow, deltaq;
//    double y, yl, yu, val, xy, initalval,eta = 20.0;
//    fl p_mut = 1.0/(dim);
//    for(int i = 0;i < popsize;i++){
////        if((i%6) == 0){
//        for (int j = 0; j < dim; j++) {
//            if (random_fl(0, 1.0, generator) <= p_mut ) {
//                if (j < 3) {
//                    y = parents[i].c.ligands[0].rigid.position[j];
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
//                        parents[i].c.ligands[0].rigid.position[j] = initalval;
//                    else
//                        parents[i].c.ligands[0].rigid.position[j] = y;
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
//                    y = parents[i].rotor_angle[j - 3];
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
//                        parents[i].rotor_angle[j - 3] = initalval;
//                    else
//                        parents[i].rotor_angle[j - 3] = y;
//                }
//                else {
//                    y = parents[i].c.ligands[0].torsions[j - 6];
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
//                        parents[i].c.ligands[0].torsions[j - 6] = initalval;
//                    else
//                        parents[i].c.ligands[0].torsions[j - 6] = y;
//                    normalize_angle(parents[i].c.ligands[0].torsions[j - 6]);
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
//        std::vector<double> quater; // 初始化四元数变量
//        quater.resize(4);
//        MyNormalize(parents[i].rotor_angle);
//        threeangles_to_quaternion(parents[i].rotor_angle, quater);
//
//        fl cos_theta1 = cos(quater[0]);
//        fl sin_theta1 = sin(quater[0]);
//        qt q(cos_theta1, quater[1] * sin_theta1, quater[2] * sin_theta1, quater[3] * sin_theta1);
//
//        parents[i].c.ligands[0].rigid.orientation = q;
////        }
//    }
//
//
//
//}
//
//void MPSOD_aux::calculate_objectives(model &m, const precalculate &p, const igrid &ig, vec v, output_type_MPSOD &indivisual,
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
//
//void MPSOD_aux::calculate_crowding_distance(std::vector<output_type_MPSOD>& pop){
//    int length = pop.size();
//    for (int i = 0; i < length; i++) {
//        pop[i].crowding_distance = 0.0;
//    }
//
//    for (int i = 0; i < 2; i++) {
//        std::sort(pop.begin(), pop.end(), [i](const output_type_MPSOD& a, const output_type_MPSOD& b) {
//            return a.objectives[i] < b.objectives[i];
//        });
//        pop[0].crowding_distance = pop[length - 1].crowding_distance = std::numeric_limits<fl>::infinity();
//        fl fmax = pop[length - 1].objectives[i];
//        fl fmin = pop[0].objectives[i];
//        if (fmax != fmin) {
//            for (int j = 1; j < length - 1; j++) {
//                pop[j].crowding_distance += (pop[j + 1].objectives[i] - pop[j - 1].objectives[i]) / (fmax - fmin);
//            }
//        }
//
//    }
//}
//
//void MPSOD_aux::initpbest(model& m,  const precalculate& p, const precalculate& prec_forcefield, const igrid& ig,
//                          const precalculate& p_widened,  const vec& corner1, const vec& corner2,
//                          rng& generator,const scoring_function& sf,const boost::optional<model>& ref){
//    vec v(10, 10, 10);
//    conf_size s = m.get_size();
//    energy_cal energy(&m, &p, &prec_forcefield,&ig, v);
////    MPSOD_aux aux(m);
//    change g(s);
//    std::vector<fl> obj(2, 0.0);
//    output_type_MPSOD tmp(s, obj);
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
//        //这里要将初始化好的三个角度转换为四元数，然后再将这个转换好的随机初始化四元数赋给tmp的orientation
//        std::vector<double> quater;
//        quater.resize(4);
//        threeangles_to_quaternion(tmp.rotor_angle, quater);
//        fl cos_theta = cos(quater[0]);
//        fl sin_theta = sin(quater[0]);
//        qt q(cos_theta, quater[1] * sin_theta, quater[2] * sin_theta, quater[3] * sin_theta);
//        quaternion_normalize(q);//归一化
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
//        pbest[i] = tmp;
//    }
//}
//
//std::vector<std::vector<fl>> MPSOD_aux::divide() {
//    std::vector<std::vector<fl>> value;
//    fl diff = 0.0;
//    for(int i = 0;i < pop.size(); i++){
//        std::vector<fl> v;
//        for (int n = 0; n < popsize; n++) {
//            fl sum = 0.0;
//            for (int m = 0; m < 2; m++) {
//                diff = fabs(pop[i].objectives[m] - z[m]);
//                if (lamda[n][m] == 0) {
//                    sum += 0.0001 * diff;
//                } else {
//                    sum += diff * lamda[n][m];
//                }
//            }
//            v.push_back(sum);
//        } // for
//        value.push_back(v);
//    }
//    return value;
//}
//
//void MPSOD_aux::initNext(model &m, const precalculate &p, const precalculate &prec_forcefield, const igrid &ig,
//                         const precalculate &p_widened, const vec &corner1, const vec &corner2, rng &generator,
//                         const scoring_function &sf, const boost::optional<model> &ref) {
//    vec v(10, 10, 10);
//    conf_size s = m.get_size();
//    energy_cal energy(&m, &p, &prec_forcefield,&ig, v);
//    change g(s);
//    std::vector<fl> obj(2, 0.0);
//    output_type_MPSOD tmp(s, obj);
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
//        next[i] = tmp;
//    }
//}
//
//void MPSOD_aux::computeSpeed(model& m,  const precalculate& p, const precalculate& prec_forcefield, const igrid& ig,
//                             const precalculate& p_widened,  const vec& corner1, const vec& corner2,
//                             rng& generator,const scoring_function& sf,const boost::optional<model>& ref){
//    for(int i = 0;i < popsize;i++){
//
//
////        std::sort(archive.begin(),archive.end(), crowding_distance_comparator);
////        gbest = archive[0];
//        r1 = random_fl(0,1.0,generator);
//        r2 = random_fl(0,1.0,generator);
//        c1 = random_fl(1.5, 2.5, generator);
//        c2 = random_fl(1.5, 2.5, generator);
//
//        for(int j = 0;j < dim;j++){
//            if (j < 3) {
//                velocity[i].c.ligands[0].rigid.position[j] =   w*velocity[i].c.ligands[0].rigid.position[j] + c1 * r1 * (pbest[i].c.ligands[0].rigid.position[j] - pop[i].c.ligands[0].rigid.position[j])
//                                                                                  + c2 * r2 * (gbest[i].c.ligands[0].rigid.position[j] - pop[i].c.ligands[0].rigid.position[j]);
////
//            }
//            else if (j >= 3 && j < 6) {
//                velocity[i].rotor_angle[j - 3] =   w*velocity[i].rotor_angle[j - 3] + c1 * r1 * (pbest[i].rotor_angle[j - 3] - pop[i].rotor_angle[j - 3])
//                                                                      + c2 * r2 * (gbest[i].rotor_angle[j - 3] - pop[i].rotor_angle[j - 3]);
//            }
//            else{
//                velocity[i].c.ligands[0].torsions[j - 6] =   w*velocity[i].c.ligands[0].torsions[j - 6] + c1 * r1 * (pbest[i].c.ligands[0].torsions[j - 6] - pop[i].c.ligands[0].torsions[j - 6])
//                                                                                + c2 * r2 * (gbest[i].c.ligands[0].torsions[j - 6] - pop[i].c.ligands[0].torsions[j - 6]);
////                normalize_angle(velocity[i].c.ligands[0].torsions[j - 6]);
//            }
//        }
//
//        std::vector<double> quater; // 初始化四元数变量
//        quater.resize(4);
////        MyNormalize(velocity[i].rotor_angle);
//        threeangles_to_quaternion(velocity[i].rotor_angle, quater);
//
//        fl cos_theta1 = cos(quater[0]);
//        fl sin_theta1 = sin(quater[0]);
//        qt q(cos_theta1, quater[1] * sin_theta1, quater[2] * sin_theta1, quater[3] * sin_theta1);
//
//        velocity[i].c.ligands[0].rigid.orientation = q;
//
//    }
//}
//
//
//void MPSOD_aux::computeNewPositions(const vec& corner1, const vec& corner2,std::vector<output_type_MPSOD>& parents) {
//    for(int i = 0;i < popsize;i++){
//        for(int j = 0;j < dim;j++){
//            if (j < 3) {
//                parents[i].c.ligands[0].rigid.position[j] = parents[i].c.ligands[0].rigid.position[j] + velocity[i].c.ligands[0].rigid.position[j];
//                if (parents[i].c.ligands[0].rigid.position[j] < corner1[j]) {
//                    parents[i].c.ligands[0].rigid.position[j] = corner1[j];
//                }
//                else if (parents[i].c.ligands[0].rigid.position[j] > corner2[j]) {
//                    parents[i].c.ligands[0].rigid.position[j] = corner2[j];
//                }
//            }
//            else if(j >= 3 && j < 6){
//                    parents[i].rotor_angle[j - 3] = parents[i].rotor_angle[j - 3] + velocity[i].rotor_angle[j - 3];
//            }
//            else{
//
//                parents[i].c.ligands[0].torsions[j - 6] = parents[i].c.ligands[0].torsions[j - 6] + velocity[i].c.ligands[0].torsions[j - 6];
//
//                normalize_angle(parents[i].c.ligands[0].torsions[j - 6]);
//            }
//        }
//        std::vector<double> quater; // 初始化四元数变量
//        quater.resize(4);
//        MyNormalize(parents[i].rotor_angle);
//        threeangles_to_quaternion(parents[i].rotor_angle, quater);
//
//        fl cos_theta1 = cos(quater[0]);
//        fl sin_theta1 = sin(quater[0]);
//        qt q(cos_theta1, quater[1] * sin_theta1, quater[2] * sin_theta1, quater[3] * sin_theta1);
//
//        parents[i].c.ligands[0].rigid.orientation = q;
//    }
//}
//
//void MPSOD_aux::DE_Correction(const vec& corner1, const vec& corner2,std::vector<output_type_MPSOD>& parents) {
//    for (int i = 0; i < popsize; i++) {
//        for (int j = 0; j < dim; j++) {
//            if (j < 3) {
//                parents[i].c.ligands[0].rigid.position[j] = parents[i].c.ligands[0].rigid.position[j] + 0.5 *
//                                                                                                (gbest[i].c.ligands[0].rigid.position[j] -
//                                                                                                 pbest[i].c.ligands[0].rigid.position[j]);
//                if (parents[i].c.ligands[0].rigid.position[j] < corner1[j]) {
//                    parents[i].c.ligands[0].rigid.position[j] = corner1[j];
//                } else if (parents[i].c.ligands[0].rigid.position[j] > corner2[j]) {
//                    parents[i].c.ligands[0].rigid.position[j] = corner2[j];
//                }
//            } else if (j >= 3 && j < 6) {
//                if (j == 3) {
//                    parents[i].rotor_angle[j - 3] = parents[i].rotor_angle[j - 3] +
//                                                0.5 * (gbest[i].rotor_angle[j - 3] - pbest[i].rotor_angle[j - 3]);
//                } else {
//
//                    parents[i].c.ligands[0].torsions[j - 6] = parents[i].c.ligands[0].torsions[j - 6] + 0.5 *
//                                                                                                (gbest[i].c.ligands[0].torsions[
//                                                                                                         j - 6] -
//                                                                                                 pbest[i].c.ligands[0].torsions[
//                                                                                                         j - 6]);
//                    normalize_angle(parents[i].c.ligands[0].torsions[j - 6]);
//                }
//            }
//            std::vector<double> quater; // 初始化四元数变量
//            quater.resize(4);
//            MyNormalize(parents[i].rotor_angle);
//            threeangles_to_quaternion(parents[i].rotor_angle, quater);
//
//            fl cos_theta1 = cos(quater[0]);
//            fl sin_theta1 = sin(quater[0]);
//            qt q(cos_theta1, quater[1] * sin_theta1, quater[2] * sin_theta1, quater[3] * sin_theta1);
//
//            parents[i].c.ligands[0].rigid.orientation = q;
//        }
//    }
//}
//
//std::unordered_map<int, std::vector<output_type_MPSOD>> MPSOD_aux::fast_nondominated_sort(std::vector<output_type_MPSOD>& pop) {
//    std::unordered_map<int, std::vector<output_type_MPSOD>> F;
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
//    std::vector<output_type_MPSOD> Q;
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
//void MPSOD_aux::initpopVector(model &m, const precalculate &p, const precalculate &prec_forcefield, const igrid &ig,
//                              const precalculate &p_widened, const vec &corner1, const vec &corner2, rng &generator,
//                              const scoring_function &sf, const boost::optional<model> &ref,
//                              std::vector<output_type_MPSOD> &popObj) {
//    vec v(10, 10, 10);
//    conf_size s = m.get_size();
//    energy_cal energy(&m, &p, &prec_forcefield,&ig, v);
//    change g(s);
//    std::vector<fl> obj(2, 0.0);
//    output_type_MPSOD tmp(s, obj);
//
//
//    for(int i = 0; i < popObj.size(); i++){
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
//        popObj[i] = tmp;
//    }
//}
//void MPSOD_aux::classification(model& m,  const precalculate& p, const precalculate& prec_forcefield, const igrid& ig,
//                               const precalculate& p_widened,  const vec& corner1, const vec& corner2,
//                               rng& generator,const scoring_function& sf,const boost::optional<model>& ref,energy_cal& energy) {
//// indiv: child solution
//    // id: the id of current subproblem
//    // type: update solutions in neighborhood (1) or whole population (otherwise)
//
//    int size;
//    //size = population_.size();
//    size = popsize;
//    next.resize(popsize);
////    initpopVector(m,p,prec_forcefield,ig,p_widened,corner1,corner2,generator,sf,ref,next);
//    std::vector<std::vector<fl>> value = divide();
//    std::vector<int> P;
//     for (const auto& row : value) {
//        auto maxElement = std::max_element(row.begin(), row.end());
//        P.push_back(std::distance(row.begin(), maxElement));
//    }
//
//    initNext(m,p,prec_forcefield,ig,p_widened,corner1,corner2,generator,sf,ref);
//
//    for(int i = 0; i < popsize; i++){
//        std::vector<output_type_MPSOD> current;
//        for(int j = 0; j < P.size(); j++){
//            if(P[j] == i){
//                current.push_back(pop[j]);
//            }
//        }
//        if(current.empty()){
//            init_individual(m, energy, generator, corner1, corner2, sf, next[i],ref);
//        }else{
//
//            std::unordered_map<int, std::vector<output_type_MPSOD>> F = fast_nondominated_sort(current);
//
////            for(int a = 0; a < F[1].size(); a++){
////                for(int b = 0; b < F[1][0].objectives.size(); b++){
////
////                }
////            }
//
//            std::vector<output_type_MPSOD> popObj;
//            popObj.resize(F[1].size());
//            initpopVector(m,p,prec_forcefield,ig,p_widened,corner1,corner2,generator,sf,ref,popObj);
//            for(int t = 0; t < F[1].size(); t++){
//                for(int v = 0; v < F[1][0].objectives.size(); v++){
//                    popObj[t].objectives[v] = F[1][t].objectives[v] - z[v];
//                }
//            }
//
//            std::vector<fl> sum_normal;
//            for(int c = 0; c < F[1].size(); c++){
//                fl normal = 0.0;
//                for(int d = 0; d < F[1][0].objectives.size(); d++){
//                    normal += F[1][c].objectives[d] * F[1][c].objectives[d];
//                }
//                sum_normal.push_back(pow(normal,0.6));
//            }
//
//
//
//            std::vector<fl> best_fl;
//            for(int e = 0; e < F[1].size(); e++){
//                fl temp = 0.0;
//                for(int f = 0; f < F[1][0].objectives.size(); f++){
//                    temp += popObj[e].objectives[f] * lamda[i][f];
//                }
//                best_fl.push_back(1.0 * temp /sum_normal[e]);
//
//            }
//
//            auto maxIterator = std::max_element(best_fl.begin(), best_fl.end());
//            int maxIndex = std::distance(best_fl.begin(), maxIterator);
//            next[i] = F[1][maxIndex];
//        }
//    }
////    int * perm = new int[size];
////
////    randomPermutation(perm, size,generator);
////
////    for (int i = 0; i < size; i++) {
////
////        int k;
////        k = perm[i];      // calculate the values of objective function regarding
////            // the current subproblem
////
////        double f1, f2;
////
////        f1 = fitnessFunction(pop[k], lamda[k]);
////
////        }
////        // the maximal number of solutions updated is not allowed to exceed 'limit'
////
////    delete[] perm;
//}
//void MPSOD::operator()(model& m, output_container_MPSOD& out, const precalculate& p, const precalculate& prec_forcefield, const igrid& ig,
//                       const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2,
//                       incrementable* increment_me, rng& generator, const scoring_function& sf,const boost::optional<model>& ref) const {
//
//    MPSOD_aux aux(m,corner1,corner2);
//    vec v(10, 10, 10);
//    conf_size s = m.get_size();
//    energy_cal energy(&m, &p, &prec_forcefield,&ig, v);
//    change g(s);
//
//    //MOPSO
//    //初始化种群
//    aux.initPopulation(m,p,prec_forcefield,ig,p_widened,corner1,corner2,generator,sf,ref);
//
//    //初始化每个粒子的速度
//    aux.initVelocity(m,p,prec_forcefield,ig,p_widened,corner1,corner2,generator,sf,ref);
//
//
//    //初始化每个粒子的pbest
//    aux.initpbest(m,p,prec_forcefield,ig,p_widened,corner1,corner2,generator,sf,ref);
//
//    //MOEAD
//    //初始化方向向量lamda
//    aux.initUniformWeight();
//
//    //初始化邻域粒子
//    aux.initNeighborhood();
//
//    //初始化关键点z
//    aux.initIdealPoint();
//
//
//    //为种群每个粒子分配拥挤距离
////    aux.calculate_crowding_distance(aux.pop);
//
//    aux.classification(m,p,prec_forcefield,ig,p_widened,corner1,corner2,generator,sf,ref,energy);
//
//
//    int g1 = 0;
//    while(g1 < num_steps){
//        if (increment_me)
//            ++(*increment_me);
////        int * permutation = new int[aux.popsize];
////        randomPermutation(permutation, aux.popsize, generator);
////        for (int i = 0; i < aux.popsize; i++) {
////
//////            int n = permutation[i]; // or int n = i;
////            //int n = i ; // or int n = i;
////            int n = i;
////            int type;
////            fl rnd = random_fl(0.0,1.0,generator);
////
////            // STEP 2.1. Mating selection based on probability
////            if (rnd < aux.delta) // if (rnd < realb)
////            {
////                type = 1;   // neighborhood
////            } else {
////                type = 2;   // whole population
////            }
////            std::vector<int> p1;
////            std::vector<output_type_MPSOD> paresnts = aux.matingSelection(p1, n, 2, type,generator);
////
////            // STEP 2.2. Reproduction
////            output_type_MPSOD child;
////            std::vector<output_type_MPSOD> parents;
////            parents.push_back(aux.pop[p1[0]]);
////            parents.push_back(aux.pop[p1[1]]);
////            parents.push_back(aux.pop[n]);
////
////
////
////
////            // Apply DE crossover
////
////
////            child = aux.DECrossover(m,energy,generator,corner1,corner2,parents,sf,ref);
////
////            // Apply mutation
////            aux.mutation(corner1,corner2,generator,child);
////            // Evaluation
////            aux.calculate_objectives(m, p, ig, v, child, sf,energy,ref,g);
////
////
////            // STEP 2.3. Repair. Not necessary
////
////            // STEP 2.4. Update z_
////            aux.updateReference(child);
////
////            // STEP 2.5. Update of solutions
////            aux.updateProblem(child, n, type,generator);
////        } // for
//
////        delete[] permutation;
//
//
//        std::vector<output_type_MPSOD> parents =  aux.matingSelection(m,p,prec_forcefield,ig,p_widened,corner1,corner2,generator,sf,ref,energy);
//
//
//        aux.computeSpeed(m,p,prec_forcefield,ig,p_widened,corner1,corner2,generator,sf,ref);
//
//        aux.computeNewPositions(corner1,corner2,parents);
//
//        aux.DE_Correction(corner1, corner2,parents);
//
//        aux.mutation(corner1, corner2, generator,parents);
////        aux.calculate_objectives(m, p, ig, v, parents, sf,energy,ref,g);
//
//        aux.updateReference(parents);
////        aux.pop.resize(aux.pop.size()+parents.size());
//        for(int i = 0; i < parents.size(); i++){
//            aux.pop.push_back(parents[i]);
//        }
//
//        aux.classification(m,p,prec_forcefield,ig,p_widened,corner1,corner2,generator,sf,ref,energy);
//        aux.pop = aux.next;
//        g1++;
//
//    }
//
//    for(auto& s : aux.pop){
//        s.e = s.objectives[0]+s.objectives[1];
//    }
//
//    for (int i = 0; i < aux.pop.size()-1; i++) {
//        for (int j = 0; j < aux.pop.size() - 1 - i; j++) {
//            if (aux.pop[j+1].e < aux.pop[j].e) {
//                output_type_MPSOD tmp = aux.pop[j];
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