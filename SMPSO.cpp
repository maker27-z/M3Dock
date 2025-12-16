////
//// Created by 91686 on 2023/7/3.
////
//
//#include "SMPSO.h"
//#include "coords.h"
//
//void threeangles_to_quaternion(const vec& three_angle, std::vector<double>& quater) {
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
//
//
//void SMPSO_aux::calculate_objectives(model& m, const precalculate& p, const igrid& ig, vec v, output_type_SMPSO& indivisual, const scoring_function& sf, const scoring_function& sf_ad4,energy_cal& energy,const boost::optional<model>& ref,change &g){
//    conf_size s = m.get_size();
//    fl intramolecular_energy;
//
//    for(int k = 0; k<indivisual.objectives.size();k++){
//        if(k==0){
////                tmp.objectives[k] = energy.energy_cal_force(tmp.c);
////            intramolecular_energy = energy.energy_cal_intra(indivisual.c,g);
////            intramolecular_energy = energy.AD4_energy_cal_intra(indivisual.c);
////            indivisual.objectives[k] = energy.energy_cal_forcefield(indivisual.c);
////            indivisual.objectives[k] = intramolecular_energy;
////            indivisual.objectives[k] = energy.energy_cal_autodock4(sf_ad4,indivisual.c);
////            indivisual.objectives[k] = energy.energy_cal_vdw(sf_ad4,indivisual.c);
//            intramolecular_energy = energy.energy_cal_intra(indivisual.c,g);
//            indivisual.objectives[k] = energy.energy_cal_vina(sf, indivisual.c, intramolecular_energy,g);
//
//        }else if(k==1){
////                tmp.objectives[k] = energy.energy_cal_hydrophobic(tmp.c);
////            indivisual.objectives[k] = energy.energy_cal_inter(sf, indivisual.c, intramolecular_energy,g)/(1+0.05846*(m.ligand_degrees_of_freedom(0)));
////            indivisual.objectives[k] = energy.AD4_energy_cal_inter(sf, indivisual.c, intramolecular_energy);
//            float bind_score = energy.energy_cal_xscore(&indivisual);
//            indivisual.objectives[k] = bind_score*(-1.364);
//        }else if(k==2){
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
////            indivisual.objectives[k] = energy.energy_cal_rmsd(ref, indivisual.c);
////            indivisual.objectives[k] = energy.energy_cal_Smog2016(indivisual);
//            indivisual.objectives[k] = energy.energy_cal_DLIGAND2(indivisual);
//        }else if(k==3){
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
////            indivisual.objectives[k] = energy.energy_cal_rmsd(ref, indivisual.c);
////            indivisual.objectives[k] = energy.energy_cal_Smog2016(indivisual);
//            indivisual.objectives[k] = energy.energy_cal_vdw(sf_ad4,indivisual.c);
//        }
//    }
//}
//
//fl SMPSO_aux::velocityConstriction(fl v, int dim_index){
//
//    fl result;
//
//    fl dmax = deltaMax[dim_index];
//    fl dmin = deltaMin[dim_index];
//
//    result = v;
//
//    if(v > dmax){
//        result = dmax;
//    }
//
//    if(v < dmin){
//        result = dmin;
//    }
//
//    return result;
//
//
//}
//
//// constriction coefficient (M. Clerc)
//fl SMPSO_aux::constrictionCoefficient(double c1, double c2) {
//    double rho = c1 + c2;
//    //rho = 1.0 ;
//    if (rho <= 4) {
//        return 1.0;
//    } else {
//        return 2 / (2 - rho - sqrt(pow(rho, 2.0) - 4.0 * rho));
//    }
//} // constrictionCoefficient
//
//
////拥挤度距离比较器，用来根据拥挤距离对同一rank的个体进行排序
//bool crowding_distance_comparator(const output_type_SMPSO& x, const output_type_SMPSO& y) {
//    return x.crowding_distance > y.crowding_distance;
//}
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
//void SMPSO_aux::computeSpeed(model& m,  const precalculate& p, const precalculate& prec_forcefield, const igrid& ig,
//                             const precalculate& p_widened,  const vec& corner1, const vec& corner2,
//                             rng& generator,const scoring_function& sf,const boost::optional<model>& ref){
//    for(int i = 0;i < popsize;i++){
//        int a1 = random_int(0, archive.size() - 1, generator);
//        int a2 = random_int(0, archive.size() - 1, generator);
//        output_type_SMPSO archive1 = archive[a1];
//        output_type_SMPSO archive2 = archive[a2];
//
//        if(crowding_distance_comparator(archive1,archive2)) gbest = archive1; else gbest = archive2;
////        std::sort(archive.begin(),archive.end(), crowding_distance_comparator);
////        gbest = archive[0];
//        r1 = random_fl(0,1.0,generator);
//        r2 = random_fl(0,1.0,generator);
//        c1 = random_fl(1.5, 2.5, generator);
//        c2 = random_fl(1.5, 2.5, generator);
//
//        for(int j = 0;j < dim;j++){
//            if (j < 3) {
//                velocity[i].c.ligands[0].rigid.position[j] = velocityConstriction(constrictionCoefficient(c1, c2) * w*velocity[i].c.ligands[0].rigid.position[j] + c1 * r1 * (pbest[i].c.ligands[0].rigid.position[j] - pop[i].c.ligands[0].rigid.position[j])
//                                                                                  + c2 * r2 * (gbest.c.ligands[0].rigid.position[j] - pop[i].c.ligands[0].rigid.position[j]),j);
////                if (velocity[i].c.ligands[0].rigid.position[j] < corner1[j]) {
////					velocity[i].c.ligands[0].rigid.position[j] = corner1[j];
////				}
////				else if (velocity[i].c.ligands[0].rigid.position[j] > corner2[j]) {
////					velocity[i].c.ligands[0].rigid.position[j] = corner2[j];
////				}
//            }
//            else if (j >= 3 && j < 6) {
//                velocity[i].rotor_angle[j - 3] = velocityConstriction(constrictionCoefficient(c1, c2) * w*velocity[i].rotor_angle[j - 3] + c1 * r1 * (pbest[i].rotor_angle[j - 3] - pop[i].rotor_angle[j - 3])
//                                                                      + c2 * r2 * (gbest.rotor_angle[j - 3] - pop[i].rotor_angle[j - 3]),j);
//            }
//            else{
//                velocity[i].c.ligands[0].torsions[j - 6] = velocityConstriction(constrictionCoefficient(c1, c2) * w*velocity[i].c.ligands[0].torsions[j - 6] + c1 * r1 * (pbest[i].c.ligands[0].torsions[j - 6] - pop[i].c.ligands[0].torsions[j - 6])
//                                                                                + c2 * r2 * (gbest.c.ligands[0].torsions[j - 6] - pop[i].c.ligands[0].torsions[j - 6]),j);
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
//void SMPSO_aux::computeNewPositions(const vec& corner1, const vec& corner2) {
//    for(int i = 0;i < popsize;i++){
//        for(int j = 0;j < dim;j++){
//            if (j < 3) {
//                pop[i].c.ligands[0].rigid.position[j] = pop[i].c.ligands[0].rigid.position[j] + velocity[i].c.ligands[0].rigid.position[j];
//                if (pop[i].c.ligands[0].rigid.position[j] < corner1[j]) {
//                    pop[i].c.ligands[0].rigid.position[j] = corner1[j];
//                    velocity[i].c.ligands[0].rigid.position[j] *= -1;
//                }
//                else if (pop[i].c.ligands[0].rigid.position[j] > corner2[j]) {
//                    pop[i].c.ligands[0].rigid.position[j] = corner2[j];
//                    velocity[i].c.ligands[0].rigid.position[j] *= -1;
//                }
//            }
//            else if(j >= 3 && j < 6){
//                if(j == 3){
//                    pop[i].rotor_angle[j - 3] = pop[i].rotor_angle[j - 3] + velocity[i].rotor_angle[j - 3];
//                    if(pop[i].rotor_angle[j - 3] < 0 || pop[i].rotor_angle[j - 3] >2*pi){
//                        velocity[i].rotor_angle[j - 3] *= -1;
//                    }
//
//                }else if(j == 4){
//                    pop[i].rotor_angle[j - 3] = pop[i].rotor_angle[j - 3] + velocity[i].rotor_angle[j - 3];
//                    if(pop[i].rotor_angle[j - 3] < 0 || pop[i].rotor_angle[j - 3] > pi){
//                        velocity[i].rotor_angle[j - 3] *= -1;
//                    }
//                }else{
//                    pop[i].rotor_angle[j - 3] = pop[i].rotor_angle[j - 3] + velocity[i].rotor_angle[j - 3];
//                    if(pop[i].rotor_angle[j - 3] < 0 || pop[i].rotor_angle[j - 3] > 2*pi){
//                        velocity[i].rotor_angle[j - 3] *= -1;
//                    }
//                }
//
//            }
//            else{
//
//                pop[i].c.ligands[0].torsions[j - 6] = pop[i].c.ligands[0].torsions[j - 6] + velocity[i].c.ligands[0].torsions[j - 6];
//                if(pop[i].c.ligands[0].torsions[j - 6] < -pi || pop[i].c.ligands[0].torsions[j - 6] > pi){
//                    velocity[i].c.ligands[0].torsions[j - 6] *= -1;
//                }
//                normalize_angle(pop[i].c.ligands[0].torsions[j - 6]);
//            }
//        }
//        std::vector<double> quater; // 初始化四元数变量
//        quater.resize(4);
//        MyNormalize(pop[i].rotor_angle);
//        threeangles_to_quaternion(pop[i].rotor_angle, quater);
//
//        fl cos_theta1 = cos(quater[0]);
//        fl sin_theta1 = sin(quater[0]);
//        qt q(cos_theta1, quater[1] * sin_theta1, quater[2] * sin_theta1, quater[3] * sin_theta1);
//
//        pop[i].c.ligands[0].rigid.orientation = q;
//    }
//}
//
//void SMPSO_aux::mutation(const vec& corner1, const vec& corner2,rng& generator) {
//    double rnd, delta1, delta2, mut_pow, deltaq;
//    double y, yl, yu, val, xy, initalval,eta = 20.0;
//    fl p_mut = 1.0/(dim);
//    for(int i = 0;i < popsize;i++){
////        if((i%6) == 0){
//        for (int j = 0; j < dim; j++) {
//            if (random_fl(0, 1.0, generator) <= p_mut ) {
//                if (j < 3) {
//                    y = pop[i].c.ligands[0].rigid.position[j];
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
//                        pop[i].c.ligands[0].rigid.position[j] = initalval;
//                    else
//                        pop[i].c.ligands[0].rigid.position[j] = y;
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
//                    y = pop[i].rotor_angle[j - 3];
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
//                        pop[i].rotor_angle[j - 3] = initalval;
//                    else
//                        pop[i].rotor_angle[j - 3] = y;
//                }
//                else {
//                    y = pop[i].c.ligands[0].torsions[j - 6];
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
//                        pop[i].c.ligands[0].torsions[j - 6] = initalval;
//                    else
//                        pop[i].c.ligands[0].torsions[j - 6] = y;
//                    normalize_angle(pop[i].c.ligands[0].torsions[j - 6]);
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
//        MyNormalize(pop[i].rotor_angle);
//        threeangles_to_quaternion(pop[i].rotor_angle, quater);
//
//        fl cos_theta1 = cos(quater[0]);
//        fl sin_theta1 = sin(quater[0]);
//        qt q(cos_theta1, quater[1] * sin_theta1, quater[2] * sin_theta1, quater[3] * sin_theta1);
//
//        pop[i].c.ligands[0].rigid.orientation = q;
//        }
////    }
//
//
//
//}
//
//int dominated_individual_archive( output_type_SMPSO x, output_type_SMPSO y){
//    if(x < y){
//        return 1;
//    }else if(y < x){
//        return -1;
//    }else{
//        return 0;
//    }
//
//}
//
//void SMPSO_aux::calculate_crowding_distance(std::vector<output_type_SMPSO>& pop){
//    int length = pop.size();
//    for (int i = 0; i < length; i++) {
//        pop[i].crowding_distance = 0.0;
//    }
//
//    for (int i = 0; i < pop[0].objectives.size(); i++) {
//        std::sort(pop.begin(), pop.end(), [i](const output_type_SMPSO& a, const output_type_SMPSO& b) {
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
////记得改archivesize!!!!!!!!!!!!!!!!!!!
//int SMPSO_aux::indexWorst(){
//    if(archive.empty()){
//        return -1;
//    }
//    int index = 0;
//    output_type_SMPSO worstKnown = archive[0];
//    output_type_SMPSO candidateSolution;
//
//    int flag;
//    for(int i = 1;i < archive.size();i++){
//        candidateSolution = archive[i];
//        flag = crowding_distance_comparator(worstKnown,candidateSolution);
//        if(flag){
//            index = i;
//            worstKnown =candidateSolution;
//        }
//    }
//    return index;
//}
//bool SMPSO_aux::add(output_type_SMPSO x){
//    int i = 0;
//    while (i < archive.size()){
//        int flag = dominated_individual_archive(x,archive[i]);
//        if(flag == -1){
//            return false;
//        }else if(flag == 1){
//            archive.erase(archive.begin()+i);
//        }else{
//            if(x == archive[i]){
//                return false;
//            }
//            i++;
//        }
//
//
//    }
//
//    bool res = true;
//    if(archive.size() > archiveMaxSize){
//        calculate_crowding_distance(archive);
//        int indexWorst_ = indexWorst();
//        if(x == archive[indexWorst_]){//??
//            res = false;
//        }else{
//            archive.erase(archive.begin()+indexWorst_);
//            archive.push_back(x);
//            return res;
//        }
//    }
//    archive.push_back(x);
//    return res;
//}
//void SMPSO_aux::updateArchive() {
//    for(int i = 0;i < popsize;i++){
//        bool isAdded = add(pop[i]);
////        if(isAdded == false){
////            pop.erase(pop.begin()+i);
////        }
//
//    }
//}
//
//void SMPSO_aux::update_particle_memory(rng& generator) {
//    for(int i = 0;i < popsize;i++){
//        int flag = dominated_individual_archive(pop[i],pbest[i]);
//        if(flag != -1){
//            if(flag == 1){
//                pbest[i] = pop[i];
//            }else if(flag == 0){
//                if(random_fl(0, 1.0, generator) <= 0.5){
//                    pbest[i] = pop[i];
//                }
//            }
//        }
//    }
//}
//
//bool areAlmostEqual(float a, float b, float tolerance) {
//    return std::fabs(a - b) < tolerance;
//}
//float floors(float number){
//    return roundf(number * 1000) / 1000.0f;
//}
//void SMPSO_aux::initpop(model& m,  model& m_ad4,  const precalculate& p, const precalculate& prec_ad4, const precalculate& prec_vdw, const precalculate& prec_forcefield, const igrid& ig,
//                        const precalculate& p_widened,  const vec& corner1, const vec& corner2,
//                        rng& generator,const scoring_function& sf,const scoring_function& sf_ad4,const boost::optional<model>& ref){
//    vec v(10, 10, 10);
//    conf_size s = m.get_size();
//    energy_cal energy(&m, &m_ad4, &p, &prec_ad4, &prec_vdw, &prec_forcefield,&ig, v);
////    SMPSO_aux aux(m);
//    change g(s);
//    std::vector<fl> obj(3, 0.0);
//    extern Ligand * liganda;
//    extern OBMolecule obl;
//
//
//    extern string lig_name;
//    string pdb_name = "../example/"+lig_name+"/"+lig_name+"_protein.pdb";
//    string mol2_name = "../example/"+ lig_name+"/"+lig_name+"_ligand.mol2";
//
//
//    Molecule_DLIGAND2 *mol = new Molecule_DLIGAND2(pdb_name);
//    mol->rdmol2(mol2_name);
//    output_type_SMPSO tmp(s, obj,*liganda,*mol);
//
//
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
//
//    int num = 0;
//    VINA_FOR(i, m_atoms.size()){
//        if(num == tmp.lig_num){
//            break;
//        }
//        for(int j =  tmp.atoms.size()-1; j > 0; j--){
//            if(areAlmostEqual(m.get_ligand_coords()[i][0], tmp.atoms[j].x[0], 0.01) && areAlmostEqual(m.get_ligand_coords()[i][1], tmp.atoms[j].x[1], 0.01) && areAlmostEqual(m.get_ligand_coords()[i][2], tmp.atoms[j].x[2], 0.01)){
//                tmp.atoms[j].ad_num = m_atoms[i].ad_number;
//                num++;
//                break;
//            }
//        }
//    }
//    init_tmp = tmp;
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
////                intra = energy.energy_cal_intra(tmp.c,g);
////                tmp.objectives[k] = energy.energy_cal_force(tmp.c);
////                intra = energy.AD4_energy_cal_intra(tmp.c);
//                //tmp.objectives[k] = intra;
////                tmp.objectives[k] = energy.energy_cal_forcefield(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_autodock4(sf_ad4,tmp.c);
////                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4,tmp.c);
//                intra = energy.energy_cal_intra(tmp.c,g);
//                tmp.objectives[k] = energy.energy_cal_vina(sf, tmp.c, intra,g);
//            }else if(k==1){
////                tmp.objectives[k] = energy.energy_cal_hydrophobic(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_inter(sf, tmp.c, intra,g)/(1+0.05846*(m.ligand_degrees_of_freedom(0)));
////                tmp.objectives[k] = energy.AD4_energy_cal_inter(sf, tmp.c, intra);
//                float bind_score = energy.energy_cal_xscore(&tmp);
//                tmp.objectives[k] = bind_score*(-1.364);
//
//            }else if(k==2){
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
////                tmp.objectives[k] = energy.energy_cal_Smog2016(tmp);
//                tmp.objectives[k] = energy.energy_cal_DLIGAND2(tmp);
//
//            }else if(k==3){
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
////                tmp.objectives[k] = energy.energy_cal_Smog2016(tmp);
//                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4,tmp.c);
//
//            }
//        }
//
//
//        pop[i] = tmp;
//
//    }
//}
//
//
//void SMPSO_aux::initVelocity(model& m,  model& m_ad4, const precalculate& p, const precalculate& prec_ad4, const precalculate& prec_vdw, const precalculate& prec_forcefield, const igrid& ig,
//                             const precalculate& p_widened,  const vec& corner1, const vec& corner2,
//                             rng& generator,const scoring_function& sf,const scoring_function& sf_ad4,const boost::optional<model>& ref) {
//    vec v(10, 10, 10);
//    conf_size s = m.get_size();
//    energy_cal energy(&m, &m_ad4, &p, &prec_ad4, &prec_vdw, &prec_forcefield,&ig, v);
////    SMPSO_aux aux(m);
//    change g(s);
//    std::vector<fl> obj(3, 0.0);
//    extern Ligand * liganda;
//    output_type_SMPSO tmp = init_tmp;
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
////                intra = energy.energy_cal_intra(tmp.c,g);
////                tmp.objectives[k] = energy.energy_cal_force(tmp.c);
////                intra = energy.AD4_energy_cal_intra(tmp.c);
//                //tmp.objectives[k] = intra;
////                tmp.objectives[k] = energy.energy_cal_forcefield(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_autodock4(sf_ad4,tmp.c);
////                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4,tmp.c);
//                intra = energy.energy_cal_intra(tmp.c,g);
//                tmp.objectives[k] = energy.energy_cal_vina(sf, tmp.c, intra,g);
//            }else if(k==1){
////                tmp.objectives[k] = energy.energy_cal_hydrophobic(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_inter(sf, tmp.c, intra,g)/(1+0.05846*(m.ligand_degrees_of_freedom(0)));
////                tmp.objectives[k] = energy.AD4_energy_cal_inter(sf, tmp.c, intra);
//                float bind_score = energy.energy_cal_xscore(&tmp);
//                tmp.objectives[k] = bind_score*(-1.364);
//
//            }else if(k==2){
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
////                tmp.objectives[k] = energy.energy_cal_Smog2016(tmp);
////                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4,tmp.c);
//                tmp.objectives[k] = energy.energy_cal_DLIGAND2(tmp);
//            }else if(k==3){
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
////                tmp.objectives[k] = energy.energy_cal_Smog2016(tmp);
//                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4,tmp.c);
//            }
//        }
////        atomv m_atoms = m.get_ligand_atoms();
////        VINA_FOR(i, m_atoms.size()){
////            VINA_FOR(j, tmp.num_atom){
////                if(m_atoms[i].coords[0] == tmp.atom[j].coor[0] && m_atoms[i].coords[1] == tmp.atom[j].coor[1] && m_atoms[i].coords[2] == tmp.atom[j].coor[2]){
////                    tmp.atom[j].ad_number = m_atoms[i].ad_number;
////                    continue;
////                }
////            }
////        }
//        velocity[i] = tmp;
//    }
//
//}
//
//void SMPSO_aux::initpbest(model& m,  model& m_ad4, const precalculate& p, const precalculate& prec_ad4, const precalculate& prec_vdw, const precalculate& prec_forcefield, const igrid& ig,
//                          const precalculate& p_widened,  const vec& corner1, const vec& corner2,
//                          rng& generator,const scoring_function& sf,const scoring_function& sf_ad4,const boost::optional<model>& ref){
//    vec v(10, 10, 10);
//    conf_size s = m.get_size();
//    energy_cal energy(&m, &m_ad4, &p, &prec_ad4, &prec_vdw, &prec_forcefield,&ig, v);
////    SMPSO_aux aux(m);
//    change g(s);
//    std::vector<fl> obj(3, 0.0);
//    extern Ligand * liganda;
//    output_type_SMPSO tmp = init_tmp;
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
////                intra = energy.energy_cal_intra(tmp.c,g);
////                tmp.objectives[k] = energy.energy_cal_force(tmp.c);
////                intra = energy.AD4_energy_cal_intra(tmp.c);
//                //tmp.objectives[k] = intra;
////                tmp.objectives[k] = energy.energy_cal_forcefield(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_autodock4(sf_ad4,tmp.c);
////                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4,tmp.c);
//                intra = energy.energy_cal_intra(tmp.c,g);
//                tmp.objectives[k] = energy.energy_cal_vina(sf, tmp.c, intra,g);
//
//            }else if(k==1){
////                tmp.objectives[k] = energy.energy_cal_hydrophobic(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_inter(sf, tmp.c, intra,g)/(1+0.05846*(m.ligand_degrees_of_freedom(0)));
////                tmp.objectives[k] = energy.AD4_energy_cal_inter(sf, tmp.c, intra);
//                float bind_score = energy.energy_cal_xscore(&tmp);
//                tmp.objectives[k] = bind_score*(-1.364);
//
//            }else if(k==2){
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
////                tmp.objectives[k] = energy.energy_cal_Smog2016(tmp);
//                tmp.objectives[k] = energy.energy_cal_DLIGAND2(tmp);
//
//            }else if(k==3){
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
////                tmp.objectives[k] = energy.energy_cal_Smog2016(tmp);
//                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4,tmp.c);
//            }
//        }
////        atomv m_atoms = m.get_ligand_atoms();
////        VINA_FOR(i, m_atoms.size()){
////            VINA_FOR(j, tmp.num_atom){
////                if(m_atoms[i].coords[0] == tmp.atom[j].coor[0] && m_atoms[i].coords[1] == tmp.atom[j].coor[1] && m_atoms[i].coords[2] == tmp.atom[j].coor[2]){
////                    tmp.atom[j].ad_number = m_atoms[i].ad_number;
////                }
////            }
////        }
//        pbest[i] = tmp;
//    }
//}
//
//void SMPSO_aux::initArchive(model &m, model& m_ad4, const precalculate &p, const precalculate& prec_ad4, const precalculate& prec_vdw, const precalculate &prec_forcefield, const igrid &ig,
//                            const precalculate &p_widened, const vec &corner1, const vec &corner2, rng &generator,
//                            const scoring_function &sf, const scoring_function& sf_ad4,const boost::optional<model> &ref) {
//    vec v(10, 10, 10);
//    conf_size s = m.get_size();
//    energy_cal energy(&m, &m_ad4, &p, &prec_ad4, &prec_vdw, &prec_forcefield,&ig, v);
////    SMPSO_aux aux(m);
//    change g(s);
//    std::vector<fl> obj(3, 0.0);
//    extern Ligand * liganda;
//    output_type_SMPSO tmp = init_tmp;
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
////                intra = energy.energy_cal_intra(tmp.c,g);
////                tmp.objectives[k] = energy.energy_cal_force(tmp.c);
////                intra = energy.AD4_energy_cal_intra(tmp.c);
//                //tmp.objectives[k] = intra;
////                tmp.objectives[k] = energy.energy_cal_forcefield(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_autodock4(sf_ad4,tmp.c);
////                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4,tmp.c);
//                intra = energy.energy_cal_intra(tmp.c,g);
//                tmp.objectives[k] = energy.energy_cal_vina(sf, tmp.c, intra,g);
//            }else if(k==1){
////                tmp.objectives[k] = energy.energy_cal_hydrophobic(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_inter(sf, tmp.c, intra,g)/(1+0.05846*(m.ligand_degrees_of_freedom(0)));
////                tmp.objectives[k] = energy.AD4_energy_cal_inter(sf, tmp.c, intra);
//                float bind_score = energy.energy_cal_xscore(&tmp);
//                tmp.objectives[k] = bind_score*(-1.364);
//
//            }else if(k==2){
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
////                tmp.objectives[k] = energy.energy_cal_Smog2016(tmp);
//                tmp.objectives[k] = energy.energy_cal_DLIGAND2(tmp);
//
//            }else if(k==3){
////                tmp.objectives[k] = energy.energy_cal_hydrogen(tmp.c);
////                tmp.objectives[k] = energy.energy_cal_rmsd(ref, tmp.c);
////                tmp.objectives[k] = energy.energy_cal_Smog2016(tmp);
//                tmp.objectives[k] = energy.energy_cal_vdw(sf_ad4,tmp.c);
//            }
//        }
////        atomv m_atoms = m.get_ligand_atoms();
////        VINA_FOR(i, m_atoms.size()){
////            VINA_FOR(j, tmp.num_atom){
////                if(m_atoms[i].coords[0] == tmp.atom[j].coor[0] && m_atoms[i].coords[1] == tmp.atom[j].coor[1] && m_atoms[i].coords[2] == tmp.atom[j].coor[2]){
////                    tmp.atom[j].ad_number = m_atoms[i].ad_number;
////                }
////            }
////        }
//        archive[i] = tmp;
//    }
//}
//void SMPSO::operator()(model& m, model& m_ad4, output_container_SMPSO& out, const precalculate& p, const precalculate& prec_ad4, const precalculate& prec_vdw, const precalculate& prec_forcefield, const igrid& ig,
//                       const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2,
//                       incrementable* increment_me, rng& generator, const scoring_function& sf,const scoring_function& sf_ad4,const boost::optional<model>& ref) const {
//
//    //创建SMPSO_aux对象
//    SMPSO_aux aux(m,corner1,corner2);
//    vec v(10, 10, 10);
//    conf_size s = m.get_size();
//    energy_cal energy(&m, &m_ad4, &p, &prec_ad4, &prec_vdw, &prec_forcefield,&ig, v);
//    change g(s);
//
////    std::vector<fl> obj(3, 0.0);
////
////    extern Ligand * liganda;
////    output_type_SMPSO tmp(s, obj, *liganda);
////
////    std::cout<<tmp.num_atom<<"个"<<std::endl;
//////    VINA_FOR(i,tmp.num_atom)
//////        std::cout<<tmp.atom[i].name<< std::endl;
//////    std::cout<<"====================="<<std::endl;
////
////    VINA_FOR(i,tmp.num_atom){
////        std::cout<<tmp.atom[i].name<< std::endl;
////        std::cout<<tmp.atom[i].type<< std::endl;
////        std::cout<<tmp.atom[i].xtype<< std::endl;
////
////    }
////    std::cout<<"====================="<<std::endl;
//////
////        std::cout<<"-----------------------"<<std::endl;
////    VINA_FOR(i,m.get_ligand_atoms().size()){
////        sz a = m.get_ligand_atoms()[i].get(m.atom_typing_used());
////        std::cout<<"atom_typing_used: "<<m.atom_typing_used()<<std::endl;
////        std::cout<<"a: "<<a<<std::endl;
////    }
////        std::cout<<"-----------------------"<<std::endl;
//
//    //创建初始化种群
//    aux.initpop(m,m_ad4,p,prec_ad4,prec_vdw,prec_forcefield,ig,p_widened,corner1,corner2,generator,sf,sf_ad4,ref);
//
//    //初始化每个粒子的速度
//    aux.initVelocity(m,m_ad4,p,prec_ad4,prec_vdw,prec_forcefield,ig,p_widened,corner1,corner2,generator,sf,sf_ad4,ref);
//
//    //初始化存档粒子
//    aux.initArchive(m,m_ad4,p,prec_ad4,prec_vdw,prec_forcefield,ig,p_widened,corner1,corner2,generator,sf,sf_ad4,ref);
//
//    //把符合要求的粒子放入存档中
//    aux.updateArchive();
//
//    //初始化每个粒子的pbest
//    aux.initpbest(m,m_ad4,p,prec_ad4,prec_vdw,prec_forcefield,ig,p_widened,corner1,corner2,generator,sf,sf_ad4,ref);
//
//    //给存档每个粒子分配拥挤距离
//    aux.calculate_crowding_distance(aux.archive);
//
//    //开始迭代
//    for (int g1 = 0; g1 < num_steps; g1++) {
//        if (increment_me)
//            ++(*increment_me);
//        aux.computeSpeed(m,p,prec_forcefield,ig,p_widened,corner1,corner2,generator,sf,ref);
//
//        aux.computeNewPositions(corner1,corner2);
//
//        aux.mutation(corner1, corner2, generator);
//
//        for(int i = 0; i < aux.popsize; i++){
//            aux.calculate_objectives(m, p, ig, v, aux.pop[i], sf,sf_ad4,energy,ref,g);
//        }
//
//        aux.updateArchive();
//
//        aux.update_particle_memory(generator);
//
//        aux.calculate_crowding_distance(aux.archive);
//
//
//    }
//
//    for(auto& s : aux.archive){
//        s.e = s.objectives[0]+s.objectives[1]+s.objectives[2];
//    }
//
//    for (int i = 0; i < aux.archive.size()-1; i++) {
//        for (int j = 0; j < aux.archive.size() - 1 - i; j++) {
//            if (aux.archive[j+1].e < aux.archive[j].e) {
//                output_type_SMPSO tmp = aux.archive[j];
//                aux.archive[j] = aux.archive[j + 1];
//                aux.archive[j + 1] = tmp;
//            }
//        }
//    }
//
//    for(auto& p : aux.archive){
//        m.set(p.c);
//        p.coords = m.get_heavy_atom_movable_coords();
//        add_to_output_container(out,p,min_rmsd,num_saved_mins);
//
//    }
//
//
//}