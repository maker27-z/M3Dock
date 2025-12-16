#ifndef VINA_ENERGY_H
#define VINA_ENERGY_H

#include "model.h"
#include "scoring_function.h"
#include "coords.h"
#include <time.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <map>
#include <unordered_map>

#include "ad4cache.h"
//Input *input = NULL;
//ForceField *ff = NULL;
//Protein *protein = NULL;
//Ligand *liganda = NULL;
#define NbrProteinAtomType 30
#define NbrLigandAtomType 14
//char timestring[256];
struct energy_cal {
	model* m;//��ʹ��ָ��
    model* m_ad4;
	const precalculate* p;
    const precalculate* prec_ad4;
    const precalculate* prec_vdw;
    const precalculate* prec_forcefield;
    const precalculate* prec_DLIGAND2;
	const igrid* ig;
    const igrid* DLIGAND_grid;
	const vec v;
    ad4cache m_ad4grid;
//    Parameters parameters;
    float KBP_1st_Shell[NbrProteinAtomType][NbrLigandAtomType];
    float KBP_2nd_Shell[NbrProteinAtomType][NbrLigandAtomType];
    float KBP_3rd_Shell[NbrProteinAtomType][NbrLigandAtomType];
    float sizeShell_1;
    float sizeShell_2;
    float sizeShell_3;
    float range;

    energy_cal(model* m_, model* m_ad4_, const precalculate* p_, const precalculate* prec_ad4_, const precalculate* prec_vdw_, const precalculate* prec_forcefield_, const precalculate* prec_DLIGAND2_, const igrid* ig_, const igrid* DLIGAND_grid_, const vec& v_) : m(m_), m_ad4(m_ad4_), p(p_), prec_ad4(prec_ad4_), prec_vdw(prec_vdw_), prec_forcefield(prec_forcefield_), prec_DLIGAND2(prec_DLIGAND2_), ig(ig_), DLIGAND_grid(DLIGAND_grid_), v(v_) {
                extern string lig_name;
        //AD4
        string map = "./"+lig_name+"/"+lig_name+"_protein";
        load_maps(map);

//        parameters.energyFile="./parameter/KBP-3.0-5.0-8.5.dat";
//        sizeShell_1=0.0;
//        sizeShell_2=0.0;
//        sizeShell_3=0.0;
//        range = 0.1;

//=================================
//        extern string lig_name;
//        //AD4
//        string map = "../coreset_pdbqt/"+lig_name+"/"+lig_name+"_protein";
//        load_maps(map);

//        parameters.energyFile="../KBP-3.0-5.0-8.5.dat";
//        sizeShell_1=0.0;
//        sizeShell_2=0.0;
//        sizeShell_3=0.0;
//        range = 0.1;


//

//        for(int i=0; i<NbrProteinAtomType; i++)  {
//            for(int j=0; j<NbrLigandAtomType; j++) { KBP_1st_Shell[i][j]=0; }
//        }
//        for(int i=0; i<NbrProteinAtomType; i++)  {
//            for(int j=0; j<NbrLigandAtomType; j++) { KBP_2nd_Shell[i][j]=0; }
//        }
//        for(int i=0; i<NbrProteinAtomType; i++)  {
//            for(int j=0; j<NbrLigandAtomType; j++) { KBP_3rd_Shell[i][j]=0; }
//        }
//
//
//        // Open the energy file.
//        ifstream KBP(parameters.energyFile.c_str(), ios::in) ;
//        if(!KBP) {
//            cerr << "ERROR1: The file " << parameters.energyFile << " is missing, or I can't read it." << endl;
//            exit (-1);
//        }
//        // Size of the shells. It is read from the first line of the energyFile. This code works for 3 shells.
//
//        KBP >> sizeShell_1 >> sizeShell_2 >> sizeShell_3;
//
//        // First shell
//        for (int i=0; i<NbrProteinAtomType; i++) {
//            for(int j=0; j<NbrLigandAtomType; j++) { KBP >> KBP_1st_Shell[i][j]; }
//        }
//        // Second shell
//        for (int i=0; i<NbrProteinAtomType; i++) {
//            for(int j=0; j<NbrLigandAtomType; j++) { KBP >> KBP_2nd_Shell[i][j]; }
//        }
//        // Third shell
//        for (int i=0; i<NbrProteinAtomType; i++) {
//            for(int j=0; j<NbrLigandAtomType; j++) { KBP >> KBP_3rd_Shell[i][j]; }
//        }
//        KBP.close();

    }
//    void Smog2016_set_conf(output_type_MODEA &outputTypeMODEA){
////        int a = outputTypeMODEA.obatom.size();
////        int b = m->get_ligand_atoms().size();
//        m_ad4->set(outputTypeMODEA.c);
//        atomv m_atoms = m_ad4->get_ligand_atoms();
//        std::unordered_map<int, int> atomIndexMap;
//        for (int i = 0;i<m_atoms.size();i++) {
//            atomIndexMap[m_atoms[i].ad_number] = i;
//        }
//
//// 在 tmp 中遍历原子并查找匹配项
//        for (int j = 0; j < outputTypeMODEA.obatom.size(); j++) {
//            int ad_number = outputTypeMODEA.atom[j].ad_number;
//            auto it = atomIndexMap.find(ad_number);
//            if (it != atomIndexMap.end()) {
//                // 找到匹配项，更新坐标
//                int i = it->second;
//                outputTypeMODEA.obatom[j].coordinates[0] = m_ad4->coords[i][0];
//                outputTypeMODEA.obatom[j].coordinates[1] = m_ad4->coords[i][1];
//                outputTypeMODEA.obatom[j].coordinates[2] = m_ad4->coords[i][2];
//            }
////        VINA_FOR(i, m_atoms.size()){
////            VINA_FOR(j, outputTypeMODEA.obatom.size()){
////                if(m_atoms[i].ad_number == outputTypeMODEA.atom[j].ad_number){
////                    outputTypeMODEA.obatom[j].coordinates[0] = m_ad4->get_ligand_coords()[i][0];
////                    outputTypeMODEA.obatom[j].coordinates[1] = m_ad4->get_ligand_coords()[i][1];
////                    outputTypeMODEA.obatom[j].coordinates[2] = m_ad4->get_ligand_coords()[i][2];
////                    break;
////                }
////            }
////        }
//        }
//    }

    // This function computes the energy between a protein and a ligand. It calls other functions defined below.
//    double KBP20161(Parameters const & parameters, OBProtein const & protein, output_type_MODEA & molecule)
//    {
//        // Initialize parameters. In SMoG2016 we use a continuous description of the contacts. For example, if one shell starts at 4.5A,
//        // if d<4.5-range there is a contact of 0, if d>4.5+range there is a contact of 1 (until the shell stops) and between the contact
//        // value varies linearly: contact=(4.5+range-distance)/(2*range).
//
//
//        // Open and read the KBP.dat file.
//        //   L I G A N D
//        // P
//        // R
//        // O
//        // T
//        // E
//        // I
//        // N
//
//
//
//        // Compute the KBP energy
//        double ScoreShell1=0;
//        double ScoreShell2=0;
//        double ScoreShell3=0;
//        for (unsigned int i=0 ; i < protein.obatom.size() ; i++) {
//            for (unsigned int j=0 ; j < molecule.obatom.size() ; j++) {
//                if   ( (protein.obatom[i].TypeNumber == 0) || (molecule.obatom[j].TypeNumber == 0) ) { }
//                else {
//                    double distance = dist(protein.obatom[i].coordinates, molecule.obatom[j].coordinates);
//                    // Shell1
//                    if      (distance <= sizeShell_1-range) {
//                        ScoreShell1 += KBP_1st_Shell[protein.obatom[i].TypeNumber-1][molecule.obatom[j].TypeNumber-1];
//                    }
//                        // Shell1-Shell2
//                    else if (distance <= sizeShell_1+range) {
//                        double c1 = (sizeShell_1+range-distance)/(2*range);
//                        double c2 = (distance-sizeShell_1+range)/(2*range);
//                        ScoreShell1 += c1*KBP_1st_Shell[protein.obatom[i].TypeNumber-1][molecule.obatom[j].TypeNumber-1] ;
//                        ScoreShell2 += c2*KBP_2nd_Shell[protein.obatom[i].TypeNumber-1][molecule.obatom[j].TypeNumber-1];
//                    }
//                        // Shell2
//                    else if (distance <= sizeShell_2-range) {
//                        ScoreShell2 += KBP_2nd_Shell[protein.obatom[i].TypeNumber-1][molecule.obatom[j].TypeNumber-1];
//                    }
//                        // Shell2-Shell3
//                    else if (distance <= sizeShell_2+range) {
//                        double c1 = (sizeShell_2+range-distance)/(2*range);
//                        double c2 = (distance-sizeShell_2+range)/(2*range);
//                        ScoreShell2 += c1*KBP_2nd_Shell[protein.obatom[i].TypeNumber-1][molecule.obatom[j].TypeNumber-1];
//                        ScoreShell3 += c2*KBP_3rd_Shell[protein.obatom[i].TypeNumber-1][molecule.obatom[j].TypeNumber-1];
//                    }
//                        // Shell3
//                    else if (distance <= sizeShell_3-range) {
//                        ScoreShell3 += KBP_3rd_Shell[protein.obatom[i].TypeNumber-1][molecule.obatom[j].TypeNumber-1];
//                    }
//                        // Shell3-EndOfShell3
//                    else if (distance <= sizeShell_3+range) {
//                        double c1 = (sizeShell_3+range-distance)/(2*range);
//                        ScoreShell3 += c1*KBP_3rd_Shell[protein.obatom[i].TypeNumber-1][molecule.obatom[j].TypeNumber-1];
//                    }
//                }
//            }
//        }
//        // We sum the contributions of the three shells. Each score is divided by the average distance in the shells.
//        double KBPScore = ScoreShell1/((0.0+sizeShell_1)/2) + ScoreShell2/((sizeShell_1+sizeShell_2)/2) + ScoreShell3/((sizeShell_2+sizeShell_3)/2);
//
//        return KBPScore;
//    }
    float energy_cal_Smog2016(output_type_MODEA &outputTypeMODEA){

//        Smog2016_set_conf(outputTypeMODEA);
//        extern OBProtein obp;
//        float energy = 0.032*KBP20161(parameters, obp, outputTypeMODEA);
//        return energy;

    }
    void load_maps(std::string maps) {
        const fl slope = 1e6; // FIXME: too large? used to be 100
        grid_dims gd;


        ad4cache grid(slope);
        grid.read(maps);
        m_ad4grid = grid;

//        // Check that all the affinity map are present for ligands/flex residues (if initialized already)
//        if (m_ligand_initialized) {
//            atom_type::t atom_typing = m_scoring_function.get_atom_typing();
//            szv atom_types = m_model.get_movable_atom_types(atom_typing);
//
//            if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
//                if(!m_grid.are_atom_types_grid_initialized(atom_types))
//                    exit(EXIT_FAILURE);
//            } else {
//                if(!m_ad4grid.are_atom_types_grid_initialized(atom_types))
//                    exit(EXIT_FAILURE);
//            }
//        }

        // Store in Vina object
//        m_map_initialized = true;
    }

    float energy_cal_autodock4(const scoring_function& f,conf& c){
        //        struct timeval t1,t2;
//        double timeuse;
//        gettimeofday(&t1,NULL);
        m_ad4->set(c);
        fl lig_grids = m_ad4grid.eval(*m_ad4,v[1]);
        fl inter_pairs = m_ad4->eval_inter(*prec_ad4,v);
        fl inter = lig_grids + inter_pairs;

        // Intra
        fl flex_grids = m_ad4grid.eval_intra(*m_ad4, v[1]); // [1] flex -- grid
        fl intra_pairs = m_ad4->evalo(*prec_ad4,v); // [1] flex_i -- flex_i and flex_i -- flex_j
        fl lig_intra = m_ad4->evali(*prec_ad4,v); // [2] ligand_i -- ligand_i
        fl intra = flex_grids + intra_pairs + lig_intra;

        // Torsion
        fl conf_independent = f.conf_independent(*m_ad4, 0); // [3] we can pass e=0 because we do not modify the energy like in vina
        // Total
        fl total = inter + conf_independent; // (+ intra - intra)
        //        gettimeofday(&t2,NULL);
//        timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
//
//        cout<<"vina_time = "<<timeuse<<endl;  //输出时间（单位：ｓ）
        return total;


    }
    float energy_cal_autodock4(const scoring_function& f,conf& c,change& ad4_g){
        m_ad4->set(c);

        fl e = m_ad4->ad4_eval_deriv(*prec_ad4, m_ad4grid, v, c, ad4_g);
        return e;


    }
    float energy_cal_autodock4_std(const scoring_function& f,conf& c,model m_){
        //        struct timeval t1,t2;
//        double timeuse;
//        gettimeofday(&t1,NULL);
        fl lig_grids = m_ad4grid.eval(m_,v[1]);
        fl inter_pairs = m_.eval_inter(*prec_ad4,v);
        fl inter = lig_grids + inter_pairs;

        // Intra
//        fl flex_grids = m_ad4grid.eval_intra(*m_ad4, v[1]); // [1] flex -- grid
//        fl intra_pairs = m_ad4->evalo(*prec_ad4,v); // [1] flex_i -- flex_i and flex_i -- flex_j
//        fl lig_intra = m_ad4->evali(*prec_ad4,v); // [2] ligand_i -- ligand_i
//        fl intra = flex_grids + intra_pairs + lig_intra;

        // Torsion
        fl conf_independent = f.conf_independent(m_, 0); // [3] we can pass e=0 because we do not modify the energy like in vina
        // Total
        fl total = inter + conf_independent; // (+ intra - intra)
        //        gettimeofday(&t2,NULL);
//        timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
//
//        cout<<"vina_time = "<<timeuse<<endl;  //输出时间（单位：ｓ）
        return total;


    }
    float energy_cal_vdw(const scoring_function& f,conf& c){
//                struct timeval t1,t2;
//        double timeuse;
//        gettimeofday(&t1,NULL);
        m_ad4->set(c);
        fl lig_grids = m_ad4grid.eval(*m_ad4,v[1]);
        fl inter_pairs = m_ad4->eval_inter(*prec_vdw,v);
        fl inter = lig_grids + inter_pairs;

        // Intra
        fl flex_grids = m_ad4grid.eval_intra(*m_ad4, v[1]); // [1] flex -- grid
        fl intra_pairs = m_ad4->evalo(*prec_vdw,v); // [1] flex_i -- flex_i and flex_i -- flex_j
        fl lig_intra = m_ad4->evali(*prec_vdw,v); // [2] ligand_i -- ligand_i
        fl intra = flex_grids + intra_pairs + lig_intra;

        // Torsion
//        fl conf_independent = f.conf_independent(*m_ad4, 0); // [3] we can pass e=0 because we do not modify the energy like in vina
        // Total
        fl total = inter + intra; // (+ intra - intra)
//                gettimeofday(&t2,NULL);
//        timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
//
//        cout<<"vdw_time = "<<timeuse<<endl;  //输出时间（单位：ｓ）
        return total;
    }

    float energy_cal_vdw(const scoring_function& f,conf& c, change &vdw_g){
////                struct timeval t1,t2;
////        double timeuse;
////        gettimeofday(&t1,NULL);
//        m_ad4->set(c);
//        fl lig_grids = m_ad4grid.eval(*m_ad4,v[1]);
//        fl inter_pairs = m_ad4->eval_inter(*prec_vdw,v);
//        fl inter = lig_grids + inter_pairs;
//
//        // Intra
//        fl flex_grids = m_ad4grid.eval_intra(*m_ad4, v[1]); // [1] flex -- grid
//        fl intra_pairs = m_ad4->evalo(*prec_vdw,v); // [1] flex_i -- flex_i and flex_i -- flex_j
//        fl lig_intra = m_ad4->evali(*prec_vdw,v); // [2] ligand_i -- ligand_i
//        fl intra = flex_grids + intra_pairs + lig_intra;
//
//        // Torsion
////        fl conf_independent = f.conf_independent(*m_ad4, 0); // [3] we can pass e=0 because we do not modify the energy like in vina
//        // Total
//        fl total = inter + intra; // (+ intra - intra)
////                gettimeofday(&t2,NULL);
////        timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
////
////        cout<<"vdw_time = "<<timeuse<<endl;  //输出时间（单位：ｓ）
//        return total;
// INTER ligand - grid
        m_ad4->set(c);
        fl e = m_ad4->vdw_eval_deriv(*prec_vdw, m_ad4grid, v, c, vdw_g);

        return e;
    }

    fl operator()(const conf& c, change& g) {
		const fl tmp = m->eval_deriv(*p, *ig, v, c, g);
		return tmp;
	}
    //计算vina分子内能量
    fl energy_cal_intra(const conf& c,change &g){
//        fl e = m->eval_intramolecular_deriv(*p,*ig,v,c,g);
        fl e = m->eval_intramolecular(*p,v,c);
        return e;
    }

    fl energy_cal_inter(const scoring_function& f,conf& c,fl intramolecular_energy,change &g){

        fl e = m->eval_inter_deriv(*p,*ig,v,c,g);
        return e;
    }

    fl energy_cal_vina(const scoring_function& f,conf& c,fl intramolecular_energy,change &g){
//        struct timeval t1,t2;
//        double timeuse;
//        gettimeofday(&t1,NULL);
        fl tmp = m->eval_adjusted(f,*p, *ig, v, c, intramolecular_energy,g);
//        gettimeofday(&t2,NULL);
//        timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
//
//        cout<<"vina_time = "<<timeuse<<endl;  //输出时间（单位：ｓ）
        return tmp;
    }
    fl energy_cal_vina_binding(const scoring_function& f,conf& c,fl intramolecular_energy,fl vina_energy,change &vina_g){
//        struct timeval t1,t2;
//        double timeuse;
//        gettimeofday(&t1,NULL);
//        fl tmp = m->eval_adjusted(f,*p, *ig, v, c, intramolecular_energy,g);
        fl binding_energy = f.conf_independent(*m, vina_energy - intramolecular_energy);
//        gettimeofday(&t2,NULL);
//        timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
//
//        cout<<"vina_time = "<<timeuse<<endl;  //输出时间（单位：ｓ）
        return binding_energy;
    }
    fl energy_cal_vina_deriv(const scoring_function& f,conf& c,fl intramolecular_energy,change &vina_g){
//        struct timeval t1,t2;
//        double timeuse;
//        gettimeofday(&t1,NULL);
//        fl tmp = m->eval_adjusted(f,*p, *ig, v, c, intramolecular_energy,g);
        fl tmp = m->eval_deriv(*p, *ig, v, c, vina_g);
        fl e = f.conf_independent(*m, tmp - intramolecular_energy);
//        gettimeofday(&t2,NULL);
//        timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
//
//        cout<<"vina_time = "<<timeuse<<endl;  //输出时间（单位：ｓ）
        return e;
    }

    mat xscore_set_orientation(const qt& q) { // does not normalize the orientation

        mat orientation_m = quaternion_to_r3(q);
        return orientation_m;
    }

     vec local_to_lab(const vec& origin, const mat& orientation_m,float * coor) const {
        vec local_coords(coor[1],coor[2],coor[3]);
        vec tmp;
        tmp = origin + orientation_m*local_coords;

         return tmp;
    }

    void set_coords(struct atom *atoms,const vec& origin, const mat& orientation_m, output_type_MODEA* tmp) const {
        VINA_FOR(i, tmp->num_atom) {
            vec vec_tmp;
            vec_tmp = local_to_lab(origin, orientation_m, atoms[i].coor);
            tmp->atom[i].coor[0] = vec_tmp.data[0];
            tmp->atom[i].coor[1] = vec_tmp.data[1];
            tmp->atom[i].coor[2] = vec_tmp.data[2];
        }

    }

    std::string vina_remark(std::vector<fl> objectives, fl lb, fl ub) {
        std::ostringstream remark;
        remark.setf(std::ios::fixed, std::ios::floatfield);
        remark.setf(std::ios::showpoint);
        remark << "REMARK VINA RESULT: "
               << std::setw(9) << std::setprecision(1) << objectives[0]<< std::setw(9) << std::setprecision(1)<<objectives[1]
               << "  " << std::setw(9) << std::setprecision(3) << lb
               << "  " << std::setw(9) << std::setprecision(3) << ub
               << '\n';
        return remark.str();
    }
    path make_path(const std::string& str) {
        return path(str);
    }
    void write_all_output(model& m, const output_container_MODEA & out, sz how_many,
                          const std::string& output_name,
                          const std::vector<std::string>& remarks) {
        if(out.size() < how_many)
            how_many = out.size();
        VINA_CHECK(how_many <= remarks.size());
        ofile f(make_path(output_name));
        VINA_FOR(i, how_many) {
            m.set(out[i].c);
            m.write_model(f, i+1, remarks[i]); // so that model numbers start with 1
        }
    }

    void xscore_set_conf(struct atom *atoms, output_type_MODEA* tmp) {

//        vec origin = tmp->c.ligands[0].rigid.position;
//        mat orientation_m = xscore_set_orientation(tmp->c.ligands[0].rigid.orientation);
//
//        set_coords(atoms, origin, orientation_m, tmp);


//        m->set(tmp->c);
//        atomv m_atoms = m->get_ligand_atoms();
//        VINA_FOR(i, tmp->num_atom){
//            VINA_FOR(j, m_atoms.size()){
//                if(tmp->atom[i].name == m_atoms[j].ad_name){
//                    tmp->atom[i].coor[0] = m_atoms[j].coords[0];
//                    tmp->atom[i].coor[1] = m_atoms[j].coords[1];
//                    tmp->atom[i].coor[2] = m_atoms[j].coords[2];
//                }
//            }
//        }



//        vec v(10, 10, 10);
//        conf_size s = m->get_size();
////    MODEA_aux aux(m);
//        change g(s);
//        scoring_function *sf;
//        tmp->objectives[0] = energy_cal_inter(*sf, tmp->c, 123,g)/(1+0.05846*(m->ligand_degrees_of_freedom(0)));
//        std::cout<<tmp->objectives[0]<<"______________"<<tmp->objectives[1]<<std::endl;
//        output_container_MODEA out;
//
//        tmp->coords = m->get_heavy_atom_movable_coords();
//        add_to_output_container(out,*tmp,2,20);
//        std::string out_name = "../coreset_pdbqt/1a30/1a30_docking_out.pdbqt";
//        std::vector<std::string> remarks;
//        remarks.push_back(vina_remark(tmp->objectives, 1.23000, 1.23000));
//        write_all_output(*m,out,1,out_name,remarks);



        extern Protein* protein;
//        int tmp_protein = protein->num_atom;
//        int m_protein = m->get_protein_atoms().size();
//        int m_num_movable_atoms = m->num_movable_atoms();
        m_ad4->set(tmp->c);
        std::unordered_map<int, int> atomIndexMap;

//        int count = 0;
        atomv m_atoms = m_ad4->get_ligand_atoms();
        for (int i = 0;i<m_atoms.size();i++) {
            atomIndexMap[m_atoms[i].ad_number] = i;
        }
        for (int j = 0; j < tmp->num_atom; ++j) {
            int ad_number = tmp->atom[j].ad_number;
            auto it = atomIndexMap.find(ad_number);
            if (it != atomIndexMap.end()) {
                int i = it->second;
                // 找到匹配项，更新坐标
                tmp->atom[j].coor[0] = m_ad4->coords[i][0];
                tmp->atom[j].coor[1] = m_ad4->coords[i][1];
                tmp->atom[j].coor[2] = m_ad4->coords[i][2];
            }
//        VINA_FOR(i, m_atoms.size()){
//            VINA_FOR(j, tmp->num_atom){
//               if(m_atoms[i].ad_number == tmp->atom[j].ad_number){
//                   tmp->atom[j].coor[0] = m_ad4->coords[i][0];
//                   tmp->atom[j].coor[1] = m_ad4->coords[i][1];
//                   tmp->atom[j].coor[2] = m_ad4->coords[i][2];
////                   count++;
//                   break;
//               }
//            }
//        }
        }

    }


    float Xtool_Score_Shortcut(output_type_MODEA* tmp,int num_molecule){


//        std::cout<<liganda->num_atom<<"��"<<std::endl;
//        std::cout<<m->get_ligand_atoms().size()<<"��"<<std::endl;
//        std::cout<<m->get_protein_atoms().size()<<"��"<<std::endl;
//        std::cout<<m->get_ligand_coords().size()<<"��"<<std::endl;
//        std::cout<<m->get_heavy_atom_movable_coords().size()<<"��"<<std::endl;
//        struct timeval t1,t2;
//        double timeuse;
//        gettimeofday(&t1,NULL);
        extern Ligand* liganda;
        xscore_set_conf(tmp->atom, tmp);
        if(num_molecule==1)	// single mol2 file
        {
//            liganda->Calculate_Binding_Score();

//            tmp.coo
            float result = tmp->Calculate_Binding_Score();
//            gettimeofday(&t2,NULL);
//            timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
//
//            cout<<"xscore_time = "<<timeuse<<endl;  //输出时间（单位：ｓ）
            return result;
//            printf("\n***********************************************\n");
//
//            printf("HPScore -log(Kd) = %-5.2f\n", liganda->pkd1);
//            printf("HMScore -log(Kd) = %-5.2f\n", liganda->pkd2);
//            printf("HSScore -log(Kd) = %-5.2f\n", liganda->pkd3);
//            printf("Predicted average -log(Kd) = %-5.2f\n", liganda->bind_score);
//
//            printf("Predicted binding energy = %-6.2f kcal/mol\n",
//                   liganda->bind_score*(-1.364));
//
//            printf("***********************************************\n\n");
//            printf("\n***********************************************\n");
//
//            printf("HPScore -log(Kd) = %-5.2f\n", tmp->pkd1);
//            printf("HMScore -log(Kd) = %-5.2f\n", tmp->pkd2);
//            printf("HSScore -log(Kd) = %-5.2f\n", tmp->pkd3);
//            printf("Predicted average -log(Kd) = %-5.2f\n", tmp->bind_score);
//
//            printf("Predicted binding energy = %-6.2f kcal/mol\n",
//                   tmp->bind_score*(-1.364));
//
//            printf("***********************************************\n\n");

        }
//    else	// multi-mol2 file
//    {
//        printf("ID: HPSCORE HMSCORE HSSCORE AVE_SCORE BIND_ENERGY NAME\n");
//
//        if((fin=fopen(input->ligand_file,"r"))==NULL) Open_File_Error(input->ligand_file);
//
//        if((fout=fopen(input->log_file,"w"))==NULL) Open_File_Error(input->log_file);
//
//        for(i=1;i<=num_molecule;i++)
//        {
//            liganda->Clear();
//            liganda->Read_From_Mol2(fin);
//            liganda->Value_Atom();
//            liganda->Calculate_Binding_Score();
//
//            printf("Molecule %6d: ", i);
//            printf("%5.2f  ", liganda->pkd1);
//            printf("%5.2f  ", liganda->pkd2);
//            printf("%5.2f  ", liganda->pkd3);
//            printf("%5.2f  ", liganda->bind_score);
//            printf("%6.2f  ", liganda->bind_score*(-1.364));
//            printf("%s\n", liganda->name);
//
//            liganda->Write_Out_Log(fout);
//        }
//
//        fclose(fin); fclose(fout);
//    }

    }

    float energy_cal_xscore(output_type_MODEA *outputTypenMODEA){
//        fl e  =
//xscore
        extern Input *input;
        int num_molecule;

        num_molecule=Check_Mol2_File(input->ligand_file);
        return Xtool_Score_Shortcut(outputTypenMODEA,num_molecule);
    }


    float energy_cal_DLIGAND2(output_type_MODEA &outputTypeMODEA){
//        struct timeval t1,t2;
//        double timeuse;
//        gettimeofday(&t1,NULL);
        DLIGAND2_set_conf(outputTypeMODEA);
//                gettimeofday(&t2,NULL);
//        timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
//
//        cout<<"DLIGAND_time = "<<timeuse<<endl;  //输出时间（单位：ｓ）
        float energy = outputTypeMODEA.score();
//        gettimeofday(&t2,NULL);
//        timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
//
//        cout<<"DLIGAND_time = "<<timeuse<<endl;  //输出时间（单位：ｓ）
        return energy;

    }

    float energy_cal_DLIGAND2(output_type_MODEA &outputTypeMODEA,change &DLIGAND2_g){
//        struct timeval t1,t2;
//        double timeuse;
//        gettimeofday(&t1,NULL);
        //跟setc()函数的作用一样，但是DLIGAND2是以167种原子类型为原型，因此要自己实现一下
        DLIGAND2_set_conf(outputTypeMODEA);
        //计算DLIGAND2能量评分函数
        float energy = DLIGAND_grid->eval_deriv_DLIGAND2(outputTypeMODEA, v[1]);
//        vecv coords;
//        coords.resize(outputTypeMODEA.lig_num);
//        for(int i = 0; i < outputTypeMODEA.lig_num; i++){
//            VINA_FOR_IN(j,coords[i]){
//                coords[i][j] = outputTypeMODEA.getx(i+outputTypeMODEA.natm-outputTypeMODEA.lig_num)[j];
//            }
//        }
        //初始化DLIGAND2的梯度
        vecv minus_forces(m_ad4->get_m_num_movable_atoms());

        //得到所有配体原子
        atomv m_atoms = m_ad4->get_ligand_atoms();
        //蛋白质和配体原子都记录在m_atoms中，前面一部分是蛋白的原子，剩余后面的一部分是配体原子，因此想遍历配体原子需要先得到配体原子的初始索引，具体为什么要这样做请详细看DLIGAND2论文中及其代码
        int start_lig = outputTypeMODEA.atoms.size() - outputTypeMODEA.lig_num;
//        int count = 0;

        //得到DLIGAND2梯度
        int flag;
        for(int i = 0; i < m_atoms.size();i++){
            flag = 1;
            for(int j = start_lig; j < outputTypeMODEA.natm; j++){
                if(m_atoms[i].ad_number == outputTypeMODEA.atoms[j].ad_num){
                    minus_forces[i] = prec_DLIGAND2->minus_forces[j - outputTypeMODEA.natm + outputTypeMODEA.lig_num];
                    flag=0;
                    break;
                }
            }
            if(flag == 1){
                minus_forces[i] = vec(0,0,0);
            }
        }
        //更新DLIGAND2梯度
        m_ad4->get_m_ligands().derivative(m_ad4->coords, minus_forces, DLIGAND2_g.ligands);
//                gettimeofday(&t2,NULL);
//        timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
//
//        cout<<"DLIGAND_time = "<<timeuse<<endl;  //输出时间（单位：ｓ）
//        float energy = outputTypeMODEA.score();
//        gettimeofday(&t2,NULL);
//        timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
//
//        cout<<"DLIGAND_time = "<<timeuse<<endl;  //输出时间（单位：ｓ）
        return energy;

    }

    void DLIGAND2_set_conf(output_type_MODEA &outputTypeMODEA){
//        int a = outputTypeMODEA.obatom.size();
//        int b = m->get_ligand_atoms().size();
//        m_ad4->set(outputTypeMODEA.c);
//        atomv m_atoms = m_ad4->get_ligand_atoms();
//        int num = 0;
//        VINA_FOR(i, m_atoms.size()){
//            if(num == outputTypeMODEA.lig_num){
//                break;
//            }
//            for(int j =  outputTypeMODEA.atoms.size()-1; j > 0; j--){
//
//                if(m_atoms[i].ad_number == outputTypeMODEA.atoms[j].ad_num){
//                    outputTypeMODEA.getx(j)[0] = m_ad4->get_ligand_coords()[i][0];
//                    outputTypeMODEA.getx(j)[1] = m_ad4->get_ligand_coords()[i][1];
//                    outputTypeMODEA.getx(j)[2] = m_ad4->get_ligand_coords()[i][2];
//                    num++;
//                    break;
//                }
//            }
//        }





        m_ad4->set(outputTypeMODEA.c);
        atomv m_atoms = m_ad4->get_ligand_atoms();
        std::unordered_map<int, int> atomIndexMap;
        int start_lig = outputTypeMODEA.atoms.size() - outputTypeMODEA.lig_num;
//        int count = 0;
        for (int i = 0;i<m_atoms.size();i++) {
            atomIndexMap[m_atoms[i].ad_number] = i;
        }
        for (int j = start_lig; j < outputTypeMODEA.atoms.size(); j++) {
            int ad_number = outputTypeMODEA.atoms[j].ad_num;
            auto it = atomIndexMap.find(ad_number);
            if (it != atomIndexMap.end()) {
                int i = it->second;
                // 找到匹配项，更新坐标
                outputTypeMODEA.getx(j)[0] = m_ad4->coords[i][0];
                outputTypeMODEA.getx(j)[1] = m_ad4->coords[i][1];
                outputTypeMODEA.getx(j)[2] = m_ad4->coords[i][2];
            }
        }
    }
    fl energy_cal_rmsd(const boost::optional<model>& r,const conf& c){
        m->set(c);
        fl rmsd = m->rmsd_upper_bound(r.get());
        return rmsd;
    }
    fl energy_cal_forcefield(conf& c){
        fl e = m->eval(*prec_forcefield,*ig,v,c);
        return e;
    }
};

#endif
