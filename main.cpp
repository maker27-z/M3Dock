/*
main
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

#include <iostream>
#include<fstream>
#include <string>
#include <exception>
#include <vector> // ligand paths
#include <cmath> // for ceila
#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/convenience.hpp> // filesystem::basename
#include <boost/thread/thread.hpp> // hardware_concurrency // FIXME rm ?
#include "parse_pdbqt.h"
#include "parallel_mc.h"
#include "file.h"
#include "cache.h"
#include "non_cache.h"
#include "naive_non_cache.h"
#include "parse_error.h"
#include "everything.h"
#include "weighted_terms.h"
#include "current_weights.h"
#include "quasi_newton.h"
#include "tee.h"
#include "coords.h" // add_to_output_container
//#include "boost/python.hpp"
//#include "python3.8
//#include "python3.8/Python.h"
#include "xtools.h"
#include "matplotlibcpp.h"
using boost::filesystem::path;

void Stringsplit(std::string str,const  char split,std::vector<std::string>& list)
{
    std::istringstream iss(str);	//
    std::string token;			//    ջ
    while (getline(iss, token, split))	//   splitΪ ָ
    {
        if(token != ""){
            list.push_back(token);

        }
    }
}
////3 obj
//void draw_pareto(const std::vector<output_container_MODEA>& out_cont, int s){
//    // Iterate through the output containers
//    VINA_FOR_IN(j, out_cont) {
//        std::vector<std::vector<double>> objectives;
//        // Collect objectives from the output container
//        for (int i = 0; i < out_cont[j].size(); i++) {
//            objectives.emplace_back(out_cont[j][i].objectives.begin(), out_cont[j][i].objectives.end());
//        }
//        // Initialize Python and import necessary modules
//        Py_Initialize();
//        PyRun_SimpleString("import matplotlib\n"
//                           "matplotlib.use('Agg')\n"
//                           "from mpl_toolkits.mplot3d import Axes3D\n");
//
//        // Get pyplot interface
//        PyObject* pyplot = PyImport_ImportModule("matplotlib.pyplot");
//        // Create a new figure and add a 3d subplot
//        PyObject* fig = PyObject_CallMethod(pyplot, "figure", NULL);
//        PyObject* ax = PyObject_CallMethod(fig, "add_subplot", "iii", 111, "projection", "3d");
//
//        // Prepare data for matplotlib
//        PyObject* x = PyList_New(objectives.size());
//        PyObject* y = PyList_New(objectives.size());
//        PyObject* z = PyList_New(objectives.size());
//        for (int i = 0; i < objectives.size(); ++i) {
//            PyList_SetItem(x, i, PyFloat_FromDouble(objectives[i][0]));
//            PyList_SetItem(y, i, PyFloat_FromDouble(objectives[i][1]));
//            PyList_SetItem(z, i, PyFloat_FromDouble(objectives[i][2]));
//        }
//        // Set plot title and labels
//        PyObject* title = PyUnicode_FromString("Pareto Front");
//        PyObject* xlabel = PyUnicode_FromString("Objective 1");
//        PyObject* ylabel = PyUnicode_FromString("Objective 2");
//        PyObject* zlabel = PyUnicode_FromString("Objective 3");
//
//        PyObject_CallMethod(ax, "set_title", "O", title);
//        PyObject_CallMethod(ax, "set_xlabel", "O", xlabel);
//        PyObject_CallMethod(ax, "set_ylabel", "O", ylabel);
//        PyObject_CallMethod(ax, "set_zlabel", "O", zlabel);
//
//        // Plot the pareto front in 3D
//        PyObject_CallMethod(ax, "scatter3D", "OOO", x, y, z);
//
//        // Save the figure to a file
//        extern std::string lig_name;
//        PyObject_CallMethod(pyplot, "savefig", "s", ("../coreset_pdbqt/"+lig_name+"/pareto_front"+"_"+std::to_string(j)+".png").c_str());
//
//        // Show the plot (optional, usually not needed when saving to a file)
//        // PyObject_CallMethod(pyplot, "show", NULL);
//
//        // Decrement reference counts if necessary
//        // Py_DECREF(x);
//        // Py_DECREF(y);
//        // Py_DECREF(z);
//        // Py_DECREF(pyplot);
//        // Py_Finalize();
//    }
//}


//void draw_pareto(const output_container_MODEA& out_cont){
//    // Iterate through the output containers
//
//        std::vector<std::vector<double>> objectives;
//        // Collect objectives from the output container
//        for (int i = 0; i < out_cont.size(); i++) {
//            objectives.emplace_back(out_cont[i].objectives.begin(), out_cont[i].objectives.end());
//        }
//
//        // Prepare data for matplotlib
//        std::vector<double> x, y, z;
//        for (const auto& obj : objectives) {
//            x.push_back(obj[0]);
//            y.push_back(obj[1]);
//            z.push_back(obj[2]);
//        }
//
//        // Create a new figure
//        matplotlibcpp::figure();
//
//        // Plot the pareto front in 3D
//        matplotlibcpp::scatter(x, y, z);
//
//        // Set plot title and labels
//        matplotlibcpp::title("Pareto Front");
//        matplotlibcpp::xlabel("Objective 1");
//        matplotlibcpp::ylabel("Objective 2");
////        matplotlibcpp::zlabel("Objective 3");
//        matplotlibcpp::set_zlabel("Z Label");
//        // Save the figure to a file
//        extern std::string lig_name;
//        matplotlibcpp::save("./" + lig_name + "/pareto_front.png");
//
//}
//void draw_pareto(const std::vector<output_container_MODEA>& out_cont, int s) {
//    // Initialize Python environment
//    Py_Initialize();
//
//    VINA_FOR_IN(j, out_cont) {
//        PyObject* matplotlib = PyImport_ImportModule("matplotlib");
//        PyObject* pyplot = PyImport_ImportModule("matplotlib.pyplot");
//
//        PyObject* fig = PyObject_CallMethod(pyplot, "figure", NULL);
//        PyObject* ax = PyObject_CallMethod(fig, "add_subplot", "iii", 1, 1, 1);  // Creates a 1x1 grid and selects the first subplot
//
//
//        PyObject* x = PyList_New(out_cont[j].size());
//        PyObject* y = PyList_New(out_cont[j].size());
//        PyObject* x_item;
//        PyObject* y_item;
//        PyObject* color;
//        PyObject* colors = PyList_New(out_cont[j].size());
//
//        for (int i = 0; i < out_cont[j].size(); i++) {
//            x_item = PyFloat_FromDouble(out_cont[j][i].objectives[0]);
//            PyList_SetItem(x, i, x_item);
//            y_item = PyFloat_FromDouble(out_cont[j][i].objectives[1]);
//            PyList_SetItem(y, i, y_item);
//
//            if (out_cont[j][i].rm <= 10.0) {
//                color = PyUnicode_FromString("r");
//            } else {
//                color = PyUnicode_FromString("b");
//            }
//            PyList_SetItem(colors, i, color);
//        }
//
//        // Set plot title and axis labels
//        PyObject* title = PyUnicode_FromString("Pareto Front");
//        PyObject* xlabel = PyUnicode_FromString("vina");
//        PyObject* ylabel = PyUnicode_FromString("xscore");
//
//        PyObject* set_title = PyObject_GetAttrString(pyplot, "title");
//        PyObject* set_xlabel = PyObject_GetAttrString(pyplot, "xlabel");
//        PyObject* set_ylabel = PyObject_GetAttrString(pyplot, "ylabel");
//        PyObject_CallFunctionObjArgs(set_title, title, NULL);
//        PyObject_CallFunctionObjArgs(set_xlabel, xlabel, NULL);
//        PyObject_CallFunctionObjArgs(set_ylabel, ylabel, NULL);
//
//        PyObject* np = PyImport_ImportModule("numpy");
//        PyObject* np_array = PyObject_CallMethod(np, "array", "O", colors);
//
//        // Use 2D scatter plot in Python
//        PyObject* scatter = PyObject_GetAttrString(ax, "scatter");
//        PyObject* scatter_args = PyTuple_Pack(3, x, y, colors);
//        PyObject_CallObject(scatter, scatter_args);
//
//        extern std::string lig_name;
//        PyObject_CallMethod(pyplot, "savefig", "s", ("../coreset_pdbqt/"+lig_name+"/pareto_front"+"_"+std::to_string(j)+".png").c_str());
//
//        PyObject_CallMethod(pyplot, "show", nullptr);
//    }
//
//    // Release Python resources, you may need to adjust this based on your needs
//    // ...
//}

//2 obj
//void draw_pareto(const std::vector<output_container_MODEA>& out_cont,int s){
//    //   Ŀ 꺯  ֵ 洢    ά
//    VINA_FOR_IN(j, out_cont) {
//        std::vector<std::vector<double>> objectives;
//        for (int i = 0; i < out_cont[j].size(); i++) {
//            objectives.emplace_back(out_cont[j][i].objectives.begin(), out_cont[j][i].objectives.end());
//
//        }
//        //   ʼ  Python
//        Py_Initialize();
//
//        //     matplotlib Ⲣ л Ϊ ǽ   ʽ
//        PyRun_SimpleString("import matplotlib\n"
//                           "matplotlib.use('Agg')\n");
//
//        //     pyplot ӿ
//        PyObject* pyplot = PyImport_ImportModule("matplotlib.pyplot");
////    Py_DECREF(matplotlib);
//
//        //   Ŀ 꺯  ֵ  Ϊ       ݸ matplotlib
//        PyObject* x = PyList_New(objectives.size());
//        PyObject* y = PyList_New(objectives.size());
//        for (int i = 0; i < objectives.size(); ++i) {
//            PyObject* x_i = PyList_New(objectives[i].size());
//            PyObject* y_i = PyList_New(objectives[i].size());
//            for (int j = 0; j < objectives[i].size(); ++j) {
//                PyList_SetItem(x_i, j, PyFloat_FromDouble(objectives[i][0]));
//                PyList_SetItem(y_i, j, PyFloat_FromDouble(objectives[i][1]));
//            }
//            PyList_SetItem(x, i, x_i);
//            PyList_SetItem(y, i, y_i);
//        }
////        // Set plot title and axis labels
//        PyObject* title = PyUnicode_FromString("Pareto Front");
//        PyObject* xlabel = PyUnicode_FromString("intra");
//        PyObject* ylabel = PyUnicode_FromString("DLIGAND2");
//
//        PyObject* set_title = PyObject_GetAttrString(pyplot, "title");
//        PyObject* set_xlabel = PyObject_GetAttrString(pyplot, "xlabel");
//        PyObject* set_ylabel = PyObject_GetAttrString(pyplot, "ylabel");
//        PyObject_CallFunctionObjArgs(set_title, title, NULL);
//        PyObject_CallFunctionObjArgs(set_xlabel, xlabel, NULL);
//        PyObject_CallFunctionObjArgs(set_ylabel, ylabel, NULL);
//        //   Python  ʹ  matplotlib    paretoǰ  ͼ
//        PyObject_CallMethod(pyplot, "scatter", "OO", x, y);
//        // save the figure to a file
//        extern std::string lig_name;
////        PyObject_CallMethod(pyplot, "savefig", "s", ("../coreset_pdbqt/"+lig_name+"/pareto_front"+"_"+std::to_string(j)+".png").c_str());
//        PyObject_CallMethod(pyplot, "savefig", "s", ("./pareto/"+lig_name+"/"+lig_name+"_"+std::to_string(j)+".png").c_str());
//
//
//        PyObject_CallMethod(pyplot, "show", nullptr);
//
//    }
//
//
////    PyObject_CallMethod(pyplot, "show", nullptr);
//
//    //  ͷ Python    ͻ
////    Py_DECREF(x);
////    Py_DECREF(y);
////    Py_DECREF(pyplot);
////    Py_Finalize();
//}
path make_path(const std::string& str) {
    return path(str);
}
bool sort_comparator(const output_type_MODEA& x,const output_type_MODEA& y);
void doing(int verbosity, const std::string& str, tee& log) {
    if(verbosity > 1) {
        log << str << std::string(" ... ");
        log.flush();
    }
}
bool not_max_multi(const output_type_MODEA& out){
    VINA_FOR_IN(i,out.objectives)
        if(out.objectives[i]>0.1 * max_fl){
            return false;
        }
    return true;
}

//bool energy_range_compare(const output_type_MODEA& a,const output_type_MODEA& b,fl energy){
//    if(a.rank != b.rank){
//        if(a.rank < b.rank){
//            return true;
//        }else{
//            return false;
//        }
//    } else if(a.crowding_distance != b.crowding_distance){
//        if(a.crowding_distance > b.crowding_distance){
//            return true;
//        }else{
//            return false;
//        }
//    } else{
//        return true;
//    }
//}
void done(int verbosity, tee& log) {
    if(verbosity > 1) {
        log << "done.";
        log.endl();
    }
}
std::string default_output(const std::string& input_name) {
    std::string tmp = input_name;
    if(tmp.size() >= 6 && tmp.substr(tmp.size()-6, 6) == ".pdbqt")
        tmp.resize(tmp.size() - 6); // FIXME?
    return tmp + "_out.pdbqt";
}

std::string default_output_data(const std::string& input_name) {
    std::string tmp = input_name;
    if(tmp.size() >= 6 && tmp.substr(tmp.size()-6, 6) == ".pdbqt")
        tmp.resize(tmp.size() - 6); // FIXME?
    return tmp + "_data.txt";
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

//void write_one_output(model& m, const conf& c, sz how_many, const std::string& output_name, const std::string remark) {
//	ofile f(make_path(output_name));
//	m.set(c);
//	m.write_model(f, 1, remark); // so that model numbers start with 1
//}

void do_randomization(model& m,
                      const std::string& out_name,
                      const vec& corner1, const vec& corner2, int seed, int verbosity, tee& log) {
    conf init_conf = m.get_initial_conf();
    rng generator(static_cast<rng::result_type>(seed));
    if(verbosity > 1) {
        log << "Using random seed: " << seed;
        log.endl();
    }
    const sz attempts = 10000;
    conf best_conf = init_conf;
    fl best_clash_penalty = 0;
    VINA_FOR(i, attempts) {
        conf c = init_conf;
        c.randomize(corner1, corner2, generator);
        m.set(c);
        fl penalty = m.clash_penalty();
        if(i == 0 || penalty < best_clash_penalty) {
            best_conf = c;
            best_clash_penalty = penalty;
        }
    }
    m.set(best_conf);
    if(verbosity > 1) {
        log << "Clash penalty: " << best_clash_penalty; // FIXME rm?
        log.endl();
    }
    m.write_structure(make_path(out_name));
}

void refine_structure(model& m, const precalculate& prec, non_cache& nc, output_type& out, const vec& cap, sz max_steps = 1000) {
    change g(m.get_size());
    quasi_newton quasi_newton_par;
    quasi_newton_par.max_steps = max_steps;
    const fl slope_orig = nc.slope;
    sz total_torsions = m.ligand_degrees_of_freedom(0); // ligands    ת   ĸ
    VINA_FOR(p, 5) {
        nc.slope = 100 * std::pow(10.0, 2.0*p);
        quasi_newton_par(m, prec, nc, out, g, cap);
        m.set(out.c); // just to be sure
        if(nc.within(m))
            break;
    }
    out.coords = m.get_heavy_atom_movable_coords();
    if(!nc.within(m))
        out.e = max_fl;
    nc.slope = slope_orig;
}

std::string vina_remark(std::vector<fl> objectives, fl lb, fl ub) {
    std::ostringstream remark;
    remark.setf(std::ios::fixed, std::ios::floatfield);
    remark.setf(std::ios::showpoint);
    remark << "REMARK VINA RESULT: "
           //                        << std::setw(9) << std::setprecision(3) << binding_energy
           << std::setw(9) << std::setprecision(3) << objectives[0]<< std::setw(9) << std::setprecision(3)<<objectives[1]
           << std::setw(9) << std::setprecision(3)<<objectives[2]
           //                        << std::setw(9) << std::setprecision(1)<<objectives[3]
           << std::setw(9) << std::setprecision(3)<<objectives[0]+objectives[1]+objectives[2]
           << "  " << std::setw(9) << std::setprecision(3) << lb
           << "  " << std::setw(9) << std::setprecision(3) << ub
           << '\n';
    return remark.str();
}

output_container_MODEA remove_redundant(const output_container_MODEA& in, fl min_rmsd) {
    output_container_MODEA tmp;
//	VINA_FOR_IN(i, in)
//		add_to_output_container(tmp, in[i], min_rmsd, in.size());
    return tmp;
}

void do_search(model& m, model& m_ad4, const boost::optional<model>& ref, const scoring_function& sf, const scoring_function& sf_ad4, const precalculate& prec, const precalculate& prec_ad4, const precalculate& prec_vdw, const precalculate& prec_forcefield, const precalculate& prec_DLIGAND2, const igrid& ig, const igrid& DLIGAND2_grid, const precalculate& prec_widened, const igrid& ig_widened, non_cache& nc, // nc.slope is changed
               const std::string& out_name,
        //std::ofstream& ofs,
               const vec& corner1, const vec& corner2,
               parallel_mc& par, fl energy_range, sz num_modes,
               int seed, int verbosity, bool score_only, bool local_only, tee& log, const terms& t, const flv& weights) {
    conf_size s = m.get_size();
    conf c = m.get_initial_conf();
    conf_size s_ad4 = m_ad4.get_size();
    conf c_ad4 = m_ad4.get_initial_conf();
    fl e = max_fl;
    const vec authentic_v(10, 10, 10);
    if(score_only) {
//		fl intramolecular_energy = m.eval_intramolecular(prec, authentic_v, c);
//		naive_non_cache nnc(&prec); // for out of grid issues
////		e = m.eval_adjusted(sf, prec, nnc, authentic_v, c, intramolecular_energy);
//		log << "Affinity: " << std::fixed << std::setprecision(5) << e << " (kcal/mol)";
//		log.endl();
//		flv term_values = t.evale_robust(m);
//		VINA_CHECK(term_values.size() == 5);
//		log << "Intermolecular contributions to the terms, before weighting:\n";
//		log << std::setprecision(5);
//		log << "    gauss 1     : " << term_values[0] << '\n';
//		log << "    gauss 2     : " << term_values[1] << '\n';
//		log << "    repulsion   : " << term_values[2] << '\n';
//		log << "    hydrophobic : " << term_values[3] << '\n';
//		log << "    Hydrogen    : " << term_values[4] << '\n';
//		VINA_CHECK(weights.size() == term_values.size() + 1);
//		fl e2 = 0;
//		VINA_FOR_IN(i, term_values)
//			e2 += term_values[i] * weights[i];
//		e2 = sf.conf_independent(m, e2);
//		if(e < 100 && std::abs(e2 - e) > 0.05) {
//			log << "WARNING: the individual terms are inconsisent with the\n";
//			log << "WARNING: affinity. Consider reporting this as a bug:\n";
//			log << "WARNING: http://vina.scripps.edu/manual.html#bugs\n";
//		}
    }
    else if(local_only) {
//		output_type out(c, e);
//		doing(verbosity, "Performing local search", log);
//	//	refine_structure(m, prec, nc, out, authentic_v, par.ls.ssd_par.evals);
//		done(verbosity, log);
//		fl intramolecular_energy = m.eval_intramolecular(prec, authentic_v, out.c);
//		e = m.eval_adjusted(sf, prec, nc, authentic_v, out.c, intramolecular_energy);
//
//		log << "Affinity: " << std::fixed << std::setprecision(5) << e << " (kcal/mol)";
//		log.endl();
//		if(!nc.within(m))
//			log << "WARNING: not all movable atoms are within the search space\n";
//
//		doing(verbosity, "Writing output", log);
//		output_container out_cont;
//		out_cont.push_back(new output_type(out));
//		std::vector<std::string> remarks(1, vina_remark(e, 0, 0));
//		write_all_output(m, out_cont, 1, out_name, remarks); // how_many == 1
//		done(verbosity, log);
    }
    else {
        int iters_ = 1;
        fl avg_distance = 0;

        //const model& r_ = ref.get();
        for (int i = 0; i < iters_; i++) {
            seed = auto_seed();
            rng generator(static_cast<rng::result_type>(seed));
            log << "Using random seed: " << seed;
            log.endl();
//			output_container_MODEA out_cont;
            output_container_MODEA out_cont;
            doing(verbosity, "Performing search", log);
//			clock_t startTime, endTime;
//			startTime = clock();
//            draw_pareto(out_cont);
            struct timeval t1,t2;
            double timeuse;
            gettimeofday(&t1,NULL);
            //主要执行该函数，用cpu的多线同时执行多个任务，进行多次搜索，执行几次取决于config文件中的cpu和exhaustiveness设置为多少
            par(m, m_ad4, out_cont, prec, prec_ad4, prec_vdw, prec_forcefield, prec_DLIGAND2, ig, DLIGAND2_grid, prec_widened, ig_widened, corner1, corner2, generator,sf,sf_ad4,ref);//????
            gettimeofday(&t2,NULL);
            timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
            log<<"run time is: "<<timeuse<<"\n";  //输出时间（单位：ｓ）
//            extern output_container_MODEA out_300;
//            extern output_container_MODEA out_500;
//            extern output_container_MODEA out_800;
//            extern output_container_MODEA out_1000;
////            extern output_container_MODEA out_process_pf;
////            if(out_process_pf.size()>0){
////                output_container_MODEA out_temp;
////                out_temp = out_process_pf;
////                std::unordered_map<int, std::vector<output_type_MODEA>> fast_nondominated_sort(output_container_MODEA& pop);
////                std::unordered_map<int, std::vector<output_type_MODEA>> F = fast_nondominated_sort(out_temp);
////////
////                output_container_MODEA out_results;
////                for(int i=0;i<F[1].size();i++){
////                    out_results.push_back(new output_type_MODEA(F[1][i]));
////                }
////                out_process_pf = out_results;
////                bool compare_e(const output_type_MODEA& a, const output_type_MODEA& b);
////                std::sort(out_process_pf.begin(),out_process_pf.end(),compare_e);
////            }
//            if(out_300.size()>0){
//                output_container_MODEA out_temp;
//                out_temp = out_300;
//                std::unordered_map<int, std::vector<output_type_MODEA>> fast_nondominated_sort(output_container_MODEA& pop);
//                std::unordered_map<int, std::vector<output_type_MODEA>> F = fast_nondominated_sort(out_temp);
//////
//                output_container_MODEA out_results;
//                for(int i=0;i<F[1].size();i++){
//                    out_results.push_back(new output_type_MODEA(F[1][i]));
//                }
//                out_300 = out_results;
//                bool compare_e(const output_type_MODEA& a, const output_type_MODEA& b);
//                std::sort(out_300.begin(),out_300.end(),compare_e);
//            }
//            if(out_500.size()>0){
//                output_container_MODEA out_temp;
//                out_temp = out_500;
//                std::unordered_map<int, std::vector<output_type_MODEA>> fast_nondominated_sort(output_container_MODEA& pop);
//                std::unordered_map<int, std::vector<output_type_MODEA>> F = fast_nondominated_sort(out_temp);
//////
//                output_container_MODEA out_results;
//                for(int i=0;i<F[1].size();i++){
//                    out_results.push_back(new output_type_MODEA(F[1][i]));
//                }
//                out_500 = out_results;
//                bool compare_e(const output_type_MODEA& a, const output_type_MODEA& b);
//                std::sort(out_500.begin(),out_500.end(),compare_e);
//            }
//            if(out_800.size()>0){
//                output_container_MODEA out_temp;
//                out_temp = out_800;
//                std::unordered_map<int, std::vector<output_type_MODEA>> fast_nondominated_sort(output_container_MODEA& pop);
//                std::unordered_map<int, std::vector<output_type_MODEA>> F = fast_nondominated_sort(out_temp);
//////
//                output_container_MODEA out_results;
//                for(int i=0;i<F[1].size();i++){
//                    out_results.push_back(new output_type_MODEA(F[1][i]));
//                }
//                out_800 = out_results;
//                bool compare_e(const output_type_MODEA& a, const output_type_MODEA& b);
//                std::sort(out_800.begin(),out_800.end(),compare_e);
//            }
//            if(out_1000.size()>0){
//                output_container_MODEA out_temp;
//                out_temp = out_1000;
//                std::unordered_map<int, std::vector<output_type_MODEA>> fast_nondominated_sort(output_container_MODEA& pop);
//                std::unordered_map<int, std::vector<output_type_MODEA>> F = fast_nondominated_sort(out_temp);
//////
//                output_container_MODEA out_results;
//                for(int i=0;i<F[1].size();i++){
//                    out_results.push_back(new output_type_MODEA(F[1][i]));
//                }
//                out_1000 = out_results;
//                bool compare_e(const output_type_MODEA& a, const output_type_MODEA& b);
//                std::sort(out_1000.begin(),out_1000.end(),compare_e);
//            }
//            std::vector<output_container_MODEA> runtimes;
//            runtimes.push_back(out_300);
//            runtimes.push_back(out_500);
//            runtimes.push_back(out_800);
//            runtimes.push_back(out_1000);
////            runtimes.push_back(out_process_pf);
//
//            std::vector<string> runtimestr;
//            extern std::string lig_name;
//
//            runtimestr.push_back("./"+lig_name+"/"+lig_name+"_300_out");
//            runtimestr.push_back("./"+lig_name+"/"+lig_name+"_500_out");
//            runtimestr.push_back("./"+lig_name+"/"+lig_name+"_800_out");
//            runtimestr.push_back("./"+lig_name+"/"+lig_name+"_1000_out");
////            runtimestr.push_back("./"+lig_name+"/"+lig_name+"_aprocess_pf_out");
//
//            VINA_FOR_IN(k,runtimes){
//
//            std::ofstream logFile(runtimestr[k]+".log");
//            logFile.setf(std::ios::fixed, std::ios::floatfield);
//            logFile.setf(std::ios::showpoint);
//            logFile<<"Reading input ... done.\n"
//                     "Setting up the scoring function ... done.\n"
//                     "Using random seed: -494562976\n"
//                     "Performing search ... done.\n"
//                     "run time is: \n"
//                     "Refining results ... done.\n";
//            logFile<<"\n";
//            logFile << "mode |    vina     |   DLIGAND2 |     vdw    |     e    | dist from best mode\n";
//            logFile << "     | (kcal/mol)  | (kcal/mol) | (kcal/mol) |  (float) | rmsd l.b.| rmsd u.b.\n";
//            logFile << "-----+------------+---------"+runtimestr[k]+"---+----------+-------------------\n";
////            logFile << "This will be written to output.log" << std::endl;
//            VINA_FOR_IN(i, runtimes[k]) {
//                logFile << std::setw(4) << i + 1
//                    << std::setw(9) << std::setprecision(3) << runtimes[k][i].objectives[0]
//                    << "    " << std::setw(9) << std::setprecision(3) << runtimes[k][i].objectives[1]
//                    << "    " << std::setw(9) << std::setprecision(3) << runtimes[k][i].objectives[2]
//                    //                        << "    " << std::setw(9) << std::setprecision(3) << out_cont[i].objectives[3]
//                    << "    " << std::setw(9) << std::setprecision(3) << runtimes[k][i].e;
//
//                m_ad4.set(runtimes[k][i].c);
//                const model& r = ref.get() ;
//                const fl lb = m_ad4.rmsd_lower_bound(r);
//                const fl ub = m_ad4.rmsd_upper_bound(r);
//                logFile << "  " << std::setw(9) << std::setprecision(3) << lb
//                    << "  " << std::setw(9) << std::setprecision(3) << ub
//                    << "  " << std::setw(9) << std::setprecision(3) << runtimes[k][i].threadd; // FIXME need user-readable error messages in case of failures
//                logFile<<std::endl;
//            }
//
//            logFile.close();
//            }

            //m.set(out_cont[0].c);
            //fl lb_ = m.rmsd_lower_bound(r_);
            //fl ub_ = m.rmsd_upper_bound(r_);
            //ofs << std::setprecision(8) << out_cont[0].e << "   " << std::setprecision(5) << lb_ << "   " << std::setprecision(5) << ub_ << std::endl;

            done(verbosity, log);
//			endTime = clock();//  ʱ
//            log<<"run time is: "<<(double)(endTime - startTime) / CLOCKS_PER_SEC<< "\n" ;  //输出时间（单位：ｓ）

            doing(verbosity, "Refining results", log);
            /*		VINA_FOR_IN(i, out_cont)
                        refine_structure(m, prec, nc, out_cont[i], authentic_v, par.mc.ssd_par.evals);
            */
//			if (!out_cont.empty()) {
//				sort(out_cont.begin(),out_cont.end(), sort_comparator);
//				const fl best_mode_intramolecular_energy = m.eval_intramolecular(prec, authentic_v, out_cont[0].c);
//				VINA_FOR_IN(i, out_cont)
//					if (not_max_multi(out_cont[i]))
//						out_cont[i].e = m.eval_adjusted(sf, prec, nc, authentic_v, out_cont[i].c, best_mode_intramolecular_energy);
//				// the order must not change because of non-decreasing g (see paper), but we'll re-sort in case g is non strictly increasing
//                sort(out_cont.begin(),out_cont.end(), sort_comparator);
//			}
//            if(!out_cont.empty()){
////                out_cont.sort();
////                sort(out_cont.begin(),out_cont.end(), sort_comparator);
//				for (int i = 0; i < out_cont.size() - 1; i++) {
//					for (int j = 0; j < out_cont.size() - 1 - i; j++) {
//						if (sort_comparator(out_cont[j+1], out_cont[j])) {
//							output_type_MODEA tmp = out_cont[j];
//							out_cont[j] = out_cont[j + 1];
//							out_cont[j + 1] = tmp;
//						}
//					}
//				}
//            }

//			const fl out_min_rmsd = 1;
//			out_cont = remove_redundant(out_cont, out_min_rmsd);
            done(verbosity, log);
//
            extern double deviceTime;
            log << "deviceTime: "<<deviceTime;
            log.setf(std::ios::fixed, std::ios::floatfield);
            log.setf(std::ios::showpoint);
//			log << '\n';
//			log << "mode |   affinity | dist from best mode\n";
//			log << "     | (kcal/mol) | rmsd l.b.| rmsd u.b.\n";
//			log << "-----+------------+----------+----------\n";
            log << '\n';
//            log << "mode |   intra    |   inter   |   rmsd   |     e    |  rank  | crowding_dist | dist from best mode\n";
//            log << "     | (kcal/mol) | (kcal/mol)| (double) | (double) |  (int) |   (double)    | rmsd l.b.| rmsd u.b.\n";
//            log << "-----+------------+----------+--------------------------------------------------------------------\n";
//            log << "mode |   intra    |   inter   |     e    |  rank  | crowding_dist | dist from best mode\n";
//            log << "     | (kcal/mol) | (kcal/mol)| (double) |  (int) |   (double)    | rmsd l.b.| rmsd u.b.\n";
//            log << "-----+------------+----------+--------------------------------------------------------------------\n";
//            log << "mode |    ad4     |     xscore    |  DLIGAND2  |     e    | dist from best mode\n";
//            log << "     | (kcal/mol) |  (kcal/mol)   | (kcal/mol) |  (float) | rmsd l.b.| rmsd u.b.\n";
//            log << "-----+------------+--------------+------------+----------+---------------------\n";
            log << "mode |  vina_bind  |  DLIGAND2  |    vdw     |     e    | dist from best mode\n";
            log << "     |  (kcal/mol) | (kcal/mol) | (kcal/mol) |  (float) | rmsd l.b.| rmsd u.b.\n";
            log << "-----+-------------+------------+------------+----------+----------+-----------\n";
            model best_mode_model = m;
            if (!out_cont.empty())
                best_mode_model.set(out_cont.front().c);
            sz how_many = 0;
            std::vector<std::string> remarks;
            fl min_ub = max_fl;
            VINA_FOR_IN(i, out_cont) {
//                fl intra = m.eval_intramolecular(prec,authentic_v,out_cont[i].c);
                // ɾȥ|| out_cont[i].e > out_cont[0].e + energy_range
//                || ((out_cont[i].objectives[0] > out_cont[0].objectives[0] + energy_range) && (out_cont[i].objectives[1] > out_cont[0].objectives[1] + energy_range))
                if ( !not_max_multi(out_cont[i]) ) break; // check energy_range sanity FIXME
                ++how_many;

                log << std::setw(4) << i + 1
                    //                    << std::setw(9) << std::setprecision(3) << out_cont[i].binding_energy
                    << std::setw(9) << std::setprecision(3) << out_cont[i].objectives[0]
                    << "    " << std::setw(9) << std::setprecision(3) << out_cont[i].objectives[1]
                    << "    " << std::setw(9) << std::setprecision(3) << out_cont[i].objectives[2]
                    //                        << "    " << std::setw(9) << std::setprecision(3) << out_cont[i].objectives[3]
                    << "    " << std::setw(9) << std::setprecision(3) << out_cont[i].e;
//                << "    " << std::setw(9) << std::setprecision(3) << out_cont[i].rank;

//                        << "    " << std::setw(9) << std::setprecision(3) << out_cont[i].objectives[2]
//                    << "    " << std::setw(9) << std::setprecision(3) << intra
//                        << "    " << std::setw(9) << std::setprecision(3) << m.eval(prec, ig, authentic_v, out_cont[i].c)  - intra
//                    << "    " << std::setw(3) << std::setprecision(3) <<out_cont[i].rank
//                    << "    " << std::setw(9) << std::setprecision(3) <<out_cont[i].crowding_distance; // intermolecular_energies[i];
                m_ad4.set(out_cont[i].c);
                model tmp_m;
                if(ref){
                    tmp_m = ref.get();
                }else{
                    tmp_m = best_mode_model;
                };
                const model& r = tmp_m;
                const fl lb = m_ad4.rmsd_lower_bound(r);
                const fl ub = m_ad4.rmsd_upper_bound(r);
                out_cont[i].rm = ub;
                if (ub < min_ub)
                    min_ub = ub;
                log << "  " << std::setw(9) << std::setprecision(3) << lb
                    << "  " << std::setw(9) << std::setprecision(3) << ub
                    << "  " << std::setw(9) << std::setprecision(3) << out_cont[i].threadd; // FIXME need user-readable error messages in case of failures

                remarks.push_back(vina_remark(out_cont[i].objectives, lb, ub));
                log.endl();
            }
//            draw_pareto(out_cont);

            avg_distance += min_ub;
            doing(verbosity, "Writing output", log);
            write_all_output(m_ad4, out_cont, how_many, out_name, remarks);
            done(verbosity, log);
//            output_container_MODEA s = out_cont;

            if (how_many < 1) {
                log << "WARNING: Could not find any conformations completely within the search space.\n"
                    << "WARNING: Check that it is large enough for all movable atoms, including those in the flexible side chains.";
                log.endl();
            }
        }
        //avg_distance = avg_distance / iters_;
        //log << "avg_distance: " << avg_distance;
        log.endl();
        extern std::vector<output_container_MODEA> out_tmp;

        int i = 0;
//        draw_pareto(out_tmp,i);
    }
}

void main_procedure(model& m, model& m_ad4, const boost::optional<model>& ref, // m is non-const (FIXME?)
                    const std::string& out_name,
                    bool score_only, bool local_only, bool randomize_only, bool no_cache,
                    const grid_dims& gd, int exhaustiveness,
                    const flv& weights, const flv& weights_ad4, const flv& weights_vdw,
                    int cpu, int seed, int verbosity, sz num_modes, fl energy_range, tee& log) {

    doing(verbosity, "Setting up the scoring function", log);

    //everything结构体，详情见everything.cpp文件中的everything::everything(const std::string term_typename)和everything::everything()函数，主要看这个t和t_vdw，其他的是废案
    everything t;

    everything t_ad4("ad4");
    everything t_forcefield("forcefield");
    everything t_xscore("xscore");

    everything t_vdw("vdw");

    //设定不同评分函数的能量函数项之前的每个系数，主要看weight_vdw和weights
    flv weights_forcefield;
    fl weight_vdw      = 1;
    fl weight_electrostatic      = 1;
    weights_forcefield.push_back(weight_vdw);
    weights_forcefield.push_back(weight_electrostatic);

//	VINA_CHECK(weights.size() == 6);

    //wt是vina的能量权重，wt_vdw是vdw的能量权重
    weighted_terms wt(&t, weights);
    weighted_terms wt_ad4(&t_ad4, weights_ad4);
    weighted_terms wt_vdw(&t_vdw, weights_vdw);
    weighted_terms wt_forcefield(&t_forcefield, weights_forcefield);
    weighted_terms_DLIGAND2 wt_Dligand2;

    //重点，需要得到三个能量评分函数的预计算对象prec_vdw、prec、prec_DLIGAND2，其余的为废案，用于迭代中快速计算能量
    precalculate prec_ad4(wt_ad4,m_ad4.atoms);
    precalculate prec_vdw(wt_vdw,m_ad4.atoms);

    precalculate prec(wt);
    precalculate prec_forcefield(wt_forcefield,m_ad4.atoms);

    precalculate prec_DLIGAND2(wt_Dligand2);
    const fl left  = 0.25;
    const fl right = 0.25;
    precalculate prec_widened(prec);
    prec_widened.widen(left, right);
//    precalculate prec_widened1(prec_ ad4); prec_widened1.widen(left, right);

    done(verbosity, log);
    //确定结合口袋的xyz的空间范围
    vec corner1(gd[0].begin, gd[1].begin, gd[2].begin);
    vec corner2(gd[0].end,   gd[1].end,   gd[2].end);

    parallel_mc par;
//	sz heuristic = m.num_movable_atoms() + 10 * m.get_size().num_degrees_of_freedom();
    par.ls.min_rmsd = 1.0;
    par.ls.num_saved_mins = 20;
    par.ns.min_rmsd = 1.0;
    par.ns.num_saved_mins = 20;
//	par.ls.hunt_cap = vec(10, 10, 10);
//	par.num_tasks = exhaustiveness;//
    par.num_tasks = exhaustiveness;//
    par.num_threads = cpu;//  ͬʱִ е
    par.display_progress = (verbosity > 1);

    const fl slope = 1e6; // FIXME: too large? used to be 100
    if(randomize_only) {
        do_randomization(m, out_name,
                         corner1, corner2, seed, verbosity, log);
    }
    else {
        extern string lig_name;
//        string pdb_name = "../example/"+lig_name+"/"+lig_name+"_protein.pdb";
//        string mol2_name = "../example/"+ lig_name+"/"+lig_name+"_ligand.mol2";
        //将蛋白质的pdb文件和配体的mol2文件供DLIGAND2读取相关分子信息
        string pdb_name = "./"+lig_name+"/"+lig_name+"_protein.pdb";
        string mol2_name = "./"+lig_name+"/"+lig_name+"_ligand.mol2";
        Molecule_DLIGAND2 *mol = new Molecule_DLIGAND2(pdb_name);
        mol->rdmol2(mol2_name);
        //用于线性插值计算
        non_cache nc        (m, gd, &prec,         slope); // if gd has 0 n's, this will not constrain anything
        non_cache nc_widened(m, gd, &prec_widened, slope); // if gd has 0 n's, this will not constrain anything

        non_cache DLIGAND2_nc        ( gd, &prec_DLIGAND2, *mol,        slope); // if gd has 0 n's, this will not constrain anything
//        DLIGAND2_nc.set_outTypeModea(*mol);
//        non_cache nc_widened(m, gd, &prec_widened, slope);
        if(1/*no_cache*/) {
            //接着进入该函数进行搜索
            do_search(m, m_ad4, ref, wt, wt_ad4, prec, prec_ad4, prec_vdw, prec_forcefield, prec_DLIGAND2, nc, DLIGAND2_nc, prec_widened, nc_widened, nc,
                      out_name,
                      corner1, corner2,
                      par, energy_range, num_modes,
                      seed, verbosity, score_only, local_only, log, t, weights);
        }
        else {
//			bool cache_needed = !(score_only || randomize_only || local_only);
//			if(cache_needed) doing(verbosity, "Analyzing the binding site", log);
//			cache c("scoring_function_version001", gd, slope, atom_type::XS);
//			if(cache_needed) c.populate(m, prec, m.get_movable_atom_types(prec.atom_typing_used()));
//			if(cache_needed) done(verbosity, log);
//			do_search(m, ref, wt, prec, prec_forcefield,c, prec, c, nc,
//					  out_name,
//					  //ofs,
//					  corner1, corner2,
//					  par, energy_range, num_modes,
//					  seed, verbosity, score_only, local_only, log, t, weights);
        }
    }
}

struct usage_error : public std::runtime_error {
    usage_error(const std::string& message) : std::runtime_error(message) {}
};

struct options_occurrence {
    bool some;
    bool all;
    options_occurrence() : some(false), all(true) {} // convenience
    options_occurrence& operator+=(const options_occurrence& x) {
        some = some || x.some;
        all  = all  && x.all;
        return *this;
    }
};

options_occurrence get_occurrence(boost::program_options::variables_map& vm, boost::program_options::options_description& d) {
    options_occurrence tmp;
    VINA_FOR_IN(i, d.options())
        if(vm.count((*d.options()[i]).long_name()))
            tmp.some = true;
        else
            tmp.all = false;
    return tmp;
}

void check_occurrence(boost::program_options::variables_map& vm, boost::program_options::options_description& d) {
    VINA_FOR_IN(i, d.options()) {
        const std::string& str = (*d.options()[i]).long_name();
        if(!vm.count(str))
            std::cerr << "Required parameter --" << str << " is missing!\n";
    }
}

model parse_bundle(const std::string& rigid_name, const boost::optional<std::string>& flex_name_opt, const std::vector<std::string>& ligand_names) {
    model tmp = (flex_name_opt) ? parse_receptor_pdbqt(make_path(rigid_name), make_path(flex_name_opt.get()))
                                : parse_receptor_pdbqt(make_path(rigid_name));
    VINA_FOR_IN(i, ligand_names)
        tmp.append(parse_ligand_pdbqt(make_path(ligand_names[i])));
    return tmp;
}

model parse_bundle(const std::vector<std::string>& ligand_names) {
    VINA_CHECK(!ligand_names.empty()); // FIXME check elsewhere
    model tmp = parse_ligand_pdbqt(make_path(ligand_names[0]));
    VINA_RANGE(i, 1, ligand_names.size())
        tmp.append(parse_ligand_pdbqt(make_path(ligand_names[i])));
    return tmp;
}

model parse_bundle(const boost::optional<std::string>& rigid_name_opt, const boost::optional<std::string>& flex_name_opt, const std::vector<std::string>& ligand_names) {
    if(rigid_name_opt)
        return parse_bundle(rigid_name_opt.get(), flex_name_opt, ligand_names);
    else
        return parse_bundle(ligand_names);
}


//int Xtool_Score_Shortcut(){
//    extern Input *input;
//    extern ForceField *ff;
//    extern Protein *protein;
//    extern Ligand *liganda;
//    int i,num_molecule;
//    FILE *fin,*fout;
//    Ligand cofactor;
//
//    num_molecule=Check_Mol2_File(input->ligand_file);
//    if(num_molecule==0)
//    {
//        puts("Error: no valid ligand molecule is given.");
//        exit(1);
//    }
//
//    liganda=new Ligand; if(liganda==NULL) Memory_Allocation_Error();
//    protein=new Protein; if(protein==NULL) Memory_Allocation_Error();
//
//    printf("Now reading parameters from '%s' ...\n",input->parameter_dir);
//
//    ff=new ForceField(input->parameter_dir);
//    if(ff==NULL) Memory_Allocation_Error();
//
//    printf("Now reading the ligand from '%s' ...\n", input->ligand_file);
//    liganda->Read_From_Mol2(input->ligand_file);
//    liganda->Value_Atom();
//    // ligand->Write_Out_Mol2("temp.mol2");
//
//    printf("Now reading the protein from '%s' ... \n",input->receptor_file);
//    protein->Read_From_PDB(input->receptor_file);
//    protein->Value_Atom();
//
//    printf("Now defining the binding pocket ...\n");
//    protein->Define_Pocket(liganda,10.0);
//
//
//    if(num_molecule==1)	// single mol2 file
//    {
//        liganda->Calculate_Binding_Score();
//
//        printf("\n***********************************************\n");
//
//        printf("HPScore -log(Kd) = %-5.2f\n", liganda->pkd1);
//        printf("HMScore -log(Kd) = %-5.2f\n", liganda->pkd2);
//        printf("HSScore -log(Kd) = %-5.2f\n", liganda->pkd3);
//        printf("Predicted average -log(Kd) = %-5.2f\n", liganda->bind_score);
//
//        printf("Predicted binding energy = %-6.2f kcal/mol\n",
//               liganda->bind_score*(-1.364));
//
//        printf("***********************************************\n\n");
//
//
//    }
////    else	// multi-mol2 file
////    {
////        printf("ID: HPSCORE HMSCORE HSSCORE AVE_SCORE BIND_ENERGY NAME\n");
////
////        if((fin=fopen(input->ligand_file,"r"))==NULL) Open_File_Error(input->ligand_file);
////
////        if((fout=fopen(input->log_file,"w"))==NULL) Open_File_Error(input->log_file);
////
////        for(i=1;i<=num_molecule;i++)
////        {
////            liganda->Clear();
////            liganda->Read_From_Mol2(fin);
////            liganda->Value_Atom();
////            liganda->Calculate_Binding_Score();
////
////            printf("Molecule %6d: ", i);
////            printf("%5.2f  ", liganda->pkd1);
////            printf("%5.2f  ", liganda->pkd2);
////            printf("%5.2f  ", liganda->pkd3);
////            printf("%5.2f  ", liganda->bind_score);
////            printf("%6.2f  ", liganda->bind_score*(-1.364));
////            printf("%s\n", liganda->name);
////
////            liganda->Write_Out_Log(fout);
////        }
////
////        fclose(fin); fclose(fout);
////    }
//
//}
#include "restype.h"
#include "dfire.h"
#include "misc.h"

// global variables for DFIRE
int etype = 1;  // 1;167*14   2; 14*14
int Atype = 0;	//0,all; 1,CA; 12,CB; 2,CA+CB; 5,main-chain
double rcut = 15.;
double ALPHA=1.61;
bool bsym=true, bag=false;
int nkind, iskind[matype], nkind2, ismol2[matype];
float edfire[mkind][mkind][mbin], Nobs0[mkind][mbin];
//
void initDFIRE(string fn){
    void rdpolar(string); rdpolar("");
    void initEobs(string); initEobs(fn);

}
int mol2Define(string s1){
/*	char mol2_type[matype_mol2][6] = {"C.2", "C.3", "C.ar", "C.cat", "N.4",
		"N.am", "N.pl3", "O.2", "O.3", "O.co2", "S.3", "P.3", "N.2", "N.ar"
	};*/
    char mol2_type[matype_mol2][6] = {"C.2", "C.3", "C.ar", "C.cat",
                                      "O.2", "O.3", "O.co2", "N.4", "N.am", "N.pl3", "S.3" //, "P.3", "N.2", "N.ar"
    };
    if(s1.substr(0,2) == "S.") s1 = "S.3";
    else if(s1=="P.3" || s1=="Cl" || s1=="Br" || s1=="Met") s1 = "S.3";

    else if(s1 == "F") s1 = "O.co2";
    else if(s1=="C.1") s1 = "C.3";

    else if(s1=="N.3" || s1=="N.1") s1 = "N.2";
    if(s1=="N.2" || s1=="N.ar") s1 = "N.pl3";

    for(int i=0; i<matype_mol2; i++){
        if(s1 == mol2_type[i]) return i;
    }
    if(DEBUG > 0) fprintf(stderr, "not known atom: %s\n", s1.c_str());
    return -1;
}
void rdpolar(string fn){
    void initRestypes_pro(); initRestypes_pro();
    nkind = matype;
    for(int i=0; i<matype; i++) iskind[i] = i;
    if(Atype == 0) return;
//
    for(int i=0; i<matype_pro; i++){
        string &an = atomtypes[i]->getname();
        if(Atype == 1){
            if(an == "CA") continue;
        } else if (Atype == 12){
            if(an == "CB") continue;
        } else if (Atype == 2){
            if(an=="CA" || an=="CB") continue;
        } else if (Atype == 5){
            if(an=="CA" || an=="CB" || an=="N" || an=="O" || an=="C") continue;
        }
        iskind[i] = -1;
    }
}
void rdmol2type(string fn){
    if(fn == "") fn = "./parameter/amino.mol2";
    FILE *fp = openfile(fn, "r");
    char str[121], ss[9][10];
    fgets(str, 120, fp);
    while(fgets(str, 120, fp) != NULL){
        if(strstr(str, "END") == str) break;
        if(strstr(str, "ATOM ") != str) continue;
        str2dat(str+4, 3, ss);
        int ia = atomDefine(ss[0], ss[1]);
        if(ia < 0) continue;
        int imol2 = mol2Define(ss[2]);
        ismol2[ia] = imol2;
    }
    nkind2 = matype_mol2;
    for(int i=0; i<matype_mol2; i++){
        ismol2[matype_pro+i] = i;// + matype_mol2;
    }
//	nkind2 += matype_mol2;
    if(DEBUG > 0){
        fprintf(stderr, "nkind2: %d\n", nkind2);
    }
}
void initEobs(string fn){
//	void initpobs1(string fn); initpobs1("gyr.dat");
    if(fn == ""){
        if(datadir=="" && getenv("DATADIR") != NULL){
            datadir = getenv("DATADIR");
            if(datadir != "") datadir += string("/");
        }
        fn = datadir + "dfire.2";
    }
    FILE *fp=openfile(fn, "r");
    char str[6001], *strs, *stre;
    if(DEBUG) cerr<<fn<<endl;
// read the Nobs libary
    int Nobs[nkind][nkind][mbin], Nobs0[nkind][mbin], dat[nkind];
    nkind = matype_pro;
    for(int i=0; i<mbin; i++){
        fgets(str,120,fp);
        if(strchr(str, ':') == NULL) die("Error lib: %s\n", str);
        for(int j=0; j<nkind; j++){
            fgets(str,6000,fp); int it = str2dat(str, dat);
            if(it != nkind) die("wrong lib: %d -- %d\n%s\n", it, nkind, str);
            for(int k=0; k<nkind; k++) Nobs[j][k][i] = dat[k];
        }
    }
    if(fgets(str,120,fp) != NULL){		// Nobs0
        for(int i=0; i<mbin; i++){
            fgets(str,6000,fp); int it = str2dat(str, dat);
            if(it != nkind) die("wrong Nobs0: %s", str);
            for(int k=0; k<nkind; k++) Nobs0[k][i] = dat[k];
        }
    }
//
    nkind = matype;
/* symtry
if(bsym){
	for(int i=0; i<mbin; i++){
		for(int j=0; j<nkind; j++)
		for(int k=j+1; k<nkind; k++){
			Nobs[k][j][i] = Nobs[j][k][i] + Nobs[k][j][i];
			Nobs[j][k][i] = Nobs[j][k][i];
		}
	}
}*/
// mol2
    rdmol2type("");
    int Nobs2[nkind2][nkind2][mbin]; bzero(Nobs2, sizeof(Nobs2));
    int Nobs3[nkind][nkind2][mbin]; bzero(Nobs3, sizeof(Nobs3));
    double pvol[3][mbin]; bzero(pvol, sizeof(pvol));
    for(int b=0; b<mbin; b++){
        for(int i=0; i<matype; i++)
            for(int j=i; j<matype; j++){
                int ik = ismol2[i], jk = ismol2[j];
                int dt = Nobs[i][j][b] + Nobs[j][i][b];
                if(j == i) dt /= 2;
                Nobs2[jk][ik][b] += dt; Nobs3[i][jk][b] += dt;
                Nobs2[ik][jk][b] += dt; Nobs3[j][ik][b] += dt;
            }
    }
/*	for(int i=0; i<nkind; i++)
	for(int j=0; j<nkind; j++){
		int ik = ismol2[i], jk = ismol2[j];
		for(int b=0; b<mbin; b++) Nobs[i][j][b] = Nobs2[ik][jk][b];
	}*/
// pobs0
    double pobs[mbin], pobs0[mbin];
    for(int m=0; m<nbin; m++) pobs0[m] = pow(bin2r(m), ALPHA);
    double ds = sum(nbin, pobs0);
    for(int m=0; m<nbin; m++) pobs0[m] /= ds;
// calculte Eobs
    double dat1[mbin];
    for(int i=0; i<matype_pro; i++)
        for(int j=matype_pro; j<matype; j++){
            int ik=ismol2[i], jk=ismol2[j]; int n0 = 75;
            if(etype == 2){
                ds = sum(nbin, Nobs2[ik][jk]) + n0;
                for(int m=0; m<nbin; m++) pobs[m] = (Nobs2[ik][jk][m]+n0*pobs0[m]) / ds;
            }else{
//		ds = sum(nbin, Nobs[i][j]) + n0;
//		for(int m=0; m<nbin; m++) pobs[m] = (Nobs[i][j][m]+n0*pobs[m]) / ds;
                ds = sum(nbin, Nobs3[i][jk]) + n0;
                for(int m=0; m<nbin; m++) pobs[m] = (Nobs3[i][jk][m]+n0*pobs0[m]) / ds;
            }
//
            for(int m=0; m<nbin; m++){
                double e, r = bin2r(m), rc=bin2r(nbin-1);
                e = -ert*log(pobs[m] / pobs0[m]);
//			e = -ert*log(pobs[m] / pow(r/rc, ALPHA) /  pobs[nbin-1]);
//			e = min(e, epen1*2);
                dat1[m] = e;
            }
            for(int m=0; m<nbin; m++) {
                edfire[j][i][m] = edfire[i][j][m] = min(epen1, dat1[m] - dat1[nbin-1]);
            }
        }
    if(DEBUG > 2)
    {
        for(int m=0; m<nbin; m++){
            printf("%d %f %f\n", m, edfire[2][167+1][m], edfire[20][167+7][m]);
        }
        exit(0);
    }
    return;
}
//#include <boost/python.hpp>

Input *input = NULL;
ForceField *ff = NULL;
Protein *protein = NULL;
Ligand *liganda = NULL;
//OBMolecule obl;
//OBProtein obp;
std::string lig_name;
char timestring[256];
double deviceTime;
std::vector<output_container_MODEA> out_tmp;
output_type_MODEA std_tmp;
model std_model;
output_container_MODEA out_300;
output_container_MODEA out_500;
output_container_MODEA out_800;
output_container_MODEA out_1000;
float min_vina = max_fl;
float min_DLIGAND2 = max_fl;
float min_vdw = max_fl;
float min_e = max_fl;
output_container_MODEA out_process_pf;
//boost::python::object find_min_norm_element;
//
//boost::python::object sys;
//boost::python::object pathss;
//boost::python::object module;
//boost::python::object solver_class;
//boost::python::object py_args;
//boost::python::object result;
//python_op po;
//boost::python::object getfunc(){
//    Py_Initialize();
//    sys = boost::python::import("sys");
//    pathss = sys.attr("path");
//    pathss.attr("append")("/home/caokun/Cplusprojects/koto_MODEA_gradient1_exp/LSHADE_Adam_final_nsga2");
//    module = boost::python::import("min_norm_solvers_numpy");
//    solver_class = module.attr("MinNormSolverNumpy");
//    find_min_norm_element = solver_class.attr("find_min_norm_element");
//}

// 函数用于从CSV文件中读取数据
bool read_csv(const std::string &filename, std::vector<double> &iterations, std::vector<double> &fitness_values) {
    std::ifstream in_file(filename);
    std::string line;
    if (!in_file.is_open()) {
        std::cerr << "Unable to open the file " << filename << std::endl;
        return false;
    }

    // 跳过标题行
    std::getline(in_file, line);

    // 读取文件中的每一行
    while (std::getline(in_file, line)) {
        std::stringstream ss(line);
        std::string value;
        double iteration, fitness_value;

        // 读取迭代次数
        std::getline(ss, value, ',');
        iteration = std::stod(value);

        // 读取适应度值
        std::getline(ss, value);
        fitness_value = std::stod(value);

        iterations.push_back(iteration);
        fitness_values.push_back(fitness_value);
    }

    in_file.close();
    return true;
}

int main(int argc, char* argv[]) {
//    // 假设我们有两个CSV文件
//    std::vector<double> iterations1, fitness_values1;
//    std::vector<double> iterations2, fitness_values2;
//    std::vector<double> iterations3, fitness_values3;
//    std::vector<double> iterations4, fitness_values4;
//    std::vector<double> iterations5, fitness_values5;
//    std::vector<double> iterations6, fitness_values6;
//
//    // 读取第一个CSV文件
//    if (!read_csv("../csv/results_1g2k_xscore.csv", iterations1, fitness_values1)) {
//        return 1; // 如果文件无法打开，则返回错误
//    }
//
//    // 读取第二个CSV文件
//    if (!read_csv("../csv/results_1qf1_xscore.csv", iterations2, fitness_values2)) {
//        return 1; // 如果文件无法打开，则返回错误
//    }
//
//    if (!read_csv("../csv/results_2zda_xscore.csv", iterations3, fitness_values3)) {
//        return 1; // 如果文件无法打开，则返回错误
//    }
//    if (!read_csv("../csv/results_3acw_xscore.csv", iterations4, fitness_values4)) {
//        return 1; // 如果文件无法打开，则返回错误
//    }
////    if (!read_csv("../csv/results_3bv9_total.csv", iterations5, fitness_values5)) {
////        return 1; // 如果文件无法打开，则返回错误
////    }
//    if (!read_csv("../csv/results_3e92_xscore.csv", iterations6, fitness_values6)) {
//        return 1; // 如果文件无法打开，则返回错误
//    }
//    std::map<std::string, std::string> keywords;
//    keywords["label"] = "1";
//    keywords["label1"] = "2";
//    // 绘制第一次运行结果
//    matplotlibcpp::plot(iterations1, fitness_values1, "r-");
//
//    // 绘制第二次运行结果
//    matplotlibcpp::plot(iterations2, fitness_values2, "b-");
//
//    matplotlibcpp::plot(iterations3, fitness_values3, "g-");
//
//    matplotlibcpp::plot(iterations4, fitness_values4, "c-");
//
////    matplotlibcpp::plot(iterations5, fitness_values5, "m-");
//
//    matplotlibcpp::plot(iterations6, fitness_values6, "y-");
//
//    // 添加图例
//    matplotlibcpp::legend();
//
//    // 设置标题和轴标签
//    matplotlibcpp::ylim(-15.0, 0.0);
//    matplotlibcpp::title("xscore Convergence Fig");
//    matplotlibcpp::xlabel("iterations");
//    matplotlibcpp::ylabel("xscore");
//
//    // 显示图表
//    matplotlibcpp::show();
//    matplotlibcpp::save("../csv/convergence_fig_xscore.png"); // 保存图像




    // 示例数据
//    std::vector<double> iterations = {1, 2, 3, 4, 5}; // 用实际的迭代次数替换
//    std::vector<double> best_fitness = {10, 5, 3, 2.5, 2}; // 用实际的最优适应度分数替换
//
//    // 绘图
//    matplotlibcpp::plot(iterations, best_fitness, "r-"); // 'r' 表示红色, '-' 表示实线
//    matplotlibcpp::xlabel("aaa");
//    matplotlibcpp::ylabel("bbb");
//    matplotlibcpp::title("pareto");
//    matplotlibcpp::save("../coreset_pdbqt/3aru/convergence_graph.png"); // 保存图像
//    matplotlibcpp::show();

    using namespace boost::program_options;

    string dfire_name = "./parameter/dfire.2";
    void initDFIRE(string fn); initDFIRE(dfire_name);
//    Py_Initialize();
//    extern boost::python::object find_min_norm_element;
//
//    extern boost::python::object sys;
//    extern boost::python::object pathss;
//    extern boost::python::object module;
//    extern boost::python::object solver_class;
//    extern python_op po;
//    sys = boost::python::import("sys");
//    pathss = sys.attr("path");
//    pathss.attr("append")("/home/caokun/Cplusprojects/koto_MODEA_gradient1_exp/LSHADE_Adam_final_nsga2");
//    module = boost::python::import("min_norm_solvers_numpy");
//    solver_class = module.attr("MinNormSolverNumpy");
//    find_min_norm_element = solver_class.attr("find_min_norm_element");
//    po = python_op();
    const std::string version_string = "AutoDock Vina 1.1.2 (May 11, 2011)";
    const std::string error_message = "\n\n\
Please contact the author, Dr. Oleg Trott <ot14@columbia.edu>, so\n\
that this problem can be resolved. The reproducibility of the\n\
error may be vital, so please remember to include the following in\n\
your problem report:\n\
* the EXACT error message,\n\
* your version of the program,\n\
* the type of computer system you are running it on,\n\
* all command line options,\n\
* configuration file (if used),\n\
* ligand file as PDBQT,\n\
* receptor file as PDBQT,\n\
* flexible side chains file as PDBQT (if used),\n\
* output file as PDBQT (if any),\n\
* input (if possible),\n\
* random seed the program used (this is printed when the program starts).\n\
\n\
Thank you!\n";
/*
	const std::string cite_message = "\
#################################################################\n\
# If you used AutoDock Vina in your work, please cite:          #\n\
#                                                               #\n\
# O. Trott, A. J. Olson,                                        #\n\
# AutoDock Vina: improving the speed and accuracy of docking    #\n\
# with a new scoring function, efficient optimization and       #\n\
# multithreading, Journal of Computational Chemistry 31 (2010)  #\n\
# 455-461                                                       #\n\
#                                                               #\n\
# DOI 10.1002/jcc.21334                                         #\n\
#                                                               #\n\
# Please see http://vina.scripps.edu for more information.      #\n\
#################################################################\n";
*/
    try {
        std::string rigid_name, ligand_name, flex_name, config_name, out_name, log_name ,target_name;
        fl center_x, center_y, center_z, size_x, size_y, size_z;
        int cpu = 0, seed, exhaustiveness, verbosity = 2, num_modes = 9;
        fl energy_range = 2.0;

        // -0.035579, -0.005156, 0.840245, -0.035069, -0.587439, 0.05846
        fl weight_gauss1      = -0.035579;
        fl weight_gauss2      = -0.005156;
        fl weight_repulsion   =  0.840245;
        fl weight_hydrophobic = -0.035069;
        fl weight_hydrogen    = -0.587439;
        fl weight_rot         =  0.05846;
        //ad4
        fl weight_ad4_vdw     = 0.1662;
        fl weight_ad4_hb      = 0.1209;
        fl weight_ad4_elec    = 0.1406;
        fl weight_ad4_dsolv   = 0.1322;
        fl weight_glue        = 50;
        fl weight_ad4_rot     = 0.2983;

        bool score_only = false, local_only = false, randomize_only = false, help = false, help_advanced = false, version = false; // FIXME

        positional_options_description positional; // remains empty

        options_description inputs("Input");
        inputs.add_options()
                ("receptor", value<std::string>(&rigid_name), "rigid part of the receptor (PDBQT)")
                ("flex", value<std::string>(&flex_name), "flexible side chains, if any (PDBQT)")
                ("ligand", value<std::string>(&ligand_name), "ligand (PDBQT)")
            //("target", value<std::string>(&target_name), "target (PDBQT)")
                ;
        //options_description search_area("Search area (required, except with --score_only)");
        options_description search_area("Search space (required)");
        search_area.add_options()
                ("center_x", value<fl>(&center_x), "X coordinate of the center")
                ("center_y", value<fl>(&center_y), "Y coordinate of the center")
                ("center_z", value<fl>(&center_z), "Z coordinate of the center")
                ("size_x", value<fl>(&size_x), "size in the X dimension (Angstroms)")
                ("size_y", value<fl>(&size_y), "size in the Y dimension (Angstroms)")
                ("size_z", value<fl>(&size_z), "size in the Z dimension (Angstroms)")
                ;
        //options_description outputs("Output prefixes (optional - by default, input names are stripped of .pdbqt\nare used as prefixes. _001.pdbqt, _002.pdbqt, etc. are appended to the prefixes to produce the output names");
        options_description outputs("Output (optional)");
        outputs.add_options()
                ("out", value<std::string>(&out_name), "output models (PDBQT), the default is chosen based on the ligand file name")
                ("log", value<std::string>(&log_name), "optionally, write log file")
                ;
        options_description advanced("Advanced options (see the manual)");
        advanced.add_options()
                ("score_only",     bool_switch(&score_only),     "score only - search space can be omitted")
                ("local_only",     bool_switch(&local_only),     "do local search only")
                ("randomize_only", bool_switch(&randomize_only), "randomize input, attempting to avoid clashes")
                ("weight_gauss1", value<fl>(&weight_gauss1)->default_value(weight_gauss1),                "gauss_1 weight")
                ("weight_gauss2", value<fl>(&weight_gauss2)->default_value(weight_gauss2),                "gauss_2 weight")
                ("weight_repulsion", value<fl>(&weight_repulsion)->default_value(weight_repulsion),       "repulsion weight")
                ("weight_hydrophobic", value<fl>(&weight_hydrophobic)->default_value(weight_hydrophobic), "hydrophobic weight")
                ("weight_hydrogen", value<fl>(&weight_hydrogen)->default_value(weight_hydrogen),          "Hydrogen bond weight")
                ("weight_rot", value<fl>(&weight_rot)->default_value(weight_rot),                         "N_rot weight")
                ;
        options_description misc("Misc (optional)");
        misc.add_options()
                ("target", value<std::string>(&target_name), "target (PDBQT)")
                ("cpu", value<int>(&cpu), "the number of CPUs to use (the default is to try to detect the number of CPUs or, failing that, use 1)")
                ("seed", value<int>(&seed), "explicit random seed")
                ("exhaustiveness", value<int>(&exhaustiveness)->default_value(8), "exhaustiveness of the global search (roughly proportional to time): 1+")
                ("num_modes", value<int>(&num_modes)->default_value(9), "maximum number of binding modes to generate")
                ("energy_range", value<fl>(&energy_range)->default_value(3.0), "maximum energy difference between the best binding mode and the worst one displayed (kcal/mol)")
                ;
        options_description config("Configuration file (optional)");
        config.add_options()
                ("config", value<std::string>(&config_name), "the above options can be put here")
                ;
        options_description info("Information (optional)");
        info.add_options()
                ("help",          bool_switch(&help), "display usage summary")
                ("help_advanced", bool_switch(&help_advanced), "display usage summary with advanced options")
                ("version",       bool_switch(&version), "display program version")
                ;
        options_description desc, desc_config, desc_simple;
        desc       .add(inputs).add(search_area).add(outputs).add(advanced).add(misc).add(config).add(info);
        desc_config.add(inputs).add(search_area).add(outputs).add(advanced).add(misc);
        desc_simple.add(inputs).add(search_area).add(outputs).add(misc).add(config).add(info);

        variables_map vm;
        try {
            //store(parse_command_line(argc, argv, desc, command_line_style::default_style ^ command_line_style::allow_guessing), vm);
            store(command_line_parser(argc, argv)
                          .options(desc)
                          .style(command_line_style::default_style ^ command_line_style::allow_guessing)
                          .positional(positional)
                          .run(),
                  vm);
            notify(vm);
        }
        catch(boost::program_options::error& e) {
            std::cerr << "Command line parse error: " << e.what() << '\n' << "\nCorrect usage:\n" << desc_simple << '\n';
            return 1;
        }
        if(vm.count("config")) {
            try {
                path name = make_path(config_name);
                ifile config_stream(name);
                store(parse_config_file(config_stream, desc_config), vm);
                notify(vm);
            }
            catch(boost::program_options::error& e) {
                std::cerr << "Configuration file parse error: " << e.what() << '\n' << "\nCorrect usage:\n" << desc_simple << '\n';
                return 1;
            }
        }
        if(help) {
            std::cout << desc_simple << '\n';
            return 0;
        }
        if(help_advanced) {
            std::cout << desc << '\n';
            return 0;
        }
        if(version) {
            //std::cout << version_string << '\n';
            return 0;
        }

        bool search_box_needed = !score_only; // randomize_only and local_only still need the search space
        bool output_produced   = !score_only;
        bool receptor_needed   = !randomize_only;

        if(receptor_needed) {
            if(vm.count("receptor") <= 0) {
                std::cerr << "Missing receptor.\n" << "\nCorrect usage:\n" << desc_simple << '\n';
                return 1;
            }
        }
        if(vm.count("ligand") <= 0) {
            std::cerr << "Missing ligand.\n" << "\nCorrect usage:\n" << desc_simple << '\n';
            return 1;
        }
        if(cpu < 1)
            cpu = 1;
        if(vm.count("seed") == 0)
            seed = auto_seed();
        if(exhaustiveness < 1)
            throw usage_error("exhaustiveness must be 1 or greater");
        if(num_modes < 1)
            throw usage_error("num_modes must be 1 or greater");
        sz max_modes_sz = static_cast<sz>(num_modes);

        boost::optional<std::string> rigid_name_opt;
        if(vm.count("receptor"))
            rigid_name_opt = rigid_name;

        boost::optional<std::string> flex_name_opt;
        if(vm.count("flex"))
            flex_name_opt = flex_name;

        if(vm.count("flex") && !vm.count("receptor"))
            throw usage_error("Flexible side chains are not allowed without the rest of the receptor"); // that's the only way parsing works, actually

        tee log;
        if(vm.count("log") > 0)
            log.init(log_name);

        if(search_box_needed) {
            options_occurrence oo = get_occurrence(vm, search_area);
            if(!oo.all) {
                check_occurrence(vm, search_area);
                std::cerr << "\nCorrect usage:\n" << desc_simple << std::endl;
                return 1;
            }
            if(size_x <= 0 || size_y <= 0 || size_z <= 0)
                throw usage_error("Search space dimensions should be positive");
        }

        //log << cite_message << '\n';

        if(search_box_needed && size_x * size_y * size_z > 27e3) {
            log << "WARNING: The search space volume > 27000 Angstrom^3 (See FAQ)\n";
        }

        if(output_produced) { // FIXME
            if(!vm.count("out")) {
                out_name = default_output(ligand_name);
                log << "Output will be " << out_name << '\n';
            }
        }

        grid_dims gd; // n's = 0 via default c'tor

        flv weights;
        weights.push_back(weight_gauss1);
        weights.push_back(weight_gauss2);
        weights.push_back(weight_repulsion);
        weights.push_back(weight_hydrophobic);
        weights.push_back(weight_hydrogen);
        weights.push_back(5 * weight_rot / 0.1 - 1); // linearly maps onto a different range, internally. see everything.cpp

        flv weights_ad4;
        weights_ad4.push_back(weight_ad4_vdw);
        weights_ad4.push_back(weight_ad4_hb);
        weights_ad4.push_back(weight_ad4_elec);
        weights_ad4.push_back(weight_ad4_dsolv);
        weights_ad4.push_back(weight_glue);
        weights_ad4.push_back(weight_ad4_rot);

        flv weights_vdw;
        weights_vdw.push_back(weight_ad4_vdw);
        //xscore




        extern Input *input;
        extern string lig_name;
        // 找到第一个斜杠的位置
        size_t first_slash = ligand_name.find('/');
        // 找到第二个斜杠的位置
        size_t second_slash = ligand_name.find('/', first_slash + 1);

        // 提取两个斜杠之间的字符串
        lig_name= ligand_name.substr(first_slash + 1, second_slash - first_slash - 1);
//        printf("\nX-Score starts to run ... %s\n\n", Get_Time());
//
//        input=new Input;
//        string receptor_pdb = "./"+lig_name+"/"+lig_name+"_receptor.pdb";
//        string ligand_li_mol2 = "./"+lig_name+"/"+lig_name+"_ligand_li.mol2";
//        strcpy(input->receptor_file,receptor_pdb.c_str());
//        strcpy(input->ligand_file,ligand_li_mol2.c_str());
//        strcpy(input->parameter_dir,"./parameter/");
//
//
////        strcpy(input->log_file,"xscore.log");
//        strcpy(input->apply_hpscore,"YES");
//        strcpy(input->apply_hmscore,"YES");
//        strcpy(input->apply_hsscore,"YES");
//        input->num_method=3;
//        strcpy(input->apply_pmfscore,"NO");
//
//
//
//        extern ForceField *ff;
//        extern Protein *protein;
//        extern Ligand *liganda;
//        int i,num_molecule;
//        FILE *fin,*fout;
//        Ligand cofactor;
//
//        num_molecule=Check_Mol2_File(input->ligand_file);
//        if(num_molecule==0)
//        {
//            puts("Error: no valid ligand molecule is given.");
//            exit(1);
//        }
//
//        liganda=new Ligand; if(liganda==NULL) Memory_Allocation_Error();
//        protein=new Protein; if(protein==NULL) Memory_Allocation_Error();
////        if(outputTypeMODEA==NULL) Memory_Allocation_Error();
//
//        printf("Now reading parameters from '%s' ...\n",input->parameter_dir);
//
//        ff=new ForceField(input->parameter_dir);
//        if(ff==NULL) Memory_Allocation_Error();
//
//        printf("Now reading the ligand from '%s' ...\n", input->ligand_file);
//        liganda->Read_From_Mol2(input->ligand_file);
//        liganda->Value_Atom();
//
////        outputTypeMODEA->Read_From_Mol2(input->ligand_file);
////        outputTypeMODEA->Value_Atom();
//
//        // ligand->Write_Out_Mol2("temp.mol2");
//
//        printf("Now reading the protein from '%s' ... \n",input->receptor_file);
//        protein->Read_From_PDB(input->receptor_file);
//        protein->Value_Atom();
//
//        printf("Now defining the binding pocket ...\n");
//        protein->Define_Pocket(liganda,10.0);


//        // Prepare the ligand.
//        OBMolecule moleculeRef;
//        OpenBabel::OBConversion obConversion;
//        OpenBabel::OBFormat *format = obConversion.FormatFromExt("./"+lig_name+"/"+lig_name+"_ligand.mol2");
//        obConversion.SetInFormat(format);
//        obConversion.ReadFile(&obl.obMol, "./"+lig_name+"/"+lig_name+"_ligand.mol2");
//        LigandTypeSMOG2016(obl);
//
//
//        // Create an OpenBabel object for the protein and open it.
//        OpenBabel::OBConversion obConversion1;
//        OpenBabel::OBFormat *format1 = obConversion1.FormatFromExt("./"+lig_name+"/"+lig_name+"_protein.pdb");
//        obConversion1.SetInFormat(format1);
//        obConversion1.ReadFile(&obp.obMol, "./"+lig_name+"/"+lig_name+"_protein.pdb");
//        // This part assign the potential atomtypes for the protein atoms.
//        for (unsigned int i=1; i<=obp.obMol.NumAtoms(); i++) {
//            OpenBabel::OBAtom *atomProt;
//            atomProt = obp.obMol.GetAtom(i);
//            OpenBabel::OBResidue *obResidue;
//            obResidue = atomProt->GetResidue();
//            string residue=obResidue->GetName();
//            string atomID=obResidue->GetAtomID(atomProt);
//            residue.erase(remove(residue.begin(), residue.end(), ' '), residue.end());
//            atomID.erase(remove(atomID.begin(), atomID.end(), ' '), atomID.end());
//            // Get the parameters.
//            Atom tempAtom;
//            string atomType="";
//            int atomTypeNumber=0;
//            float LJe=0.0;
//            float LJr=0.0;
//            ProteinTypeSMOG2016(residue, atomID, atomTypeNumber, LJr, LJe);
//            tempAtom.coordinates[0] = atomProt->GetX();
//            tempAtom.coordinates[1] = atomProt->GetY();
//            tempAtom.coordinates[2] = atomProt->GetZ();
//            tempAtom.Type = atomType;
//            tempAtom.TypeNumber = atomTypeNumber;
//            tempAtom.LJr = LJr;
//            tempAtom.LJe = LJe;
//            tempAtom.index = i;
//            obp.obatom.push_back(tempAtom);
//        }
//        // Get the PDB ID.
//        string PDBID="./"+lig_name+"/"+lig_name+"_protein.pdb";
//        unsigned int positionToErase = PDBID.find_last_of("_");
//        PDBID.erase(positionToErase);


























//===================
//        extern Input *input;
//        extern string lig_name;
//        lig_name = ligand_name.substr(17, 4);
//        printf("\nX-Score starts to run ... %s\n\n", Get_Time());
//
//        input=new Input;
//        strcpy(input->receptor_file,("../coreset_pdbqt/"+lig_name+"/"+lig_name+"_receptor.pdb").c_str());
//        strcpy(input->ligand_file,("../coreset_pdbqt/"+lig_name+"/"+lig_name+"_ligand_li.mol2").c_str());
//        strcpy(input->parameter_dir,"../parameter/");
//
//
////        strcpy(input->log_file,"xscore.log");
//        strcpy(input->apply_hpscore,"YES");
//        strcpy(input->apply_hmscore,"YES");
//        strcpy(input->apply_hsscore,"YES");
//        input->num_method=3;
//        strcpy(input->apply_pmfscore,"NO");
//
//
//
//        extern ForceField *ff;
//        extern Protein *protein;
//        extern Ligand *liganda;
//        int i,num_molecule;
//        FILE *fin,*fout;
//        Ligand cofactor;
//
//        num_molecule=Check_Mol2_File(input->ligand_file);
//        if(num_molecule==0)
//        {
//            puts("Error: no valid ligand molecule is given.");
//            exit(1);
//        }
//
//        liganda=new Ligand; if(liganda==NULL) Memory_Allocation_Error();
//        protein=new Protein; if(protein==NULL) Memory_Allocation_Error();
////        if(outputTypeNSGA2==NULL) Memory_Allocation_Error();
//
//        printf("Now reading parameters from '%s' ...\n",input->parameter_dir);
//
//        ff=new ForceField(input->parameter_dir);
//        if(ff==NULL) Memory_Allocation_Error();
//
//        printf("Now reading the ligand from '%s' ...\n", input->ligand_file);
//        liganda->Read_From_Mol2(input->ligand_file);
//        liganda->Value_Atom();
//
////        outputTypeNSGA2->Read_From_Mol2(input->ligand_file);
////        outputTypeNSGA2->Value_Atom();
//
//        // ligand->Write_Out_Mol2("temp.mol2");
//
//        printf("Now reading the protein from '%s' ... \n",input->receptor_file);
//        protein->Read_From_PDB(input->receptor_file);
//        protein->Value_Atom();
//
//        printf("Now defining the binding pocket ...\n");
//        protein->Define_Pocket(liganda,10.0);


//        // Prepare the ligand.
//        OBMolecule moleculeRef;
//        OpenBabel::OBConversion obConversion;
//        OpenBabel::OBFormat *format = obConversion.FormatFromExt("../coreset_pdbqt/3aru/3aru_ligand.mol2");
//        obConversion.SetInFormat(format);
//        obConversion.ReadFile(&obl.obMol, "../coreset_pdbqt/3aru/3aru_ligand.mol2");
//        LigandTypeSMOG2016(obl);
//
//
//        // Create an OpenBabel object for the protein and open it.
//        OpenBabel::OBConversion obConversion1;
//        OpenBabel::OBFormat *format1 = obConversion1.FormatFromExt("../coreset_pdbqt/3aru/3aru_protein.pdb");
//        obConversion1.SetInFormat(format1);
//        obConversion1.ReadFile(&obp.obMol, "../coreset_pdbqt/3aru/3aru_protein.pdb");
//        // This part assign the potential atomtypes for the protein atoms.
//        for (unsigned int i=1; i<=obp.obMol.NumAtoms(); i++) {
//            OpenBabel::OBAtom *atomProt;
//            atomProt = obp.obMol.GetAtom(i);
//            OpenBabel::OBResidue *obResidue;
//            obResidue = atomProt->GetResidue();
//            string residue=obResidue->GetName();
//            string atomID=obResidue->GetAtomID(atomProt);
//            residue.erase(remove(residue.begin(), residue.end(), ' '), residue.end());
//            atomID.erase(remove(atomID.begin(), atomID.end(), ' '), atomID.end());
//            // Get the parameters.
//            Atom tempAtom;
//            string atomType="";
//            int atomTypeNumber=0;
//            float LJe=0.0;
//            float LJr=0.0;
//            ProteinTypeSMOG2016(residue, atomID, atomTypeNumber, LJr, LJe);
//            tempAtom.coordinates[0] = atomProt->GetX();
//            tempAtom.coordinates[1] = atomProt->GetY();
//            tempAtom.coordinates[2] = atomProt->GetZ();
//            tempAtom.Type = atomType;
//            tempAtom.TypeNumber = atomTypeNumber;
//            tempAtom.LJr = LJr;
//            tempAtom.LJe = LJe;
//            tempAtom.index = i;
//            obp.obatom.push_back(tempAtom);
//        }
//        // Get the PDB ID.
//        string PDBID="../coreset_pdbqt/3aru/3aru_protein.pdb";
//        unsigned int positionToErase = PDBID.find_last_of("_");
//        PDBID.erase(positionToErase);

        if(search_box_needed) {
            const fl granularity = 0.375;
            vec span(size_x,   size_y,   size_z);
            vec center(center_x, center_y, center_z);
            VINA_FOR_IN(i, gd) {
                gd[i].n = sz(std::ceil(span[i] / granularity));
                fl real_span = granularity * gd[i].n;
                gd[i].begin = center[i] - real_span/2;
                gd[i].end = gd[i].begin + real_span;
            }
        }
        if(vm.count("cpu") == 0) {
            unsigned num_cpus = boost::thread::hardware_concurrency();
            if(verbosity > 1) {
                if(num_cpus > 0)
                    log << "Detected " << num_cpus << " CPU" << ((num_cpus > 1) ? "s" : "") << '\n';
                else
                    log << "Could not detect the number of CPUs, using 1\n";
            }
            if(num_cpus > 0)
                cpu = num_cpus;
            else
                cpu = 1;
        }
        if(cpu < 1)
            cpu = 1;
        if(verbosity > 1 && exhaustiveness < cpu)
            log << "WARNING: at low exhaustiveness, it may be impossible to utilize all CPUs\n";

        doing(verbosity, "Reading input", log);

        model m       = parse_bundle(rigid_name_opt, flex_name_opt, std::vector<std::string>(1, ligand_name));
        model m_ad4       = parse_bundle(rigid_name_opt, flex_name_opt, std::vector<std::string>(1, ligand_name));
        extern model std_model;
        std_model = parse_bundle(rigid_name_opt, flex_name_opt, std::vector<std::string>(1, ligand_name));
        m_ad4.set_m_atom_typing_used();

        //std::string result_name = default_output_data(ligand_name);
        //std::ofstream ofs(result_name.c_str());

        boost::optional<model> ref;
        if (target_name.empty() != true)
            ref = parse_bundle(std::vector<std::string>(1, target_name));

        done(verbosity, log);

        main_procedure(m, m_ad4, ref,
                       out_name,
                       score_only, local_only, randomize_only, false, // no_cache == false
                       gd, exhaustiveness,
                       weights,weights_ad4,weights_vdw,
                       cpu, seed, verbosity, max_modes_sz, energy_range, log);
    }
    catch(file_error& e) {
        std::cerr << "\n\nError: could not open \"" << e.name.filename() << "\" for " << (e.in ? "reading" : "writing") << ".\n";
        return 1;
    }
    catch(boost::filesystem::filesystem_error& e) {
        std::cerr << "\n\nFile system error: " << e.what() << '\n';
        return 1;
    }
    catch(usage_error& e) {
        std::cerr << "\n\nUsage error: " << e.what() << ".\n";
        return 1;
    }
    catch(parse_error& e) {
        std::cerr << "\n\nParse error on line " << e.line << " in file \"" << e.file.filename() << "\": " << e.reason << '\n';
        return 1;
    }
    catch(std::bad_alloc&) {
        std::cerr << "\n\nError: insufficient memory!\n";
        return 1;
    }

        // Errors that shouldn't happen:

    catch(std::exception& e) {
        std::cerr << "\n\nAn error occurred: " << e.what() << ". " << error_message;
        return 1;
    }
    catch(internal_error& e) {
        std::cerr << "\n\nAn internal error occurred in " << e.file << "(" << e.line << "). " << error_message;
        return 1;
    }
    catch(...) {
        std::cerr << "\n\nAn unknown error occurred. " << error_message;
        return 1;
    }
}
