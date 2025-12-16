/*
parallel_mc
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
// add
#include <boost/thread/condition_variable.hpp>
#include <boost/throw_exception.hpp>


#include "parallel.h"
#include "parallel_mc.h"
#include "coords.h"
#include "parallel_progress.h"
//#include "python3.8/Python.h"
#include "matplotlibcpp.h"
bool sort_comparator(const output_type_MODEA& x,const output_type_MODEA& y);
struct parallel_mc_task {
    model m;
    model m_ad4;
    output_container_MODEA out;
    rng generator;
    parallel_mc_task(const model& m_, const model& m_ad4_, int seed) : m(m_), m_ad4(m_ad4_),generator(static_cast<rng::result_type>(seed)) {}
};

typedef boost::ptr_vector<parallel_mc_task> parallel_mc_task_container;

struct parallel_mc_aux {
    const lshade* ls;
    const MODEA* modea;
    const precalculate* p;
    const precalculate* prec_ad4;
    const precalculate* prec_vdw;
    const precalculate* prec_forcefield;
    const precalculate* prec_DLIGAND2;
    const igrid* ig;
    const igrid* DLIGAND2_grid;
    const precalculate* p_widened;
    const igrid* ig_widened;
    const vec* corner1;
    const vec* corner2;
    parallel_progress* pg;
    const scoring_function* sf;
    const scoring_function* sf_ad4;
    const boost::optional<model>& ref;
    parallel_mc_aux(const MODEA* modea_, const precalculate* p_, const precalculate* prec_ad4_, const precalculate* prec_vdw_, const precalculate* prec_forcefield_, const precalculate* prec_DLIGAND2_, const igrid* ig_, const igrid* DLIGAND2_grid_, const precalculate* p_widened_, const igrid* ig_widened_, const vec* corner1_, const vec* corner2_, parallel_progress* pg_,const scoring_function* sf_,const scoring_function* sf_ad4_,const boost::optional<model>& ref_)
            : modea(modea_), p(p_), prec_ad4(prec_ad4_), prec_vdw(prec_vdw_), prec_forcefield(prec_forcefield_), prec_DLIGAND2(prec_DLIGAND2_), ig(ig_), DLIGAND2_grid(DLIGAND2_grid_), p_widened(p_widened_), ig_widened(ig_widened_), corner1(corner1_), corner2(corner2_), pg(pg_),sf(sf_),sf_ad4(sf_ad4_),ref(ref_) {}
    void operator()(parallel_mc_task& t) const {
        (*modea)(t.m, t.m_ad4, t.out, *p, *prec_ad4, *prec_vdw, *prec_forcefield, *prec_DLIGAND2, *ig, *DLIGAND2_grid, *p_widened, *ig_widened, *corner1, *corner2, pg, t.generator,*sf,*sf_ad4,ref);
    }
};

void Stringsplit_(std::string str,const  char split,std::vector<std::string>& list)
{
    std::istringstream iss(str);	// ??????
    std::string token;			// ?????????
    while (getline(iss, token, split))	// ??split??????
    {
        if(token != ""){
            list.push_back(token);

        }
    }
}
//void draw_pareto_(const output_container_MODEA& out_cont,int i){
//    // ?????????��???????????
//    std::vector<std::vector<double>> objectives;
//    for (auto& output : out_cont) {
//        objectives.emplace_back(output.objectives.begin(), output.objectives.end());
//    }
//
//    // ?????Python????
//    Py_Initialize();
//
//    // ????matplotlib???��???????????
//    PyRun_SimpleString("import matplotlib\n"
//                       "matplotlib.use('Agg')\n");
//    // ????pyplot???
//    PyObject* pyplot = PyImport_ImportModule("matplotlib.pyplot");
//    PyObject* mpl_toolkits = PyImport_ImportModule("mpl_toolkits");
//    PyObject* mplot3d = PyObject_GetAttrString(mpl_toolkits, "mplot3d");
//    PyObject* fig = PyObject_CallMethod(pyplot, "figure", NULL);
//    PyObject* ax = PyObject_CallMethod(mplot3d, "Axes3D", "O", fig);
//
//
////    Py_DECREF(matplotlib);
//
//    PyObject* x = PyList_New(out_cont.size());
//    PyObject* y = PyList_New(out_cont.size());
//    PyObject* z = PyList_New(out_cont.size());
//    PyObject* x_item;
//    PyObject* y_item;
//    PyObject* z_item;
//    PyObject* color;
//    PyObject* colors = PyList_New(out_cont.size());
//    for (int i = 0; i < out_cont.size(); i++) {
//        x_item = PyFloat_FromDouble(out_cont[i].objectives[0]);
//        PyList_SetItem(x, i, x_item);
//        y_item = PyFloat_FromDouble(out_cont[i].objectives[1]);
//        PyList_SetItem(y, i, y_item);
//        z_item = PyFloat_FromDouble(out_cont[i].objectives[2]);
//        PyList_SetItem(z, i, z_item);
//
//        if (out_cont[i].rm <= 10.0) {
//            color = PyUnicode_FromString("r");
//        } else {
//            color = PyUnicode_FromString("b");
//        }
//        PyList_SetItem(colors, i, color);
//    }
//
//
//// ??????��?????????????????
//    PyObject* title = PyUnicode_FromString("Pareto Front");
//    PyObject* xlabel = PyUnicode_FromString("vina");
//    PyObject* ylabel = PyUnicode_FromString("xscore");
//    PyObject* zlabel = PyUnicode_FromString("DLIGAND2");
//
////    PyObject_SetAttrString(pyplot, "title", title);
//
//
//    // ??????????????
//    PyObject* set_title = PyObject_GetAttrString(pyplot, "title");
//    PyObject* set_xlabel = PyObject_GetAttrString(pyplot, "xlabel");
//    PyObject* set_ylabel = PyObject_GetAttrString(pyplot, "ylabel");
//    PyObject* set_zlabel = PyObject_GetAttrString(pyplot, "zlabel");
//    PyObject_CallFunctionObjArgs(set_title, title, NULL);
//    PyObject_CallFunctionObjArgs(set_xlabel, xlabel, NULL);
//    PyObject_CallFunctionObjArgs(set_ylabel, ylabel, NULL);
//    PyObject_CallFunctionObjArgs(set_zlabel, zlabel, NULL);
//// ????numpy????
//    PyObject* np = PyImport_ImportModule("numpy");
//    PyObject* np_array = PyObject_CallMethod(np, "array", "O", colors);
//
//    // ??Python?????matplotlib????3D????
//    PyObject* scatter = PyObject_GetAttrString(ax, "scatter");
//    PyObject* scatter_args = PyTuple_Pack(4, x, y, z, colors);
//    PyObject_CallObject(scatter, scatter_args);
//
////    PyObject* draw = PyObject_GetAttrString(pyplot, "draw");
////    PyObject_CallObject(draw, NULL);
//    // save the figure to a file
////    PyObject_CallMethod(fig, "savefig", "s", "../coreset_pdbqt/1a30/pareto_front.png");
////    PyObject_CallMethod(pyplot, "show", nullptr);
//    extern std::string lig_name;
//    std::vector<std::string> out_list;
//    Stringsplit_(lig_name, '_',out_list);
//
////    std::string path = "./pareto/";
////    std::string out = out_list[0];
////    std::string index = out_list[2];
////    std::string filename = path + out + "/pareto_front_" + index + ".png";
////    PyObject_CallMethod(pyplot, "savefig", "s", filename.c_str());
//    extern string lig_name;
//    PyObject_CallMethod(pyplot, "savefig", "s", ("../coreset_pdbqt/"+lig_name+"/pareto_front"+"_"+std::to_string(i)+".png").c_str());
////    PyObject_CallMethod(pyplot, "savefig", "s", "/home/caokun/Cplusprojects/koto_MODEA_exp/LSHADE_Adam_final_nsga2/coreset_pdbqt/3aru");
//
//    PyObject_CallMethod(pyplot, "show", nullptr);
//    // ???Python????????
////    Py_DECREF(draw);
//    Py_DECREF(colors);
//    Py_DECREF(x);
//    Py_DECREF(y);
//    Py_DECREF(z);
//    Py_DECREF(colors);
//    Py_DECREF(fig);
//    Py_DECREF(ax);
//    Py_DECREF(pyplot);
//    Py_DECREF(mplot3d);
//    Py_DECREF(mpl_toolkits);
//    Py_Finalize();
//}

void merge_output_containers(const output_container_MODEA & in, output_container_MODEA& out, fl min_rmsd, sz max_size) {
//	VINA_FOR_IN(i, in)
//		add_to_output_container(out, in[i], min_rmsd, max_size);
}

void merge_output_containers(const parallel_mc_task_container& many, output_container_MODEA & out, fl min_rmsd, sz max_size) {
    min_rmsd = 2; // FIXME? perhaps it's necessary to separate min_rmsd during search and during output?
    VINA_FOR_IN(i, many)
        merge_output_containers(many[i].out, out, min_rmsd, max_size);
//    out.sort();
//    sort(out.begin(),out.end(), sort_comparator);
//	for (int i = 0; i < out.size() - 1; i++) {
//		for (int j = 0; j < out.size() - 1 - i; j++) {
//			if (sort_comparator(out[j+1], out[j])) {
//				output_type_MODEA tmp = out[j];
//				out[j] = out[j + 1];
//				out[j + 1] = tmp;
//			}
//		}
//	}
}
std::unordered_map<int, std::vector<output_type_MODEA>> fast_nondominated_sort(output_container_MODEA& pop) {
    std::unordered_map<int, std::vector<output_type_MODEA>> F;
    for (auto& p : pop) {
        //p.Slaves.clear();
        p.Slaves.clear();

        p.dominated_count = 0;
        for (auto& q : pop) {
            if (p < q) {
                // �����ǰ����p֧��q,��ô����֧��ĸ���q����Sp��
                p.Slaves.push_back(&q);
            }
            else if (q < p) {
                //�������q֧��p,��ô��np��ֵ��һ
                p.dominated_count += 1;
            }
        }
        //ѡ�����ʼ��Ⱥ��rank1�ĸ��壬���ҽ�rank1�еĸ������F[1]��
        if (p.dominated_count == 0) {
            p.rank = 1;
            F[1].push_back(p);
        }

    }

    int i = 1;
    std::vector<output_type_MODEA> Q;
    while (!F[i].empty()) {
        Q.clear();
        //ѭ����ǰF[i]�����������Ž⼯����rankn�еĽ⼯
        for (auto& p : F[i]) {
            //ѭ����ǰ���������Ž⼯�е�Sp���ϣ�����ǰ�����֧��Ľ�
            for (auto& q : p.Slaves) {
                //��Sp�е�ÿ�������np��-1����֧�䵱ǰ����Sp��ĸ�������-1
                q->dominated_count = q->dominated_count - 1;
                // �ж�Sp�����np�Ƿ�Ϊ0,���Ϊ0,���ʾ��ǰSp�еĸø����岻���κ�����������֧��
                // ����ǰSp��npΪ0�Ľ����rank+1�У�������Q�У���Ϊ��֧���
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
// �ҵ���������������ĵ�����
int getNearestCenterIndex(const output_type_MODEA& point, const std::vector<output_type_MODEA>& centers) {
    float minDist = std::numeric_limits<double>::max();
    int index = -1;
    for (size_t i = 0; i < centers.size(); ++i) {
        float dist = output_type_MODEA::distance(point, centers[i]);
        if (dist < minDist) {
            minDist = dist;
            index = i;
        }
    }
    return index;
}

// �����µ����ĵ�
output_type_MODEA* computeCenter(const std::vector<output_type_MODEA>& points) {
    std::vector<fl> objectives(3, 0.0);
    float sumX = 0, sumY = 0, sumZ = 0;
    for (const auto& point : points) {
        sumX += point.objectives[0];
        sumY += point.objectives[1];
        sumZ += point.objectives[2];

    }
    objectives[0] = sumX/points.size();
    objectives[1] = sumY/points.size();
    objectives[2] = sumZ/points.size();
    return new output_type_MODEA(objectives);
}

void KMeans(output_container_MODEA & out, int k){
    std::vector<output_type_MODEA> centers(k);
    std::vector<std::vector<output_type_MODEA>> clusters(k);
    bool areAlmostEqual(float a, float b, float tolerance);
    //�����ʼ�����ĵ�
    for(auto& center: centers){
        center = out[rand()%out.size()];
    }

    bool changed;
    do{
        for(auto& cluster: clusters){
            clusters.clear();
        }

        for(const auto& point: out){
            int centerIndex = getNearestCenterIndex(point,centers);
            clusters[centerIndex].push_back(point);
        }

        changed = false;
        for (size_t i = 0; i < k; ++i) {
            if (!clusters[i].empty()) {
                output_type_MODEA *newCenter = computeCenter(clusters[i]);
                if (!areAlmostEqual((*newCenter).objectives[0], centers[i].objectives[0],0.01) || !areAlmostEqual((*newCenter).objectives[1], centers[i].objectives[1],0.01) || !areAlmostEqual((*newCenter).objectives[2], centers[i].objectives[2],0.01)) {
                    centers[i] = *newCenter;
                    changed = true;
                }
            }
        }
    } while (changed);

    //��άɢ��ͼ
    // Ԥ����һЩ��ɫ
//    std::vector<std::string> colors = {"red", "blue", "green", "yellow", "cyan", "magenta"};
//    // ��ÿ���ؽ��е���
//    for (size_t i = 0; i < centers.size(); ++i) {
//        std::vector<double> x, y;
//
//        // ��ȡÿ���ص�x��y����
//        for (const auto& point : clusters[i]) {
//            x.push_back(point.objectives[0]);
//            y.push_back(point.objectives[1]);
//        }
//
//        // Ϊÿ����ѡ��һ����ɫ��������
//        matplotlibcpp::scatter(x, y, 50, {{"color", colors[i % colors.size()]}}); // ѡ����ɫ
//    }
//
//    // ��ʾͼ��
//    matplotlibcpp::legend();
//
//    // ��ʾ���Ƶ�ͼ
//    matplotlibcpp::show();
//    matplotlibcpp::save("../kmeans/clusters_2d.png"); // ����ͼ��
//    // ��ӡ���
//    for (int i = 0; i < centers.size(); ++i) {
//        std::cout << "Cluster " << i << ":\n";
//        for (const auto& point : clusters[i]) {
//            std::cout << "(" << point.objectives[0] << ", " << point.objectives[1] << ", " << point.objectives[2]<<")\n";
//        }
//    }


//������άɢ���ļ�
    std::ofstream file("../kmeans/clusters.csv");
    for (size_t i = 0; i < centers.size(); ++i) {
        for (const auto& point : clusters[i]) {
            file << point.objectives[0] << "," << point.objectives[1] << "," << point.objectives[2] << "," << point.e << "," << i << "\n";
        }
    }

    file.close();
}
// �ȽϺ���������out[i].e��ֵ��С��������
bool compare_e(const output_type_MODEA& a, const output_type_MODEA& b) {
    return a.e < b.e;
}
void merge_output_results( parallel_mc_task_container& many, output_container_MODEA & out, fl min_rmsd, sz max_size){
    min_rmsd = 2;
    output_container_MODEA out_temp;
    extern std::vector<output_container_MODEA> out_tmp;
//    float get_objective_distance(output_type_MODEA outt);
    VINA_FOR_IN(i,many)
    {
//        draw_pareto_(many[i].out,i);
        out_tmp.push_back(many[i].out);
        VINA_FOR_IN(j,many[i].out){
//            if(get_objective_distance(many[i].out[j]) < 2.0){
            many[i].out[j].threadd = i;
                add_to_output_container(out, many[i].out[j], min_rmsd, max_size);

//            }

        }
    }
//    KMeans(out,3);

//            out_temp.push_back(new output_type_MODEA(many[i].out[j]));
    out_temp = out;
    std::unordered_map<int, std::vector<output_type_MODEA>> F = fast_nondominated_sort(out_temp);
////
    output_container_MODEA out_results;
    for(int i=0;i<F[1].size();i++){
        out_results.push_back(new output_type_MODEA(F[1][i]));
    }

//    for(int i = 0; i < out_temp.size(); i++){
//        const fl ub = ;
//        m.set(out_temp[i].c);
//        if(m.rmsd_upper_bound(r))
//    }


//    for(int j = 0;j<out_results.size();j++){
//        std::pair<sz, fl> closest_rmsd = find_closest(out_results[j].coords, out);
//        if(closest_rmsd.first < out.size() && closest_rmsd.second < min_rmsd) { // have a very similar one
//            if(out_results[j].e < out[closest_rmsd.first].e) { // the new one is better, apparently
//                out[closest_rmsd.first] = out_results[j]; // FIXME? slow
//            }
//        }
//        else { // nothing similar
//            if(out.size() < max_size)
//                out.push_back(new output_type_MODEA(out_results[j])); // the last one had the worst energy - replacing
//            else
//            if(!out.empty() && out_results[j].e < out.back().e) // FIXME? - just changed
//                out.back() = out_results[j]; // FIXME? slow
//        }
//    }
    out = out_results;

//    VINA_FOR_IN(m,out)
//        out[m].e = out[m].objectives[

    std::sort(out.begin(),out.end(),compare_e);

//    for (int i = 0; i < out.size() - 1; i++) {
//        for (int j = 0; j < out.size() - 1 - i; j++) {
//            if (sort_comparator(out[j+1], out[j])) {
//                output_type_MODEA tmp = out[j];
//                out[j] = out[j + 1];
//                out[j + 1] = tmp;
//            }
//        }
//    }


}


void parallel_mc::operator()(const model& m, const model& m_ad4, output_container_MODEA& out, const precalculate& p, const precalculate& prec_ad4, const precalculate& prec_vdw, const precalculate& prec_forcefield, const precalculate& prec_DLIGAND2, const igrid& ig, const igrid& DLIGAND2_grid, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, rng& generator,const scoring_function& sf,const scoring_function& sf_ad4,const boost::optional<model>& ref) const {
    parallel_progress pp;
    parallel_mc_aux parallel_mc_aux_instance(&modea, &p, &prec_ad4, &prec_vdw, &prec_forcefield, &prec_DLIGAND2, &ig, &DLIGAND2_grid, &p_widened, &ig_widened, &corner1, &corner2, (display_progress ? (&pp) : NULL),&sf,&sf_ad4,ref);
    parallel_mc_task_container task_container;
    VINA_FOR(i, num_tasks)
        task_container.push_back(new parallel_mc_task(m, m_ad4, random_int(0, 1000000, generator)));
    if(display_progress)
        pp.init(num_tasks * modea.num_steps);
    parallel_iter<parallel_mc_aux, parallel_mc_task_container, parallel_mc_task, true> parallel_iter_instance(&parallel_mc_aux_instance, num_threads);
    parallel_iter_instance.run(task_container);

    merge_output_results(task_container, out, modea.min_rmsd, modea.num_saved_mins);
}





//SMPSO
/*bool sort_comparator(const output_type_SMPSO& x,const output_type_SMPSO& y);
struct parallel_mc_task {
	model m;
	output_container_SMPSO out;
	rng generator;
	parallel_mc_task(const model& m_, int seed) : m(m_), generator(static_cast<rng::result_type>(seed)) {}
};

typedef boost::ptr_vector<parallel_mc_task> parallel_mc_task_container;

struct parallel_mc_aux {
	const lshade* ls;
    const SMPSO* smpso;
	const precalculate* p;
    const precalculate* prec_forcefield;
	const igrid* ig;
	const precalculate* p_widened;
	const igrid* ig_widened;
	const vec* corner1;
	const vec* corner2;
	parallel_progress* pg;
    const scoring_function* sf;
    const boost::optional<model>& ref;
	parallel_mc_aux(const SMPSO* smpso_, const precalculate* p_, const precalculate* prec_forcefield_, const igrid* ig_, const precalculate* p_widened_, const igrid* ig_widened_, const vec* corner1_, const vec* corner2_, parallel_progress* pg_,const scoring_function* sf_,const boost::optional<model>& ref_)
		: smpso(smpso_), p(p_), prec_forcefield(prec_forcefield_), ig(ig_), p_widened(p_widened_), ig_widened(ig_widened_), corner1(corner1_), corner2(corner2_), pg(pg_),sf(sf_),ref(ref_) {}
	void operator()(parallel_mc_task& t) const {
		(*smpso)(t.m, t.out, *p, *prec_forcefield,*ig, *p_widened, *ig_widened, *corner1, *corner2, pg, t.generator,*sf,ref);
	}
};

void merge_output_containers(const output_container_SMPSO & in, output_container_SMPSO& out, fl min_rmsd, sz max_size) {
//	VINA_FOR_IN(i, in)
//		add_to_output_container(out, in[i], min_rmsd, max_size);
}

void merge_output_containers(const parallel_mc_task_container& many, output_container_SMPSO & out, fl min_rmsd, sz max_size) {
	min_rmsd = 2; // FIXME? perhaps it's necessary to separate min_rmsd during search and during output?
	VINA_FOR_IN(i, many)
		merge_output_containers(many[i].out, out, min_rmsd, max_size);
//    out.sort();
//    sort(out.begin(),out.end(), sort_comparator);
//	for (int i = 0; i < out.size() - 1; i++) {
//		for (int j = 0; j < out.size() - 1 - i; j++) {
//			if (sort_comparator(out[j+1], out[j])) {
//				output_type_SMPSO tmp = out[j];
//				out[j] = out[j + 1];
//				out[j + 1] = tmp;
//			}
//		}
//	}
}
//std::unordered_map<int, std::vector<output_type_SMPSO>> fast_nondominated_sort1(output_container_SMPSO& pop) {
//    std::unordered_map<int, std::vector<output_type_SMPSO>> F;
//    for (auto& p : pop) {
//        //p.Slaves.clear();
//        p.Slaves.clear();
//
//        p.dominated_count = 0;
//        for (auto& q : pop) {
//            if (p < q) {
//                // �����ǰ����p֧��q,��ô����֧��ĸ���q����Sp��
//                p.Slaves.push_back(&q);
//            }
//            else if (q < p) {
//                //�������q֧��p,��ô��np��ֵ��һ
//                p.dominated_count += 1;
//            }
//        }
//        //ѡ�����ʼ��Ⱥ��rank1�ĸ��壬���ҽ�rank1�еĸ������F[1]��
//        if (p.dominated_count == 0) {
//            p.rank = 1;
//            F[1].push_back(p);
//        }
//
//    }
//
//    int i = 1;
//    std::vector<output_type_SMPSO> Q;
//    while (!F[i].empty()) {
//        Q.clear();
//        //ѭ����ǰF[i]�����������Ž⼯����rankn�еĽ⼯
//        for (auto& p : F[i]) {
//            //ѭ����ǰ���������Ž⼯�е�Sp���ϣ�����ǰ�����֧��Ľ�
//            for (auto& q : p.Slaves) {
//                //��Sp�е�ÿ�������np��-1����֧�䵱ǰ����Sp��ĸ�������-1
//                q->dominated_count = q->dominated_count - 1;
//                // �ж�Sp�����np�Ƿ�Ϊ0,���Ϊ0,���ʾ��ǰSp�еĸø����岻���κ�����������֧��
//                // ����ǰSp��npΪ0�Ľ����rank+1�У�������Q�У���Ϊ��֧���
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
void draw_pareto(const output_container_SMPSO& out_cont);
// �ȽϺ���������out[i].e��ֵ��С��������
bool compare_e(const output_type_SMPSO& a, const output_type_SMPSO& b) {
    return a.e < b.e;
}
void merge_output_results(const parallel_mc_task_container& many, output_container_SMPSO & out, fl min_rmsd, sz max_size){
    min_rmsd = 2;
    output_container_SMPSO out_temp;
    VINA_FOR_IN(i,many)
        VINA_FOR_IN(j,many[i].out)
            out_temp.push_back(new output_type_SMPSO(many[i].out[j]));

//    std::unordered_map<int, std::vector<output_type_SMPSO>> F = fast_nondominated_sort1(out_temp);
//
//    output_container_SMPSO out_results;
//    for(int i=0;i<F[1].size();i++){
//        out_results.push_back(new output_type_SMPSO(F[1][i]));
//    }
//    for(int j = 0;j<out_results.size();j++){
//        std::pair<sz, fl> closest_rmsd = find_closest(out_results[j].coords, out);
//        if(closest_rmsd.first < out.size() && closest_rmsd.second < min_rmsd) { // have a very similar one
//            if(out_results[j].e < out[closest_rmsd.first].e) { // the new one is better, apparently
//                out[closest_rmsd.first] = out_results[j]; // FIXME? slow
//            }
//        }
//        else { // nothing similar
//            if(out.size() < max_size)
//                out.push_back(new output_type_SMPSO(out_results[j])); // the last one had the worst energy - replacing
//            else
//            if(!out.empty() && out_results[j].e < out.back().e) // FIXME? - just changed
//                out.back() = out_results[j]; // FIXME? slow
//        }
//    }
    out = out_temp;
//    VINA_FOR_IN(m,out)
//        out[m].e = out[m].objectives[
    std::sort(out.begin(),out.end(),compare_e);

//    for (int i = 0; i < out.size() - 1; i++) {
//        for (int j = 0; j < out.size() - 1 - i; j++) {
//            if (sort_comparator(out[j+1], out[j])) {
//                output_type_SMPSO tmp = out[j];
//                out[j] = out[j + 1];
//                out[j + 1] = tmp;
//            }
//        }
//    }


}

void parallel_mc::operator()(const model& m, output_container_SMPSO& out, const precalculate& p, const precalculate& prec_forcefield,const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, rng& generator,const scoring_function& sf,const boost::optional<model>& ref) const {
	parallel_progress pp;
	parallel_mc_aux parallel_mc_aux_instance(&smpso, &p, &prec_forcefield, &ig, &p_widened, &ig_widened, &corner1, &corner2, (display_progress ? (&pp) : NULL),&sf,ref);
	parallel_mc_task_container task_container;
	VINA_FOR(i, num_tasks)
		task_container.push_back(new parallel_mc_task(m, random_int(0, 1000000, generator)));
	if(display_progress)
		pp.init(num_tasks * smpso.num_steps);
	parallel_iter<parallel_mc_aux, parallel_mc_task_container, parallel_mc_task, true> parallel_iter_instance(&parallel_mc_aux_instance, num_threads);
	parallel_iter_instance.run(task_container);
    merge_output_results(task_container, out, smpso.min_rmsd, smpso.num_saved_mins);
}*/








//nsga2
/*
#include "parallel.h"
#include "parallel_mc.h"
#include "coords.h"
#include "parallel_progress.h"
bool sort_comparator(const output_type_nsga2& x,const output_type_nsga2& y);
struct parallel_mc_task {
    model m;
    output_container_nsga2 out;
    rng generator;
    parallel_mc_task(const model& m_, int seed) : m(m_), generator(static_cast<rng::result_type>(seed)) {}
};

typedef boost::ptr_vector<parallel_mc_task> parallel_mc_task_container;

struct parallel_mc_aux {
    const lshade* ls;
    const nsga2* ns;
    const precalculate* p;
    const precalculate* prec_forcefield;
    const igrid* ig;
    const precalculate* p_widened;
    const igrid* ig_widened;
    const vec* corner1;
    const vec* corner2;
    parallel_progress* pg;
    const scoring_function* sf;
    const boost::optional<model>& ref;
    parallel_mc_aux(const nsga2* ns_, const precalculate* p_, const precalculate* prec_forcefield_, const igrid* ig_, const precalculate* p_widened_, const igrid* ig_widened_, const vec* corner1_, const vec* corner2_, parallel_progress* pg_,const scoring_function* sf_,const boost::optional<model>& ref_)
            : ns(ns_), p(p_), prec_forcefield(prec_forcefield_), ig(ig_), p_widened(p_widened_), ig_widened(ig_widened_), corner1(corner1_), corner2(corner2_), pg(pg_),sf(sf_),ref(ref_) {}
    void operator()(parallel_mc_task& t) const {
        (*ns)(t.m, t.out, *p, *prec_forcefield,*ig, *p_widened, *ig_widened, *corner1, *corner2, pg, t.generator,*sf,ref);
    }
};

void merge_output_containers(const output_container_nsga2 & in, output_container_nsga2& out, fl min_rmsd, sz max_size) {
//	VINA_FOR_IN(i, in)
//		add_to_output_container(out, in[i], min_rmsd, max_size);
}

void merge_output_containers(const parallel_mc_task_container& many, output_container_nsga2 & out, fl min_rmsd, sz max_size) {
    min_rmsd = 2; // FIXME? perhaps it's necessary to separate min_rmsd during search and during output?
    VINA_FOR_IN(i, many)
        merge_output_containers(many[i].out, out, min_rmsd, max_size);
//    out.sort();
//    sort(out.begin(),out.end(), sort_comparator);
//	for (int i = 0; i < out.size() - 1; i++) {
//		for (int j = 0; j < out.size() - 1 - i; j++) {
//			if (sort_comparator(out[j+1], out[j])) {
//				output_type_nsga2 tmp = out[j];
//				out[j] = out[j + 1];
//				out[j + 1] = tmp;
//			}
//		}
//	}
}
std::unordered_map<int, std::vector<output_type_nsga2>> fast_nondominated_sort1(output_container_nsga2& pop) {
    std::unordered_map<int, std::vector<output_type_nsga2>> F;
    for (auto& p : pop) {
        //p.Slaves.clear();
        p.Slaves.clear();

        p.dominated_count = 0;
        for (auto& q : pop) {
            if (p < q) {
                // �����ǰ����p֧��q,��ô����֧��ĸ���q����Sp��
                p.Slaves.push_back(&q);
            }
            else if (q < p) {
                //�������q֧��p,��ô��np��ֵ��һ
                p.dominated_count += 1;
            }
        }
        //ѡ�����ʼ��Ⱥ��rank1�ĸ��壬���ҽ�rank1�еĸ������F[1]��
        if (p.dominated_count == 0) {
            p.rank = 1;
            F[1].push_back(p);
        }

    }

    int i = 1;
    std::vector<output_type_nsga2> Q;
    while (!F[i].empty()) {
        Q.clear();
        //ѭ����ǰF[i]�����������Ž⼯����rankn�еĽ⼯
        for (auto& p : F[i]) {
            //ѭ����ǰ���������Ž⼯�е�Sp���ϣ�����ǰ�����֧��Ľ�
            for (auto& q : p.Slaves) {
                //��Sp�е�ÿ�������np��-1����֧�䵱ǰ����Sp��ĸ�������-1
                q->dominated_count = q->dominated_count - 1;
                // �ж�Sp�����np�Ƿ�Ϊ0,���Ϊ0,���ʾ��ǰSp�еĸø����岻���κ�����������֧��
                // ����ǰSp��npΪ0�Ľ����rank+1�У�������Q�У���Ϊ��֧���
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
void draw_pareto(const output_container_nsga2& out_cont);
// �ȽϺ���������out[i].e��ֵ��С��������
bool compare_e(const output_type_nsga2& a, const output_type_nsga2& b) {
    return a.e < b.e;
}
void merge_output_results(const parallel_mc_task_container& many, output_container_nsga2 & out, fl min_rmsd, sz max_size){
    min_rmsd = 2;
    output_container_nsga2 out_temp;
    VINA_FOR_IN(i,many)
        VINA_FOR_IN(j,many[i].out)
            out_temp.push_back(new output_type_nsga2(many[i].out[j]));

    std::unordered_map<int, std::vector<output_type_nsga2>> F = fast_nondominated_sort1(out_temp);

    output_container_nsga2 out_results;
    for(int i=0;i<F[1].size();i++){
        out_results.push_back(new output_type_nsga2(F[1][i]));
    }
//    for(int j = 0;j<out_results.size();j++){
//        std::pair<sz, fl> closest_rmsd = find_closest(out_results[j].coords, out);
//        if(closest_rmsd.first < out.size() && closest_rmsd.second < min_rmsd) { // have a very similar one
//            if(out_results[j].e < out[closest_rmsd.first].e) { // the new one is better, apparently
//                out[closest_rmsd.first] = out_results[j]; // FIXME? slow
//            }
//        }
//        else { // nothing similar
//            if(out.size() < max_size)
//                out.push_back(new output_type_nsga2(out_results[j])); // the last one had the worst energy - replacing
//            else
//            if(!out.empty() && out_results[j].e < out.back().e) // FIXME? - just changed
//                out.back() = out_results[j]; // FIXME? slow
//        }
//    }
    out = out_results;
//    VINA_FOR_IN(m,out)
//        out[m].e = out[m].objectives[
    std::sort(out.begin(),out.end(),compare_e);

//    for (int i = 0; i < out.size() - 1; i++) {
//        for (int j = 0; j < out.size() - 1 - i; j++) {
//            if (sort_comparator(out[j+1], out[j])) {
//                output_type_nsga2 tmp = out[j];
//                out[j] = out[j + 1];
//                out[j + 1] = tmp;
//            }
//        }
//    }


}

void parallel_mc::operator()(const model& m, output_container_nsga2& out, const precalculate& p, const precalculate& prec_forcefield,const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, rng& generator,const scoring_function& sf,const boost::optional<model>& ref) const {
    parallel_progress pp;
    parallel_mc_aux parallel_mc_aux_instance(&ns, &p, &prec_forcefield, &ig, &p_widened, &ig_widened, &corner1, &corner2, (display_progress ? (&pp) : NULL),&sf,ref);
    parallel_mc_task_container task_container;
    VINA_FOR(i, num_tasks)
        task_container.push_back(new parallel_mc_task(m, random_int(0, 1000000, generator)));
    if(display_progress)
        pp.init(num_tasks * ns.num_steps);
    parallel_iter<parallel_mc_aux, parallel_mc_task_container, parallel_mc_task, true> parallel_iter_instance(&parallel_mc_aux_instance, num_threads);
    parallel_iter_instance.run(task_container);
    merge_output_results(task_container, out, ns.min_rmsd, ns.num_saved_mins);
}*/

