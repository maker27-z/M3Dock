//
// Created by 91686 on 2023/10/16.
//
#include "model.h"
#include "conf.h"
#include "incrementable.h"
#include <algorithm>
#include "random.h"
#include "energy.h"
#include "unordered_map"
//#include <boost/python.hpp>
#ifndef LSHADE_ADAM_FINAL_MODEA_H
#define LSHADE_ADAM_FINAL_MODEA_H



struct MODEA {
    sz num_saved_mins;//���out������
    int num_steps;//��������
    fl min_rmsd;//�������Ľ��
    MODEA() : num_steps(800), num_saved_mins(20), min_rmsd(1.0) { }
    void operator()(model& m, model& m_ad4,output_container_MODEA& out, const precalculate& p, const precalculate& prec_ad4, const precalculate& prec_vdw, const precalculate& prec_forcefield, const precalculate& prec_DLIGAND2, const igrid& ig, const igrid& DLIGAND2_grid, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator, const scoring_function& sf,const scoring_function& sf_ad4,const boost::optional<model>& ref) const;
    void init_individual1(model& m, energy_cal& energy, rng& generator, const vec& corner1, const vec& corner2, const scoring_function& sf, output_type_MODEA& tmp) const;
};

struct MODEA_aux {
public:
    std::vector<output_type_MODEA> archive_gradient;
    int popsize;//��Ⱥ��ģ
    int dim;//�����ά�ȣ�7+n
    int gen;//��Ⱥ��������
//    std::vector<output_type_MODEA> P;//��һ����Ⱥ����P
//    std::vector<output_type_MODEA> Q;//���ɵ�һ����Ⱥ���Ӵ�Q
    std::vector<output_type_MODEA> R;//������Ⱥ�ϲ����R
    output_type_MODEA init_tmp;
    fl F;//????????
    fl CR;//???????

    int A_size;//�ⲿ�浵������
    std::vector<output_type_MODEA> oppositepopulation;//��ʼ�෴��Ⱥ
    std::vector<output_type_MODEA> A;//�ⲿ�浵
    std::vector<output_type_MODEA> nowpopulation;//��ʼ��Ⱥ
    std::vector<output_type_MODEA> mutatepopulation;//���콻������Ⱥ
    int best_index; //��¼��ǰ��Ⱥ�������ŵĸ���
    MODEA_aux() :popsize(100), dim(7) {}
    MODEA_aux(const model& m) {
        archive_gradient.resize(20);
        dim = m.ligand_degrees_of_freedom(0) + 6;//ֻ����ֻ��һ����������
        popsize = 100;
        nowpopulation.resize(popsize);
        oppositepopulation.resize(popsize);
        mutatepopulation.resize(popsize);
        F = 0.4;
        CR = 0.9;
//        P.resize(popsize);
        //Q.resize(popsize);
        R.resize(popsize * 2);
    }
    std::vector<output_type_MODEA> make_new_pop(model& m, const precalculate& p, const igrid& ig, vec v, std::vector<output_type_MODEA> pop, energy_cal& energy, const vec& corner1, const vec& corner2, rng& generator, const scoring_function& sf,const boost::optional<model>& ref,change &g,const scoring_function& sf_ad4);//��������Ⱥ
    std::vector<output_type_MODEA> crossover_mutation(model& m, const precalculate& p, const igrid& ig, vec v, energy_cal& energy, rng& generator, const vec& corner1, const vec& corner2, output_type_MODEA parent1, output_type_MODEA parent2, const scoring_function& sf,const boost::optional<model>& ref,change &g,const scoring_function& sf_ad4);//�����������
    output_type_MODEA binary_tournament(output_type_MODEA indivisual1, output_type_MODEA indivisual2);//��Ԫ������ѡ��
    void init_individual(model& m, energy_cal& energy, rng& generator, const vec& corner1, const vec& corner2, const scoring_function& sf, output_type_MODEA& tmp,const boost::optional<model>& ref,const scoring_function& sf_ad4);



    //MODEA
    void mutate_cross(model& m, energy_cal& energy, rng& generator, const vec& corner1, const vec& corner2);
    void select(model& m, energy_cal& energy, rng& generator);
    void calculate_objectives(model& m, const precalculate& p, const igrid& ig, vec v, output_type_MODEA& indivisual, const scoring_function& sf, const scoring_function& sf_ad4, energy_cal& energy,const boost::optional<model>& ref,change &g, change &vdw_g, change &vina_g, change &DLIGAND2_g);
    std::unordered_map<int, std::vector<output_type_MODEA>> fast_nondominated_sort(std::vector<output_type_MODEA>& pop);//���ٷ�֧������
    void calculate_crowding_distance(std::vector<output_type_MODEA>& rank_i_vector);//����ӵ������
    void nd_cd_sort(std::vector<output_type_MODEA>& pop);
    int dominated_individual_archive( output_type_MODEA x, output_type_MODEA y){
        if(x < y){
            return 1;
        }else if(y < x){
            return -1;
        }else{
            return 0;
        }

    }
//    bool crowding_distance_comparator(const output_type_MODEA& x, const output_type_MODEA& y) {
//        return x.crowding_distance > y.crowding_distance;
//    }
    int indexWorst(){
        if(archive_gradient.empty()){
            return -1;
        }
        int index = 0;
        output_type_MODEA worstKnown = archive_gradient[0];
        output_type_MODEA candidateSolution;

        int flag;
        for(int i = 1;i < archive_gradient.size();i++){
            candidateSolution = archive_gradient[i];
            bool crowding_distance_comparator(const output_type_MODEA& x, const output_type_MODEA& y);
            flag = crowding_distance_comparator(worstKnown,candidateSolution);
            if(flag){
                index = i;
                worstKnown =candidateSolution;
            }
        }
        return index;
    }
    int add(output_type_MODEA outt){
        int i = 0;
        while (i < archive_gradient.size()){
            int flag = dominated_individual_archive(outt,archive_gradient[i]);
            if(flag == -1){
                return -1;
            }else if(flag == 1){
                archive_gradient.erase(archive_gradient.begin()+i);
            }else{
                if(outt == archive_gradient[i]){
                    return -1;
                }
                i++;
            }


        }

        bool res = 1;
        if(archive_gradient.size() > 20){
            calculate_crowding_distance(archive_gradient);
            int indexWorst_ = indexWorst();
            if(outt == archive_gradient[indexWorst_]){//??
                res = -1;
            }else{
                archive_gradient.erase(archive_gradient.begin()+indexWorst_);
                archive_gradient.push_back(outt);
                return res;
            }
        }
        archive_gradient.push_back(outt);
        return res;
    }
};

//struct python_op{
//    boost::python::object find_min_norm_element;
//
//    boost::python::object sys;
//    boost::python::object pathss;
//    boost::python::object module;
//    boost::python::object solver_class;
//    boost::python::object py_args;
//    boost::python::object result;
//    python_op(){
//        Py_Initialize();
//
//        sys = boost::python::import("sys");
//        pathss = sys.attr("path");
//        pathss.attr("append")("/home/caokun/Cplusprojects/koto_MODEA_gradient1_exp/LSHADE_Adam_final_nsga2");
//        module = boost::python::import("min_norm_solvers_numpy");
//        solver_class = module.attr("MinNormSolverNumpy");
//        find_min_norm_element = solver_class.attr("find_min_norm_element");
//    }
//};

struct Solution {
    std::pair<int, int> indexPair;
    float c;
    float d;
};
class MinNormSolverNumpy {
    int MAX_ITER = 250;
    int STOP_CRIT = 1e-6;


public: float dot(const std::vector<float>& a, const std::vector<float>& b) {
        float sum = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            sum += a[i] * b[i];
        }
        return sum;
    }
    std::vector<float> calculate_grad_dir(const std::vector<std::vector<float>>& grad_mat,
        const std::vector<float>& sol_vec)//sol_vec里面元素和为1
    {

        std::vector<float> grad_dir(sol_vec.size(), 0.0);

        for (size_t i = 0; i < grad_mat.size(); ++i) {
            grad_dir[i] = -1.0 * std::inner_product(grad_mat[i].begin(), grad_mat[i].end(), sol_vec.begin(), 0.0);
        }

        return grad_dir;
    }

    std::vector<float> projection2simplex(std::vector<float> y){
        size_t m = y.size();
        // 对y进行降序排序
        std::sort(y.begin(), y.end(), std::greater<float>());
        float tmpsum = 0.0;
        float tmax_f = std::accumulate(y.begin(), y.end(), 0.0) / m - 1.0 / m;

        for (size_t i = 0; i < m - 1; ++i) {
            tmpsum += y[i];
            float tmax = (tmpsum - 1.0) / (i + 1.0);
            if (tmax > y[i + 1]) {
                tmax_f = tmax;
                break;
            }
        }
        // 计算投影
        std::vector<float> projection(m);
        for (size_t i = 0; i < m; ++i) {
            projection[i] = std::max(y[i] - tmax_f, 0.0f);
        }
        return projection;
    }

    std::vector<float> next_point(const std::vector<float>& cur_val, const std::vector<float>& grad, int n){
        std::vector<float> proj_grad(n);
        float grad_sum = std::accumulate(grad.begin(), grad.end(), 0.0);
        for (int i = 0; i < n; ++i) {
            proj_grad[i] = grad[i] - (grad_sum / n);
        }
        std::vector<float> tm1;
        std::vector<float> tm2;
        for (int i = 0; i < n; ++i) {
            if (proj_grad[i] < 0) {
                tm1.push_back(-1.0 * cur_val[i] / proj_grad[i]);
            } else if (proj_grad[i] > 0) {
                tm2.push_back((1.0 - cur_val[i]) / proj_grad[i]);
            }
        }
        float t = 1.0;
        // 寻找 tm1 中大于 1e-7 的最小值
        float vmin = std::numeric_limits<float>::max();
        for (float val : tm1) {
            if (val > 1e-7) {
                if(val<vmin){
                    vmin = val;
                    t =  vmin;
                }
            }
        }

// 寻找 tm2 中大于 1e-7 的最小值，并与当前 t 的值比较
        for (float val : tm2) {
            if (val > 1e-7) {
                t = std::min(t, val);
            }
        }

        std::vector<float> next_point(n);
        for (int i = 0; i < n; ++i) {
            next_point[i] = proj_grad[i] * t + cur_val[i];
        }
        next_point = projection2simplex(next_point);
        return next_point;
//        if (!tm1.empty()) {
//            t = *std::min_element(tm1.begin(), tm1.end());
//        }
//        if (!tm2.empty()) {
//            t = std::min(t, *std::min_element(tm2.begin(), tm2.end()));
//        }
    }

    std::pair<float, float> min_norm_element_from2(float v1v1, float v1v2, float v2v2){
        float gamma, cost;
        if (v1v2 >= v1v1) {
            gamma = 0.999;
            cost = v1v1;
        } else if (v1v2 >= v2v2) {
            gamma = 0.001;
            cost = v2v2;
        } else {
            gamma = -1.0 * ((v1v2 - v2v2) / (v1v1 + v2v2 - 2 * v1v2));
            cost = v2v2 + gamma * (v1v2 - v2v2);
        }
        return {gamma, cost};
    }


    std::pair<Solution, std::map<std::pair<int, int>, float>> min_norm_2d(const std::vector<std::vector<float>>& vecs, std::map<std::pair<int, int>, float>& dps){
        float dmin = std::numeric_limits<float>::max();
        Solution sol;

        int n = vecs.size();
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {

                if (!dps.count({i, j})) {
                    dps[{i, j}] = dot(vecs[i], vecs[j]);
                    dps[{j, i}] = dps[{i, j}];
                }

                if (!dps.count({i, i})) {
                    dps[{i, i}] = dot(vecs[i], vecs[i]);
                }

                if (!dps.count({j, j})) {
                    dps[{j, j}] = dot(vecs[j], vecs[j]);
                }
                //c为两个向量间的权重α，d暂时看不出来
                auto [c, d] =
                    min_norm_element_from2(dps[{i, i}], dps[{i, j}], dps[{j, j}]);
                if (d < dmin) {
                    dmin = d;
                    sol = {{i, j}, c, d};
                }

            }
        }
        return {sol, dps};
    }
    std::pair<std::vector<float>, float> find_min_norm_element(const std::vector<std::vector<float>>& vecs) {

        std::map<std::pair<int, int>, float> dps;

        int n = vecs.size();
        std::vector<float> sol_vec(n, 0.0);

        //dpss的大小与dps一样，
        //dps里面是存vecs元素两两内积，对称矩阵
        auto [init_sol, dpss] =
            min_norm_2d(vecs,dps);

        //c为两个向量间的权重α
        sol_vec[init_sol.indexPair.first] = init_sol.c;
        sol_vec[init_sol.indexPair.second] = 1 - init_sol.c;

        //n=3,所以不执行
        if (n < 3) {
            return {sol_vec, init_sol.d};
        }

        int iter_count = 0;
        std::vector<std::vector<float>> grad_mat(n, std::vector<float>(n, 0.0));

        //grad_mat是二维vector数组，dpss是二维矩阵，把dpss里面的值一一复制到grad_mat里面
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                grad_mat[i][j] = dpss[{i, j}];
            }
        }

        while (iter_count < MAX_ITER){
            //grad_dir为grad_mat中各个向量与向量sol_vec的内积，sol_vec各元素之和为1即α权重
            std::vector<float> grad_dir = calculate_grad_dir(grad_mat, sol_vec);
            std::vector<float> new_point = next_point(sol_vec, grad_dir, n);
            float v1v1 = 0.0, v1v2 = 0.0, v2v2 = 0.0;
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    v1v1 += sol_vec[i] * sol_vec[j] * dps[{i, j}];
                    v1v2 += sol_vec[i] * new_point[j] * dps[{i, j}];
                    v2v2 += new_point[i] * new_point[j] * dps[{i, j}];
                }
            }
            auto [nc, nd] = min_norm_element_from2(v1v1, v1v2, v2v2);
            std::vector<float> new_sol_vec(sol_vec.size());
            for (int i = 0; i < n; ++i) {
                new_sol_vec[i] = nc * sol_vec[i] + (1 - nc) * new_point[i];
            }
            std::vector<float> change(sol_vec.size());
            for (size_t i = 0; i < sol_vec.size(); ++i) {
                change[i] = std::abs(new_sol_vec[i] - sol_vec[i]);
            }
            float sum_change = std::accumulate(change.begin(), change.end(), 0.0);
            if (sum_change < MinNormSolverNumpy::STOP_CRIT) {
                return {sol_vec, nd};
            }
            sol_vec = new_sol_vec; // 更新 sol_vec 为新的解向量
            return {sol_vec, nd};
        }
    }
};

class PCGrad {
public:
    float dot(const std::vector<float>& a, const std::vector<float>& b) {
        float sum = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            sum += a[i] * b[i];
        }
        return sum;
    }
    std::vector<float> co_vector(float coefficient,std::vector<float>& gradient){
        for(int i = 0; i < gradient.size(); i++){
            gradient[i]*=coefficient;
        }
        return gradient;

    }
    void substrct_vector(std::vector<float>& gradient1,std::vector<float> gradient2){
        VINA_FOR(i,gradient1.size()){
            gradient1[i]-=gradient2[i];
        }
//        return gradient1;
    }
    void plus_vector(std::vector<float>& gradient1,std::vector<float> gradient2){
        VINA_FOR(i,gradient1.size()){
            gradient1[i]+=gradient2[i];
        }
//        return gradient1;
    }

    void project_gradients(std::vector<std::vector<float>> gradients,std::vector<std::vector<float>>& gradients_PC){
        for(int i = 0; i < gradients_PC.size();i++){
            for(int j = 0; j < gradients.size();j++){
                if(dot(gradients_PC[i],gradients[j])<0.0){
                    substrct_vector(gradients_PC[i] ,co_vector(
                            dot(gradients_PC[i],gradients[j])
                            /
                            dot(gradients[j],gradients[j])
                            ,gradients[j]));
                }
            }
        }
    }


};

#endif //LSHADE_ADAM_FINAL_MODEA_H
