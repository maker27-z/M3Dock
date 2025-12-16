/*
everything
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

#include "everything.h"
#include "int_pow.h"

inline fl gaussian(fl x, fl width) {
    return std::exp(-sqr(x/width));
}

// distance_additive terms

template<unsigned i>
struct electrostatic : public distance_additive {
    fl cap;
    electrostatic(fl cap_, fl cutoff_) : distance_additive(cutoff_), cap(cap_) {
        name = std::string("electrostatic(i=") + to_string(i) + ", ^=" + to_string(cap) + ", c=" + to_string(cutoff) + ")";
    }
    fl eval(const atom_base& a, const atom_base& b, fl r) const {
        fl tmp = int_pow<i>(r);
        fl q1q2 = a.charge * b.charge;
        if(tmp < epsilon_fl) return q1q2 * cap;
        else                 return q1q2 * (std::min)(cap, 1/int_pow<i>(r));
    }
};

fl solvation_parameter(const atom_type& a) {
    if(a.ad < AD_TYPE_SIZE) return ad_type_property(a.ad).solvation;
    else if(a.xs == XS_TYPE_Met_D) return metal_solvation_parameter;
    VINA_CHECK(false);
    return 0; // placating the compiler
}

fl      volume(const atom_type& a) {
    if(a.ad < AD_TYPE_SIZE) return ad_type_property(a.ad).volume;
    else if(a.xs < XS_TYPE_SIZE) return 4*pi / 3 * int_pow<3>(xs_radius(a.xs));
    VINA_CHECK(false);
    return 0; // placating the compiler
}


struct forcefiled_solvation : public distance_additive {
    fl desolvation_sigma;
    fl solvation_q;
    bool charge_dependent;
    forcefiled_solvation(fl desolvation_sigma_, fl solvation_q_, bool charge_dependent_, fl cutoff_) : distance_additive(cutoff_), solvation_q(solvation_q_), charge_dependent(charge_dependent_), desolvation_sigma(desolvation_sigma_) {
        name = std::string("forcefiled_solvation(d-sigma=") + to_string(desolvation_sigma) + ", s/q=" + to_string(solvation_q) + ", q=" + to_string(charge_dependent) + ", c=" + to_string(cutoff) + ")";
    }
    fl eval(const atom_base& a, const atom_base& b, fl r) const {
        fl q1 = a.charge;
        fl q2 = b.charge;

        VINA_CHECK(not_max(q1));
        VINA_CHECK(not_max(q2));

        sz t1 = a.ad;
        sz t2 = b.ad;

        fl solv1 = solvation_parameter(a);
        fl solv2 = solvation_parameter(b);

        fl volume1 = volume(a);
        fl volume2 = volume(b);

        fl my_solv = charge_dependent ? solvation_q : 0;

        fl tmp = ((solv1 + my_solv * std::abs(q1)) * volume2 +
                  (solv2 + my_solv * std::abs(q2)) * volume1) * std::exp(-sqr(r/(2*desolvation_sigma)));

        VINA_CHECK(not_max(tmp));
        return tmp;
    }
};
// forcefiled2 common functions
fl smoothen(fl r, fl rij, fl smoothing) {
    fl out;
    smoothing *= 0.5;

    if (r > rij + smoothing)
        out = r - smoothing;
    else if(r < rij - smoothing)
        out = r + smoothing;
    else
        out = rij;

    return out;
}
//ad4 score function
fl ad4_hb_eps(sz& a) {
    if (a < AD_TYPE_SIZE) return ad_type_property(a).hb_depth;
    VINA_CHECK(false);
    return 0; // placating the compiler
}

fl ad4_hb_radius(sz& t) {
    if (t < AD_TYPE_SIZE) return ad_type_property(t).hb_radius;
    VINA_CHECK(false);
    return 0; // placating the compiler
}
fl ad4_vdw_eps(sz& a) {
    if(a < AD_TYPE_SIZE) return ad_type_property(a).depth;
    VINA_CHECK(false);
    return 0; // placating the compiler
}

fl ad4_vdw_radius(sz& t) {
    if(t < AD_TYPE_SIZE) return ad_type_property(t).radius;
    VINA_CHECK(false);
    return 0; // placating the compiler
}
template<unsigned i, unsigned j>
struct forcefiled_vdw : public usable {
    fl smoothing;
    fl cap;
    forcefiled_vdw(fl smoothing_, fl cap_, fl cutoff_)
            : usable(cutoff_), smoothing(smoothing_), cap(cap_) {
        name = "vdw(i=" + to_string(i) + ", j=" + to_string(j) + ", s=" + to_string(smoothing) + ", ^=" + to_string(cap) + ", c=" + to_string(cutoff) + ")";
    }
    fl eval(const atom& a, const atom& b, fl r) const {
//        if (r >= cutoff)
//            return 0.0;
//        sz t1 = a.ad;
//        sz t2 = b.ad;
//        fl hb_depth = forcefiled_hb_eps(t1) * forcefiled_hb_eps(t2);
//        fl vdw_rij = forcefiled_vdw_radius(t1) + forcefiled_vdw_radius(t2);
//        fl vdw_depth = std::sqrt(forcefiled_vdw_eps(t1) * forcefiled_vdw_eps(t2));
//        if (hb_depth < 0) return 0.0; // interaction is hb, not vdw.
//        r = smoothen(r, vdw_rij, smoothing);
//        fl c_12 = int_pow<12>(vdw_rij) * vdw_depth;
//        fl c_6  = int_pow<6>(vdw_rij)  * vdw_depth * 2.0;
//        fl r6   = int_pow<6>(r);
//        fl r12  = int_pow<12>(r);
//        if(r12 > epsilon_fl && r6 > epsilon_fl)
//            return (std::min)(cap, c_12 / r12 - c_6 / r6);
//        else
//            return cap;
//        VINA_CHECK(false);
//        return 0.0; // placating the compiler
    }
};

template<unsigned i>
struct forcefiled_electrostatic : public usable {
    fl smoothing;
    fl cap;
    forcefiled_electrostatic(fl smoothing_, fl cap_, fl cutoff_)
            : usable(cutoff_), smoothing(smoothing_), cap(cap_) {
        name = std::string("electrostatic(i=") + to_string(i) + ", ^=" + to_string(cap) + ", c=" + to_string(cutoff) + ")";
    }
    fl eval(const atom& a, const atom& b, fl r) const {
        if (r >= cutoff)
            return 0.0;
        fl q1q2 = a.charge * b.charge * 332.0;
        fl B = 78.4 + 8.5525;
        fl lB = -B * 0.003627;
        fl diel = -8.5525 + (B / (1 + 7.7839 * std::exp(lB * r)));
        if (r < epsilon_fl)
            return q1q2 * cap / diel;
        else {
            return q1q2 * (std::min)(cap, 1.0 / (r * diel));
        }
    }
};
inline fl optimal_distance(sz xs_t1, sz xs_t2) {
    return xs_radius(xs_t1) + xs_radius(xs_t2);
}

struct gauss : public usable {
    fl offset; // added to optimal distance
    fl width;
    gauss(fl offset_, fl width_, fl cutoff_) : usable(cutoff_), offset(offset_), width(width_) {
        name = std::string("gauss(o=") + to_string(offset) + ", w=" + to_string(width) + ", c=" + to_string(cutoff) + ")";
    }
    fl eval(sz t1, sz t2, fl r) const {
        return gaussian(r - (optimal_distance(t1, t2) + offset), width);
    }
};

struct repulsion : public usable {
    fl offset; // added to vdw
    repulsion(fl offset_, fl cutoff_) : usable(cutoff_), offset(offset_) {
        name = std::string("repulsion(o=") + to_string(offset) + ")";
    }
    fl eval(sz t1, sz t2, fl r) const {
        fl d = r - (optimal_distance(t1, t2) + offset);
        if(d > 0)
            return 0;
        return d*d;
    }
};

inline fl slope_step(fl x_bad, fl x_good, fl x) {
    if(x_bad < x_good) {
        if(x <= x_bad) return 0;
        if(x >= x_good) return 1;
    }
    else {
        if(x >= x_bad) return 0;
        if(x <= x_good) return 1;
    }
    return (x - x_bad) / (x_good - x_bad);
}

struct hydrophobic : public usable {
    fl good;
    fl bad;
    hydrophobic(fl good_, fl bad_, fl cutoff_) : usable(cutoff_), good(good_), bad(bad_) {
        name = "hydrophobic(g=" + to_string(good) + ", b=" + to_string(bad) + ", c=" + to_string(cutoff) + ")";
    }
    fl eval(sz t1, sz t2, fl r) const {
        if(xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2))
            return slope_step(bad, good, r - optimal_distance(t1, t2));
        else return 0;
    }
};

struct non_hydrophobic : public usable {
    fl good;
    fl bad;
    non_hydrophobic(fl good_, fl bad_, fl cutoff_) : usable(cutoff_), good(good_), bad(bad_) {
        name = "non_hydrophobic(g=" + to_string(good) + ", b=" + to_string(bad) + ", c=" + to_string(cutoff) + ")";
    }
    fl eval(sz t1, sz t2, fl r) const {
        if(!xs_is_hydrophobic(t1) && !xs_is_hydrophobic(t2))
            return slope_step(bad, good, r - optimal_distance(t1, t2));
        else return 0;
    }
};

template<unsigned n, unsigned m>
void find_vdw_coefficients(fl position, fl depth, fl& c_n, fl& c_m) {
    BOOST_STATIC_ASSERT(n != m);
    c_n = int_pow<n>(position) * depth * m / (fl(n)-fl(m));
    c_m = int_pow<m>(position) * depth * n / (fl(m)-fl(n));
}


template<unsigned i, unsigned j>
struct vdw : public usable {
    fl smoothing;
    fl cap;
    vdw(fl smoothing_, fl cap_, fl cutoff_)
            : usable(cutoff_), smoothing(smoothing_), cap(cap_) {
        name = "vdw(i=" + to_string(i) + ", j=" + to_string(j) + ", s=" + to_string(smoothing) + ", ^=" + to_string(cap) + ", c=" + to_string(cutoff) + ")";
    }
    fl eval(sz t1, sz t2, fl r) const {
        fl d0 = optimal_distance(t1, t2);
        fl depth = 1;
        fl c_i = 0;
        fl c_j = 0;
        find_vdw_coefficients<i, j>(d0, depth, c_i, c_j);
        if     (r > d0 + smoothing) r -= smoothing;
        else if(r < d0 - smoothing) r += smoothing;
        else r = d0;

        fl r_i = int_pow<i>(r);
        fl r_j = int_pow<j>(r);
        if(r_i > epsilon_fl && r_j > epsilon_fl)
            return (std::min)(cap, c_i / r_i + c_j / r_j);
        else
            return cap;
    }
};

struct non_dir_h_bond : public usable {
    fl good;
    fl bad;
    non_dir_h_bond(fl good_, fl bad_, fl cutoff_) : usable(cutoff_), good(good_), bad(bad_) {
        name = std::string("non_dir_h_bond(g=") + to_string(good) + ", b=" + to_string(bad) + ")";
    }
    fl eval(sz t1, sz t2, fl r) const {
        if(xs_h_bond_possible(t1, t2))
            return slope_step(bad, good, r - optimal_distance(t1, t2));
        return 0;
    }
};

inline fl read_iterator(flv::const_iterator& i) {
    fl x = *i;
    ++i;
    return x;
}

fl smooth_div(fl x, fl y) {
    if(std::abs(x) < epsilon_fl) return 0;
    if(std::abs(y) < epsilon_fl) return ((x*y > 0) ? max_fl : -max_fl); // FIXME I hope -max_fl does not become NaN
    return x / y;
}

struct num_tors_add : public conf_independent {
    num_tors_add() { name = "num_tors_add"; }
    sz size() const { return 1; }
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
        //fl w = 0.1 * read_iterator(i); // [-1 .. 1]
        fl w = read_iterator(i); // FIXME?
        return x + w * in.num_tors;

    }
};

struct num_tors_sqr : public conf_independent {
    num_tors_sqr() { name = "num_tors_sqr"; }
    sz size() const { return 1; }
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
        fl w = 0.1 * read_iterator(i); // [-1 .. 1]
        return x + w * sqr(fl(in.num_tors)) / 5;
    }
};

struct num_tors_sqrt : public conf_independent {
    num_tors_sqrt() { name = "num_tors_sqrt"; }
    sz size() const { return 1; }
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
        fl w = 0.1 * read_iterator(i); // [-1 .. 1]
        return x + w * std::sqrt(fl(in.num_tors)) / sqrt(5.0);
    }
};

struct num_tors_div : public conf_independent {
    num_tors_div() { name = "num_tors_div"; }
    sz size() const { return 1; }
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
        fl w = 0.1 * (read_iterator(i) + 1); // w is in [0..0.2]
        return smooth_div(x, 1 + w * in.num_tors/5.0);
    }
};

struct ligand_length : public conf_independent {
    ligand_length() { name = "ligand_length"; }
    sz size() const { return 1; }
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
        fl w = read_iterator(i);
        return x + w * in.ligand_lengths_sum;
    }
};

struct num_ligands : public conf_independent {
    num_ligands() { name = "num_ligands"; }
    sz size() const { return 1; }
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
        fl w = 1 * read_iterator(i); // w is in [-1.. 1]
        return x + w * in.num_ligands;
    }
};

struct num_heavy_atoms_div : public conf_independent {
    num_heavy_atoms_div() { name = "num_heavy_atoms_div"; }
    sz size() const { return 1; }
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
        fl w = 0.05 * read_iterator(i);
        return smooth_div(x, 1 + w * in.num_heavy_atoms);
    }
};

struct num_heavy_atoms : public conf_independent {
    num_heavy_atoms() { name = "num_heavy_atoms"; }
    sz size() const { return 1; }
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
        fl w = 0.05 * read_iterator(i);
        return x + w * in.num_heavy_atoms;
    }
};

struct num_hydrophobic_atoms : public conf_independent {
    num_hydrophobic_atoms() { name = "num_hydrophobic_atoms"; }
    sz size() const { return 1; }
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
        fl w = 0.05 * read_iterator(i);
        return x + w * in.num_hydrophobic_atoms;
    }
};


//template<unsigned i, unsigned j>
//struct xscore_vdw : public usable {
//    fl smoothing;
//    fl cap;
//    xscore_vdw(fl smoothing_, fl cap_, fl cutoff_)
//            : usable(cutoff_), smoothing(smoothing_), cap(cap_) {
//        name = "xscore_vdw(i=" + to_string(i) + ", j=" + to_string(j) + ", s=" + to_string(smoothing) + ", ^=" + to_string(cap) + ", c=" + to_string(cutoff) + ")";
//    }
//    // now calculate the P-L vdw interaction
//    fl eval(const atom& a, const atom& b, fl r) const {
//
//        if (r >= cutoff)
//            return 0.0;
//        sz t1 = a.xs;
//        sz t2 = b.xs;
//        fl vdw_rij = forcefiled_vdw_radius(t1) + forcefiled_vdw_radius(t2);
//        fl distance_rij = sqrt(vec_distance_sqr(a.coords, b.coords));
//
//        // Lennard-Jones 8-4 potential
//
//        fl tmp1=vdw_rij/distance_rij;
//        tmp1=tmp1*tmp1*tmp1*tmp1; fl tmp2=tmp1*tmp1;
//        fl tmp=tmp2-2.00*tmp1;
//
//
//        return tmp; // placating the compiler
//    }
//};
//
//
//
//struct xscore_RT : public usable {
//    fl smoothing;
//    fl cap;
//    xscore_RT(fl smoothing_, fl cap_, fl cutoff_)
//            : usable(cutoff_), smoothing(smoothing_), cap(cap_) {
//    }
//    // now calculate the P-L vdw interaction
//    fl eval(model m,const atom& a, const atom& b, fl r) const {
//        fl tmp = 0;
//        if (r >= cutoff)
//            return 0.0;
//
////        sz t1 = a.ad;
////        sz t2 = b.ad;
////        fl vdw_rij = forcefiled_vdw_radius(t1) + forcefiled_vdw_radius(t2);
//        fl distance_rij = sqrt(vec_distance_sqr(a.coords, b.coords));
//        int mark = 0;
//        for(int i = 0;i<a.bonds.size();i++){
//            if(a.bonds[i].rotatable && m.atoms[a.bonds[i].connected_atom_index.i] == b){
//                mark++;
//            }
//        }
//
//        if(mark == 1) tmp += 0.50;
//        else if(mark == 2) tmp += 0.10;
//        else tmp += 0.50;
//
//        return tmp; // placating the compiler
//    }
//};
//
//
//struct xscore_HP : public usable {
//    fl smoothing;
//    fl cap;
//    xscore_HP(fl smoothing_, fl cap_, fl cutoff_)
//            : usable(cutoff_), smoothing(smoothing_), cap(cap_) {
//    }
//    // now calculate the P-L vdw interaction
//    fl eval(model m,const atom& a, const atom& b, fl r) const {
//        if (r >= cutoff)
//            return 0.0;
//
//        sz t1 = a.xs;
//        sz t2 = b.xs;
//        fl vdw_rij = forcefiled_vdw_radius(t1) + forcefiled_vdw_radius(t2);
//        fl distance_rij = sqrt(vec_distance_sqr(a.coords, b.coords));
//
//        fl d1 = vdw_rij + 0.50;
//        fl d2 = vdw_rij + 2.20;
//
//        fl tmp = 0.0;
//        if(distance_rij < d1) tmp = 1.0;
//        else if(distance_rij < d2) tmp=(1/(d1-d2))*(distance_rij-d2);
//        else tmp = 0.0;
//        return tmp;
//    }
//};
//
//
//
//
//
//
//
//class Dot
//{
//public:
//    int valid;		// status indicator
//    sz type;		// type
//    float coor[3];		// coordinates
//    float unit;		// contribution, in either A^2 or A^3
//    float score;		// score on this dot
//
//    Dot(); ~Dot();
//
//    void Clear();
//};
//
//class DotSet{
//public:
//    int num_dot;
//    std::vector <Dot> dot;
//    float r;
//    sz type;
//    float unit;		// default contribution of each dot to total
//    float total;		// total volume or surface
//
//    DotSet(); ~DotSet();
//
//    DotSet(const DotSet &original);
//    DotSet& operator = (const DotSet &original);
//
//    void Clear();
//    void Show_Contents() const;
//    void Show_Dots(char *filename, char *header, char *show="unit") const;
//};
//class ForceField{
//    int num_sdot_type;
//    int num_xatomtype;
//    DotSet *sdot;
//    std::vector <Dot> sur_dot;
//
//
//    void Read_SURFACE_DEF(char *filename)
//    {
//        FILE *fp;
//        int i,count;
//        char buf[256],head[256];
//        Dot tmp_dot;
//
//        if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);
//
//        // first, determine num_sdot_type;
//
//        count=0;
//
//        for(;;)
//        {
//            if(fgets(buf,256,fp)==NULL) break;
//            else if(buf[0]=='#') continue;
//            else if(Blank_Line_Check(buf)==true) continue;
//            else
//            {
//                sscanf(buf,"%s", head);
//                if(strcmp(head,"DOTSET")) continue;
//                else count++;
//            }
//        }
//
//        rewind(fp);
//
//        num_sdot_type=count;
//
//        sdot=new DotSet[num_sdot_type];
//        if(sdot==NULL) Memory_Allocation_Error();
//
//        // now read pre-calculated volume dots
//
//        count=0;
//
//        for(;;)
//        {
//            if(fgets(buf,256,fp)==NULL) break;
//            else if(Blank_Line_Check(buf)==true) continue;
//            else if(buf[0]=='#') continue;
//
//            sscanf(buf,"%s",head);
//
//            if(!strcmp(head,"END")) break;
//            else if(!strcmp(head,"DOTSET"))
//            {
//                sscanf(buf,"%*s%f%d%f%f",
//                       &sdot[count].r,
//                       &sdot[count].num_dot,
//                       &sdot[count].unit,
//                       &sdot[count].total);
//
////                strcpy(sdot[count].type,"Un");
//
//                sdot[count].dot.clear();
//
//                for(i=0;i<sdot[count].num_dot;i++)
//                {
//                    fgets(buf,256,fp);
//                    sscanf(buf,"%f%f%f",
//                           &tmp_dot.coor[0],
//                           &tmp_dot.coor[1],
//                           &tmp_dot.coor[2]);
//
//                    tmp_dot.valid=1;
//                    tmp_dot.unit=sdot[count].unit;
//                    strcpy(tmp_dot.type,"Un");
//
//                    sdot[count].dot.push_back(tmp_dot);
//                }
//
//                if(sdot[count].dot.size()!=sdot[count].num_dot)
//                    Read_File_Error(filename);
//
//                count++;
//            }
//        }
//
//        fclose(fp);
//
//        if(count!=num_sdot_type) Read_File_Error(filename);
//
//        return;
//    }
//
//    DotSet Get_Surface_Dot(float R, float x, float y, float z) const
//    {
//        int i,num;
//        bool mark;
//        float tmp,theta,phi,theta_step,phi_step;
//        float r,d,total;
//        DotSet tmp_set;
//        Dot tmp_dot;
//
//        // check the pre-calculated surface dot sets
//
//        mark=false;
//
//        for(i=0;i<num_sdot_type;i++)
//        {
//            if(fabs(R-sdot[i].r)>0.025) continue;
//            else {tmp_set=sdot[i]; mark=true; break;}
//        }
//
//        if(mark==true)
//        {
//            for(i=0;i<tmp_set.num_dot;i++)
//            {
//                tmp_set.dot[i].coor[0]+=x;
//                tmp_set.dot[i].coor[1]+=y;
//                tmp_set.dot[i].coor[2]+=z;
//            }
//            return tmp_set;
//        }
//
//        // if it is not pre-calculated, calculate it now
//
//        total=(4.000*pi*R*R);
//
//        // d=0.500;		// spacing between two dots
//
//        d=sqrt(total/300);
//        if(d<0.500) d=0.500;
//
//        tmp=(int)(pi*R/d+0.500); theta_step=pi/tmp;
//
//        num=0; tmp_set.dot.clear();
//
//        for(theta=0.00;theta<pi;theta+=theta_step)
//        {
//            r=R*sin(theta);
//            tmp=(int)(2*pi*r/d+0.500); phi_step=2*pi/tmp;
//
//            for(phi=0.00;phi<(2*pi);phi+=phi_step)
//            {
//                tmp_dot.coor[0]=r*cos(phi);
//                tmp_dot.coor[1]=r*sin(phi);
//                tmp_dot.coor[2]=R*cos(theta);
//                tmp_dot.valid=1;
//                strcpy(tmp_dot.type,"Un");
//                tmp_set.dot.push_back(tmp_dot);
//                num++;
//            }
//        }
//
////        strcpy(tmp_set.type, "Un");
//        tmp_set.r=R;
//        tmp_set.num_dot=num; tmp_set.unit=total/num;
//        tmp_set.total=total;
//
//        for(i=0;i<tmp_set.num_dot;i++)
//        {
//            tmp_set.dot[i].unit=tmp_set.unit;
//            tmp_set.dot[i].coor[0]+=x;
//            tmp_set.dot[i].coor[1]+=y;
//            tmp_set.dot[i].coor[2]+=z;
//        }
//
//        return tmp_set;
//    }
//
//
//    DotSet Get_Surface_Dot(const atom& atom,float r){
//        int i;
//        float R;
//        DotSet tmp_set;
//
//        sz t1 = atom.xs;
//        fl vdw_rij = forcefiled_vdw_radius(t1);
//
//        R = vdw_rij+r;
//
//        tmp_set=Get_Surface_Dot(R,atom.coords[0],atom.coords[1],atom.coords[2]);
//
//        tmp_set.type = atom.xs;
//
//        for(i=0;i<tmp_set.num_dot;i++)
//        {
//            tmp_set.dot[i].type = atom.xs;
//        }
//
//        return tmp_set;
//    }
//
//    float Atom_Buried_Surface(int id, float &total, float &buried,const atom& a, const atom& b){
//        total=buried=0.000;
//        int check;
//        float ratio;
//        sz t1 = a.xs;
//        sz t2 = b.xs;
//        fl vdw_rij = forcefiled_vdw_radius(t1) + forcefiled_vdw_radius(t2);
//        fl distance_rij = sqrt(vec_distance_sqr(a.coords, b.coords));
//
//        fl WATER_R = 1.40;
//        if(distance_rij > forcefiled_vdw_radius(t1) + forcefiled_vdw_radius(t2) + 2 * WATER_R){
//            return 0 ;
//        }else{
//            check = 1;
//        }
//
//        int num = sur_dot.size();
//        for(int j = 0; j < num; j++){
//            bool mark = false;
//
//            for(int k = 0; k < protein.num; k++){
//                if(atom_check_list[k]==0) continue;
//
//                float d = sqrt(vec_distance_sqr(sur_dot[j].coor, b.coords));
//
//                if(d>(protein->atom[k].R+WATER_R)) continue;
//                else {mark=true; break;}
//            }
//
//
//        total+=sur_dot[j].unit;
//
//        if(mark==false) continue;
//        else buried+=this->sur_dot[j].unit;
//    }
//
//    if(atom_check_list) delete [] atom_check_list;
//
//    if(total<=0.00) ratio=0.00;
//    else ratio=buried/total;
//
//    return ratio;
//
//
//    }
//
//
//
//    struct Xatom_Def
//    {
//        char type[20];
//        float logp;
//    };
//    Xatom_Def *xatom;
//
//    void Read_XATOM_DEF(char *filename)
//    {
//        FILE *fp;
//        int count;
//        char buf[256];
//
//        if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);
//
//        // first, determine num_xatomtype
//
//        count=0;
//
//        for(;;)
//        {
//            if(fgets(buf,256,fp)==NULL) break;
//            else if(buf[0]=='#') continue;
//            else if(Blank_Line_Check(buf)==true) continue;
//            else count++;
//        }
//
//        rewind(fp);
//
//        num_xatomtype=count;
//
//        xatom=new Xatom_Def[num_xatomtype];
////        if(xatom==NULL) Memory_Allocation_Error();
//
//        // second, read parameters for each atom type
//
//        count=0;
//
//        for(;;)
//        {
//            if(fgets(buf,256,fp)==NULL) break;
//            else if(buf[0]=='#') continue;
//            else if(Blank_Line_Check(buf)==true) continue;
//            else
//            {
//                sscanf(buf,"%*d%s%*s%f",
//                       xatom[count].type,
//                       &xatom[count].logp);
//                count++;
//            }
//        }
//
//        fclose(fp);
//
//        if(num_xatomtype!=count) Read_File_Error(filename);
//
//        return;
//    }
//
//
//
//    float Get_Atom_LogP(char *type) const
//    {
//        int i;
//        bool mark;
//        float logp;
//
//        mark=false;
//
//        for(i=0;i<num_xatomtype;i++)
//        {
//            if(strcmp(type,xatom[i].type)) continue;
//            else {logp=xatom[i].logp; mark=true; break;}
//        }
//
//        if(mark==true) return logp;
//        else
//        {
//            printf("Warning: no LogP parameter for atom type %s ... ", type);
//            printf("zero value assigned\n");
//            return 0.000;
//        }
//    }
//};	// full version
//
//
//
///*-------------------HM-----------------------*/
//
//fl Calculate_LogP(){
//    fl * logp_factor;
//
//    logp_factor[0] = 0.211;
//    logp_factor[1] = 0.429;
//    logp_factor[2] = 0.137;
//    logp_factor[3] = 0.485;
//    logp_factor[4] =-0.268;
//    logp_factor[5] = 0.580;
//    logp_factor[6] =-0.423;
//    logp_factor[7] =-2.166;
//    logp_factor[8] = 0.554;
//    logp_factor[9] =-0.501;
//
//    fl xlogp=0.000;
//
//    ForceField ff = new ForceField();
//    for(int i = 0; i < num_atom; i++){
//        atom[i].logp = ff.Get_Atom_LogP(type);
//        xlogp+=atom[i].logp;
//    }
//
//    for(int j = 0; j < 10; j++){
//        logp_factor[j].num = 0;
//    }
//
//    for(int i=0;i<10;i++) xlogp+=(logp_factor[i].num*logp_factor[i].coeff);
//
//    return xlogp;
//
//}
//struct xscore_HM : public usable {
//    fl smoothing;
//    fl cap;
//    xscore_HM(fl smoothing_, fl cap_, fl cutoff_)
//            : usable(cutoff_), smoothing(smoothing_), cap(cap_) {
//    }
//    // now calculate the P-L vdw interaction
//    fl eval(model m,const atom& a, const atom& b, fl r) const {
//        if (r >= cutoff)
//            return 0.0;
//
//        sz t1 = a.xs;
//        sz t2 = b.xs;
//        fl vdw_rij = forcefiled_vdw_radius(t1) + forcefiled_vdw_radius(t2);
//        fl distance_rij = sqrt(vec_distance_sqr(a.coords, b.coords));
//
//
//
//
//    }
//};
//
//struct xscore_HS : public usable {
//    fl smoothing;
//    fl cap;
//    float sum,total,buried;
//
//    xscore_HS(fl smoothing_, fl cap_, fl cutoff_)
//            : usable(cutoff_), smoothing(smoothing_), cap(cap_) {
//    }
//    // now calculate the P-L vdw interaction
//    fl eval(model m,const atom& a, const atom& b, fl r) const {
//        if (r >= cutoff)
//            return 0.0;
//
//        sz t1 = a.xs;
//        sz t2 = b.xs;
//        fl vdw_rij = forcefiled_vdw_radius(t1) + forcefiled_vdw_radius(t2);
//
//        Atom_Buried_Surface(i,total,buried);
//        fl sum = 0.0;
//
//        sum+=buried; this->atom[i].score=buried;
//
//        return sum;
//
//
//    }
//};

//这个函数是设定vina的term项
everything::everything() { // enabled according to design.out227
    const unsigned d = 0; // default
    const fl cutoff = 8; //6;

    // FIXME? enable some?
    //// distance_additive
    //add(d, new forcefiled_solvation(3.6, 0.01097,  true, cutoff)); // desolvation_sigma, solvation_q, charge_dependent, cutoff
    //add(d, new forcefiled_solvation(3.6, 0.01097, false, cutoff)); // desolvation_sigma, solvation_q, charge_dependent, cutoff

    //add(d, new electrostatic<1>(100, cutoff)); // cap, cutoff
    //add(d, new electrostatic<2>(100, cutoff)); // cap, cutoff

    //add(d, new gauss(0,   0.3, cutoff)); // offset, width, cutoff
    //add(d, new gauss(0.5, 0.3, cutoff)); // offset, width, cutoff
    //add(d, new gauss(1,   0.3, cutoff)); // offset, width, cutoff
    //add(d, new gauss(1.5, 0.3, cutoff)); // offset, width, cutoff
    //add(d, new gauss(2,   0.3, cutoff)); // offset, width, cutoff
    //add(d, new gauss(2.5, 0.3, cutoff)); // offset, width, cutoff

    add(1, new gauss(0, 0.5, cutoff)); // offset, width, cutoff // WEIGHT: -0.035579
    //add(d, new gauss(1, 0.5, cutoff)); // offset, width, cutoff
    //add(d, new gauss(2, 0.5, cutoff)); // offset, width, cutoff

    //add(d, new gauss(0, 0.7, cutoff)); // offset, width, cutoff
    //add(d, new gauss(1, 0.7, cutoff)); // offset, width, cutoff
    //add(d, new gauss(2, 0.7, cutoff)); // offset, width, cutoff

    //add(d, new gauss(0, 0.9, cutoff)); // offset, width, cutoff
    //add(d, new gauss(1, 0.9, cutoff)); // offset, width, cutoff
    //add(d, new gauss(2, 0.9, cutoff)); // offset, width, cutoff
    //add(d, new gauss(3, 0.9, cutoff)); // offset, width, cutoff

    //add(d, new gauss(0, 1.5, cutoff)); // offset, width, cutoff
    //add(d, new gauss(1, 1.5, cutoff)); // offset, width, cutoff
    //add(d, new gauss(2, 1.5, cutoff)); // offset, width, cutoff
    //add(d, new gauss(3, 1.5, cutoff)); // offset, width, cutoff
    //add(d, new gauss(4, 1.5, cutoff)); // offset, width, cutoff

    //add(d, new gauss(0, 2.0, cutoff)); // offset, width, cutoff
    //add(d, new gauss(1, 2.0, cutoff)); // offset, width, cutoff
    //add(d, new gauss(2, 2.0, cutoff)); // offset, width, cutoff
    add(1, new gauss(3, 2.0, cutoff)); // offset, width, cutoff // WEIGHT: -0.005156
    //add(d, new gauss(4, 2.0, cutoff)); // offset, width, cutoff

    //add(d, new gauss(0, 3.0, cutoff)); // offset, width, cutoff
    //add(d, new gauss(1, 3.0, cutoff)); // offset, width, cutoff
    //add(d, new gauss(2, 3.0, cutoff)); // offset, width, cutoff
    //add(d, new gauss(3, 3.0, cutoff)); // offset, width, cutoff
    //add(d, new gauss(4, 3.0, cutoff)); // offset, width, cutoff

    //add(d, new repulsion( 0.4, cutoff)); // offset, cutoff
    //add(d, new repulsion( 0.2, cutoff)); // offset, cutoff
    add(1, new repulsion( 0.0, cutoff)); // offset, cutoff // WEIGHT:  0.840245
    //add(d, new repulsion(-0.2, cutoff)); // offset, cutoff
    //add(d, new repulsion(-0.4, cutoff)); // offset, cutoff
    //add(d, new repulsion(-0.6, cutoff)); // offset, cutoff
    //add(d, new repulsion(-0.8, cutoff)); // offset, cutoff
    //add(d, new repulsion(-1.0, cutoff)); // offset, cutoff

    //add(d, new hydrophobic(0.5, 1, cutoff)); // good, bad, cutoff
    add(1, new hydrophobic(0.5, 1.5, cutoff)); // good, bad, cutoff // WEIGHT:  -0.035069
    //add(d, new hydrophobic(0.5, 2, cutoff)); // good, bad, cutoff
    //add(d, new hydrophobic(0.5, 3, cutoff)); // good, bad, cutoff

    //add(1, new non_hydrophobic(0.5, 1.5, cutoff));

    //add(d, new vdw<4,  8>(   0, 100, cutoff)); // smoothing, cap, cutoff

    add(1, new non_dir_h_bond(-0.7, 0, cutoff)); // good, bad, cutoff // WEIGHT:  -0.587439
    //add(d, new non_dir_h_bond(-0.7, 0, cutoff)); // good, bad, cutoff
    //add(d, new non_dir_h_bond(-0.7, 0.2, cutoff)); // good, bad, cutoff
    //add(d, new non_dir_h_bond(-0.7, 0.4, cutoff)); // good, bad, cutoff
    // additive

    // conf-independent
    //add(d, new num_ligands());

    add(1, new num_tors_div()); // WEIGHT: 1.923 -- FIXME too close to limit?
    //add(d, new num_heavy_atoms_div());
    //add(d, new num_heavy_atoms());
    //add(1, new num_tors_add());
    //add(d, new num_tors_sqr());
    //add(d, new num_tors_sqrt());
    //add(d, new num_hydrophobic_atoms());
    ///add(1, new ligand_length());

    //add(d, new num_tors(100, 100, false)); // cap, past_cap, heavy_only
    //add(1, new num_tors(100, 100,  true)); // cap, past_cap, heavy_only
    //add(d, new num_tors(  2,   1,  true)); // cap, past_cap, heavy_only
    //add(d, new num_heavy_atoms());
    //add(d, new ligand_max_num_h_bonds());
    //add(1, new num_ligands());
}


//ad4
//template<unsigned i, unsigned j>
struct ad4_vdw : public usable {
    fl smoothing;
    fl cap;
    ad4_vdw(fl smoothing_, fl cap_, fl cutoff_)
            : usable(cutoff_), smoothing(smoothing_), cap(cap_) {
        atom_typing_used = atom_type::AD;
//        name = "vdw(i=" + to_string(i) + ", j=" + to_string(j) + ", s=" + to_string(smoothing) + ", ^=" + to_string(cap) + ", c=" + to_string(cutoff) + ")";
    }
    fl eval(const atom& a, const atom& b, fl r) const {
        if (r >= cutoff)
            return 0.0;
        sz t1 = a.ad;
        sz t2 = b.ad;
        fl hb_depth = ad4_hb_eps(t1) * ad4_hb_eps(t2);
        fl vdw_rij = ad4_vdw_radius(t1) + ad4_vdw_radius(t2);
        fl vdw_depth = std::sqrt(ad4_vdw_eps(t1) * ad4_vdw_eps(t2));
        if (hb_depth < 0) return 0.0; // interaction is hb, not vdw.
        r = smoothen(r, vdw_rij, smoothing);
        fl c_12 = int_pow<12>(vdw_rij) * vdw_depth;
        fl c_6  = int_pow<6>(vdw_rij)  * vdw_depth * 2.0;
        fl r6   = int_pow<6>(r);
        fl r12  = int_pow<12>(r);
        if(r12 > epsilon_fl && r6 > epsilon_fl)
            return (std::min)(cap, c_12 / r12 - c_6 / r6);
        else
            return cap;
        VINA_CHECK(false);
        return 0.0; // placating the compiler
    }
};

//template<unsigned i, unsigned j>
struct ad4_hb : public usable {
    fl smoothing;
    fl cap;
    ad4_hb(fl smoothing_, fl cap_, fl cutoff_)
            : usable(cutoff_), smoothing(smoothing_), cap(cap_) {
        atom_typing_used = atom_type::AD;
//        name = "hb(i=" + to_string(i) + ", j=" + to_string(j) + ", s=" + to_string(smoothing) + ", ^=" + to_string(cap) + ", c=" + to_string(cutoff) + ")";
    }
    fl eval(const atom& a, const atom& b, fl r) const {
        if (r >= cutoff)
            return 0.0;
        sz t1 = a.ad;
        sz t2 = b.ad;
        fl hb_rij = ad4_hb_radius(t1) + ad4_hb_radius(t2);
        fl hb_depth = ad4_hb_eps(t1) * ad4_hb_eps(t2);
        fl vdw_rij = ad4_vdw_radius(t1) + ad4_vdw_radius(t2);
        if (hb_depth >= 0)
            return 0.0; // interaction is vdw, not hb.
        r = smoothen(r, hb_rij, smoothing);
        fl c_12 = int_pow<12>(hb_rij) * -hb_depth * 10 / 2.0;
        fl c_10 = int_pow<10>(hb_rij) * -hb_depth * 12 / 2.0;
        fl r10  = int_pow<10>(r);
        fl r12  = int_pow<12>(r);
        if (r12 > epsilon_fl && r10 > epsilon_fl)
            return (std::min)(cap, c_12 / r12 - c_10 / r10);
        else
            return cap;
        VINA_CHECK(false);
        return 0.0; // placating the compiler
    }
};

struct ad4_electrostatic : public usable {
    fl smoothing;
    fl cap;
    ad4_electrostatic( fl cap_, fl cutoff_)
            : usable(cutoff_), cap(cap_) {
        atom_typing_used = atom_type::AD;
    }
    fl eval(const atom& a, const atom& b, fl r) const {
        if (r >= cutoff)
            return 0.0;
        fl q1q2 = a.charge * b.charge * 332.0;
        fl B = 78.4 + 8.5525;
        fl lB = -B * 0.003627;
        fl diel = -8.5525 + (B / (1 + 7.7839 * std::exp(lB * r)));
        if (r < epsilon_fl)
            return q1q2 * cap / diel;
        else {
            return q1q2 * (std::min)(cap, 1.0 / (r * diel));
        }
    }
};

struct ad4_solvation : public usable {
    fl desolvation_sigma;
    fl solvation_q;
    bool charge_dependent;
    ad4_solvation(fl desolvation_sigma_, fl solvation_q_, bool charge_dependent_, fl cutoff_) : usable(cutoff_), solvation_q(solvation_q_), charge_dependent(charge_dependent_), desolvation_sigma(desolvation_sigma_) {
        name = std::string("forcefiled_solvation(d-sigma=") + to_string(desolvation_sigma) + ", s/q=" + to_string(solvation_q) + ", q=" + to_string(charge_dependent) + ", c=" + to_string(cutoff) + ")";
        atom_typing_used = atom_type::AD;
    }
    fl eval(const atom& a, const atom& b, fl r) const {
        if (r >= cutoff)
            return 0.0;
        fl q1 = a.charge;
        fl q2 = b.charge;
        VINA_CHECK(not_max(q1));
        VINA_CHECK(not_max(q2));
        fl solv1 = solvation_parameter(a);
        fl solv2 = solvation_parameter(b);
        fl volume1 = volume(a);
        fl volume2 = volume(b);
        fl my_solv = charge_dependent ? solvation_q : 0;
        fl tmp = ((solv1 + my_solv * std::abs(q1)) * volume2 +
                  (solv2 + my_solv * std::abs(q2)) * volume1) * std::exp(-0.5 * sqr(r / desolvation_sigma));
        VINA_CHECK(not_max(tmp));
        return tmp;
    }
};

// Macrocycle - Vina and AD42
inline bool is_glued(sz xs_t1, sz xs_t2) {
    return (xs_t1 == XS_TYPE_G0 && xs_t2 == XS_TYPE_C_H_CG0) ||
           (xs_t1 == XS_TYPE_G0 && xs_t2 == XS_TYPE_C_P_CG0) ||
           (xs_t2 == XS_TYPE_G0 && xs_t1 == XS_TYPE_C_H_CG0) ||
           (xs_t2 == XS_TYPE_G0 && xs_t1 == XS_TYPE_C_P_CG0) ||

           (xs_t1 == XS_TYPE_G1 && xs_t2 == XS_TYPE_C_H_CG1) ||
           (xs_t1 == XS_TYPE_G1 && xs_t2 == XS_TYPE_C_P_CG1) ||
           (xs_t2 == XS_TYPE_G1 && xs_t1 == XS_TYPE_C_H_CG1) ||
           (xs_t2 == XS_TYPE_G1 && xs_t1 == XS_TYPE_C_P_CG1) ||

           (xs_t1 == XS_TYPE_G2 && xs_t2 == XS_TYPE_C_H_CG2) ||
           (xs_t1 == XS_TYPE_G2 && xs_t2 == XS_TYPE_C_P_CG2) ||
           (xs_t2 == XS_TYPE_G2 && xs_t1 == XS_TYPE_C_H_CG2) ||
           (xs_t2 == XS_TYPE_G2 && xs_t1 == XS_TYPE_C_P_CG2) ||

           (xs_t1 == XS_TYPE_G3 && xs_t2 == XS_TYPE_C_H_CG3) ||
           (xs_t1 == XS_TYPE_G3 && xs_t2 == XS_TYPE_C_P_CG3) ||
           (xs_t2 == XS_TYPE_G3 && xs_t1 == XS_TYPE_C_H_CG3) ||
           (xs_t2 == XS_TYPE_G3 && xs_t1 == XS_TYPE_C_P_CG3);
}

// Macrocycle - Vina and AD42
struct linearattraction : public usable {
    fl smoothing;
    fl cap;
    linearattraction(fl cutoff_)
            : usable(cutoff_){
        atom_typing_used = atom_type::AD;
    }
    fl eval(const atom& a, const atom& b, fl r) const {
        if (r >= cutoff)
            return 0.0;
        if (is_glued(a.xs, b.xs))
            return r;
        else
            return 0.0;
    }

    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff)
            return 0.0;
        if (is_glued(t1, t2))
            return r;
        else
            return 0.0;
    }
};

struct ad4_tors_add : public conf_independent {
    ad4_tors_add() { name = "num_tors_add"; }
    sz size() const { return 1; }
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
        //fl w = 0.1 * read_iterator(i); // [-1 .. 1]
        fl weight = read_iterator(i);
        return x + weight * in.torsdof;

    }
};

everything::everything(const std::string term_typename) {
    const unsigned d = 0; // default
    const fl cutoff = 8; //6;

    if((term_typename.compare("forcefield")) == 0){//废案
        add(1,new forcefiled_vdw<12,  6>(0, 100, cutoff));
        add(1, new forcefiled_electrostatic<2>(0,100, cutoff));
    }else if((term_typename.compare("xscore")) == 0){//废案
        add(1, new gauss(0, 0.5, cutoff)); // offset, width, cutoff // WEIGHT: -0.035579
        add(1, new gauss(3, 2.0, cutoff)); // offset, width, cutoff // WEIGHT: -0.005156
        add(1, new repulsion( 0.0, cutoff)); // offset, cutoff // WEIGHT:  0.840245
        add(1, new hydrophobic(0.5, 1.5, cutoff)); // good, bad, cutoff // WEIGHT:  -0.035069
        add(1, new non_dir_h_bond(-0.7, 0, cutoff)); // good, bad, cutoff // WEIGHT:  -0.587439
        add(1, new num_tors_div()); // WEIGHT: 1.923 -- FIXME too close to limit?

    }else if((term_typename.compare("ad4")) == 0){//废案
        add(1, new ad4_vdw(0.5, 100000, 8.0)); // offset, width, cutoff // WEIGHT: -0.035579
        add(1, new ad4_hb(0.5, 100000, 8.0)); // offset, width, cutoff // WEIGHT: -0.005156
        add(1, new ad4_electrostatic(100, 20.48)); // offset, cutoff // WEIGHT:  0.840245
        add(1, new ad4_solvation(3.6, 0.01097, true, 20.48)); // good, bad, cutoff // WEIGHT:  -0.035069
        add(1, new linearattraction(20.0)); // good, bad, cutoff // WEIGHT:  -0.587439
        add(1, new ad4_tors_add()); // WEIGHT: 1.923 -- FIXME too close to limit?
    }else if((term_typename.compare("vdw")) == 0){//主要是得到vdw的term项
        add(1, new ad4_vdw(0.5, 100000, 8.0));
    }
}