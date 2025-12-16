/*
atom
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

#ifndef VINA_ATOM_H
#define VINA_ATOM_H

#include "atom_base.h"

struct atom_index {
	sz i;
	bool in_grid;
	atom_index() : i(max_sz), in_grid(false) {}
	atom_index(sz i_, bool in_grid_) : i(i_), in_grid(in_grid_) {}
private:
	friend class boost::serialization::access;
	template<class Archive> 
	void serialize(Archive& ar, const unsigned version) {
		ar & i;
		ar & in_grid;
	}
};

inline bool operator==(const atom_index& i, const atom_index& j) {
	return i.i == j.i && i.in_grid == j.in_grid;
}

struct bond {
	atom_index connected_atom_index;
	fl length;
	bool rotatable;
	bond() : length(0), rotatable(false) {}
	bond(const atom_index& connected_atom_index_, fl length_, bool rotatable_) : connected_atom_index(connected_atom_index_), length(length_), rotatable(rotatable_) {}
private:
	friend class boost::serialization::access;
	template<class Archive> 
	void serialize(Archive& ar, const unsigned version) {
		ar & connected_atom_index;
		ar & length;
		ar & rotatable;
	}
};

struct atom : public atom_base {
    vec coords;
	std::vector<struct bond> bonds;
    char ad_name[10];
    unsigned ad_number;
	atom()  : coords(max_vec){Clear();}

	~atom() {}

	void Clear()
	{
		id=0; valid=0; part=1; origin=1;

		strcpy(name,"Un"); strcpy(type,"Un");
		strcpy(xtype,"Un"); strcpy(type2,"Un");
		strcpy(residue,"Un"); strcpy(res_id,"0"); chain=' ';

		coor[0]=coor[1]=coor[2]=0.000;
		root[0]=root[1]=root[2]=0.000;

		r=0.000; eps=0.000; q=0.000; weight=0.000; R=0.000;
		strcpy(hb,"N"); logp=0.000; solv=0.000; score=0.000;
		ring=0;

		occupancy=1.000; bfactor=0.000;

		num_neib=num_nonh=0;

		int i;

		for(i=0;i<MAX_ATOM_NEIB;i++) neib[i]=0;
		for(i=0;i<MAX_ATOM_NEIB;i++) bond[i]=0;

		temp=0;
	}

	void Show_Contents() const
	{
		printf("Atom: ");
		printf("%-4d ",id);
		printf("%-5s ",name);
		printf("%8.3f ",coor[0]);
		printf("%8.3f ",coor[1]);
		printf("%8.3f ",coor[2]);
		printf("%-5s ",type);
		printf("\n");

		return;
	}

//	atom() : coords(max_vec) {}
    bool operator==(atom b){
        if(this->coords == b.coords){
            return true;
        }
    }
    int Get_Donor_Type() const;
    int Get_Acceptor_Type() const;
private:
	friend class boost::serialization::access;
	template<class Archive> 
	void serialize(Archive& ar, const unsigned version) {
		ar & boost::serialization::base_object<atom_base>(*this);
		ar & coords;
		ar & bonds;
	}
};

typedef std::vector<atom> atomv;

#endif
