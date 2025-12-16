/*
atom_base
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

#ifndef VINA_ATOM_BASE_H
#define VINA_ATOM_BASE_H

# define MAX_ATOM_NEIB 6
# define MAX_BOND_NEIB 10

#include "atom_type.h"

struct atom_base : public atom_type {

    //xscore
    int id;  	        // atom id
    short int valid; 	// valid indicator
    char name[10];          // atom name
    char type[10];          // basic atom type
    char xtype[10];         // xtool atom type
    char type2[10];         // another atom type when necessary
    char residue[10];       // residue name
    char res_id[10];       	// residue id: must be a string for PDB
    char chain;             // chain label
    float coor[3];          // coordinates
    float root[3];          // HB root's coordinates
    float weight;       	// atomic weight
    float r;                // vdw radius
    float eps;              // vdw epsilon value
    float q;		// partial atomic charge
    float R;		// X radius
    float logp;             // atomic hydrophobic scale
    float solv;		// atomic solvation parameter
    char hb[3];             // HB property
    float occupancy;	// occupancy probability, for protein atoms
    float bfactor;		// B-factor, for protein atoms
    float score;            // atomic binding score
    short int ring;         // ring indicator: 1=normal; 2=aromatic
    short int origin;	// atom origin indicator: 1=ligand; 2=protein
    short int part;         // component indicator:
    // for protein atoms: 1=ATOM; 2=HETATM
    short int num_neib;     // number of neighboring atoms
    short int num_nonh;     // number of non-H neighboring atoms
    int neib[MAX_ATOM_NEIB];  // ID of neighboring atoms
    int bond[MAX_ATOM_NEIB];  // ID of neighboring bonds
    short int temp;		// for misc uses

	fl charge;
	atom_base() : charge(0) {}
private:
	friend class boost::serialization::access;
	template<class Archive> 
	void serialize(Archive& ar, const unsigned version) {
		ar & boost::serialization::base_object<atom_type>(*this);
		ar & charge;
	}
};

#endif
