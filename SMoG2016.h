////
//// Created by 91686 on 2023/8/31.
////
//#include <openbabel/mol.h>
//#include <openbabel/obconversion.h>
//#include <openbabel/forcefield.h>
//#include <openbabel/atom.h>
//#include <iomanip>
//#ifndef KOTO_NSGA2_MPSOD_XSCORE_SMOG2016_H
//#define KOTO_NSGA2_MPSOD_XSCORE_SMOG2016_H
//
//using namespace std;
//struct Parameters{
//    string       energyFile;       // To store the file with the knowledge-based potential.
//};
//struct Atom{
//    double       coordinates[3];   // Coordinates of the atom.
//    string       Type;             // Atom type for the scoring function when it is assigned as a string.
//    int          TypeNumber;       // Atom type for the scoring function when it is assigned as a number.
//    float        LJr;              // r0 value for the Lennard-Jones potential.
//    float        LJe;              // Epsilon value for the Lennard-Jones potential.
//    int          index;            // Atom number in the protein (we need it for StericClash).
//};
//struct OBMolecule{
//    OpenBabel::OBMol     obMol;    // Molecule from OpenBabel.
//    vector<Atom>         obatom;
////    OBMolecule& operator=(const OBMolecule& tmp){
////        obMol = tmp.obMol;
////        obatom = tmp.obatom;
////    }// To store the potential types of the atoms and their coordinates.
//};
//struct OBProtein{
//    OpenBabel::OBMol obMol;        // Molecule from OpenBabel.
//    vector<Atom>     obatom;         // To store the potential types of the protein atoms and their coordinates.
//};
//
//double dist(double const atom1[3], double const atom2[3]);
//bool   IsThereABond(OBMolecule const & molecule, unsigned int const i, unsigned int const j);
//double KBP2016(Parameters const & parameters, OBProtein const & protein, OBMolecule & molecule);
//double LJP(Parameters const & parameters, OBProtein const & protein, OBMolecule & molecule);
//double RotorXScore(OpenBabel::OBMol obMol);
//void   LigandTypeSMOG2016(OBMolecule & molecule);
//void   ProteinTypeSMOG2016(string const residue, string const atomID, int & atomTypeNumber, float & LJr, float & LJe);
//#endif //KOTO_NSGA2_MPSOD_XSCORE_SMOG2016_H
