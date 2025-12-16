#include "protein.h"
#include "dfire.h"
#include "array1.h"

class Molecule_DLIGAND2:public Protein_DLIGAND2{
public:
    int lig_num = 0;
    Molecule_DLIGAND2(){}
    Molecule_DLIGAND2(string fn):Protein_DLIGAND2(fn){};
	void rdmol2(string fn);
	void stats();
	double score();
};
