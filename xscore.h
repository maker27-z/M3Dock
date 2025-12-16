////
//// Created by 91686 on 2023/8/14.
////
//
//#ifndef LSHADE_ADAM_FINAL_XSCORE_H
//#define LSHADE_ADAM_FINAL_XSCORE_H
//#include "int_pow.h"
//#include "model.h"
//# define TRUE 1
//# define FALSE 0
//# define TINY 1.0e-6
//# define LARGE 1.0e+6
//# define POCKET_DEPTH 4.00
//# define LAYER_DEPTH 3.00
//# define MAX_ATOM_NEIB 6
//# define MAX_BOND_NEIB 10
//# define SWAP(a,b) temp=(a); (a)=(b); (b)=temp
//# define MAX(a,b) (a)>=(b) ? (a):(b)
//# define MIN(a,b) (a)<=(b) ? (a):(b)
//# define SQR(a) (a)==0.0 ? 0.0:(a)*(a)
//# define SHIFT(a,b,c,d) (a)=(b); (b)=(c); (c)=(d)
//# define SIGN(a,b) ((b)>=0.0 ? fabs(a):-fabs(a))
//const float LOGP_HYDROPHOBIC_CARBON = 0.211;
//const float LOGP_INTERNAL_HBOND = 0.429;
//const float LOGP_HALOGEN_PAIR = 0.137;
//const float LOGP_NAR_PAIR = 0.485;
//const float LOGP_O3_PAIR = -0.268;
//const float LOGP_ACCEPTOR_PAIR = 0.580;
//const float LOGP_AMINO_ACID = -2.166;
//const float LOGP_SALICYLIC_ACID = 0.554;
//const float LOGP_SULFONIC_ACID = -0.501;
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
//
////HM
//struct LogP_Factor
//{
//    char symbol[30];
//    float num;
//    float coeff;
//} logp_factor[10];
//
//struct Xatom_Def
//{
//    char type[20];
//    float logp;
//};
//
////HB
//class HBond
//{
//public:
//    int valid;
//    int type;	// 1: donor=latom, acceptor=patom
//    // 2: donor=patom, acceptor=latom
//    // 3: donor=metal, acceptor=latom/patom
//    // 4: donor=latom, acceptor=latom
//    // 5: donor=patom, acceptor=patom
//    // 0: invalid, no H-bond
//    int d_type;	// donor type: 1=straight; 2=angled
//    int a_type;	// acceptor type: 1=straight; 2=angled
//    bool sb;	// flag for neutral HB or SB
//
//    int latom;	// id of the ligand atom
//    int patom;	// id of the protein atom
//
//    atom H,D,A;
//
//    float d;	// D-A distance
//    float a0;	// D-H-A angle
//    float a1;	// DR-D-A angle
//    float a2;	// D-A-AR angle
//
//    float score;	// strength of this HBond
//
//    HBond(); ~HBond();
//
//    void Clear();
//    void Show_Contents() const;
//
//    int Value_HBond();
//    void Determine_DA_Type();
//    int Value_HBond_2();
//    int Value_SBond();
//};
//
//class Group	// full version
//{
//public:
//    short int valid;
//    short int num_neib;    	// number of neighboring atoms
//    short int num_nonh;    	// number of neighboring non-h atoms
//    short int num_h;       	// number of neighboring hydrogen atoms
//    short int num_hetero;  	// number of neighboring heteroatoms
//    short int num_pi;      	// number of neighboring pi atoms
//    short int num_car;     	// number of neighboring C.ar
//    short int num_nar;     	// number of neighboring N.ar
//    short int num_db;	// number of double bonds
//    short int num_tb;	// number of triple bonds
//    short int db_type;      // double bond type
//    short int amide;	// indicator for amide group
//
//    atom center;            // center atom of the group
//    atom neib[6];   // neighboring atoms
//    bond bond[6];  /////////////////////////// // bonds
//
//    Group(); ~Group();
//
//    void Clear();
//    void Show_Contents() const;
//};
//
//
//
//
//struct xscore{
////    ForceField ff;
//    model m;
//    float cutoff;
//    atomv ligand_atoms;
//    atomv protein_atoms;
//    float WATER_R;
//
//    //HS
//    std::vector <Dot> sur_dot;
//    int num_sdot_type;
//    DotSet *sdot;
//
//    //HM
//    LogP_Factor logp_factor[10];
//    int num_xatomtype;
//    Xatom_Def *xatom;
//
//    xscore(const model& m_):m(m_){
//        cutoff = 8.0;
//        WATER_R = 1.40;
//        ligand_atoms = m.get_ligand_atoms();
//        protein_atoms = m.get_protein_atoms();
//
//        Read_XATOM_DEF("ATOM_DEF_XLOGP");
//        Read_SURFACE_DEF("SURFACE_DEF_XTOOL");
//        Generate_Surface_Dots(WATER_R);
//    }
//
//    float xscore_vdw();
//
//    float xscore_RT();
//
//    float calculate_HP();
//
//    float calculate_HM();
//
//    //HS
//    float calculate_HS();
//    float Atom_Buried_Surface(int id, float &total, float &buried);
//    void Generate_Surface_Dots(float probe_r);
//    DotSet Get_Surface_Dot(const atom& atom,float r);
//    DotSet Get_Surface_Dot(float R, float x, float y, float z);
//    void Read_SURFACE_DEF(char *filename);
//
//
//    //HM
//    float Calculate_LogP();
//    float Calculate_LogP_protein();
//
//    float Get_Atom_LogP(char *type);
//    void Read_XATOM_DEF(char *filename);
//    int Get_XLOGP_Type(int atom_id, char *type);
//    Group Find_X_Group(int atom_id);
//
//    int Connection_1_2_Check(int id1, int id2);
//    int Connection_1_3_Check(int id1, int id2);
//    int Connection_1_4_Check(int id1, int id2);
//    int Connection_1_5_Check(int id1, int id2);
//    int Connection_1_6_Check(int id1, int id2);
//
//    float Count_Hydrophobic_Carbon(int flag);
//    int Hydrophobic_Neighbor_Check(int id);
//
//    int Adjacent_Ring_Check(int id);
//    float Count_Internal_HBond(int flag);
//
//
//    float Count_Halogen_1_3_Pair(int flag);
//
//    float Count_Nar_1_4_Pair(int flag);
//
//    float Count_O3_1_4_Pair(int flag);
//
//    float Count_Acceptor_1_5_Pair(int flag);
//
//    // float Count_Remote_Donor_Pair(int flag);
//    float Count_Amino_Acid(int flag);
//
//    int Adjacent_Aromatic_Check(int id);
//    float Count_Salicylic_Acid(int flag);
//
//    float Count_Sulfonic_Acid(int flag);
//
//    //HB
//    float calculate_HB();
//    int Get_HBond_Pair_PL(HBond candidate[],bool sb_flag=false);
//    void Sum_HBonds(int num_candidate, HBond candidate[]);
//    float Angle_of_Two_Vectors(const float v1[3], const float v2[3]);
//    Group Find_A_Group(int atom_id);
//};
//
//
//
//#endif //LSHADE_ADAM_FINAL_XSCORE_H
