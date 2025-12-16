//
// Created by 91686 on 2023/8/17.
//


#include <sys/time.h>
#include "xtools.h"

Ligand::Ligand()
{
    Clear();
}

Ligand::~Ligand()
{
//    if(abs_inf) delete [] abs_inf;
}

Ligand::Ligand(int max_atom_num, int max_bond_num)
{
    Clear();

    atom=new struct atom[max_atom_num];
    if(atom==NULL) Memory_Allocation_Error();

    bond=new Bond[max_bond_num];
    if(bond==NULL) Memory_Allocation_Error();
}

Ligand::Ligand(const Ligand &original)
{
    // this part is from Molecule

    xtool_format=original.xtool_format;
    id=original.id;
    valid=original.valid;
    strcpy(name,original.name);
    weight=original.weight;
    strcpy(formula,original.formula);
    logp=original.logp;
    surface=original.surface;
    nsur=original.nsur;
    psur=original.psur;
    volume=original.volume;
    num_hb_atom=original.num_hb_atom;
    num_rotor=original.num_rotor;

    num_subst=original.num_subst;
    num_feature=original.num_feature;
    num_set=original.num_set;
    strcpy(mol_type,original.mol_type);
    strcpy(charge_type,original.charge_type);

    int i,n;

    num_atom=original.num_atom;

    if(original.atom!=NULL)
    {
        atom=new struct atom[num_atom]; if(atom==NULL) Memory_Allocation_Error();
        for(i=0;i<num_atom;i++) atom[i]=original.atom[i];
    }

    num_bond=original.num_bond;

    if(original.bond!=NULL)
    {
        bond=new Bond[num_bond]; if(bond==NULL) Memory_Allocation_Error();
        for(i=0;i<num_bond;i++) bond[i]=original.bond[i];
    }

    num_ring=original.num_ring;
    ring.clear(); n=original.ring.size();
    for(i=0;i<n;i++) ring.push_back(original.ring[i]);

    vol_dot.clear(); n=original.vol_dot.size();
    for(i=0;i<n;i++) vol_dot.push_back(original.vol_dot[i]);

    sur_dot.clear(); n=original.sur_dot.size();
    for(i=0;i<n;i++) sur_dot.push_back(original.sur_dot[i]);

    // this part is for Ligand

//    bind_score=original.bind_score;
    chem_score=original.chem_score;
    vdw=original.vdw;
    sb=original.sb;
    hb=original.hb;
    hp=original.hp;
    hm=original.hm;
    hs=original.hs;
    ar=original.ar;
    rt=original.rt;
    pmf=original.pmf;
    uhb=original.uhb;
    bnsur=original.bnsur;
    bpsur=original.bpsur;
    pkd1=original.pkd1;
    pkd2=original.pkd2;
    pkd3=original.pkd3;

//    if(original.abs_inf!=NULL)
//    {
//        abs_inf=new ABS[num_atom]; if(abs_inf==NULL) Memory_Allocation_Error();
//        for(i=0;i<num_atom;i++) abs_inf[i]=original.abs_inf[i];
//    }
}

Ligand& Ligand::operator = (const Ligand &original)
{
    if(this==&original) return *this;

    // this part is from Molecule

    xtool_format=original.xtool_format;
    id=original.id;
    valid=original.valid;
    strcpy(name,original.name);
    weight=original.weight;
    strcpy(formula,original.formula);
    logp=original.logp;
    surface=original.surface;
    nsur=original.nsur;
    psur=original.psur;
    volume=original.volume;
    num_hb_atom=original.num_hb_atom;
    num_rotor=original.num_rotor;

    num_subst=original.num_subst;
    num_feature=original.num_feature;
    num_set=original.num_set;
    strcpy(mol_type,original.mol_type);
    strcpy(charge_type,original.charge_type);

    int i,n;

    // you need to delete atom and bond first before re-create them
    // this is slightly different from the copy constructor

    num_atom=original.num_atom;
    if(atom) delete [] atom; atom=NULL;

    if(original.atom!=NULL)
    {
        atom=new struct atom[num_atom]; if(atom==NULL) Memory_Allocation_Error();
        for(i=0;i<num_atom;i++) atom[i]=original.atom[i];
    }

    num_bond=original.num_bond;
    if(bond) delete [] bond; bond=NULL;

    if(original.bond!=NULL)
    {
        bond=new Bond[num_bond]; if(bond==NULL) Memory_Allocation_Error();
        for(i=0;i<num_bond;i++) bond[i]=original.bond[i];
    }

    num_ring=original.num_ring;
    ring.clear(); n=original.ring.size();
    for(i=0;i<n;i++) ring.push_back(original.ring[i]);

    vol_dot.clear(); n=original.vol_dot.size();
    for(i=0;i<n;i++) vol_dot.push_back(original.vol_dot[i]);

    sur_dot.clear(); n=original.sur_dot.size();
    for(i=0;i<n;i++) sur_dot.push_back(original.sur_dot[i]);

    // this part is for Ligand

    bind_score=original.bind_score;
    chem_score=original.chem_score;
    vdw=original.vdw;
    sb=original.sb;
    hb=original.hb;
    hp=original.hp;
    hm=original.hm;
    hs=original.hs;
    ar=original.ar;
    rt=original.rt;
    pmf=original.pmf;
    uhb=original.uhb;
    bnsur=original.bnsur;
    bpsur=original.bpsur;
    pkd1=original.pkd1;
    pkd2=original.pkd2;
    pkd3=original.pkd3;

//    if(abs_inf) delete [] abs_inf; abs_inf=NULL;
//
//    if(original.abs_inf!=NULL)
//    {
//        abs_inf=new ABS[num_atom]; if(abs_inf==NULL) Memory_Allocation_Error();
//        for(i=0;i<num_atom;i++) abs_inf[i]=original.abs_inf[i];
//    }

    return *this;
}

void Ligand::Clear()
{
    Molecule::Clear();

    bind_score=chem_score=0.000;
    vdw=0.000;
    sb=0.000;
    hb=0.000;
    hp=0.000;
    hm=0.000;
    hs=0.000;
    ar=0.000;
    rt=0.000;
    pmf=0.000;
    uhb=0.000;
    bnsur=bpsur=0.000;
    pkd1=pkd2=pkd3=0.000;

    abs_inf=NULL;

    return;
}
void Ligand::Calculate_HB_Root()
{
    int i,j,tmp,num_nonh;
    float tmpx,tmpy,tmpz;

    for(i=0;i<num_atom;i++)
    {
        if(atom[i].valid==0) continue;
        else if(!strcmp(atom[i].hb,"N")) continue;
        else if(!strcmp(atom[i].hb,"H")) continue;
        else if(!strcmp(atom[i].hb,"P")) continue;

        tmpx=tmpy=tmpz=0.000; num_nonh=0;

        for(j=0;j<atom[i].num_neib;j++)
        {
            tmp=atom[i].neib[j]-1;
            if(atom[tmp].type[0]=='H') continue;
            else
            {
                tmpx+=atom[tmp].coor[0];
                tmpy+=atom[tmp].coor[1];
                tmpz+=atom[tmp].coor[2];
                num_nonh++;
            }
        }

        if(num_nonh==0) strcpy(atom[i].hb,"P");
        else
        {
            tmpx/=num_nonh;
            tmpy/=num_nonh;
            tmpz/=num_nonh;
            atom[i].root[0]=tmpx;
            atom[i].root[1]=tmpy;
            atom[i].root[2]=tmpz;
        }
    }

    return;
}


int Ligand::Value_Atom()
{
    int mark;

    mark=Molecule::Value_Atom();

    // now compute the H-bond root for each HB atom

    Calculate_HB_Root();

    // assign atomic logp values, which is needed by binding score
    // atomic logp is not read from ATOM_DEF_XTOOL!

    this->logp=Calculate_LogP();

    // some other molecular properties

    this->num_hb_atom=this->Get_Num_HB_Atom();
    this->num_rotor=(int)this->Count_Rotor();

    return mark;
}
//float Ligand::Calculate_VDW_ij(atom) {
//    extern Protein *protein;
//    float cutoff,d0,d,tmp,tmp1,tmp2,sum,asum;
//
//    // clear the variables
//
//    cutoff=DIST_CUTOFF; sum=0.000;
//
//}
float Ligand::Calculate_VDW()
{
    extern Protein *protein;
    int i,j;
    float cutoff,d0,d,tmp,tmp1,tmp2,sum,asum;

    // clear the variables

    cutoff=DIST_CUTOFF; sum=0.000;
    for(i=0;i<this->num_atom;i++) this->atom[i].score=0.000;

    // now calculate the P-L vdw interaction
    int count = 0;
    for(i=0;i<this->num_atom;i++)
    {
        if(this->atom[i].valid<=0) continue;
        else if(this->atom[i].type[0]=='H') continue;

        asum=0.000;

        for(j=0;j<protein->num_atom;j++)
        {
            if(protein->atom[j].valid!=2) continue;
            else if(protein->atom[j].type[0]=='H') continue;
            else if(!strcmp(protein->atom[j].type,"O.w")) continue;

            d0=this->atom[i].R+protein->atom[j].R;
            d=Distance(this->atom[i].coor,protein->atom[j].coor);

            if(d>cutoff) continue;

            // Lennard-Jones 8-4 potential

            tmp1=d0/d;
            tmp1=tmp1*tmp1*tmp1*tmp1; tmp2=tmp1*tmp1;
            tmp=tmp2-2.00*tmp1;

            asum+=tmp;
            count++;
        }

        asum*=(-1.00); // change the sign so that positive values are favorable

        // if this atom has unfavorable contribution then neglect it

        if(asum<0.00) continue;
        else {sum+=asum; this->atom[i].score=asum;}
    }

    return sum;
}

float Ligand::Calculate_HB(char *pdb_entry)
{
    int i,j;
    float sum;
    HBond candidate[1000],temp;
    int num_candidate;

    // clear the variables

    sum=0.000; num_candidate=0;
    for(i=0;i<this->num_atom;i++) this->atom[i].score=0.000;

    // first, get the HB pairs between protein and ligand

    num_candidate=this->Get_HBond_Pair_PL(candidate);
    this->Sum_HBonds(num_candidate,candidate);

    // then sum their contributions

    for(i=0;i<num_candidate;i++)
    {
        if(fabs(candidate[i].score)<0.01) continue;

        sum+=candidate[i].score;
        this->atom[candidate[i].latom-1].score+=candidate[i].score;
    }

    // output H-bonds for analysis

//    if(strcmp(pdb_entry,"unknown"))
//    {
//        // arrange them according to ligand atom id
//
//        for(i=0;i<num_candidate-1;i++)
//            for(j=i+1;j<num_candidate;j++)
//            {
//                if(candidate[i].latom<=candidate[j].latom) continue;
//                else {SWAP(candidate[i],candidate[j]);}
//            }
//
//        for(i=0;i<num_candidate;i++)
//        {
//            if(fabs(candidate[i].score)<=0.000) continue;
//            else candidate[i].Show_Contents();
//        }
//    }

    return sum;
}


float Ligand::Calculate_RT()
{
    int i,j,mark,tmp;
    int id1,id2;
    float sum;

    // if a single bond is normal, bond.valid=1;
    // if a single bond is a rotor, bond.valid=2;

    // clean the variables

    for(i=0;i<this->num_atom;i++) this->atom[i].score=0.000;

    // eliminate all the bonds in_ring and all the non-single bonds

    for(i=0;i<num_bond;i++)
    {
        if(bond[i].ring>0) continue;
        else if(strcmp(bond[i].type,"1")) continue;
        else bond[i].valid=2;
    }

    // eliminate all the R-H, R-X, R-OH, R-NH2, R-CH3 bonds

    for(i=0;i<num_bond;i++)
    {
        if(bond[i].valid!=2) continue;

        id1=bond[i].atom_1; id2=bond[i].atom_2;

        if(atom[id1-1].num_nonh==1||
           atom[id2-1].num_nonh==1) bond[i].valid=1;
        else continue;
    }

    // sp2-sp2 rotors

    for(i=0;i<num_bond;i++)
    {
        if(bond[i].valid!=2) continue;

        id1=bond[i].atom_1; id2=bond[i].atom_2;

        mark=0;

        tmp=Get_Atom_Hybridizing_Type(atom[id1-1].type);
        if(tmp==1||tmp==2) mark++;

        tmp=Get_Atom_Hybridizing_Type(atom[id2-1].type);
        if(tmp==1||tmp==2) mark++;

        if(mark==2) {bond[i].valid=1; continue;}
    }

    // eliminate terminal rotors, e.g. -PO3, -CF3, -CMe3, -NMe3

    for(i=0;i<num_bond;i++)
    {
        if(bond[i].valid!=2) continue;

        id1=bond[i].atom_1; id2=bond[i].atom_2;

        if(Judge_Terminal_Atom(atom[id1-1],id2)==TRUE)
        {
            bond[i].valid=1;
        }
        else if(Judge_Terminal_Atom(atom[id2-1],id1)==TRUE)
        {
            bond[i].valid=1;
        }
        else continue;
    }

    // eliminate abnormal rotors

    for(i=0;i<num_bond;i++)
    {
        if(bond[i].valid!=2) continue;

        id1=bond[i].atom_1; id2=bond[i].atom_2;

        if(atom[id1-1].valid<=0) bond[i].valid=1;
        else if(atom[id2-1].valid<=0) bond[i].valid=1;
        else continue;
    }

    // now count the frozen rotors, all the rotors have been labeled as 2

    // notice: the following part is different from Molecule::Count_Rotor()

    sum=0.000;

    for(i=0;i<num_atom;i++)
    {
        if(atom[i].valid<=0) continue;

        mark=0;

        for(j=0;j<atom[i].num_neib;j++)
        {
            tmp=atom[i].bond[j];
            if(tmp==0) continue;
            else if(bond[tmp-1].valid!=2) continue;
            else mark++;
        }

        if(mark==1) this->atom[i].score+=0.50;
        else if(mark==2) this->atom[i].score+=1.00;
        else if(mark>=3) this->atom[i].score+=0.50;

        sum+=this->atom[i].score;
    }

    return sum;
}


float Ligand::Calculate_HP()
{
    extern Protein *protein;
    int i,j;
    float d,d1,d2,cutoff;
    float tmp,asum,sum;

    sum=0.000;
    for(i=0;i<this->num_atom;i++) this->atom[i].score=0.000;

    cutoff=DIST_CUTOFF;

    for(i=0;i<this->num_atom;i++)
    {
        if(this->atom[i].valid<=0) continue;
        else if(this->atom[i].type[0]=='H') continue;
        else if(strcmp(this->atom[i].hb,"H")) continue;

        asum=0.000;

        for(j=0;j<protein->num_atom;j++)
        {
            if(protein->atom[j].valid!=2) continue;
            else if(protein->atom[j].type[0]=='H') continue;
            else if(!strcmp(protein->atom[j].type,"O.w")) continue;
            else if(strcmp(protein->atom[j].hb,"H")) continue;

            d=Distance(this->atom[i].coor, protein->atom[j].coor);
            if(d>=cutoff) continue;

            d1=this->atom[i].R+protein->atom[j].R+0.500;
            d2=this->atom[i].R+protein->atom[j].R+2.200;

            if(d<d1) tmp=1.000;
            else if(d<d2) tmp=(1/(d1-d2))*(d-d2);
            else tmp=0.000;

            asum+=tmp;
        }

        this->atom[i].score=asum; sum+=asum;
    }

    return sum;
}


float Ligand::Calculate_HM()
{
    extern Protein *protein;
    int i,j;
    float asum,sum,tmp,d,d1,d2,total,cutoff;

    for(i=0;i<this->num_atom;i++) this->atom[i].score=0.000;

    sum=0.000; cutoff=DIST_CUTOFF;

    for(i=0;i<this->num_atom;i++)
    {
        if(this->atom[i].valid<=0) continue;
        else if(this->atom[i].type[0]=='H') continue;
        else if(strcmp(this->atom[i].hb,"H")) continue;
        else if(this->atom[i].logp<=0.00) continue;

        total=0.000; asum=0.000;

        for(j=0;j<protein->num_atom;j++)
        {
            if(protein->atom[j].valid!=2) continue;
            else if(protein->atom[j].type[0]=='H') continue;
            else if(!strcmp(protein->atom[j].type,"O.w")) continue;

            d=Distance(this->atom[i].coor, protein->atom[j].coor);
            if(d>cutoff) continue;

            d1=this->atom[i].R+protein->atom[j].R+0.50;
            d2=this->atom[i].R+protein->atom[j].R+2.20;

            if(d<d1) tmp=1.000;
            else if(d<d2) tmp=(1/(d1-d2))*(d-d2);
            else tmp=0.000;

            total+=(protein->atom[j].logp*tmp);
        }

        if(this->atom[i].logp>=0.50) asum=this->atom[i].logp;
        else if(total>-0.50) asum=this->atom[i].logp;
        else asum=0.000;

        this->atom[i].score=asum; sum+=asum;
    }

    return sum;
}


float Ligand::Calculate_HS()
{
    int i;
    float sum,total,buried;

    // clear the variables first

    sum=0.000;
    for(i=0;i<this->num_atom;i++) this->atom[i].score=0.000;

    // then get the buried surface atom by atom

    for(i=0;i<this->num_atom;i++)
    {
        if(this->atom[i].valid<=0) continue;
        else if(this->atom[i].type[0]=='H') continue;
        else if(!strcmp(this->atom[i].hb,"DH")) continue;
        else if(!strcmp(this->atom[i].hb,"D")) continue;
        else if(!strcmp(this->atom[i].hb,"A")) continue;
        else if(!strcmp(this->atom[i].hb,"DA")) continue;
        else if(!strcmp(this->atom[i].hb,"P")) continue;
        else if(!strcmp(this->atom[i].hb,"N")) continue;

        this->Atom_Buried_Surface(i+1,total,buried);

        sum+=buried; this->atom[i].score=buried;

        // sum+=(buried*this->atom[i].logp);
        // this option works *slightly* better
    }

    return sum;
}



float Ligand::Calculate_Binding_Score()
{
    extern Input *input;
    int i;
    Ligand bak;

    if(input->num_method==0) return 0.000;
    else bak=(*this);  // backup the atomic charge information

    // these scoring functions may need formal charges

    this->Assign_Apparent_Charge();

    // prepare for atomic binding score

//    if(this->abs_inf) delete [] this->abs_inf;
//    this->abs_inf=new ABS[this->num_atom];
//    if(abs_inf==NULL) Memory_Allocation_Error();
//
//    for(i=0;i<this->num_atom;i++)
//    {
//        abs_inf[i].pkd1=0.000;
//        abs_inf[i].pkd2=0.000;
//        abs_inf[i].pkd3=0.000;
//        abs_inf[i].vdw=0.000;
//        abs_inf[i].hb=0.000;
//        abs_inf[i].hm=0.000;
//        abs_inf[i].hp=0.000;
//        abs_inf[i].hs=0.000;
//        abs_inf[i].rt=0.000;
//        abs_inf[i].score=0.000;
//    }

    // clean other variables

    this->pkd1=this->pkd2=this->pkd3=0.000;
    this->vdw=this->hb=this->hm=this->hp=this->hs=this->rt=0.000;

    // first, generate dot surface

    this->Generate_Surface_Dots(WATER_R);

    // now calculate the van der Waals term
//    struct timeval t1,t2;
//    double timeuse;
//    gettimeofday(&t1,NULL);

    this->vdw=this->Calculate_VDW();
//    gettimeofday(&t2,NULL);
//    timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
//    std::cout<<"xscore_vdw_time = "<<timeuse<<std::endl;  //���ʱ�䣨��λ����
//    for(i=0;i<this->num_atom;i++)
//    {
//        abs_inf[i].vdw=this->atom[i].score;
//    }

    // now calculate the H-bond term
//    gettimeofday(&t1,NULL);
    this->hb=this->Calculate_HB();
//    gettimeofday(&t2,NULL);
//    timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
//    std::cout<<"xscore_HB_time = "<<timeuse<<std::endl;  //���ʱ�䣨��λ����
//    for(i=0;i<this->num_atom;i++)
//    {
//        abs_inf[i].hb=this->atom[i].score;
//    }

    // now calculate the rotor term
//    gettimeofday(&t1,NULL);
    this->rt=this->Calculate_RT();
//    gettimeofday(&t2,NULL);
//    timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
//    std::cout<<"xscore_RT_time = "<<timeuse<<std::endl;  //���ʱ�䣨��λ����
//    for(i=0;i<this->num_atom;i++)
//    {
//        abs_inf[i].rt=this->atom[i].score;
//    }

    // now calculate the hydrophobic term

    if(!strcmp(input->apply_hpscore,"YES"))
    {
//        gettimeofday(&t1,NULL);
        this->hp=this->Calculate_HP();
//        gettimeofday(&t2,NULL);
//        timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
//        std::cout<<"xscore_HP_time = "<<timeuse<<std::endl;  //���ʱ�䣨��λ����
//        for(i=0;i<this->num_atom;i++)
//        {
//            abs_inf[i].hp=this->atom[i].score;
//        }

        this->pkd1=input->hpscore_c0+
                   input->hpscore_cvdw*this->vdw+
                   input->hpscore_chb*this->hb+
                   input->hpscore_chp*this->hp+
                   input->hpscore_crt*this->rt;
    }

    if(!strcmp(input->apply_hmscore,"YES"))
    {
//        gettimeofday(&t1,NULL);
        this->hm=this->Calculate_HM();
//        gettimeofday(&t2,NULL);
//        timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
//        std::cout<<"xscore_HM_time = "<<timeuse<<std::endl;  //���ʱ�䣨��λ����
//        for(i=0;i<this->num_atom;i++)
//        {
//            abs_inf[i].hm=this->atom[i].score;
//        }

        this->pkd2=input->hmscore_c0+
                   input->hmscore_cvdw*this->vdw+
                   input->hmscore_chb*this->hb+
                   input->hmscore_chm*this->hm+
                   input->hmscore_crt*this->rt;
    }

    if(!strcmp(input->apply_hsscore,"YES"))
    {
//        gettimeofday(&t1,NULL);
        this->hs=this->Calculate_HS();
//        gettimeofday(&t2,NULL);
//        timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
//        std::cout<<"xscore_HS_time = "<<timeuse<<std::endl;  //���ʱ�䣨��λ����
//        for(i=0;i<this->num_atom;i++)
//        {
//            abs_inf[i].hs=this->atom[i].score;
//        }

        this->pkd3=input->hsscore_c0+
                   input->hsscore_cvdw*this->vdw+
                   input->hsscore_chb*this->hb+
                   input->hsscore_chs*this->hs+
                   input->hsscore_crt*this->rt;
    }

    this->bind_score=(pkd1+pkd2+pkd3)/(input->num_method);
//    std::cout<<""<<std::endl;
    // now compute atomic binding scores
//
//    int num_nonh=0; num_nonh=this->Get_Num_Heavy_Atom();
//
//    if(!strcmp(input->apply_hpscore,"YES"))
//    {
//        for(i=0;i<this->num_atom;i++)
//        {
//            if(this->atom[i].valid<=0) continue;
//            else if(!strcmp(this->atom[i].type,"H")) continue;
//
//            abs_inf[i].pkd1=input->hpscore_c0/num_nonh;
//            abs_inf[i].pkd1+=input->hpscore_cvdw*abs_inf[i].vdw;
//            abs_inf[i].pkd1+=input->hpscore_chb*abs_inf[i].hb;
//            abs_inf[i].pkd1+=input->hpscore_chp*abs_inf[i].hp;
//            abs_inf[i].pkd1+=input->hpscore_crt*abs_inf[i].rt;
//        }
//    }
//
//    if(!strcmp(input->apply_hmscore,"YES"))
//    {
//        for(i=0;i<this->num_atom;i++)
//        {
//            if(this->atom[i].valid<=0) continue;
//            else if(!strcmp(this->atom[i].type,"H")) continue;
//
//            abs_inf[i].pkd2=input->hmscore_c0/num_nonh;
//            abs_inf[i].pkd2+=input->hmscore_cvdw*abs_inf[i].vdw;
//            abs_inf[i].pkd2+=input->hmscore_chb*abs_inf[i].hb;
//            abs_inf[i].pkd2+=input->hmscore_chm*abs_inf[i].hm;
//            abs_inf[i].pkd2+=input->hmscore_crt*abs_inf[i].rt;
//        }
//    }
//
//    if(!strcmp(input->apply_hsscore,"YES"))
//    {
//        for(i=0;i<this->num_atom;i++)
//        {
//            if(this->atom[i].valid<=0) continue;
//            else if(!strcmp(this->atom[i].type,"H")) continue;
//
//            abs_inf[i].pkd3=input->hsscore_c0/num_nonh;
//            abs_inf[i].pkd3+=input->hsscore_cvdw*abs_inf[i].vdw;
//            abs_inf[i].pkd3+=input->hsscore_chb*abs_inf[i].hb;
//            abs_inf[i].pkd3+=input->hsscore_chs*abs_inf[i].hs;
//            abs_inf[i].pkd3+=input->hsscore_crt*abs_inf[i].rt;
//        }
//    }
//
//    float sum=0.000;
//
//    for(i=0;i<this->num_atom;i++)
//    {
//        abs_inf[i].score=0.000;
//        abs_inf[i].score+=abs_inf[i].pkd1;
//        abs_inf[i].score+=abs_inf[i].pkd2;
//        abs_inf[i].score+=abs_inf[i].pkd3;
//        abs_inf[i].score/=input->num_method;
//        this->atom[i].score=abs_inf[i].score;
//        sum+=this->atom[i].score;
//    }

//    if(!strcmp(input->show_abs,"YES"))  // use abs as charges
//    {
//        strcpy(this->charge_type,"USER_CHARGES");
//
//        for(i=0;i<this->num_atom;i++)
//        {
//            this->atom[i].q=abs_inf[i].score;
//        }
//    }
//    else		// restore original charges
//    {
//        strcpy(this->charge_type, bak.charge_type);
//
//        for(i=0;i<this->num_atom;i++)
//        {
//            this->atom[i].q=bak.atom[i].q;
//        }
//    }

    return this->bind_score;
}


int Ligand::Get_HBond_Pair_PL(HBond candidate[], bool sb_flag) const
{
    extern Protein *protein;
    int i,j,num,type;
    bool sb;
    float d,cutoff;
    HBond tmp_candidate;

    num=0; cutoff=5.00;

    // note that the following H-bond algorithm is not based on any
    // explicit hydrogen atom at all, and this is the right thing to do.

    for(i=0;i<this->num_atom;i++)
    {
        if(this->atom[i].valid<=0) continue;
        else if(!strcmp(this->atom[i].type,"H")) continue;
        else if(!strcmp(this->atom[i].hb,"N")) continue;
        else if(!strcmp(this->atom[i].hb,"H")) continue;
        else if(!strcmp(this->atom[i].hb,"P")) continue;
        else if(!strcmp(this->atom[i].hb,"DH")) continue;

        for(j=0;j<protein->num_atom;j++)
        {
            if(protein->atom[j].valid!=2) continue;
            else if(!strcmp(protein->atom[j].type,"H")) continue;
            else if(!strcmp(protein->atom[j].type,"O.w")) continue;
            else if(!strcmp(protein->atom[j].hb,"H")) continue;
            else if(!strcmp(protein->atom[j].hb,"P")) continue;
            else if(!strcmp(protein->atom[j].hb,"N")) continue;

            // determine the type of this H-bond first
            // type=0, no H-bond
            // type=1, latom is the donor, patom is the acceptor;
            // type=2, latom is the acceptor, patom is the donor;
            // type=3, latom bound with metal ion;

            if(!strcmp(this->atom[i].hb,"D"))
            {
                if(!strcmp(protein->atom[j].hb,"D")) type=0;
                else if(!strcmp(protein->atom[j].hb,"A")) type=1;
                else if(!strcmp(protein->atom[j].hb,"DA")) type=1;
                else if(!strcmp(protein->atom[j].hb,"M")) type=0;
                else type=0;
            }
            else if(!strcmp(this->atom[i].hb,"A"))
            {
                if(!strcmp(protein->atom[j].hb,"D")) type=2;
                else if(!strcmp(protein->atom[j].hb,"A")) type=0;
                else if(!strcmp(protein->atom[j].hb,"DA")) type=2;
                else if(!strcmp(protein->atom[j].hb,"M")) type=3;
                else type=0;
            }
            else if(!strcmp(this->atom[i].hb,"DA"))
            {
                if(!strcmp(protein->atom[j].hb,"A")) type=1;
                else if(!strcmp(protein->atom[j].hb,"D")) type=2;
                else if(!strcmp(protein->atom[j].hb,"DA")) type=1;
                else if(!strcmp(protein->atom[j].hb,"M")) type=3;
                else type=0;
            }
            else type=0;

            if(type==0) continue;  // no H-bond

            // a crude distance check

            d=Distance(this->atom[i].coor,protein->atom[j].coor);
            if(d>cutoff) continue;

            // this section is used to differentiate HB and SB

            if(!strcmp(this->atom[i].type,"O.co2")&&
               !strcmp(protein->atom[j].type,"O.co2"))
            {
                sb=false;
            }
            else if((fabs(this->atom[i].q)>0.01)&&
                    (fabs(protein->atom[j].q)>0.01))
            {
                sb=true;
            }
            else
            {
                sb=false;
            }

            // now handle this h-bond

            tmp_candidate.Clear();
            tmp_candidate.type=type;
            tmp_candidate.sb=sb;
            tmp_candidate.latom=i+1;
            tmp_candidate.patom=j+1;

            if(type==1)
            {
                tmp_candidate.D=this->atom[i];
                tmp_candidate.A=protein->atom[j];
            }
            else if(type==2)
            {
                tmp_candidate.A=this->atom[i];
                tmp_candidate.D=protein->atom[j];
            }
            else if(type==3)
            {
                tmp_candidate.A=this->atom[i];
                tmp_candidate.D=protein->atom[j];
            }

            // if need to treat charged H-bonds differently

            if(sb_flag==true&&sb==true) tmp_candidate.Value_SBond();
            else tmp_candidate.Value_HBond_2();

            if(fabs(tmp_candidate.score)>0.000)
            {
                candidate[num]=tmp_candidate; num++;
            }
            else continue;
        }
    }

    return num;
}

void Ligand::Sum_HBonds(int num_candidate, HBond candidate[]) const
{
    extern Protein *protein;
    int i,j,k,count,limit,latom,patom;
    float angle,v1[3],v2[3];
    HBond temp;
    Group tmp_group;

    // Step 1: rank candidates according their strength in decreasing order
    // note "fabs" is applied because SB strength could be negative

    for(i=0;i<num_candidate-1;i++)
        for(j=i+1;j<num_candidate;j++)
        {
            if(fabs(candidate[i].score)>=fabs(candidate[j].score)) continue;
            else {SWAP(candidate[i],candidate[j]);}
        }

    // Step 2: check the angular limit: the angle between any two H-bonds
    // on the same atom must be larger than 45 degrees
    // note this this filter could be applied to both ligand and protein

    for(i=0;i<num_candidate-1;i++)
        for(j=i+1;j<num_candidate;j++)
        {
            if(candidate[i].latom!=candidate[j].latom) continue;

            for(k=0;k<3;k++)
            {
                latom=candidate[i].latom; patom=candidate[i].patom;
                v1[k]=protein->atom[patom-1].coor[k]-this->atom[latom-1].coor[k];
                latom=candidate[j].latom; patom=candidate[j].patom;
                v2[k]=protein->atom[patom-1].coor[k]-this->atom[latom-1].coor[k];
            }

            angle=fabs(Angle_of_Two_Vectors(v1,v2));

            if(angle<45.0) candidate[j].score=0.000;
            else continue;
        }

/*
 for(i=0;i<num_candidate-1;i++)
 for(j=i+1;j<num_candidate;j++)
        {
         if(candidate[i].patom!=candidate[j].patom) continue;

         for(k=0;k<3;k++)
        {
	 latom=candidate[i].latom; patom=candidate[i].patom;
	 v1[k]=this->atom[latom-1].coor[k]-protein->atom[patom-1].coor[k];
	 latom=candidate[j].latom; patom=candidate[j].patom;
	 v2[k]=this->atom[latom-1].coor[k]-protein->atom[patom-1].coor[k];
        }

         angle=fabs(Angle_of_Two_Vectors(v1,v2));

         if(angle<45.0) candidate[j].score=0.000;
         else continue;
        }
*/

    // Step 3: an donor atom shall not form more H-bonds than the H atoms it has.
    // note that this filter is applied only to the ligand side

    for(i=0;i<num_candidate-1;i++)
    {
        if(candidate[i].type!=1) continue;

        count=1; latom=candidate[i].latom;

        if(!strcmp(this->atom[latom-1].hb,"DA"))
        {
            limit=1;
        }
        else
        {
            tmp_group=this->Find_A_Group(latom);
            limit=tmp_group.num_h;
        }

        for(j=i+1;j<num_candidate;j++)
        {
            if(candidate[j].type!=1) continue;
            else if(candidate[i].latom!=candidate[j].latom) continue;

            count++;

            if(count<=limit) continue;
            else candidate[j].score=0.000;
        }
    }

    // Step 4: an acceptor atom shall not form more H-bonds than its LPs
    // note that this filter is applied only to the ligand side

    for(i=0;i<num_candidate-1;i++)
    {
        if(candidate[i].type!=2&&candidate[i].type!=3) continue;

        count=1; latom=candidate[i].latom;

        if(this->atom[latom-1].type[0]=='O') limit=2;
        else if(this->atom[latom-1].type[0]=='N') limit=1;
        else if(this->atom[latom-1].type[0]=='S') limit=2;

        for(j=i+1;j<num_candidate;j++)
        {
            if(candidate[j].type!=2&&candidate[j].type!=3) continue;
            else if(candidate[i].latom!=candidate[j].latom) continue;

            count++;

            if(count<=limit) continue;
            else candidate[j].score=0.000;
        }
    }

    return;
}


float Ligand::Atom_Buried_Surface(int id, float &total, float &buried) const
{
    extern Protein *protein;
    int j,k,num;
    bool mark;
    float d,ratio;
    int *atom_check_list;

    atom_check_list=new int[protein->num_atom];
    if(atom_check_list==NULL) Memory_Allocation_Error();

    total=buried=0.000;

    for(j=0;j<protein->num_atom;j++) atom_check_list[j]=0;

    for(j=0;j<protein->num_atom;j++)
    {
        if(protein->atom[j].valid!=2) continue;
        else if(protein->atom[j].xtype[0]=='H') continue;
        else if(!strcmp(protein->atom[j].xtype,"O.w")) continue;

        d=Distance(this->atom[id-1].coor,protein->atom[j].coor);

        if(d>(this->atom[id-1].R+protein->atom[j].R+2*WATER_R)) continue;
        else atom_check_list[j]=1;
    }

    num=this->sur_dot.size();

    for(j=0;j<num;j++)
    {
        if(this->sur_dot[j].valid!=id) continue;

        // check if this dot is buried

        mark=false;

        for(k=0;k<protein->num_atom;k++)
        {
            if(atom_check_list[k]==0) continue;

            d=Distance(this->sur_dot[j].coor, protein->atom[k].coor);

            if(d>(protein->atom[k].R+WATER_R)) continue;
            else {mark=true; break;}
        }

        total+=this->sur_dot[j].unit;

        if(mark==false) continue;
        else buried+=this->sur_dot[j].unit;
    }

    if(atom_check_list) delete [] atom_check_list;

    if(total<=0.00) ratio=0.00;
    else ratio=buried/total;

    return ratio;
}
