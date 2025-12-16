////
//// Created by 91686 on 2023/8/14.
////
//
//#include "xscore.h"
//#include "file_error.h"
///*
// * id?
// * Count()?
// * xtype?
// * Distance
// * Group Neib
// *atom.score
// * r和R不一样
// * Calculate_Log_protein的protein_atoms
// * #define
// *
// * */
//float xscore::xscore_vdw() {
//    fl sum = 0.0;
//    fl asum;
//    float d0,d,tmp,tmp1,tmp2;
//    for(int i = 0; i < ligand_atoms.size();i++){
//        if(ligand_atoms[i].valid<=0) continue;
//        else if(ligand_atoms[i].type[0] == 'H') continue;
//
//        asum =0.0;
//
//        for(int j = 0; j < protein_atoms.size(); j++){
//            if(protein_atoms[j].valid!=2){ continue;}
//            else if(protein_atoms[j].type[0] == 'H') {continue;}
//            else if(!strcmp(protein_atoms[j].type,"O.w")){ continue;}
//
//            d0 = ligand_atoms[i].r + protein_atoms[i].r;
//            d = sqrt(vec_distance_sqr(ligand_atoms[i].coords, protein_atoms[j].coords));
//            if(d>cutoff) continue;
//
//            // Lennard-Jones 8-4 potential
//
//            tmp1=d0/d;
//            tmp1=tmp1*tmp1*tmp1*tmp1; tmp2=tmp1*tmp1;
//            tmp=tmp2-2.00*tmp1;
//
//            asum+=tmp;
//        }
//    asum*=(-1.00); // change the sign so that positive values are favorable
//
//    // if this atom has unfavorable contribution then neglect it
//
//    if(asum<0.00) continue;
//    else {sum+=asum; ligand_atoms[i].score=asum;}
//    }
//
//    return sum;
//}
////问题？修改的是ligand_atoms变量而不是model里真正的变量
//float xscore::xscore_RT() {
//    float sum = 0.0;
//    for(int i = 0; i < ligand_atoms.size(); i++){
//        ligand_atoms[i].score = 0.00;
//    }
//
//
//    for(int i = 0; i < ligand_atoms.size(); i++){
//        if(ligand_atoms[i].valid<=0) continue;
//        int mark = 0;
//
//        for(int j = 0; j < ligand_atoms[i].bonds.size(); j++){
//            if(ligand_atoms[i].bonds[j].rotatable){
//                mark++;
//            }
//        }
//        if(mark == 1) ligand_atoms[i].score += 0.50;
//        else if(mark == 2) ligand_atoms[i].score += 0.10;
//        else if(mark >= 3) ligand_atoms[i].score += 0.50;
//
//        sum += ligand_atoms[i].score;
//    }
////    for(int i = 0; i < ligand_atoms[i].bonds.size(); i++){
////        if(ligand_atoms[i].bonds[i].rotatable){/////////////////////////
////            mark++;
////        }
////        if(mark == 1) tmp += 0.50;
////        else if(mark == 2) tmp += 0.10;
////        else tmp += 0.50;
////    }
//
//
//    return sum; // placating the compiler
//}
//
//float xscore::calculate_HP() {
//
//    float asum,tmp,d,d1,d2,total,sum = 0.00;
//    for(int i = 0; i < ligand_atoms.size(); i++){
//        ligand_atoms[i].score = 0.00;
//    }
//
//    for(int i = 0; i < ligand_atoms.size(); i++){
//        if(ligand_atoms[i].valid<=0) continue;
//        else if(ligand_atoms[i].type[0]=='H') continue;
//        else if(strcmp(ligand_atoms[i].hb,"H")) continue;
//
//        asum=0.000;
//
//        for(int j = 0; j < protein_atoms.size(); j++){
//            if(protein_atoms[j].valid!=2) continue;
//            else if(protein_atoms[j].type[0]=='H') continue;
//            else if(!strcmp(protein_atoms[j].type,"O.w")) continue;
//            else if(strcmp(protein_atoms[j].hb,"H")) continue;
//
//
//            d = sqrt(vec_distance_sqr(ligand_atoms[i].coords, protein_atoms[j].coords));
//            if(d>cutoff) continue;
//
//            d1=ligand_atoms[i].r+protein_atoms[j].r+0.50;
//            d2=ligand_atoms[i].r+protein_atoms[j].r+2.20;
//
//            if(d<d1) tmp=1.000;
//            else if(d<d2) tmp=(1/(d1-d2))*(d-d2);
//            else tmp=0.000;
//
//            asum+=tmp;
//        }
//
//        ligand_atoms[i].score=asum; sum+=asum;
//    }
//
//    return sum;
//
//}
//
////HM
//void xscore::Read_XATOM_DEF(char *filename)
//{
//    FILE *fp;
//    int count;
//    char buf[256];
//
//    if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);
//
//    // first, determine num_xatomtype
//
//    count=0;
//
//    for(;;)
//    {
//        if(fgets(buf,256,fp)==NULL) break;
//        else if(buf[0]=='#') continue;
//        else if(Blank_Line_Check(buf)==true) continue;
//        else count++;
//    }
//
//    rewind(fp);
//
//    num_xatomtype=count;
//
//    xatom=new Xatom_Def[num_xatomtype];
////        if(xatom==NULL) Memory_Allocation_Error();
//
//    // second, read parameters for each atom type
//
//    count=0;
//
//    for(;;)
//    {
//        if(fgets(buf,256,fp)==NULL) break;
//        else if(buf[0]=='#') continue;
//        else if(Blank_Line_Check(buf)==true) continue;
//        else
//        {
//            sscanf(buf,"%*d%s%*s%f",
//                   xatom[count].type,
//                   &xatom[count].logp);
//            count++;
//        }
//    }
//
//    fclose(fp);
//
//    if(num_xatomtype!=count) Read_File_Error(filename);
//
//    return;
//}
//
//
//float xscore::Get_Atom_LogP(char *type) {
//    bool mark;
//    float logp;
//
//    mark=false;
//
//    for(int i=0;i<num_xatomtype;i++)
//    {
//        if(strcmp(type,xatom[i].type)) continue;
//        else {logp=xatom[i].logp; mark=true; break;}
//    }
//
//    if(mark==true) return logp;
//    else
//    {
//        printf("Warning: no LogP parameter for atom type %s ... ", type);
//        printf("zero value assigned\n");
//        return 0.000;
//    }
//}
//
//Group xscore::Find_X_Group(int atom_id) {
//    int i,num;
//    Group group;
//
//    // define the center
//
//    group.center=atom[atom_id-1];
//
//    // find the center's neighbours
//
//    num=group.center.num_neib;
//
//    for(i=0;i<num;i++)
//    {
//        group.neib[i]=this->atom[group.center.neib[i]-1];
//        group.bond[i]=this->bond[group.center.bond[i]-1];
//    }
//
//    // count the necessary parameters
//
//    group.num_neib=num;
//
//    group.num_h=0; group.num_nonh=0;
//
//    for(i=0;i<num;i++)
//    {
//        if(group.neib[i].type[0]=='H') group.num_h++;
//        else group.num_nonh++;
//    }
//
//    group.num_hetero=0;
//
//    for(i=0;i<group.num_nonh;i++)
//    {
//        if(!strcmp(group.bond[i].type,"2")) continue;
//        else if(!strcmp(group.bond[i].type,"3")) continue;
//        else if(!strcmp(group.bond[i].type,"ar")) continue;
//        else if(group.neib[i].type[0]=='N') group.num_hetero++;
//        else if(group.neib[i].type[0]=='O') group.num_hetero++;
//        else continue;
//    }
//
//    group.num_pi=0;
//
//    for(i=0;i<group.num_nonh;i++)
//    {
//        if(!strcmp(group.bond[i].type,"2")) continue;
//        else if(!strcmp(group.bond[i].type,"3")) continue;
//        else if(!strcmp(group.neib[i].type,"C.ar")) group.num_pi++;
//        else if(!strcmp(group.neib[i].type,"C.2")) group.num_pi++;
//        else if(!strcmp(group.neib[i].type,"C.1"))group.num_pi++;
//        else if(!strcmp(group.neib[i].type,"C.cat")) group.num_pi++;
//        else if(!strcmp(group.neib[i].type,"N.2")) group.num_pi++;
//        else if(!strcmp(group.neib[i].type,"N.1")) group.num_pi++;
//        else if(!strcmp(group.neib[i].type,"N.ar")) group.num_pi++;
//        else continue;
//    }
//
//    group.num_nar=group.num_car=0;
//
//    for(i=0;i<group.num_nonh;i++)
//    {
//        if(strcmp(group.bond[i].type,"ar")) continue;
//        else if(!strcmp(group.neib[i].type,"N.ar")) group.num_nar++;
//        else if(!strcmp(group.neib[i].type,"C.ar")) group.num_car++;
//        else if(group.neib[i].type[0]=='C') group.num_car++;
//        else group.num_nar++;
//    }
//
//    group.db_type=0;
//
//    for(i=0;i<group.num_nonh;i++)
//    {
//        if(!strcmp(group.center.type,"O.co2")||
//           !strcmp(group.center.type,"O.2")||
//           !strcmp(group.center.type,"S.2"))
//        {
//            if(group.neib[i].type[0]=='C') group.db_type=1;
//            else if(group.neib[i].type[0]=='N') group.db_type=2;
//            else if(group.neib[i].type[0]=='O') group.db_type=3;
//            else if(group.neib[i].type[0]=='S') group.db_type=4;
//            else if(group.neib[i].type[0]=='P') group.db_type=5;
//            else continue;
//        }
//        else if(!strcmp(group.bond[i].type,"2")||
//                !strcmp(group.neib[i].type,"O.co2")||
//                !strcmp(group.neib[i].type,"O.2")||
//                !strcmp(group.neib[i].type,"S.2"))
//        {
//            if(group.neib[i].type[0]=='C') group.db_type=1;
//            else if(group.neib[i].type[0]=='N') group.db_type=2;
//            else if(group.neib[i].type[0]=='O') group.db_type=3;
//            else if(group.neib[i].type[0]=='S') group.db_type=4;
//            else if(group.neib[i].type[0]=='P') group.db_type=5;
//            else continue;
//        }
//        else continue;
//    }
//
//    group.valid=1; return group;
//}
//int xscore::Get_XLOGP_Type(int atom_id, char *type) {
//    Group group;
//    atom atom;
//
//    atom=ligand_atoms[atom_id-1]; group=Find_X_Group(atom.id);
//
//    strcpy(type,"Un");
//
//    if(!strcmp(atom.type,"H")||!strcmp(atom.type,"H.spc"))
//    {
//        if(group.neib[0].type[0]=='O') strcpy(type,"H.hb");
//        else if(group.neib[0].type[0]=='N') strcpy(type,"H.hb");
//            // else if(!strcmp(group.neib[0].type,"S.3")) strcpy(type,"H.hb");
//        else strcpy(type,"H");
//    }
//
//    if(!strcmp(atom.type,"C.3"))
//    {
//        if(group.num_nonh==1)
//        {
//            if(group.num_hetero==0)
//            {
//                if(group.num_pi==0) strcpy(type,"C.3.h3.pi=0");
//                else strcpy(type,"C.3.h3.pi=1");
//            }
//            else
//            {
//                strcpy(type,"C.3.h3.x");
//            }
//        }
//        else if(group.num_nonh==2)
//        {
//            if(group.num_hetero==0)
//            {
//                if(group.num_pi==0)
//                {
//                    strcpy(type,"C.3.h2.pi=0");
//                }
//                else if(group.num_pi==1)
//                {
//                    strcpy(type,"C.3.h2.pi=1");
//                }
//                else
//                {
//                    strcpy(type,"C.3.h2.pi=2");
//                }
//            }
//            else
//            {
//                if(group.num_pi==0)
//                {
//                    strcpy(type,"C.3.h2.x.pi=0");
//                }
//                else if(group.num_pi==1)
//                {
//                    strcpy(type,"C.3.h2.x.pi=1");
//                }
//                else
//                {
//                    strcpy(type,"C.3.h2.x.pi=2");
//                }
//            }
//        }
//        else if(group.num_nonh==3)
//        {
//            if(group.num_hetero==0)
//            {
//                if(group.num_pi==0)
//                {
//                    strcpy(type,"C.3.h.pi=0");
//                }
//                else if(group.num_pi==1)
//                {
//                    strcpy(type,"C.3.h.pi=1");
//                }
//                else
//                {
//                    strcpy(type,"C.3.h.pi>1");
//                }
//            }
//            else
//            {
//                if(group.num_pi==0)
//                {
//                    strcpy(type,"C.3.h.x.pi=0");
//                }
//                else if(group.num_pi==1)
//                {
//                    strcpy(type,"C.3.h.x.pi=1");
//                }
//                else
//                {
//                    strcpy(type,"C.3.h.x.pi>1");
//                }
//            }
//        }
//        else if(group.num_nonh==4)
//        {
//            if(group.num_hetero==0)
//            {
//                if(group.num_pi==0)
//                {
//                    strcpy(type,"C.3.pi=0");
//                }
//                else if(group.num_pi==1)
//                {
//                    strcpy(type,"C.3.pi=1");
//                }
//                else
//                {
//                    strcpy(type,"C.3.pi>1");
//                }
//            }
//            else
//            {
//                if(group.num_pi==0) strcpy(type,"C.3.x.pi=0");
//                else strcpy(type,"C.3.x.pi>0");
//            }
//        }
//        else
//        {
//            strcpy(type,"C.3.unknown");
//        }
//    }
//
//    if(!strcmp(atom.type,"C.2"))
//    {
//        if(group.num_nonh==1)
//        {
//            strcpy(type,"C.2.h2");
//        }
//        else if(group.num_nonh==2)
//        {
//            if(group.num_hetero==0)
//            {
//                if(group.num_pi==0) strcpy(type,"C.2.h.pi=0");
//                else strcpy(type,"C.2.h.pi=1");
//            }
//            else
//            {
//                if(group.num_pi==0)
//                {
//                    strcpy(type,"C.2.h.x.pi=0");
//                }
//                else strcpy(type,"C.2.h.x.pi=1");
//            }
//        }
//        else if(group.num_nonh==3)
//        {
//            if(group.num_hetero==0)
//            {
//                if(group.num_pi==0) strcpy(type,"C.2.pi=0");
//                else strcpy(type,"C.2.pi>0");
//            }
//            else if(group.num_hetero==1)
//            {
//                if(group.num_pi==0) strcpy(type,"C.2.x.pi=0");
//                else strcpy(type,"C.2.x.pi>0");
//            }
//            else
//            {
//                if(group.num_pi==0) strcpy(type,"C.2.x2.pi=0");
//                else strcpy(type,"C.2.x2.pi>0");
//            }
//        }
//        else
//        {
//            strcpy(type,"C.2.unknown");
//        }
//    }
//
//    if(!strcmp(atom.type,"C.cat"))
//    {
//        strcpy(type,"C.2.x2.pi>0");
//    }
//
//    if(!strcmp(atom.type,"C.ar"))
//    {
//        if(group.num_nonh==2)
//        {
//            if(group.num_nar==0) strcpy(type,"C.ar.h");
//            else strcpy(type,"C.ar.h.(X)");
//        }
//        else if(group.num_nonh==3)
//        {
//            if(group.num_nar==0)
//            {
//                if(group.num_hetero==0) strcpy(type,"C.ar");
//                else strcpy(type,"C.ar.x");
//            }
//            else
//            {
//                if(group.num_hetero==0)
//                {
//                    strcpy(type,"C.ar.(X)");
//                }
//                else
//                {
//                    strcpy(type,"C.ar.(X).x");
//                }
//            }
//        }
//        else
//        {
//            strcpy(type,"C.ar.unknown");
//        }
//    }
//
//    if(!strcmp(atom.type,"C.1"))
//    {
//        if(group.db_type!=0) strcpy(type,"C.1.==");
//        else if(group.num_nonh==1) strcpy(type,"C.1.h");
//        else if(group.num_nonh==2) strcpy(type,"C.1");
//        else strcpy(type,"C.1.unknown");
//    }
//
//    if(!strcmp(atom.type,"N.4")||
//       !strcmp(atom.type,"N.3")||
//       !strcmp(atom.type,"N.pl3"))
//    {
//        if(group.num_nonh==1)
//        {
//            if(group.num_hetero==0)
//            {
//                if(group.num_pi==0)
//                {
//                    strcpy(type,"N.3.h2.pi=0");
//                }
//                else
//                {
//                    strcpy(type,"N.3.h2.pi=1");
//                }
//            }
//            else
//            {
//                strcpy(type,"N.3.h2.x");
//            }
//        }
//        else if(group.num_nonh==2)
//        {
//            if(group.num_hetero==0)
//            {
//                if(group.num_pi==0)
//                {
//                    strcpy(type,"N.3.h.pi=0");
//                }
//                else if(group.num_pi==1)
//                {
//                    strcpy(type,"N.3.h.pi>0");
//                }
//                else if(atom.ring==0)
//                {
//                    strcpy(type,"N.3.h.pi>0");
//                }
//                else strcpy(type,"N.3.h.ring");
//            }
//            else
//            {
//                if(group.num_pi<=1)
//                {
//                    strcpy(type,"N.3.h.x");
//                }
//                else if(atom.ring==0)
//                {
//                    strcpy(type,"N.3.h.x");
//                }
//                else
//                {
//                    strcpy(type,"N.3.h.x.ring");
//                }
//            }
//        }
//        else if(group.num_nonh==3)
//        {
//            if(group.num_hetero==0)
//            {
//                if(group.num_pi==0)
//                {
//                    strcpy(type,"N.3.pi=0");
//                }
//                else if(group.num_pi==1)
//                {
//                    strcpy(type,"N.3.pi>0");
//                }
//                else if(atom.ring==0)
//                {
//                    strcpy(type,"N.3.pi>0");
//                }
//                else
//                {
//                    strcpy(type,"N.3.ring");
//                }
//            }
//            else
//            {
//                if(group.num_pi<=1)
//                {
//                    strcpy(type,"N.3.x");
//                }
//                else if(atom.ring==0)
//                {
//                    strcpy(type,"N.3.x");
//                }
//                else
//                {
//                    strcpy(type,"N.3.x.ring");
//                }
//            }
//        }
//        else
//        {
//            strcpy(type,"N.3.unknown");
//        }
//    }
//
//    if(!strcmp(atom.type,"N.am"))
//    {
//        if(group.num_nonh==1) strcpy(type,"N.am.h2");
//        else if(group.num_nonh==2)
//        {
//            if(group.num_hetero==0) strcpy(type,"N.am.h");
//            else strcpy(type,"N.am.h.x");
//        }
//        else if(group.num_nonh==3)
//        {
//            if(group.num_hetero==0) strcpy(type,"N.am");
//            else strcpy(type,"N.am.x");
//        }
//        else strcpy(type,"N.am.unknown");
//    }
//
//    if(!strcmp(atom.type,"N.2"))
//    {
//        // N=C, N=S, N=P
//        if(group.db_type==1||group.db_type==4||group.db_type==5)
//        {
//            if(group.num_hetero==0)
//            {
//                if(group.num_pi==0)
//                {
//                    strcpy(type,"N.2.(=C).pi=0");
//                }
//                else
//                {
//                    strcpy(type,"N.2.(=C).pi=1");
//                }
//            }
//            else
//            {
//                if(group.num_pi==0)
//                {
//                    strcpy(type,"N.2.(=C).x.pi=0");
//                }
//                else
//                {
//                    strcpy(type,"N.2.(=C).x.pi=1");
//                }
//            }
//        }
//        else if(group.db_type==2)
//        {
//            if(group.num_hetero==0) strcpy(type,"N.2.(=N)");
//            else strcpy(type,"N.2.(=N).x");
//        }
//        else if(group.db_type==3)
//        {
//            if(group.num_nonh==2) strcpy(type,"N.2.o");
//            else if(group.num_nonh==3) strcpy(type,"N.2.o2");
//            else strcpy(type,"N.2.o");
//        }
//        else
//        {
//            strcpy(type,"N.2.unknown");
//        }
//    }
//
//    if(!strcmp(atom.type,"N.ar")) strcpy(type,"N.ar");
//
//    if(!strcmp(atom.type,"N.1")) strcpy(type,"N.1");
//
//    if(!strcmp(atom.type,"O.3"))
//    {
//        if(group.num_nonh==1)
//        {
//            if(group.num_hetero==0)
//            {
//                if(group.num_pi==0) strcpy(type,"O.3.h.pi=0");
//                else strcpy(type,"O.3.h.pi=1");
//            }
//            else
//            {
//                strcpy(type,"O.3.h.x");
//            }
//        }
//        else if(group.num_nonh==2)
//        {
//            if(group.num_hetero==0)
//            {
//                if(group.num_pi==0) strcpy(type,"O.3.pi=0");
//                else strcpy(type,"O.3.pi>0");
//            }
//            else
//            {
//                strcpy(type,"O.3.x");
//            }
//        }
//        else
//        {
//            strcpy(type,"O.3.unknown");
//        }
//    }
//
//    if(!strcmp(atom.type,"O.2")) strcpy(type,"O.2");
//
//    if(!strcmp(atom.type,"O.co2")) strcpy(type,"O.co2");
//
//    if(!strcmp(atom.type,"S.3"))
//    {
//        if(group.num_nonh==1) strcpy(type,"S.3.h");
//        else if(group.num_nonh==2) strcpy(type,"S.3");
//        else strcpy(type,"S.3.unknown");
//    }
//
//    if(!strcmp(atom.type,"S.2")) strcpy(type,"S.2");
//
//    if(!strcmp(atom.type,"S.o")) strcpy(type,"S.o");
//
//    if(!strcmp(atom.type,"S.o2")) strcpy(type,"S.o2");
//
//    if(!strcmp(atom.type,"P.3"))
//    {
//        if(group.db_type==3) strcpy(type,"P.3.(=O)");
//        else if(group.db_type==4) strcpy(type,"P.3.(=S)");
//        else strcpy(type,"P.3.unknown");
//    }
//
//    if(!strcmp(atom.type,"F"))
//    {
//        if(group.num_pi==0) strcpy(type,"F.pi=0");
//        else if(group.num_pi==1) strcpy(type,"F.pi=1");
//        else strcpy(type,"F.unknown");
//    }
//
//    if(!strcmp(atom.type,"Cl"))
//    {
//        if(group.num_pi==0) strcpy(type,"Cl.pi=0");
//        else if(group.num_pi==1) strcpy(type,"Cl.pi=1");
//        else strcpy(type,"Cl.unknown");
//    }
//
//    if(!strcmp(atom.type,"Br"))
//    {
//        if(group.num_pi==0) strcpy(type,"Br.pi=0");
//        else if(group.num_pi==1) strcpy(type,"Br.pi=1");
//        else strcpy(type,"Br.unknown");
//    }
//
//    if(!strcmp(atom.type,"I"))
//    {
//        if(group.num_pi==0) strcpy(type,"I.pi=0");
//        else if(group.num_pi==1) strcpy(type,"I.pi=1");
//        else strcpy(type,"I.unknown");
//    }
//
//    if(!strcmp(atom.type,"Si")) strcpy(type,"Si");
//
//    if(!strcmp(type,"Un")) return FALSE;
//    return TRUE;
//}
//
//int xscore::Connection_1_2_Check(int id1, int id2) {
//    int i;
//
//    if(id1==id2) return FALSE;
//
//    for(i=0;i<ligand_atoms[id2-1].num_neib;i++)
//    {
//        if(id1==ligand_atoms[id2-1].neib[i]) return TRUE;
//        else continue;
//    }
//
//    return FALSE;
//}
//
//int xscore::Connection_1_3_Check(int id1, int id2) {
//    int i,j;
//
//    if(id1==id2) return FALSE;
//    if(Connection_1_2_Check(id1,id2)==TRUE) return FALSE;
//
//    id1--; id2--;
//
//    for(i=0;i<ligand_atoms[id1].num_neib;i++)
//        for(j=0;j<ligand_atoms[id2].num_neib;j++)
//        {
//            if(ligand_atoms[id1].neib[i]==atom[id2].neib[j]) return TRUE;
//            else continue;
//        }
//
//    return FALSE;
//}
//
//int xscore::Connection_1_4_Check(int id1, int id2) {
//    int i,j;
//
//    if(id1==id2) return FALSE;
//    if(Connection_1_2_Check(id1,id2)==TRUE) return FALSE;
//    if(Connection_1_3_Check(id1,id2)==TRUE) return FALSE;
//
//    id1--; id2--;
//
//    for(i=0;i<ligand_atoms[id1].num_neib;i++)
//        for(j=0;j<ligand_atoms[id2].num_neib;j++)
//        {
//            if(Connection_1_2_Check(ligand_atoms[id1].neib[i],
//                                    ligand_atoms[id2].neib[j])==TRUE) return TRUE;
//            else continue;
//        }
//
//    return FALSE;
//}
//
//int xscore::Connection_1_5_Check(int id1, int id2) {
//    int i,j;
//
//    if(id1==id2) return FALSE;
//    if(Connection_1_2_Check(id1,id2)==TRUE) return FALSE;
//    if(Connection_1_3_Check(id1,id2)==TRUE) return FALSE;
//    if(Connection_1_4_Check(id1,id2)==TRUE) return FALSE;
//
//    id1--; id2--;
//
//    for(i=0;i<ligand_atoms[id1].num_neib;i++)
//        for(j=0;j<ligand_atoms[id2].num_neib;j++)
//        {
//            if(Connection_1_3_Check(ligand_atoms[id1].neib[i],
//                                    ligand_atoms[id2].neib[j])==TRUE) return TRUE;
//            else continue;
//        }
//
//    return FALSE;
//}
//
//int xscore::Connection_1_6_Check(int id1, int id2) {
//    int i,j;
//
//    if(id1==id2) return FALSE;
//    if(Connection_1_2_Check(id1,id2)==TRUE) return FALSE;
//    if(Connection_1_3_Check(id1,id2)==TRUE) return FALSE;
//    if(Connection_1_4_Check(id1,id2)==TRUE) return FALSE;
//    if(Connection_1_5_Check(id1,id2)==TRUE) return FALSE;
//
//    id1--; id2--;
//
//    for(i=0;i<ligand_atoms[id1].num_neib;i++)
//        for(j=0;j<ligand_atoms[id2].num_neib;j++)
//        {
//            if(Connection_1_4_Check(ligand_atoms[id1].neib[i],
//                                    ligand_atoms[id2].neib[j])==TRUE) return TRUE;
//            else continue;
//        }
//
//    return FALSE;
//}
//
//int xscore::Hydrophobic_Neighbor_Check(int id) {
//    int i,mark;
//
//    id--;
//
//    mark=TRUE;
//
//    for(i=0;i<ligand_atoms.size();i++)
//    {
//        if(i==id) continue;
//        else if(!strcmp(ligand_atoms[i].type,"F")||
//                !strcmp(ligand_atoms[i].type,"Cl")||
//                !strcmp(ligand_atoms[i].type,"Br")||
//                !strcmp(ligand_atoms[i].type,"I")||
//                !strcmp(ligand_atoms[i].type,"Si")||
//                ligand_atoms[i].type[0]=='N'||
//                ligand_atoms[i].type[0]=='O'||
//                ligand_atoms[i].type[0]=='S'||
//                ligand_atoms[i].type[0]=='P')
//        {
//            if(Connection_1_2_Check(ligand_atoms[id].id,ligand_atoms[i].id)==TRUE)
//            {mark=FALSE; break;}
//            else if(Connection_1_3_Check(ligand_atoms[id].id,ligand_atoms[i].id)==TRUE)
//            {mark=FALSE; break;}
//            else if(Connection_1_4_Check(ligand_atoms[id].id,ligand_atoms[i].id)==TRUE)
//            {mark=FALSE; break;}
//            else continue;
//        }
//        else continue;
//    }
//
//    return mark;
//}
//
//float xscore::Count_Hydrophobic_Carbon(int flag) {
//    int i,num;
//
//    num=0;
//
//    for(i=0;i<ligand_atoms.size();i++)
//    {
//        if(strcmp(ligand_atoms[i].type,"C.3")&&
//           strcmp(ligand_atoms[i].type,"C.2")) continue;
//        else if(Hydrophobic_Neighbor_Check(ligand_atoms[i].id)
//                ==FALSE) continue;
//
//        num++;
//
//        if(flag) ligand_atoms[i].logp+=(LOGP_HYDROPHOBIC_CARBON);
//        else continue;
//    }
//
//    if(num>=10) num/=2;
//
//    return (float)num;
//}
//
//
//int xscore::Adjacent_Ring_Check(int id) {
//    int i,num,tmp,mark;
//
//    id--;
//
//    mark=FALSE; num=ligand_atoms[id].num_neib;
//
//    for(i=0;i<num;i++)
//    {
//        tmp=ligand_atoms[id].neib[i];
//        if(atom[tmp-1].ring!=0) {mark=tmp; break;}
//        else continue;
//    }
//
//    return mark;
//}
//
//float xscore::Count_Internal_HBond(int flag) {
//    int i,j,mark1,mark2;
//    float num;
//    int *record;
//
//    record=new int[ligand_atoms.size()];
//    if(record==NULL) Memory_Allocation_Error();
//
//    for(i=0;i<ligand_atoms.size();i++) record[i]=0;
//
//    num=0;
//
//    for(i=0;i<ligand_atoms.size();i++)
//    {
//        if(strcmp(ligand_atoms[i].hb,"D")&&strcmp(ligand_atoms[i].hb,"DA")) continue;
//        else if(ligand_atoms[i].ring!=0) continue; // not allowed in ring
//
//        if(Adjacent_Ring_Check(ligand_atoms[i].id)==FALSE) mark1=FALSE;
//        else mark1=TRUE;
//
//        for(j=0;j<ligand_atoms.size();j++)
//        {
//            if(i==j) continue;
//            else if(strcmp(ligand_atoms[j].hb,"A")&&strcmp(ligand_atoms[j].hb,"DA")) continue;
//            else if(!strcmp(ligand_atoms[j].type,"O.3")&&!strcmp(ligand_atoms[j].hb,"A")) continue;
//            else if(!strcmp(ligand_atoms[j].type,"N.2")) continue;
//            else if(!strcmp(ligand_atoms[j].type,"N.ar")) continue;
//            else if(atom[j].ring!=0) continue; // not in ring
//
//            if(Adjacent_Ring_Check(ligand_atoms[j].id)==FALSE) mark2=FALSE;
//            else mark2=TRUE;
//
//            if(mark1==TRUE&&mark2==TRUE)
//            {
//                if(Connection_1_4_Check(ligand_atoms[i].id,
//                                        ligand_atoms[j].id)==FALSE) continue;
//                else
//                {
//                    record[i]++; record[j]++;
//                }
//            }
//            else if(mark1==TRUE&&mark2==FALSE)
//            {
//                if(Connection_1_5_Check(ligand_atoms[i].id,
//                                        ligand_atoms[j].id)==FALSE) continue;
//                else
//                {
//                    record[i]++; record[j]++;
//                }
//            }
//            else if(mark2==FALSE&&mark2==TRUE)
//            {
//                if(Connection_1_5_Check(ligand_atoms[i].id,
//                                        ligand_atoms[j].id)==FALSE) continue;
//                else
//                {
//                    record[i]++; record[j]++;
//                }
//            }
//            else continue;
//        }
//    }
//
//    for(i=0;i<ligand_atoms.size();i++)
//    {
//        if(record[i]==0) continue;
//
//        num+=0.500;
//
//        if(flag) ligand_atoms[i].logp+=(LOGP_INTERNAL_HBOND/2.0);
//        else continue;
//    }
//
//    if(record) delete [] record;
//
//    return num;
//}
//
//float xscore::Count_Halogen_1_3_Pair(int flag) {
//    int i,j;
//    int num1,num2;
//
//    num1=num2=0;
//
//    for(i=0;i<ligand_atoms.size()-1;i++)
//    {
//        if(strcmp(ligand_atoms[i].type,"F")) continue;
//
//        for(j=i+1;j<ligand_atoms.size();j++)
//        {
//            if(strcmp(ligand_atoms[j].type,"F")) continue;
//            else if(Connection_1_3_Check(ligand_atoms[i].id,ligand_atoms[j].id)
//                    ==FALSE) continue;
//
//            num1++;
//
//            if(flag)
//            {
//                ligand_atoms[i].logp+=((LOGP_HALOGEN_PAIR)/2.0);
//                ligand_atoms[j].logp+=((LOGP_HALOGEN_PAIR)/2.0);
//            }
//            else continue;
//        }
//    }
//
//    for(i=0;i<ligand_atoms.size()-1;i++)
//    {
//        if(strcmp(ligand_atoms[i].type,"Cl")&&
//           strcmp(ligand_atoms[i].type,"Br")&&
//           strcmp(ligand_atoms[i].type,"I")) continue;
//
//        for(j=i+1;j<ligand_atoms.size();j++)
//        {
//            if(strcmp(ligand_atoms[j].type,"Cl")&&
//               strcmp(ligand_atoms[j].type,"Br")&&
//               strcmp(ligand_atoms[j].type,"I")) continue;
//            else if(Connection_1_3_Check(ligand_atoms[i].id,ligand_atoms[j].id)
//                    ==FALSE) continue;
//
//            num2++;
//
//            if(flag)
//            {
//                ligand_atoms[i].logp+=((LOGP_HALOGEN_PAIR)/2.0);
//                ligand_atoms[j].logp+=((LOGP_HALOGEN_PAIR)/2.0);
//            }
//            else continue;
//        }
//    }
//
//    return (float)(num1+num2);
//}
//
//float xscore::Count_Nar_1_4_Pair(int flag) {
//    int i,j,num,tmp1,tmp2,tmp3,tmp4;
//
//    num=0;
//
//    for(i=0;i<ligand_atoms.size()-1;i++)
//    {
//        if(strcmp(ligand_atoms[i].type,"N.ar")) continue;
//
//        tmp1=ligand_atoms[i].neib[0]; tmp2=ligand_atoms[i].neib[1];
//
//        for(j=i+1;j<ligand_atoms.size();j++)
//        {
//            if(strcmp(ligand_atoms[j].type,"N.ar")) continue;
//            else if(Connection_1_4_Check(ligand_atoms[i].id, ligand_atoms[j].id)
//                    ==FALSE) continue;
//            else
//            {
//                tmp3=ligand_atoms[j].neib[0]; tmp4=ligand_atoms[j].neib[1];
//
//                if(Connection_1_2_Check(tmp1,tmp3)==TRUE)
//                {
//                    if(Connection_1_2_Check(tmp2,tmp4)==TRUE)
//                    {
//                        num++;
//                        if(flag)
//                        {
//                            ligand_atoms[i].logp+=((LOGP_NAR_PAIR)/2.0);
//                            ligand_atoms[j].logp+=((LOGP_NAR_PAIR)/2.0);
//                        }
//                    }
//                    else continue;
//                }
//                else if(Connection_1_2_Check(tmp1,tmp4)==TRUE)
//                {
//                    if(Connection_1_2_Check(tmp2,tmp3)
//                       ==TRUE)
//                    {
//                        num++;
//                        if(flag)
//                        {
//                            ligand_atoms[i].logp+=((LOGP_NAR_PAIR)/2.0);
//                            ligand_atoms[j].logp+=((LOGP_NAR_PAIR)/2.0);
//                        }
//                    }
//                    else continue;
//                }
//                else continue;
//            }
//        }
//    }
//
//    return (float)num;
//}
//
//float xscore::Count_O3_1_4_Pair(int flag) {
//    int i,j;
//    float num;
//    int *record;
//
//    record=new int[ligand_atoms.size()];
//    if(record==NULL) Memory_Allocation_Error();
//
//    for(i=0;i<ligand_atoms.size();i++) record[i]=0;
//
//    num=0;
//
//    for(i=0;i<ligand_atoms.size()-1;i++)
//    {
//        if(strcmp(ligand_atoms[i].type,"O.3")) continue;
//        else if(!strcmp(ligand_atoms[i].hb,"DA")) continue;
//        else if(ligand_atoms[i].ring!=0) continue;
//        else if(Adjacent_Aromatic_Check(ligand_atoms[i].id)==FALSE) continue;
//
//        for(j=i+1;j<ligand_atoms.size();j++)
//        {
//            if(strcmp(ligand_atoms[j].type,"O.3")) continue;
//            else if(!strcmp(ligand_atoms[j].hb,"DA")) continue;
//            else if(ligand_atoms[j].ring!=0) continue;
//            else if(Adjacent_Aromatic_Check(ligand_atoms[j].id)==FALSE)
//                continue;
//            else if(Connection_1_4_Check(ligand_atoms[i].id,ligand_atoms[j].id)
//                    ==FALSE) continue;
//            else
//            {
//                record[i]++; record[j]++;
//            }
//        }
//    }
//
//    for(i=0;i<ligand_atoms.size();i++)
//    {
//        if(record[i]==0) continue;
//
//        num+=0.500;
//
//        if(flag) ligand_atoms[i].logp+=(LOGP_O3_PAIR/2.0);
//        else continue;
//    }
//
//    if(record) delete [] record;
//
//    return (float)num;
//}
//
//float xscore::Count_Acceptor_1_5_Pair(int flag) {
//    int i,j,num,tmp1,tmp2;
//
//    num=0;
//
//    for(i=0;i<ligand_atoms.size()-1;i++)
//    {
//        if(strcmp(ligand_atoms[i].hb,"A")) continue;
//        else if(!strcmp(ligand_atoms[i].type,"O.3")) continue;
//        else if(!strcmp(ligand_atoms[i].type,"N.2")) continue;
//        else if(!strcmp(ligand_atoms[i].type,"N.ar")) continue;
//
//        tmp1=ligand_atoms[i].neib[0]-1;
//        if(ligand_atoms[tmp1].type[0]=='S') continue;
//        else if(ligand_atoms[tmp1].type[0]=='P') continue;
//        else if(ligand_atoms[tmp1].ring!=0) continue;
//
//        for(j=i+1;j<ligand_atoms.size();j++)
//        {
//            if(strcmp(ligand_atoms[j].hb,"A")) continue;
//            else if(!strcmp(ligand_atoms[j].type,"O.3")) continue;
//            else if(!strcmp(ligand_atoms[j].type,"N.2")) continue;
//            else if(!strcmp(ligand_atoms[j].type,"N.ar")) continue;
//
//            tmp2=ligand_atoms[j].neib[0]-1;
//            if(ligand_atoms[tmp2].type[0]=='S') continue;
//            else if(ligand_atoms[tmp2].type[0]=='P') continue;
//            else if(ligand_atoms[tmp2].ring!=0) continue;
//
//            if(Connection_1_5_Check(ligand_atoms[i].id, ligand_atoms[j].id)
//               ==FALSE) continue;
//
//            num++;
//
//            if(flag)
//            {
//                ligand_atoms[i].logp+=((LOGP_ACCEPTOR_PAIR)/2.0);
//                ligand_atoms[j].logp+=((LOGP_ACCEPTOR_PAIR)/2.0);
//            }
//            else continue;
//        }
//    }
//
//    return (float)num;
//}
//
//int xscore::Adjacent_Aromatic_Check(int id) {
//    int i,num,tmp,mark;
//
//    id--;
//
//    mark=FALSE; num=ligand_atoms[id].num_neib;
//
//    for(i=0;i<num;i++)
//    {
//        tmp=ligand_atoms[id].neib[i];
//        if(!strcmp(ligand_atoms[tmp-1].type,"C.ar")) {mark=tmp; break;}
//        else continue;
//    }
//
//    return mark;
//}
//
//
//float xscore::Count_Salicylic_Acid(int flag) {
//    int i,j,num,tmp,mark;
//
//    num=0;
//
//    for(i=0;i<ligand_atoms.size();i++)
//    {
//        if(strcmp(ligand_atoms[i].type,"O.2")) continue;
//
//        tmp=ligand_atoms[i].neib[0]-1;
//        if(ligand_atoms[tmp].type[0]!='C') continue;
//        else if(ligand_atoms[tmp].ring!=0) continue;
//        else if(Adjacent_Aromatic_Check(ligand_atoms[tmp].id)==FALSE) continue;
//        mark=FALSE;
//
//        for(j=0;j<ligand_atoms.size();j++)
//        {
//            if(i==j) continue;
//            else if(strcmp(ligand_atoms[j].type,"O.3")) continue;
//            else if(ligand_atoms[j].ring!=0) continue;
//            else if(Connection_1_3_Check(ligand_atoms[i].id, ligand_atoms[j].id)
//                    ==FALSE) continue;
//            else {mark=TRUE; break;}
//        }
//
//        if(mark==FALSE) continue;
//
//        mark=FALSE;
//
//        for(j=0;j<ligand_atoms.size();j++)
//        {
//            if(i==j) continue;
//            else if(strcmp(ligand_atoms[j].type,"O.3")) continue;
//            else if(strcmp(ligand_atoms[j].hb,"DA")) continue;
//            else if(Adjacent_Aromatic_Check(ligand_atoms[j].id)==FALSE)
//                continue;
//            else if(Connection_1_5_Check(ligand_atoms[i].id, ligand_atoms[j].id)
//                    ==FALSE) continue;
//
//            num++;
//
//            if(flag)
//            {
//                ligand_atoms[i].logp+=((LOGP_SALICYLIC_ACID)/2.0);
//                ligand_atoms[j].logp+=((LOGP_SALICYLIC_ACID)/2.0);
//            }
//
//            break;
//        }
//
//        if(num!=0) break;
//    }
//
//    return (float)num;
//}
//
//float xscore::Count_Amino_Acid(int flag) {
//    int i,j,tmp,num,mark;
//
//    num=0;
//
//    for(i=0;i<ligand_atoms.size();i++)
//    {
//        if(strcmp(ligand_atoms[i].type,"O.2")) continue;
//
//        tmp=ligand_atoms[i].neib[0]-1;
//        if(ligand_atoms[tmp].type[0]!='C') continue;
//        else if(ligand_atoms[tmp].ring!=0) continue;
//
//        mark=FALSE;
//
//        for(j=0;j<ligand_atoms.size();j++)
//        {
//            if(i==j) continue;
//            else if(strcmp(ligand_atoms[j].type,"O.3")) continue;
//            else if(strcmp(ligand_atoms[j].hb,"DA")) continue;
//            else if(ligand_atoms[j].ring!=0) continue;
//            else if(Connection_1_3_Check(ligand_atoms[i].id, ligand_atoms[j].id)
//                    ==FALSE) continue;
//            else {mark=TRUE; break;}
//        }
//
//        if(mark==FALSE) continue;
//
//        for(j=0;j<ligand_atoms.size();j++)
//        {
//            if(strcmp(ligand_atoms[j].xtype,"N.3.h2.pi=0")) continue;
//            else if(Connection_1_4_Check(ligand_atoms[i].id, ligand_atoms[j].id)
//                    ==FALSE) continue;
//
//            num++;
//
//            if(flag)
//            {
//                ligand_atoms[i].logp+=((LOGP_AMINO_ACID)/2.0);
//                ligand_atoms[j].logp+=((LOGP_AMINO_ACID)/2.0);
//            }
//
//            break;
//        }
//
//        if(num!=0) break;
//
//        for(j=0;j<ligand_atoms.size();j++)
//        {
//            if(strcmp(ligand_atoms[j].type,"N.ar")) continue;
//            else if(Connection_1_4_Check(ligand_atoms[i].id, ligand_atoms[j].id)
//                    ==FALSE) continue;
//
//            num++;
//
//            if(flag)
//            {
//                ligand_atoms[i].logp+=((LOGP_AMINO_ACID)/2.0);
//                ligand_atoms[j].logp+=((LOGP_AMINO_ACID)/2.0);
//            }
//
//            break;
//        }
//
//        if(num!=0) break;
//    }
//
//    return (float)num;
//}
//
//float xscore::Count_Sulfonic_Acid(int flag) {
//    int i,j,num,tmp1,tmp2;
//
//    num=0;
//
//    for(i=0;i<ligand_atoms.size();i++)
//    {
//        if(strcmp(ligand_atoms[i].type,"S.o2")) continue;
//        else if(ligand_atoms[i].ring!=0) continue;
//
//        tmp1=Adjacent_Aromatic_Check(ligand_atoms[i].id);
//        if(tmp1==FALSE) continue;
//
//        for(j=0;j<ligand_atoms.size();j++)
//        {
//            if(strcmp(ligand_atoms[j].xtype,"N.3.h2.pi=1")) continue;
//
//            tmp2=Adjacent_Aromatic_Check(ligand_atoms[j].id);
//            if(tmp2==FALSE) continue;
//
//            if(Connection_1_6_Check(ligand_atoms[i].id,ligand_atoms[j].id)==FALSE) continue;
//            else if(Connection_1_4_Check(tmp1,tmp2)==FALSE) continue;
//
//            num++;
//
//            if(flag)
//            {
//                ligand_atoms[i].logp+=(LOGP_SULFONIC_ACID/2.0);
//                ligand_atoms[j].logp+=(LOGP_SULFONIC_ACID/2.0);
//            }
//
//            break;
//        }
//    }
//
//    return (float)num;
//}
//
//
//float xscore::Calculate_LogP() {
//    int flag;
//    float xlogp;
//    char type[20];
//    LogP_Factor logp_factor[10];
//
//    strcpy(logp_factor[0].symbol,"Hydrophobic carbon");
//    strcpy(logp_factor[1].symbol,"Internal H-bond");
//    strcpy(logp_factor[2].symbol,"Halogen 1-3 pair");
//    strcpy(logp_factor[3].symbol,"Aromatic nitrogen 1-4 pair");
//    strcpy(logp_factor[4].symbol,"Ortho sp3 oxygen pair");
//    strcpy(logp_factor[5].symbol,"Acceptor 1-5 pair");
//    strcpy(logp_factor[6].symbol,"Paralleled donor pair");
//    strcpy(logp_factor[7].symbol,"Alpha amino acid");
//    strcpy(logp_factor[8].symbol,"Salicylic acid");
//    strcpy(logp_factor[9].symbol,"P-amino sulfonic acid");
//
//    logp_factor[0].coeff = 0.211;
//    logp_factor[1].coeff = 0.429;
//    logp_factor[2].coeff = 0.137;
//    logp_factor[3].coeff = 0.485;
//    logp_factor[4].coeff =-0.268;
//    logp_factor[5].coeff = 0.580;
//    logp_factor[6].coeff =-0.423;
//    logp_factor[7].coeff =-2.166;
//    logp_factor[8].coeff = 0.554;
//    logp_factor[9].coeff =-0.501;
//
//    xlogp=0.000;
//
//    for(int i=0;i<ligand_atoms.size();i++)
//    {
//        if(ligand_atoms[i].valid==0) continue;
//
//        Get_XLOGP_Type(ligand_atoms[i].id,type);
//        ligand_atoms[i].logp=Get_Atom_LogP(type);
//
//        xlogp+=ligand_atoms[i].logp;
//    }
//
//    for(int i=0;i<10;i++) logp_factor[i].num=0;
//
//    // this flag determines if correction factors are also distributed
//    // onto each atoms or not
//
//    flag=0;
//
//    logp_factor[0].num=Count_Hydrophobic_Carbon(flag);
//    logp_factor[1].num=Count_Internal_HBond(flag);
//    logp_factor[2].num=Count_Halogen_1_3_Pair(flag);
//    logp_factor[3].num=Count_Nar_1_4_Pair(flag);
//    logp_factor[4].num=Count_O3_1_4_Pair(flag);
//    logp_factor[5].num=Count_Acceptor_1_5_Pair(flag);
//    logp_factor[6].num=0.000;
//    logp_factor[7].num=Count_Amino_Acid(flag);
//    logp_factor[8].num=Count_Salicylic_Acid(flag);
//    logp_factor[9].num=Count_Sulfonic_Acid(flag);
//
//    for(int i=0;i<10;i++) xlogp+=(logp_factor[i].num*logp_factor[i].coeff);
//
//    return xlogp;
//}
//
//float xscore::Calculate_LogP_protein() {
//    int flag;
//    float xlogp;
//    char type[20];
//    LogP_Factor logp_factor[10];
//
//    strcpy(logp_factor[0].symbol,"Hydrophobic carbon");
//    strcpy(logp_factor[1].symbol,"Internal H-bond");
//    strcpy(logp_factor[2].symbol,"Halogen 1-3 pair");
//    strcpy(logp_factor[3].symbol,"Aromatic nitrogen 1-4 pair");
//    strcpy(logp_factor[4].symbol,"Ortho sp3 oxygen pair");
//    strcpy(logp_factor[5].symbol,"Acceptor 1-5 pair");
//    strcpy(logp_factor[6].symbol,"Paralleled donor pair");
//    strcpy(logp_factor[7].symbol,"Alpha amino acid");
//    strcpy(logp_factor[8].symbol,"Salicylic acid");
//    strcpy(logp_factor[9].symbol,"P-amino sulfonic acid");
//
//    logp_factor[0].coeff = 0.211;
//    logp_factor[1].coeff = 0.429;
//    logp_factor[2].coeff = 0.137;
//    logp_factor[3].coeff = 0.485;
//    logp_factor[4].coeff =-0.268;
//    logp_factor[5].coeff = 0.580;
//    logp_factor[6].coeff =-0.423;
//    logp_factor[7].coeff =-2.166;
//    logp_factor[8].coeff = 0.554;
//    logp_factor[9].coeff =-0.501;
//
//    xlogp=0.000;
//
//    for(int i=0;i<protein_atoms.size();i++)
//    {
//        if(protein_atoms[i].valid==0) continue;
//
//        Get_XLOGP_Type(protein_atoms[i].id,type);
//        protein_atoms[i].logp=Get_Atom_LogP(type);
//
//        xlogp+=protein_atoms[i].logp;
//    }
//
//    for(int i=0;i<10;i++) logp_factor[i].num=0;
//
//    // this flag determines if correction factors are also distributed
//    // onto each atoms or not
//
//    flag=0;
///////////////////要改为protein_atoms
//    logp_factor[0].num=Count_Hydrophobic_Carbon(flag);
//    logp_factor[1].num=Count_Internal_HBond(flag);
//    logp_factor[2].num=Count_Halogen_1_3_Pair(flag);
//    logp_factor[3].num=Count_Nar_1_4_Pair(flag);
//    logp_factor[4].num=Count_O3_1_4_Pair(flag);
//    logp_factor[5].num=Count_Acceptor_1_5_Pair(flag);
//    logp_factor[6].num=0.000;
//    logp_factor[7].num=Count_Amino_Acid(flag);
//    logp_factor[8].num=Count_Salicylic_Acid(flag);
//    logp_factor[9].num=Count_Sulfonic_Acid(flag);
//
//    for(int i=0;i<10;i++) xlogp+=(logp_factor[i].num*logp_factor[i].coeff);
//
//    return xlogp;
//}
//float xscore::calculate_HM() {
//    float asum,sum,tmp,d,d1,d2,total,cutoff;
//    Calculate_LogP();
//    Calculate_LogP_protein();
//    for(int i=0;i<ligand_atoms.size();i++) ligand_atoms[i].score=0.000;
//
//    sum=0.000;
//
//    for(int i=0;i<ligand_atoms.size();i++)
//    {
//        if(ligand_atoms[i].valid<=0) continue;
//        else if(ligand_atoms[i].type[0]=='H') continue;
//        else if(strcmp(ligand_atoms[i].hb,"H")) continue;
//        else if(ligand_atoms[i].logp<=0.00) continue;
//
//        total=0.000; asum=0.000;
//
//        for(int j=0;j<protein_atoms.size();j++)
//        {
//            if(protein_atoms[j].valid!=2) continue;
//            else if(protein_atoms[j].type[0]=='H') continue;
//            else if(!strcmp(protein_atoms[j].type,"O.w")) continue;
//
//            d = sqrt(vec_distance_sqr(ligand_atoms[i].coords, protein_atoms[j].coords));
//            if(d>cutoff) continue;
//
//            d1=ligand_atoms[i].r+protein_atoms[j].r+0.50;
//            d2=ligand_atoms[i].r+protein_atoms[j].r+2.20;
//
//            if(d<d1) tmp=1.000;
//            else if(d<d2) tmp=(1/(d1-d2))*(d-d2);
//            else tmp=0.000;
//
//            total+=(protein_atoms[j].logp*tmp);
//        }
//
//        if(ligand_atoms[i].logp>=0.50) asum=ligand_atoms[i].logp;
//        else if(total>-0.50) asum=ligand_atoms[i].logp;
//        else asum=0.000;
//
//        ligand_atoms[i].score=asum; sum+=asum;
//    }
//
//    return sum;
//}
//
//
//
////HS
//void xscore::Read_SURFACE_DEF(char *filename)
//{
//    FILE *fp;
//    int i,count;
//    char buf[256],head[256];
//    Dot tmp_dot;
//
//    if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);
//
//    // first, determine num_sdot_type;
//
//    count=0;
//
//    for(;;)
//    {
//        if(fgets(buf,256,fp)==NULL) break;
//        else if(buf[0]=='#') continue;
//        else if(Blank_Line_Check(buf)==true) continue;
//        else
//        {
//            sscanf(buf,"%s", head);
//            if(strcmp(head,"DOTSET")) continue;
//            else count++;
//        }
//    }
//
//    rewind(fp);
//
//    num_sdot_type=count;
//
//    sdot=new DotSet[num_sdot_type];
//    if(sdot==NULL) Memory_Allocation_Error();
//
//    // now read pre-calculated volume dots
//
//    count=0;
//
//    for(;;)
//    {
//        if(fgets(buf,256,fp)==NULL) break;
//        else if(Blank_Line_Check(buf)==true) continue;
//        else if(buf[0]=='#') continue;
//
//        sscanf(buf,"%s",head);
//
//        if(!strcmp(head,"END")) break;
//        else if(!strcmp(head,"DOTSET"))
//        {
//            sscanf(buf,"%*s%f%d%f%f",
//                   &sdot[count].r,
//                   &sdot[count].num_dot,
//                   &sdot[count].unit,
//                   &sdot[count].total);
//
//            strcpy(sdot[count].type,"Un");
//
//            sdot[count].dot.clear();
//
//            for(i=0;i<sdot[count].num_dot;i++)
//            {
//                fgets(buf,256,fp);
//                sscanf(buf,"%f%f%f",
//                       &tmp_dot.coor[0],
//                       &tmp_dot.coor[1],
//                       &tmp_dot.coor[2]);
//
//                tmp_dot.valid=1;
//                tmp_dot.unit=sdot[count].unit;
//                strcpy(tmp_dot.type,"Un");
//
//                sdot[count].dot.push_back(tmp_dot);
//            }
//
//            if(sdot[count].dot.size()!=sdot[count].num_dot)
//                Read_File_Error(filename);
//
//            count++;
//        }
//    }
//
//    fclose(fp);
//
//    if(count!=num_sdot_type) Read_File_Error(filename);
//
//    return;
//}
//
//
//void xscore::Generate_Surface_Dots(float probe_r) {
//    int i,j,k;
//    bool mark;
//    float d,dd,dmin;
//    DotSet tmp_set;
//    Dot tmp_dot;
//
//    sur_dot.clear();
//
//    int *ligand_check_list;
//
//    ligand_check_list=new int[ligand_atoms.size()];
//    if(ligand_check_list==NULL) Memory_Allocation_Error();
//
//    for(i=0;i<ligand_atoms.size();i++)
//    {
//        if(ligand_atoms[i].valid<=0) continue;
//        else if(!strcmp(ligand_atoms[i].type,"H")) continue;////////////////////
//
//        for(j=0;j<ligand_atoms.size();j++)
//        {
//            if(i==j||ligand_atoms[j].valid<=0)
//            {
//                ligand_check_list[j]=0; continue;
//            }
//            else if(!strcmp(ligand_atoms[j].type,"H"))/////////////////
//            {
//                ligand_check_list[j]=0; continue;
//            }
//
//            d = sqrt(vec_distance_sqr(ligand_atoms[i].coords, ligand_atoms[j].coords));
//
//            if(d>(ligand_atoms[i].r+ligand_atoms[j].r+2*probe_r))
//            {
//                ligand_check_list[j]=0;
//            }
//            else
//            {
//                ligand_check_list[j]=1;
//            }
//        }
//
//        tmp_set=Get_Surface_Dot(ligand_atoms[i],probe_r);
//
//        for(j=0;j<tmp_set.num_dot;j++)
//        {
//            // check whether this dot is on the surface or not
//
//            mark=true; dmin=1000.0;
//
//            for(k=0;k<ligand_atoms.size();k++)
//            {
//                if(ligand_check_list[k]==0) continue;
//
//                d = sqrt(vec_distance_sqr(tmp_set.dot[j].coor, ligand_atoms[k].coords));
//                dd=d/(ligand_atoms[k].r+probe_r)-1.00;
//
//                if(dd<0.000) {mark=false; break;}
//                else if(dd<dmin) {dmin=dd; continue;}
//                else continue;
//            }
//
//            if(mark==false) tmp_set.dot[j].valid=0;
//            else if(dmin>=0.10) tmp_set.dot[j].valid=1; // regular
//            else tmp_set.dot[j].valid=2;	// dots at edge
//        }
//
//        // now record the surface dots of the current atom
//
//        for(j=0;j<tmp_set.num_dot;j++)
//        {
//            if(tmp_set.dot[j].valid==0) continue;
//
//            tmp_dot=tmp_set.dot[j]; tmp_dot.valid=i+1;
//
//            if(tmp_set.dot[j].valid==1)
//            {
//                tmp_dot.unit=tmp_set.unit;
//            }
//            else  // correct the overlapping of edging dots
//            {
//                tmp_dot.unit=tmp_set.unit*0.500;
//            }
//
//            sur_dot.push_back(tmp_dot);
//        }
//    }
//
//    if(ligand_check_list) delete [] ligand_check_list;
//
//    return;
//}
//
//DotSet xscore::Get_Surface_Dot(const atom &atom, float r) {
//
//    float R;
//    DotSet tmp_set;
//
//    R=atom.R+r;/////////////////////
//
//    tmp_set=Get_Surface_Dot(R,atom.coords[0],atom.coords[1],atom.coords[2]);
//
//    strcpy(tmp_set.type, atom.type);
//
//    for(int i=0;i<tmp_set.num_dot;i++)
//    {
//        tmp_set.dot[i].valid=atom.id;
//        strcpy(tmp_set.dot[i].type, atom.type);
//    }
//
//    return tmp_set;
//}
//
//DotSet xscore::Get_Surface_Dot(float R, float x, float y, float z) {
//    int i,num;
//    bool mark;
//    float tmp,theta,phi,theta_step,phi_step;
//    float r,d,total;
//    DotSet tmp_set;
//    Dot tmp_dot;
//
//    // check the pre-calculated surface dot sets
//
//    mark=false;
//
//    for(i=0;i<num_sdot_type;i++)
//    {
//        if(fabs(R-sdot[i].r)>0.025) continue;
//        else {tmp_set=sdot[i]; mark=true; break;}
//    }
//
//    if(mark==true)
//    {
//        for(i=0;i<tmp_set.num_dot;i++)
//        {
//            tmp_set.dot[i].coor[0]+=x;
//            tmp_set.dot[i].coor[1]+=y;
//            tmp_set.dot[i].coor[2]+=z;
//        }
//        return tmp_set;
//    }
//
//    // if it is not pre-calculated, calculate it now
//
//    total=(4.000*pi*R*R);
//
//    // d=0.500;		// spacing between two dots
//
//    d=sqrt(total/300);
//    if(d<0.500) d=0.500;
//
//    tmp=(int)(pi*R/d+0.500); theta_step=pi/tmp;
//
//    num=0; tmp_set.dot.clear();
//
//    for(theta=0.00;theta<pi;theta+=theta_step)
//    {
//        r=R*sin(theta);
//        tmp=(int)(2*pi*r/d+0.500); phi_step=2*pi/tmp;
//
//        for(phi=0.00;phi<(2*pi);phi+=phi_step)
//        {
//            tmp_dot.coor[0]=r*cos(phi);
//            tmp_dot.coor[1]=r*sin(phi);
//            tmp_dot.coor[2]=R*cos(theta);
//            tmp_dot.valid=1;
//            strcpy(tmp_dot.type,"Un");
//            tmp_set.dot.push_back(tmp_dot);
//            num++;
//        }
//    }
//
//    strcpy(tmp_set.type, "Un"); tmp_set.r=R;
//    tmp_set.num_dot=num; tmp_set.unit=total/num;
//    tmp_set.total=total;
//
//    for(i=0;i<tmp_set.num_dot;i++)
//    {
//        tmp_set.dot[i].unit=tmp_set.unit;
//        tmp_set.dot[i].coor[0]+=x;
//        tmp_set.dot[i].coor[1]+=y;
//        tmp_set.dot[i].coor[2]+=z;
//    }
//
//    return tmp_set;
//}
//float xscore::Atom_Buried_Surface(int id, float &total, float &buried){
//    int num;
//    bool mark;
//    float d,ratio;
//    int *atom_check_list;
//
//    atom_check_list=new int[protein_atoms.size()];
//    if(atom_check_list==NULL) Memory_Allocation_Error();
//
//    total=buried=0.000;
//
//    for(int j=0;j<protein_atoms.size();j++) atom_check_list[j]=0;
//
//    for(int j=0;j<protein_atoms.size();j++)
//    {
//        if(protein_atoms[j].valid!=2) continue;
//        else if(protein_atoms[j].type[0]=='H') continue;////////////////
//        else if(!strcmp(protein_atoms[j].type,"O.w")) continue;///////////////////////
//
//        d = sqrt(vec_distance_sqr(ligand_atoms[id].coords, protein_atoms[j].coords));
//
//        if(d>(ligand_atoms[id].r+protein_atoms[j].r+2*WATER_R)) continue;
//        else atom_check_list[j]=1;
//    }
//
//    num=this->sur_dot.size();//////////////
//
//    for(int j=0;j<num;j++)
//    {
//        if(this->sur_dot[j].valid!=id) continue;
//
//        // check if this dot is buried
//
//        mark=false;
//
//        for(int k=0;k<protein_atoms.size();k++)
//        {
//            if(atom_check_list[k]==0) continue;
//
//
//            d = sqrt(vec_distance_sqr(sur_dot[j].coor, protein_atoms[k].coords));
//            if(d>(protein_atoms[k].r+WATER_R)) continue;
//            else {mark=true; break;}
//        }
//
//        total+=this->sur_dot[j].unit;
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
//}
//
//
//float xscore::calculate_HS() {
//    float total,buried,sum = 0.00;
//
//    for(int i = 0; i < ligand_atoms.size(); i++){
//        ligand_atoms[i].score = 0.00;
//    }
//
//    // then get the buried surface atom by atom
//
//    for(int i=0;i<ligand_atoms.size();i++)
//    {
//        if(ligand_atoms[i].valid<=0) continue;
//        else if(ligand_atoms[i].type[0]=='H') continue;
//        else if(!strcmp(ligand_atoms[i].hb,"DH")) continue;
//        else if(!strcmp(ligand_atoms[i].hb,"D")) continue;
//        else if(!strcmp(ligand_atoms[i].hb,"A")) continue;
//        else if(!strcmp(ligand_atoms[i].hb,"DA")) continue;
//        else if(!strcmp(ligand_atoms[i].hb,"P")) continue;
//        else if(!strcmp(ligand_atoms[i].hb,"N")) continue;
//
//        this->Atom_Buried_Surface(i,total,buried);///////////////
//
//        sum+=buried; ligand_atoms[i].score=buried;
//
//        // sum+=(buried*this->atom[i].logp);
//        // this option works *slightly* better
//    }
//
//    return sum;
//
//}
//
//
//int xscore::Get_HBond_Pair_PL(HBond *candidate,bool sb_flag) {
//    int num,type;
//    bool sb;
//    float d,cutoff;
//    HBond tmp_candidate;
//
//    num=0; cutoff=5.00;
//
//    // note that the following H-bond algorithm is not based on any
//    // explicit hydrogen atom at all, and this is the right thing to do.
//
//    for(int i=0;i<ligand_atoms.size();i++)
//    {
//        if(ligand_atoms[i].valid<=0) continue;
//        else if(!strcmp(ligand_atoms[i].type,"H")) continue;
//        else if(!strcmp(ligand_atoms[i].hb,"N")) continue;
//        else if(!strcmp(ligand_atoms[i].hb,"H")) continue;
//        else if(!strcmp(ligand_atoms[i].hb,"P")) continue;
//        else if(!strcmp(ligand_atoms[i].hb,"DH")) continue;
//
//        for(int j=0;j<protein_atoms.size();j++)
//        {
//            if(protein_atoms[j].valid!=2) continue;
//            else if(!strcmp(protein_atoms[j].type,"H")) continue;
//            else if(!strcmp(protein_atoms[j].type,"O.w")) continue;
//            else if(!strcmp(protein_atoms[j].hb,"H")) continue;
//            else if(!strcmp(protein_atoms[j].hb,"P")) continue;
//            else if(!strcmp(protein_atoms[j].hb,"N")) continue;
//
//            // determine the type of this H-bond first
//            // type=0, no H-bond
//            // type=1, latom is the donor, patom is the acceptor;
//            // type=2, latom is the acceptor, patom is the donor;
//            // type=3, latom bound with metal ion;
//
//            if(!strcmp(ligand_atoms[i].hb,"D"))
//            {
//                if(!strcmp(protein_atoms[j].hb,"D")) type=0;
//                else if(!strcmp(protein_atoms[j].hb,"A")) type=1;
//                else if(!strcmp(protein_atoms[j].hb,"DA")) type=1;
//                else if(!strcmp(protein_atoms[j].hb,"M")) type=0;
//                else type=0;
//            }
//            else if(!strcmp(ligand_atoms[i].hb,"A"))
//            {
//                if(!strcmp(protein_atoms[j].hb,"D")) type=2;
//                else if(!strcmp(protein_atoms[j].hb,"A")) type=0;
//                else if(!strcmp(protein_atoms[j].hb,"DA")) type=2;
//                else if(!strcmp(protein_atoms[j].hb,"M")) type=3;
//                else type=0;
//            }
//            else if(!strcmp(ligand_atoms[i].hb,"DA"))
//            {
//                if(!strcmp(protein_atoms[j].hb,"A")) type=1;
//                else if(!strcmp(protein_atoms[j].hb,"D")) type=2;
//                else if(!strcmp(protein_atoms[j].hb,"DA")) type=1;
//                else if(!strcmp(protein_atoms[j].hb,"M")) type=3;
//                else type=0;
//            }
//            else type=0;
//
//            if(type==0) continue;  // no H-bond
//
//            // a crude distance check
//
//            d = sqrt(vec_distance_sqr(ligand_atoms[i].coords, protein_atoms[j].coords));
//            if(d>cutoff) continue;
//
//            // this section is used to differentiate HB and SB
//
//            if(!strcmp(ligand_atoms[i].type,"O.co2")&&
//               !strcmp(protein_atoms[j].type,"O.co2"))
//            {
//                sb=false;
//            }
//            else if((fabs(ligand_atoms[i].q)>0.01)&&
//                    (fabs(protein_atoms[j].q)>0.01))
//            {
//                sb=true;
//            }
//            else
//            {
//                sb=false;
//            }
//
//            // now handle this h-bond
//
//            tmp_candidate.Clear();
//            tmp_candidate.type=type;
//            tmp_candidate.sb=sb;
//            tmp_candidate.latom=i+1;
//            tmp_candidate.patom=j+1;
//
//            if(type==1)
//            {
//                tmp_candidate.D=ligand_atoms[i];
//                tmp_candidate.A=protein_atoms[j];
//            }
//            else if(type==2)
//            {
//                tmp_candidate.A=ligand_atoms[i];
//                tmp_candidate.D=protein_atoms[j];
//            }
//            else if(type==3)
//            {
//                tmp_candidate.A=ligand_atoms[i];
//                tmp_candidate.D=protein_atoms[j];
//            }
//
//            // if need to treat charged H-bonds differently
//
//            if(sb_flag==true&&sb==true) tmp_candidate.Value_SBond();
//            else tmp_candidate.Value_HBond_2();
//
//            if(fabs(tmp_candidate.score)>0.000)
//            {
//                candidate[num]=tmp_candidate; num++;
//            }
//            else continue;
//        }
//    }
//
//    return num;
//}
//
//
//float Angle_of_Two_Vectors(const float v1[3], const float v2[3])
//{
//    double angle;
//    double l1,l2,tmp1,tmp2;
//
//    l1=sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
//    l2=sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);
//
//    tmp1=v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
//    tmp2=l1*l2;
//
//    angle=acos(tmp1/tmp2); angle=angle/(pi)*180.0;
//
//    return (float)angle;  // return angle in degree, 0-180
//}
//
//Group xscore::Find_A_Group(int atom_id) {
//    int i,j,id,num;
//    bool mark;
//    Group group;
//
//    // define the center
//
//    group.center=ligand_atoms[atom_id-1];
//
//    // find the center's neighbours
//
//    num=group.center.num_neib;
//
//    for(i=0;i<num;i++)
//    {
//        group.neib[i]=ligand_atoms[group.center.neib[i]-1];
//        group.bond[i]=this->bond[group.center.bond[i]-1];
//    }
//
//    // count the necessary parameters
//
//    group.num_neib=group.center.num_neib;
//    group.num_h=0; group.num_nonh=0;
//
//    for(i=0;i<num;i++)
//    {
//        if(group.neib[i].type[0]=='H') group.num_h++;
//        else group.num_nonh++;
//    }
//
//    group.num_hetero=0;
//
//    for(i=0;i<group.num_nonh;i++)
//    {
//        if(!strcmp(group.neib[i].type,"F")) group.num_hetero++;
//        else if(!strcmp(group.neib[i].type,"Cl")) group.num_hetero++;
//        else if(!strcmp(group.neib[i].type,"Br")) group.num_hetero++;
//        else if(!strcmp(group.neib[i].type,"I")) group.num_hetero++;
//        else if(!strcmp(group.neib[i].type,"Si")) continue;
//        else if(group.neib[i].type[0]=='N') group.num_hetero++;
//        else if(group.neib[i].type[0]=='O') group.num_hetero++;
//        else if(group.neib[i].type[0]=='P') group.num_hetero++;
//        else if(group.neib[i].type[0]=='S') group.num_hetero++;
//        else continue;
//    }
//
//    group.db_type=0; group.num_db=0; group.num_tb=0;
//
//    for(i=0;i<group.num_nonh;i++)
//    {
//        if(!strcmp(group.bond[i].type,"2"))
//        {
//            group.num_db++;
//            if(group.neib[i].type[0]=='C') group.db_type=1;
//            else if(group.neib[i].type[0]=='N') group.db_type=2;
//            else if(group.neib[i].type[0]=='O') group.db_type=3;
//            else if(group.neib[i].type[0]=='S') group.db_type=4;
//            else if(group.neib[i].type[0]=='P') group.db_type=5;
//            else continue;
//        }
//        else if(!strcmp(group.bond[i].type,"1")&&
//                !strcmp(group.neib[i].type,"O.co2"))
//        {
//            group.db_type=3; group.num_db++;
//        }
//        else if(!strcmp(group.bond[i].type,"ar")&&
//                !strcmp(group.neib[i].type,"O.co2"))
//        {
//            group.db_type=3; group.num_db++;
//        }
//        else if(!strcmp(group.bond[i].type,"3"))
//        {
//            group.num_tb++;
//        }
//        else continue;
//    }
//
//    group.num_pi=0;
//
//    for(i=0;i<group.num_nonh;i++)
//    {
//        if(!strcmp(group.bond[i].type,"2")) continue;
//        else if(!strcmp(group.bond[i].type,"3")) continue;
//        else if(!strcmp(group.neib[i].type,"C.ar")) group.num_pi++;
//        else if(!strcmp(group.neib[i].type,"C.2")) group.num_pi++;
//        else if(!strcmp(group.neib[i].type,"C.1"))group.num_pi++;
//        else if(!strcmp(group.neib[i].type,"C.cat")) group.num_pi++;
//        else if(!strcmp(group.neib[i].type,"N.2")) group.num_pi++;
//        else if(!strcmp(group.neib[i].type,"N.1")) group.num_pi++;
//        else if(!strcmp(group.neib[i].type,"N.ar")) group.num_pi++;
//        else continue;
//    }
//
//    // check if the central atom is adjacent to any -SO-, -PO-, or -CO-
//
//    mark=false;
//
//    for(i=0;i<group.num_nonh;i++)
//    {
//        if(!strcmp(group.neib[i].type,"P.3")||
//           !strcmp(group.neib[i].type,"S.o")||
//           !strcmp(group.neib[i].type,"S.o2")||
//           !strcmp(group.neib[i].type,"C.2"))
//        {
//            num=group.neib[i].num_nonh;
//
//            for(j=0;j<num;j++)
//            {
//                id=this->atom[group.neib[i].id-1].neib[j];
//                if(id==group.center.id) continue;
//                else if(!strcmp(this->atom[id-1].type,"O.2"))
//                {
//                    mark=true; break;
//                }
//                else if(!strcmp(this->atom[id-1].type,"O.co2"))
//                {
//                    mark=true; break;
//                }
//                else continue;
//            }
//
//            if(mark==true) break;
//            else continue;
//        }
//        else continue;
//    }
//
//    if(mark==false) group.amide=0;
//    else group.amide=j+1;   // assign the value of the amide bond
//
//    group.valid=1; return group;
//}
//
//void xscore::Sum_HBonds(int num_candidate, HBond *candidate) {
//    int count,limit,latom,patom;
//    float angle,v1[3],v2[3];
//    HBond temp;
//    Group tmp_group;
//
//// Step 1: rank candidates according their strength in decreasing order
//    // note "fabs" is applied because SB strength could be negative
//
//    for(int i=0;i<num_candidate-1;i++)
//        for(int j=i+1;j<num_candidate;j++)
//        {
//            if(fabs(candidate[i].score)>=fabs(candidate[j].score)) continue;
//            else {SWAP(candidate[i],candidate[j]);}
//        }
//
//    // Step 2: check the angular limit: the angle between any two H-bonds
//    // on the same atom must be larger than 45 degrees
//    // note this this filter could be applied to both ligand and protein
//
//    for(int i=0;i<num_candidate-1;i++)
//        for(int j=i+1;j<num_candidate;j++)
//        {
//            if(candidate[i].latom!=candidate[j].latom) continue;
//
//            for(int k=0;k<3;k++)
//            {
//                latom=candidate[i].latom; patom=candidate[i].patom;
//                v1[k]=protein_atoms[patom-1].coords[k]-ligand_atoms[latom-1].coords[k];
//                latom=candidate[j].latom; patom=candidate[j].patom;
//                v2[k]=protein_atoms[patom-1].coords[k]-ligand_atoms[latom-1].coords[k];
//            }
//
//            angle=fabs(Angle_of_Two_Vectors(v1,v2));
//
//            if(angle<45.0) candidate[j].score=0.000;
//            else continue;
//        }
//
//
//
//    // Step 3: an donor atom shall not form more H-bonds than the H atoms it has.
//    // note that this filter is applied only to the ligand side
//
//    for(int i=0;i<num_candidate-1;i++)
//    {
//        if(candidate[i].type!=1) continue;
//
//        count=1; latom=candidate[i].latom;
//
//        if(!strcmp(ligand_atoms[latom-1].hb,"DA"))
//        {
//            limit=1;
//        }
//        else
//        {
//            tmp_group=this->Find_A_Group(latom);
//            limit=tmp_group.num_h;
//        }
//
//        for(int j=i+1;j<num_candidate;j++)
//        {
//            if(candidate[j].type!=1) continue;
//            else if(candidate[i].latom!=candidate[j].latom) continue;
//
//            count++;
//
//            if(count<=limit) continue;
//            else candidate[j].score=0.000;
//        }
//    }
//
//    // Step 4: an acceptor atom shall not form more H-bonds than its LPs
//    // note that this filter is applied only to the ligand side
//
//    for(int i=0;i<num_candidate-1;i++)
//    {
//        if(candidate[i].type!=2&&candidate[i].type!=3) continue;
//
//        count=1; latom=candidate[i].latom;
//
//        if(ligand_atoms[latom-1].type[0]=='O') limit=2;
//        else if(ligand_atoms[latom-1].type[0]=='N') limit=1;
//        else if(ligand_atoms[latom-1].type[0]=='S') limit=2;
//
//        for(int j=i+1;j<num_candidate;j++)
//        {
//            if(candidate[j].type!=2&&candidate[j].type!=3) continue;
//            else if(candidate[i].latom!=candidate[j].latom) continue;
//
//            count++;
//
//            if(count<=limit) continue;
//            else candidate[j].score=0.000;
//        }
//    }
//
//    return;
//}
//float xscore::calculate_HB() {
//    float sum;
//    HBond candidate[1000],temp;
//    int num_candidate;
//
//// clear the variables
//
//    sum=0.000; num_candidate=0;
//    for(int i=0;i<ligand_atoms.size();i++) ligand_atoms[i].score=0.000;
//
//    // first, get the HB pairs between protein and ligand
//
//    num_candidate=this->Get_HBond_Pair_PL(candidate);
//    this->Sum_HBonds(num_candidate,candidate);
//
//    // then sum their contributions
//
//    for(int i=0;i<num_candidate;i++)
//    {
//        if(fabs(candidate[i].score)<0.01) continue;
//
//        sum+=candidate[i].score;
//        ligand_atoms[candidate[i].latom-1].score+=candidate[i].score;
//    }
//
//
//
//    return sum;
//}