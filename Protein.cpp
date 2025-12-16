//
// Created by 91686 on 2023/8/17.
//
#include "xtools.h"
Protein::Protein()
{
    Clear();
}

Protein::~Protein()
{
    atom.clear(); chain.clear(); ring.clear();
    vol_dot.clear(); sur_dot.clear();
}

Protein::Protein(const Protein &original)
{
    this->xtool_format=original.xtool_format;
    strcpy(this->name,original.name);
    this->surface=original.surface;
    this->bnsur=original.bnsur;
    this->bpsur=original.bpsur;

    int i,n;

    this->num_atom=original.num_atom;
    this->atom.clear(); n=original.num_atom;
    for(i=0;i<n;i++) this->atom.push_back(original.atom[i]);

    this->num_chain=original.num_chain;
    this->chain.clear(); n=original.num_chain;
    for(i=0;i<n;i++) this->chain.push_back(original.chain[i]);

    this->num_ring=original.num_ring;
    this->ring.clear(); n=original.num_ring;
    for(i=0;i<n;i++) this->ring.push_back(original.ring[i]);

    this->vol_dot.clear(); n=original.vol_dot.size();
    for(i=0;i<n;i++) this->vol_dot.push_back(original.vol_dot[i]);

    this->sur_dot.clear(); n=original.sur_dot.size();
    for(i=0;i<n;i++) this->sur_dot.push_back(original.sur_dot[i]);
}

Protein& Protein::operator = (const Protein &original)
{
    if(this==&original) return *this;

    this->xtool_format=original.xtool_format;
    strcpy(this->name,original.name);
    this->surface=original.surface;
    this->bnsur=original.bnsur;
    this->bpsur=original.bpsur;

    int i,n;

    this->num_atom=original.num_atom;
    this->atom.clear(); n=original.num_atom;
    for(i=0;i<n;i++) this->atom.push_back(original.atom[i]);

    this->num_chain=original.num_chain;
    this->chain.clear(); n=original.num_chain;
    for(i=0;i<n;i++) this->chain.push_back(original.chain[i]);

    this->num_ring=original.num_ring;
    this->ring.clear(); n=original.num_ring;
    for(i=0;i<n;i++) this->ring.push_back(original.ring[i]);

    this->vol_dot.clear(); n=original.vol_dot.size();
    for(i=0;i<n;i++) this->vol_dot.push_back(original.vol_dot[i]);

    this->sur_dot.clear(); n=original.sur_dot.size();
    for(i=0;i<n;i++) this->sur_dot.push_back(original.sur_dot[i]);

    return *this;
}

void Protein::Clear()
{
    strcpy(name,"");
    surface=bnsur=bpsur=0.000;
    num_atom=0; atom.clear();
    num_chain=0; chain.clear();
    num_ring=0; ring.clear();
    vol_dot.clear(); sur_dot.clear();
}

void Protein::Read_SEQRES(FILE *fp)
{
    char line[256],head[256],temp[256],label,residue[10];
    int i,id,len;
    char *p;
    Residue tmp_residue;
    Chain tmp_chain;

    // clear the record first

    this->chain.clear(); this->num_chain=0;

    while(fgets(line,256,fp))
    {
        strcpy(temp,line); temp[6]='\0';
        strcpy(head,""); sscanf(temp,"%s",head);

        if(strcmp(head,"SEQRES")) continue;

        label=line[11]; line[11]=' ';

        // now check if this line starts a new chain

        id=-1;

        for(i=0;i<this->num_chain;i++)
        {
            if(label!=this->chain[i].label) continue;
            else {id=i; break;}
        }

        if(id>=0)  // an existing chain
        {
            for(i=0;i<=17;i++) line[i]=' '; line[70]='\0';

            for(;;)
            {
                strcpy(residue,"XXX"); sscanf(line,"%s", residue);
                if(!strcmp(residue,"XXX")) break;
                else
                {
                    tmp_residue.Clear();
                    tmp_residue.valid=1;
                    strcpy(tmp_residue.name,residue);
                    tmp_residue.chain=label;
                    strcpy(tmp_residue.id,"0"); // not known yet

                    this->chain[id].residue.push_back(tmp_residue);
                    this->chain[id].length++;

                    // erase this residue from the line

                    len=strlen(residue);
                    p=strstr(line,residue);
                    for(i=1;i<=len;i++) {*p=' '; p++;}
                }
            }
        }
        else  // a new chain
        {
            tmp_chain.Clear();
            tmp_chain.label=label;
            tmp_chain.valid=1;

            this->chain.push_back(tmp_chain);
            this->num_chain=this->chain.size();
            id=this->chain.size()-1;

            for(i=0;i<=17;i++) line[i]=' '; line[70]='\0';

            for(;;)
            {
                strcpy(residue,"XXX"); sscanf(line,"%s", residue);
                if(!strcmp(residue,"XXX")) break;
                else
                {
                    tmp_residue.Clear();
                    tmp_residue.valid=1;
                    strcpy(tmp_residue.name,residue);
                    tmp_residue.chain=label;
                    strcpy(tmp_residue.id,"0");  // not known yet

                    this->chain[id].residue.push_back(tmp_residue);
                    this->chain[id].length++;

                    // erase this residue from the line

                    len=strlen(residue);
                    p=strstr(line,residue);
                    for(i=1;i<=len;i++) {*p=' '; p++;}
                }
            }
        }
    }

    rewind(fp); return;
}

void Protein::Read_ATOM(FILE *fp)
{
    int i,len;
    char line[256],head[80],temp[256];
    char name[80],id[80],residue[80],chain,res_id[80];
    float x,y,z,occupancy,bfactor;
    struct atom tmp_atom;

    while(fgets(line,256,fp))
    {
        // get the heading and erase that section

        strcpy(temp,line); temp[6]='\0';
        strcpy(head,""); sscanf(temp,"%s",head);
        for(i=0;i<6;i++) line[i]=' ';

        if(!strcmp(head,"END")) break;
        else if(!strcmp(head,"ATOM")||!strcmp(head,"HETATM"))
        {
            // get the atom id and erase that section

            strcpy(temp,line); temp[11]='\0';
            sscanf(temp,"%s",id);
            for(i=0;i<11;i++) line[i]=' ';

            // get the atom name and erase that section

            strcpy(temp,line); temp[17]='\0';
            sscanf(temp,"%s",name);
            for(i=0;i<17;i++) line[i]=' ';

            if(strstr(name,"LP")) continue;  // skip lone pairs

            // get the residue name and erase that section

            strcpy(temp,line); temp[20]='\0';
            sscanf(temp,"%s",residue);
            for(i=0;i<20;i++) line[i]=' ';

            // get the chain label and erase that section

            chain=line[21];
            for(i=0;i<22;i++) line[i]=' ';

            // get the residue id and erase that section

            sscanf(line,"%s", res_id);
            for(i=0;i<27;i++) line[i]=' ';

            len=strlen(res_id);
            if(isdigit(res_id[len-1])) strcat(res_id," ");

            // get coordinate x and erase that section

            strcpy(temp,line); temp[38]='\0';
            sscanf(temp,"%f",&x);
            for(i=0;i<38;i++) line[i]=' ';

            // get coordinate y and erase that section

            strcpy(temp,line); temp[46]='\0';
            sscanf(temp,"%f",&y);
            for(i=0;i<46;i++) line[i]=' ';

            // get coordinate z and erase that section

            strcpy(temp,line); temp[54]='\0';
            sscanf(temp,"%f",&z);
            for(i=0;i<54;i++) line[i]=' ';

            // get occupancy probability and erase that section
            // note that not all PDB files have this field

            len=strlen(line);

            if(len>=60)
            {
                strcpy(temp,line); temp[60]='\0';
                sscanf(temp,"%f",&occupancy);
                for(i=0;i<60;i++) line[i]=' ';
            }

            // get B-factor and erase that section
            // note that not all PDB files have this field

            len=strlen(line);

            if(len>=66)
            {
                strcpy(temp,line); temp[66]='\0';
                sscanf(temp,"%f",&bfactor);
                for(i=0;i<66;i++) line[i]=' ';
            }

            // now summarize

            sscanf(id,"%d", &tmp_atom.id);
            sscanf(name,"%s", tmp_atom.name);
            sscanf(residue,"%s", tmp_atom.residue);
            tmp_atom.chain=chain;
            strcpy(tmp_atom.res_id,res_id);
            tmp_atom.coor[0]=x;
            tmp_atom.coor[1]=y;
            tmp_atom.coor[2]=z;
            tmp_atom.occupancy=occupancy;
            tmp_atom.bfactor=bfactor;

            if(!strcmp(head,"ATOM")) tmp_atom.part=1;  // regular atoms
            else tmp_atom.part=2;  // hetero-atoms

            tmp_atom.valid=1; tmp_atom.origin=2; // protein atom

            atom.push_back(tmp_atom);
        }
        else continue;
    }

    num_atom=atom.size();

    rewind(fp); return;
}

void Protein::Analyze_Sequence()
{
    int i,j,id,total;
    bool mark;
    Residue tmp_residue;
    Chain tmp_chain;

    this->chain.clear(); this->num_chain=0;

    for(i=0;i<this->num_atom;i++)
    {
        if(this->atom[i].part>1) continue;   // hetero atom

        // first, check if this atom belongs to a known chain

        mark=false; id=-1; total=this->chain.size();

        for(j=0;j<total;j++)
        {
            if(this->atom[i].chain!=this->chain[j].label) continue;
            else {id=j; mark=true; break;}
        }

        if(mark==true)
        {
            // check if this atom belongs to a known residue

            mark=false; total=this->chain[id].residue.size();

            for(j=0;j<total;j++)
            {
                if(strcmp(atom[i].res_id, chain[id].residue[j].id)) continue;
                else if(strcmp(atom[i].residue,chain[id].residue[j].name)) continue;
                else {mark=true; break;}
            }

            if(mark==true) continue;
            else
            {
                // add a new residue to this chain

                tmp_residue.Clear();
                tmp_residue.valid=1;
                strcpy(tmp_residue.name,atom[i].residue);
                strcpy(tmp_residue.id,atom[i].res_id);
                tmp_residue.chain=chain[id].label;

                chain[id].residue.push_back(tmp_residue);
                chain[id].length++;
            }
        }
        else
        {
            id=this->chain.size();

            // add a new chain
            tmp_chain.Clear();
            tmp_chain.valid=1;
            tmp_chain.label=this->atom[i].chain;
            tmp_chain.residue.clear();

            this->chain.push_back(tmp_chain);

            // add the first residue to this chain
            tmp_residue.Clear();
            tmp_residue.valid=1;
            strcpy(tmp_residue.name,atom[i].residue);
            strcpy(tmp_residue.id,atom[i].res_id);
            tmp_residue.chain=chain[id].label;

            chain[id].residue.push_back(tmp_residue);
            chain[id].length++;
        }
    }

    this->num_chain=this->chain.size();
    for(i=0;i<num_chain;i++) chain[i].length=chain[i].residue.size();

    return;
}

void Protein::Read_CONECT(FILE *fp)
{
    char line[256],head[256],temp[256];
    int i,j,id,atom_id,len,count;
    bool mark;

    while(fgets(line,256,fp))
    {
        strcpy(temp,line); temp[6]='\0';
        strcpy(head,""); sscanf(temp,"%s",head);

        if(strcmp(head,"CONECT")) continue;

        for(i=0;i<6;i++) line[i]=' '; // remove the header

        // read the atom id first

        strcpy(temp,line); temp[11]='\0';
        sscanf(temp,"%d",&atom_id);
        for(i=0;i<11;i++) line[i]=' ';

        mark=false;

        for(i=0;i<this->num_atom;i++)
        {
            if(this->atom[i].valid==0) continue;
            else if(this->atom[i].part==1) continue;
            else if(this->atom[i].id!=atom_id) continue;
            else {id=i; mark=true; break;}
        }

        if(mark==false) continue;	// no such atom in hetero atom list

        for(i=1;i<=MAX_ATOM_NEIB;i++)
        {
            len=11+i*5; strcpy(temp,line); temp[len]='\0';
            if(!isdigit(temp[len-1])) break;  // no more numbers
            else
            {
                sscanf(temp,"%d",&atom_id);
                for(j=0;j<len;j++) line[j]=' ';

                count=this->atom[id].num_neib;
                this->atom[id].neib[count]=atom_id;
                this->atom[id].num_neib++;
            }
        }

/*
	 printf("CONECT read: %d ", this->atom[id].id);
	 for(i=0;i<this->atom[id].num_neib;i++)
		printf("%d ", this->atom[id].neib[i]);
	 printf("\n");
*/
    }

    rewind(fp); return;
}


void Protein::Read_From_PDB(char *filename)
{
    FILE *fp;
    int i,j,n,count,len;
    char line[256],head[80],temp[256];
    bool mark,mark1;

    if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

    Clear();	// clear this protein

    // first, check the PDB file

    count=0; this->xtool_format=false;

    mark=false;	// SEQRES information indicator
    mark1=false;	// COMPND information indicator

    for(;;)
    {
        if(fgets(line,256,fp)==NULL) break;
        else if(line[0]=='#') continue;
        else if(Blank_Line_Check(line)==TRUE) continue;
        else
        {
            sscanf(line,"%s", head);
            if(!strcmp(head,"ATOM")) count++;
            else if(!strcmp(head,"HETATM")) count++;
            else if(!strcmp(head,"SEQRES")) mark=true;
            else if(!strcmp(head,"COMPND")&&mark1==false)
            {
                len=strlen(line);

                if(len<70) n=len;
                else n=70;

                j=0;
                for(i=10;i<n;i++)
                {
                    if(line[i]=='\n') continue;
                    else {temp[j]=line[i];j++;}
                }
                temp[j]='\0';

                strcpy(this->name,temp); mark1=true;

                // note that sometimes COMPND spans more than
                // one line, but we only take the first line.
                // why the hell the name should be that long?
            }
            else if(!strcmp(head,"REMARK"))
            {
                if(strstr(line,"X-Tool")||
                   strstr(line,"X-TOOL")||
                   strstr(line,"XTool")||
                   strstr(line,"XTOOL"))
                {
                    this->xtool_format=true;
                }
            }
            else continue;
        }
    }

    rewind(fp);

    if(count<5) PDB_Format_Error(filename);

    // read the ATOM/HETATM section

    this->Read_ATOM(fp);

    // read the input sequence information

    if(mark==true) this->Read_SEQRES(fp);
    else this->Analyze_Sequence();

    // read the CONECT section if it is generated by X-Tool itself

    if(this->xtool_format==true) this->Read_CONECT(fp);

    fclose(fp); return;
}

void Protein::Check_Atom_Type()
{
    int i,j,n,len;
    struct atom atm;
    char tmp_name[10];

    for(i=0;i<this->num_atom;i++)
    {
        atm=this->atom[i];

        // first, correct residue names -------------------------------------

        if(!strcmp(atm.residue,"WAT")) strcpy(atm.residue,"HOH");
        if(!strcmp(atm.residue,"WTR")) strcpy(atm.residue,"HOH");
        if(!strcmp(atm.residue,"CYX")) strcpy(atm.residue,"CYS");
        if(!strcmp(atm.residue,"HID")) strcpy(atm.residue,"HIS");
        if(!strcmp(atm.residue,"HIE")) strcpy(atm.residue,"HIS");
        if(!strcmp(atm.residue,"HIP")) strcpy(atm.residue,"HIS");

        // the following three mistakes are often seen in SYBYL-generated
        // PDB files

        if(!strcmp(atm.residue,"SO")&&
           (!strcmp(atm.name,"S")|| !strcmp(atm.name,"O1")||
            !strcmp(atm.name,"O2")|| !strcmp(atm.name,"O3")||
            !strcmp(atm.name,"O4")))
        {
            strcpy(atm.residue,"SO4");
        }

        if(!strcmp(atm.residue,"PO")&&
           (!strcmp(atm.name,"P")|| !strcmp(atm.name,"O1")||
            !strcmp(atm.name,"O2")|| !strcmp(atm.name,"O3")||
            !strcmp(atm.name,"O4")))
        {
            strcpy(atm.residue,"PO4");
        }

        if(!strcmp(atm.residue,"CA")&&!strcmp(atm.name,"O"))
        {
            strcpy(atm.residue,"HOH");
        }

        // now correct atom names -------------------------------------------

        if(isdigit(atm.name[0]))
        {
            len=strlen(atm.name); n=0;
            for(j=1;j<len;j++)
            {
                tmp_name[n]=atm.name[j]; n++;
            }
            tmp_name[n]=atm.name[0]; n++; tmp_name[n]='\0';
            strcpy(atm.name,tmp_name);
        }

        if(!strcmp(atm.name,"OCT")) strcpy(atm.name,"OXT");  // all residue
        if(!strcmp(atm.name,"HN")) strcpy(atm.name,"H");     // all residue

        if(!strcmp(atm.residue,"ACE"))
        {
            if(!strcmp(atm.name,"CH3")) strcpy(atm.name,"CA");
        }
        else if(!strcmp(atm.residue,"ASN"))
        {
            if(!strcmp(atm.name,"AD1")) strcpy(atm.name,"OD1");
            else if(!strcmp(atm.name,"AD2")) strcpy(atm.name,"ND2");
        }
        else if(!strcmp(atm.residue,"GLN"))
        {
            if(!strcmp(atm.name,"AE1")) strcpy(atm.name,"OE1");
            else if(!strcmp(atm.name,"AE2")) strcpy(atm.name,"NE2");
        }
        else if(!strcmp(atm.residue,"LEU"))
        {
            if(!strcmp(atm.name,"CD")) strcpy(atm.name,"CD1");
            else if(!strcmp(atm.name,"CE")) strcpy(atm.name,"CD2");
        }

        this->atom[i]=atm;
    }

    return;
}


void Protein::Detect_Connections()
{
    extern ForceField *ff;
    int i,j,k,count,id;
    bool mark;
    float cutoff,d;

    for(i=0;i<this->num_atom;i++)
    {
        this->atom[i].num_neib=this->atom[i].num_nonh=0;
        for(j=0;j<MAX_ATOM_NEIB;j++) this->atom[i].neib[j]=0;
    }

    // build the connections for regular atoms first

    cutoff=2.00;  // covalent bond distance cutoff

    for(i=0;i<this->num_atom;i++)
    {
        if(this->atom[i].valid==0) continue;
        else if(this->atom[i].part!=1) continue;
        else if(this->atom[i].num_neib>0) continue; // already done

        count=0;

        // find internal connections inside the same residue first

        for(j=0;j<this->num_atom;j++)
        {
            if(i==j) continue;
            else if(this->atom[j].valid==0) continue;
            else if(this->atom[j].part!=1) continue;
            else if(atom[i].chain!=atom[j].chain) continue;
            else if(strcmp(atom[i].residue,atom[j].residue)) continue;
            else if(strcmp(atom[i].res_id,atom[j].res_id)) continue;

            // there is no bond between hydrogen atoms

            if(!strcmp(atom[i].type,"H")&&
               !strcmp(atom[j].type,"H")) continue;

            // now check whether atom i and atom j are covalently bound

            if(ff->Patom_Connection_Test(atom[i],atom[j])==true)
            {
                if(count>=MAX_ATOM_NEIB) continue;  // already full
                else
                {
                    this->atom[i].neib[count]=j+1; count++;
                }
            }
            else continue;
        }

        // now detect peptide amide bonds

        if(!strcmp(atom[i].name,"N")&&!strcmp(atom[i].type,"N.am"))
        {
            for(j=0;j<this->num_atom;j++)
            {
                if(this->atom[j].valid==0) continue;
                else if(this->atom[j].part!=1) continue;
                else if(atom[i].chain!=atom[j].chain) continue;
                else if(!strcmp(atom[i].res_id,atom[j].res_id)) continue;
                else if(strcmp(atom[j].type,"C.2")) continue;
                else if(strcmp(atom[j].name,"C")) continue;

                d=Distance(atom[i].coor,atom[j].coor);

                if(d>cutoff) continue;
                else if(count>=MAX_ATOM_NEIB) continue;
                else
                {
                    this->atom[i].neib[count]=j+1;
                    count++;
                }
            }
        }
        else if(!strcmp(atom[i].name,"C")&&!strcmp(atom[i].type,"C.2"))
        {
            for(j=0;j<this->num_atom;j++)
            {
                if(this->atom[j].valid==0) continue;
                else if(this->atom[j].part!=1) continue;
                else if(atom[i].chain!=atom[j].chain) continue;
                else if(!strcmp(atom[i].res_id,atom[j].res_id)) continue;
                else if(strcmp(atom[j].type,"N.am")) continue;
                else if(strcmp(atom[j].name,"N")) continue;

                d=Distance(atom[i].coor,atom[j].coor);

                if(d>cutoff) continue;
                else if(count>=MAX_ATOM_NEIB) continue;
                else
                {
                    this->atom[i].neib[count]=j+1;
                    count++;
                }
            }
        }

        this->atom[i].num_neib=count;
    }

    // now find the patoms bound to the metal ions

    cutoff=3.00;  // M-bond distance cutoff

    for(i=0;i<this->num_atom;i++)
    {
        if(this->atom[i].valid==0) continue;
        else if(this->atom[i].part==1) continue;
        else if(strcmp(this->atom[i].xtype,"M+")) continue;
        else if(this->atom[i].num_neib>0) continue;  // already done

        count=0;

        for(j=0;j<this->num_atom;j++)
        {
            if(this->atom[j].valid==0) continue;
            else if(this->atom[j].part!=1) continue;
            else if(strcmp(this->atom[j].hb,"A")&&
                    strcmp(this->atom[j].hb,"DA")) continue;

            d=Distance(this->atom[i].coor,this->atom[j].coor);

            if(d>cutoff) continue;
            else if(count>=MAX_ATOM_NEIB) continue;
            else
            {
                this->atom[i].neib[count]=j+1;  // note this
                count++;
                // here we do not add connections to atom[j]
            }
        }

        this->atom[i].num_neib=count;
    }

    // now build the connections for SO4s and PO4s

    for(i=0;i<this->num_atom;i++)
    {
        if(this->atom[i].valid==0) continue;
        else if(this->atom[i].part==1) continue;
        else if(strcmp(this->atom[i].residue,"SO4")&&
                strcmp(this->atom[i].residue,"PO4")) continue;
        else if(strcmp(this->atom[i].name,"S")&&
                strcmp(this->atom[i].name,"P")) continue;

        if(this->atom[i].num_neib>0) continue;  // done

        // now check the satellite atoms for this P or S atom

        for(j=0;j<this->num_atom;j++)
        {
            if(this->atom[j].valid==0) continue;
            else if(this->atom[j].part==1) continue;
            else if(strcmp(this->atom[j].residue,"SO4")&&
                    strcmp(this->atom[j].residue,"PO4")) continue;
            else if(strcmp(this->atom[j].res_id,
                           this->atom[i].res_id)) continue;
            else if(strstr(this->atom[j].name,"S")) continue;
            else if(strstr(this->atom[j].name,"P")) continue;

            count=this->atom[i].num_neib; mark=false;

            for(k=0;k<count;k++)
            {
                if(this->atom[i].neib[k]==(j+1))
                {
                    mark=true; break;
                }
                else continue;
            }

            if(mark==false)
            {
                this->atom[i].neib[count]=j+1;  // note this
                this->atom[i].num_neib++;
            }

            count=this->atom[j].num_neib; mark=false;

            for(k=0;k<count;k++)
            {
                if(this->atom[j].neib[k]==(i+1))
                {
                    mark=true; break;
                }
                else continue;
            }

            if(mark==false)
            {
                this->atom[j].neib[count]=i+1;  // note this
                this->atom[j].num_neib++;
            }
        }
    }

/*
 // now list what we get

 for(i=0;i<this->num_atom;i++)
	{
	 if(this->atom[i].valid<=0) continue;

	 printf("ATOM %d %s %s %s ",
		this->atom[i].id,
		this->atom[i].name,
		this->atom[i].residue,
		this->atom[i].res_id);

	 for(j=0;j<this->atom[i].num_neib;j++)
		{
		 id=this->atom[i].neib[j]; printf("%d ", id);
		}

	 printf("\n");
	}
*/

    // now determine number of heavy atoms for each valid atom

    for(i=0;i<this->num_atom;i++)
    {
        if(this->atom[i].valid<=0) continue;

        count=0;

        for(j=0;j<this->atom[i].num_neib;j++)
        {
            id=this->atom[i].neib[j];
            if(!strcmp(this->atom[id-1].type,"H")) continue;
            else count++;
        }

        this->atom[i].num_nonh=count;
    }

    return;
}

void Protein::Calculate_HB_Root()
{
    int i,j,id,count;
    float tmpx,tmpy,tmpz;

    // (1) we do calculate root for polar hydrogens
    // (2) we do NOT calculate root for non-HB atoms
    // (3) we do NOT calculate root for metal ions
    // (4) we do NOT calculate root for waters
    // (5) we do calculate root for the oxygen atoms on SO4 and PO4

    for(i=0;i<num_atom;i++)
    {
        if(atom[i].valid<=0) continue;
        else if(!strcmp(atom[i].xtype,"H")) continue;
        else if(!strcmp(atom[i].xtype,"O.w")) continue;
        else if(!strcmp(atom[i].hb,"N")) continue;
        else if(!strcmp(atom[i].hb,"H")) continue;
        else if(!strcmp(atom[i].hb,"P")) continue;
        else if(!strcmp(atom[i].hb,"M")) continue;

        if(atom[i].num_neib==0) continue;

        tmpx=tmpy=tmpz=0.000; count=0;

        for(j=0;j<atom[i].num_neib;j++)
        {
            id=atom[i].neib[j]-1;
            if(!strcmp(atom[id].type,"H")) continue;
            else
            {
                tmpx+=(atom[id].coor[0]);
                tmpy+=(atom[id].coor[1]);
                tmpz+=(atom[id].coor[2]);
                count++;
            }
        }

        if(count==0) continue;

        atom[i].root[0]=tmpx/count;
        atom[i].root[1]=tmpy/count;
        atom[i].root[2]=tmpz/count;
    }

    return;
}

void Protein::Value_Atom(int flag)
{
    extern ForceField *ff;
    int i;

    // note that within this subroutine the protein is considered to be
    // all-atom model. this is because this subroutine is shared by
    // all of the functions related to protein

    // correct atom names and residue names first.

    this->Check_Atom_Type();

    // assign the parameters for all atoms

    for(i=0;i<num_atom;i++) ff->Assign_Patom_Parameters(atom[i]);

    // build connection tables, must assign parameters first

    this->Detect_Connections();

    // calculate H-bond root

    if(flag!=0) Calculate_HB_Root();

    return;
}

void Protein::Define_Pocket(float origin[], float radius)
{
    // need to be finished
    return;
}

//int Protein::Check_Buried_Ratio(float coor[]) const
//{
//    extern ForceField *ff;
//    extern Ligand *ligand;
//    int j,k,count;
//    float d,buried_ratio;
//    bool mark;
//    DotSet tmp_set;
//
//    // notice that only protein atoms are considered in this check
//    // current criterion is 90%
//
//    buried_ratio=0.90;
//
//    tmp_set=ff->Get_Surface_Dot(2*WATER_R,coor[0],coor[1],coor[2]);
//
//    count=0;
//
//    for(j=0;j<tmp_set.num_dot;j++)
//    {
//        mark=false;	// check the protein
//
//        for(k=0;k<this->num_atom;k++)
//        {
//            if(this->atom[k].valid<2) continue;
//            else if(!strcmp(this->atom[k].type,"H")) continue;
//            else if(!strcmp(this->atom[k].type,"O.w")) continue;
//
//            d=Distance(tmp_set.dot[j].coor,this->atom[k].coor);
//
//            if(d>(WATER_R+this->atom[k].R)) continue;
//            else {mark=true; break;}
//        }
//
//        if(mark==true) {count++; tmp_set.dot[j].valid=0; continue;}
//
//        mark=false;	// check the ligand
//
//        for(k=0;k<ligand->num_atom;k++)
//        {
//            if(ligand->atom[k].valid<=0) continue;
//            else if(!strcmp(ligand->atom[k].type,"H")) continue;
//
//            d=Distance(tmp_set.dot[j].coor,ligand->atom[k].coor);
//
//            if(d>(WATER_R+ligand->atom[k].R)) continue;
//            else {mark=true; break;}
//        }
//
//        if(mark==true) {count++; tmp_set.dot[j].valid=0;}
//        else continue;
//    }
//
//    // now determine this water is buried or not
//
//    if(((float)count/tmp_set.num_dot)<buried_ratio) return FALSE;
//    else return TRUE;
//}
//
//// ***************************************************************************
//// determine which water to keep (valid=2) and which to disregard (valid=1).
//// latest update 03/04/03
//// ***************************************************************************
//void Protein::Define_Water(const Ligand *ligand, bool w_flag)
//{
//    int i,j,mark;
//    float d;
//
//    if(w_flag==false)	// water molecules are not considered at all
//    {
//        for(i=0;i<this->num_atom;i++)
//        {
//            if(this->atom[i].valid<2) continue;
//            else if(strcmp(this->atom[i].type,"O.w")) continue;
//            else this->atom[i].valid=1;
//        }
//        return;
//    }
//
//    for(i=0;i<this->num_atom;i++)
//    {
//        if(this->atom[i].valid<2) continue;
//        else if(strcmp(this->atom[i].type,"O.w")) continue;
//
//        // first, this water should be close enough to the ligand
//        // only keep the first solvation shell
//
//        mark=false;
//
//        for(j=0;j<ligand->num_atom;j++)
//        {
//            if(ligand->atom[j].valid<=0) continue;
//            else if(!strcmp(ligand->atom[j].type,"H")) continue;
//
//            d=Distance(this->atom[i].coor,ligand->atom[j].coor);
//
//            if(d<(ligand->atom[j].R+2*WATER_R))
//            {
//                mark=true; break;
//            }
//            else continue;
//        }
//
//        if(mark==false) {this->atom[i].valid=1; continue;}
//
//        // second, buried status check
//
//        mark=Check_Buried_Ratio(this->atom[i].coor);
//
//        if(mark==true) this->atom[i].valid=2;
//        else this->atom[i].valid=1;
//    }
//
//    return;
//}

void Protein::Define_Pocket(const Ligand *ligand, float cutoff)
{
    extern ForceField *ff;
    int i,j,mark,count;
    float d;
    int num_pocket_res;
    Residue *pocket_res=NULL;

    num_pocket_res=0;

    pocket_res=new Residue[300];
    if(pocket_res==NULL) Memory_Allocation_Error();

    // find pocket atoms on the protein and label their valid as '2',
    // including the waters close to the ligand

    for(i=0;i<num_atom;i++)
    {
        if(atom[i].valid<=0) continue;

        // filter out all of the non-polar hydrogen atoms
        // since X-Score uses unit-atom model

        if(!strcmp(atom[i].xtype,"H")&&strcmp(atom[i].residue,"COF"))
        {
            atom[i].valid=0; continue;
        }

        // check if this atom is close to the ligand

        mark=FALSE;

        for(j=0;j<ligand->num_atom;j++)
        {
            if(ligand->atom[j].valid<=0) continue;
            else if(!strcmp(ligand->atom[j].type,"H")) continue;
            else
            {
                d=Distance(atom[i].coor,ligand->atom[j].coor);
                if(d>cutoff) continue;
                else {mark=TRUE; break;}
            }
        }

        if(mark==FALSE) continue;	// not a pocket atom

        atom[i].valid=2;

        // if this atom is a water or a metal ion, do not check its residue

        if(!strcmp(atom[i].xtype,"O.w")) continue;
        if(!strcmp(atom[i].xtype,"M+")) continue;

        // for regular atoms, check if it has been included in a
        // known pocket residue; if not, add that newly found pocket residue

        mark=FALSE;

        for(j=0;j<num_pocket_res;j++)
        {
            if(!strcmp(atom[i].residue,pocket_res[j].name)&&
               !strcmp(atom[i].res_id,pocket_res[j].id)&&
               (atom[i].chain==pocket_res[j].chain)) {mark=TRUE;break;}
            else continue;
        }

        if(mark==FALSE)  // new pocket residue found
        {
            strcpy(pocket_res[num_pocket_res].name,atom[i].residue);
            strcpy(pocket_res[num_pocket_res].id,atom[i].res_id);
            pocket_res[num_pocket_res].chain=atom[i].chain;
            pocket_res[num_pocket_res].valid=1;
            num_pocket_res++;
        }
    }

    // now check if the binding pocket is well defined

    count=0;

    for(i=0;i<num_pocket_res;i++)
    {
        if(!strcmp(pocket_res[i].name,"HET")) continue;
        else if(!strcmp(pocket_res[i].name,"COF")) continue;
        else if(!strcmp(pocket_res[i].name,"WAT")) continue;
        else if(!strcmp(pocket_res[i].name,"HOH")) continue;
        else count++;
    }

    if(count<3) // there should be at least 3 residues as pocket!
    {
        puts("Error: cannot find binding pocket residues on the protein.");
        puts("Probably the ligand has not been docked with the protein.");
        exit(1);
    }

    // define all the left atoms in pocket residues as pocket atoms

    for(i=0;i<num_atom;i++)
    {
        if(atom[i].valid<=0) continue;
        else if(atom[i].valid==2) continue;
        else if(!strcmp(atom[i].type,"O.w")) continue;
        else if(!strcmp(atom[i].hb,"M")) continue;

        for(j=0;j<num_pocket_res;j++)
        {
            if(strcmp(atom[i].residue,pocket_res[j].name)) continue;
            else if(strcmp(atom[i].res_id,pocket_res[j].id)) continue;
            else if(atom[i].chain!=pocket_res[j].chain) continue;
            else {atom[i].valid=2; break;}
        }
    }

    // now detect all of the aromatic rings within pocket residues
    // this information is need for later scoring purposes

    Ring tmp_ring;

    this->ring.clear(); this->num_ring=0;

    for(i=0;i<num_pocket_res;i++)
    {
        if(strcmp(pocket_res[i].name,"PHE")&&
           strcmp(pocket_res[i].name,"TYR")&&
           strcmp(pocket_res[i].name,"HIS")&&
           strcmp(pocket_res[i].name,"TRP")) continue;

        tmp_ring.Clear();

        for(j=0;j<num_atom;j++)
        {
            if(atom[j].valid!=2) continue;
            else if(atom[j].part!=1) continue;
            else if(atom[j].ring!=2) continue;
            else if(atom[j].chain!=pocket_res[i].chain) continue;
            else if(strcmp(atom[j].residue,pocket_res[i].name)) continue;
            else if(strcmp(atom[j].res_id,pocket_res[i].id)) continue;

            tmp_ring.atom_id.push_back(j+1);
            tmp_ring.centroid[0]+=atom[j].coor[0];
            tmp_ring.centroid[1]+=atom[j].coor[1];
            tmp_ring.centroid[2]+=atom[j].coor[2];
        }

        tmp_ring.num_member=tmp_ring.atom_id.size();

        if(tmp_ring.num_member>0)
        {
            tmp_ring.centroid[0]/=tmp_ring.num_member;
            tmp_ring.centroid[1]/=tmp_ring.num_member;
            tmp_ring.centroid[2]/=tmp_ring.num_member;

            tmp_ring.valid=1; tmp_ring.type=2;
            this->ring.push_back(tmp_ring);
        }
    }

    this->num_ring=this->ring.size();

    if(pocket_res) if(pocket_res) delete [] pocket_res; return;
}