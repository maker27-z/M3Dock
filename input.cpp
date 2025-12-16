//
// Created by 91686 on 2023/8/18.
//
# include "xtools.h"

Input::Input()
{
    strcpy(function,"none");

    strcpy(input_file,"none");
    strcpy(output_file,"none");
    strcpy(log_file,"none");

    strcpy(receptor_file,"none");
    strcpy(cofactor_file,"none");
    strcpy(reference_file,"none");
    strcpy(ligand_file,"none");

    num_method=3;
    strcpy(apply_hpscore,"YES");
    strcpy(apply_hmscore,"YES");
    strcpy(apply_hsscore,"YES");
    strcpy(apply_pmfscore,"NO");
    strcpy(show_abs,"NO");

    hpscore_cvdw=0.004;
    hpscore_chb=0.054;
    hpscore_chp=0.009;
    hpscore_crt=-0.061;
    hpscore_c0=3.441;

    hmscore_cvdw=0.004;
    hmscore_chb=0.101;
    hmscore_chm=0.387;
    hmscore_crt=-0.097;
    hmscore_c0=3.567;

    hsscore_cvdw=0.004;
    hsscore_chb=0.073;
    hsscore_chs=0.004;
    hsscore_crt=-0.090;
    hsscore_c0=3.328;

    if(getenv("XSCORE_PARAMETER")==NULL)
    {
        puts("Warning: XSCORE_PARAMETER is not set ... use default setting");
        strcpy(parameter_dir,"./parameter/");
    }
    else strcpy(parameter_dir, getenv("XSCORE_PARAMETER"));

    Check_Directory(parameter_dir);
}

Input::~Input()
{
    // destructor
}

void Input::Read_Inputs(char *filename)
{
    FILE *fp;
    char line[256],head[256];

    if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

    for(;;)
    {
        if(fgets(line,256,fp)==NULL) break;
        else if(line[0]=='#') continue;
        else if(Blank_Line_Check(line)==TRUE) continue;
        else sscanf(line,"%s",head);

        if(!strcmp(head,"FUNCTION"))
        {
            sscanf(line,"%*s%s",function);
        }
        else if(!strcmp(head,"INPUT_FILE"))
        {
            sscanf(line,"%*s%s",input_file);
        }
        else if(!strcmp(head,"OUTPUT_FILE"))
        {
            sscanf(line,"%*s%s",output_file);
        }
        else continue;
    }

    fclose(fp); return;
}

void Input::Show_Contents() const
{
    printf("FUNCTION = %s\n", this->function);
    printf("INPUT_FILE = %s\n", this->input_file);
    printf("OUTPUT_FILE = %s\n", this->output_file);

    return;
}

void Input::Missing_Parameter_Error(const char *name) const
{
    printf("\n");
    printf("Error: You have probably forgot telling me %s.\n", name);
    return;
}

void Input::Invalid_Parameter_Error(const char *name) const
{
    printf("\n");
    printf("Error: %s has an invalid value.\n", name);
    return;
}

// ***********************************************************************
// add '/' to the directory name automatically
// if the directory is set for output, the 'flag' needs to be non-zero
// latest update: 11/21/2003
// ***********************************************************************
void Input::Check_Directory(char *dirname, int flag)
{
    int len;
    FILE *fp;
    char filename[256],command[256];

    len=strlen(dirname);

    if(dirname[len-1]!='/') strcat(dirname,"/");

    if(flag==FALSE) return;

    // check if the directory exists

    strcpy(filename,dirname); strcat(filename,".");

    if((fp=fopen(filename,"r"))==NULL)
    {
        strcpy(command, "mkdir ");
        strcat(command, dirname);
        system(command);
        printf("Directory '%s' does not exist ... ", dirname);
        printf("it is created.\n");
        return;
    }
    else
    {
        fclose(fp); return;
    }
}

SCORE_Input::SCORE_Input()
{
    strcpy(function,"SCORE");
    num_hits=0; strcpy(hits_dir,"none");
}

SCORE_Input::~SCORE_Input()
{
}

void SCORE_Input::Read_Inputs(char *filename)
{
    FILE *fp;
    char line[256],head[256];
    int num_error;

    if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

    for(;;)
    {
        if(fgets(line,256,fp)==NULL) break;
        else if(line[0]=='#') continue;
        else if(Blank_Line_Check(line)==TRUE) continue;
        else sscanf(line,"%s",head);

        if(!strcmp(head,"RECEPTOR_PDB_FILE"))
        {
            sscanf(line,"%*s%s",receptor_file);
        }
        else if(!strcmp(head,"COFACTOR_MOL2_FILE"))
        {
            sscanf(line,"%*s%s",cofactor_file);
        }
        else if(!strcmp(head,"REFERENCE_MOL2_FILE"))
        {
            sscanf(line,"%*s%s",reference_file);
        }
        else if(!strcmp(head,"LIGAND_MOL2_FILE"))
        {
            sscanf(line,"%*s%s",ligand_file);
        }
        else if(!strcmp(head,"OUTPUT_TABLE_FILE"))
        {
            sscanf(line,"%*s%s",output_file);
        }
        else if(!strcmp(head,"OUTPUT_LOG_FILE"))
        {
            sscanf(line,"%*s%s",log_file);
        }
        else if(!strcmp(head,"NUMBER_OF_HITS"))
        {
            sscanf(line,"%*s%d",&num_hits);
        }
        else if(!strcmp(head,"HITS_DIRECTORY"))
        {
            sscanf(line,"%*s%s",hits_dir);
        }

        if(!strcmp(head,"APPLY_HPSCORE"))
        {
            sscanf(line,"%*s%s",apply_hpscore);
        }
        else if(!strcmp(head,"HPSCORE_CVDW"))
        {
            sscanf(line,"%*s%f", &hpscore_cvdw);
        }
        else if(!strcmp(head,"HPSCORE_CHB"))
        {
            sscanf(line,"%*s%f", &hpscore_chb);
        }
        else if(!strcmp(head,"HPSCORE_CHP"))
        {
            sscanf(line,"%*s%f", &hpscore_chp);
        }
        else if(!strcmp(head,"HPSCORE_CRT"))
        {
            sscanf(line,"%*s%f", &hpscore_crt);
        }
        else if(!strcmp(head,"HPSCORE_C0"))
        {
            sscanf(line,"%*s%f", &hpscore_c0);
        }
        if(!strcmp(head,"APPLY_HMSCORE"))
        {
            sscanf(line,"%*s%s",apply_hmscore);
        }
        else if(!strcmp(head,"HMSCORE_CVDW"))
        {
            sscanf(line,"%*s%f", &hmscore_cvdw);
        }
        else if(!strcmp(head,"HMSCORE_CHB"))
        {
            sscanf(line,"%*s%f", &hmscore_chb);
        }
        else if(!strcmp(head,"HMSCORE_CHM"))
        {
            sscanf(line,"%*s%f", &hmscore_chm);
        }
        else if(!strcmp(head,"HMSCORE_CRT"))
        {
            sscanf(line,"%*s%f", &hmscore_crt);
        }
        else if(!strcmp(head,"HMSCORE_C0"))
        {
            sscanf(line,"%*s%f", &hmscore_c0);
        }
        if(!strcmp(head,"APPLY_HSSCORE"))
        {
            sscanf(line,"%*s%s",apply_hsscore);
        }
        else if(!strcmp(head,"HSSCORE_CVDW"))
        {
            sscanf(line,"%*s%f", &hsscore_cvdw);
        }
        else if(!strcmp(head,"HSSCORE_CHB"))
        {
            sscanf(line,"%*s%f", &hsscore_chb);
        }
        else if(!strcmp(head,"HSSCORE_CHS"))
        {
            sscanf(line,"%*s%f", &hsscore_chs);
        }
        else if(!strcmp(head,"HSSCORE_CRT"))
        {
            sscanf(line,"%*s%f", &hsscore_crt);
        }
        else if(!strcmp(head,"HSSCORE_C0"))
        {
            sscanf(line,"%*s%f", &hsscore_c0);
        }

        if(!strcmp(head,"SHOW_ATOM_BIND_SCORE"))
        {
            sscanf(line,"%*s%s",show_abs);
        }
    }

    num_error=0;

    if(!strcmp(receptor_file,"none"))
    {
        Missing_Parameter_Error("RECEPTOR_PDB_FILE"); num_error++;
    }
    if(!strcmp(ligand_file,"none"))
    {
        Missing_Parameter_Error("LIGAND_MOL2_FILE"); num_error++;
    }
    if(!strcmp(output_file,"none"))
    {
        Missing_Parameter_Error("OUTPUT_TABLE_FILE"); num_error++;
    }
    if(!strcmp(parameter_dir,"none"))
    {
        Missing_Parameter_Error("PARAMETER_DIRECTORY"); num_error++;
    }

    if(!strncasecmp(show_abs,"Y",1)) strcpy(show_abs,"YES");
    else if(!strncasecmp(show_abs,"N",1)) strcpy(show_abs,"NO");
    else {Invalid_Parameter_Error("SHOW_ATOM_BIND_SCORE"); num_error++;}

    if((num_hits>0)&&(!strcmp(hits_dir,"none")))
    {
        Missing_Parameter_Error("HITS_DIRECTORY"); num_error++;
    }

    if(num_hits>0) Check_Directory(hits_dir,TRUE);
    else Check_Directory(hits_dir,FALSE);

    num_method=0;

    if(!strncasecmp(apply_hpscore,"Y",1))
    {
        strcpy(apply_hpscore,"YES"); num_method++;
    }
    else if(!strncasecmp(apply_hpscore,"N",1))
    {
        strcpy(apply_hpscore,"NO");
    }
    else
    {
        Invalid_Parameter_Error("APPLY_HPSCORE");
        num_error++;
    }

    if(!strncasecmp(apply_hmscore,"Y",1))
    {
        strcpy(apply_hmscore,"YES");
        num_method++;
    }
    else if(!strncasecmp(apply_hmscore,"N",1))
    {
        strcpy(apply_hmscore,"NO");
    }
    else
    {
        Invalid_Parameter_Error("APPLY_HMSCORE");
        num_error++;
    }

    if(!strncasecmp(apply_hsscore,"Y",1))
    {
        strcpy(apply_hsscore,"YES");
        num_method++;
    }
    else if(!strncasecmp(apply_hsscore,"N",1))
    {
        strcpy(apply_hsscore,"NO");
    }
    else
    {
        Invalid_Parameter_Error("APPLY_HSSCORE");
        num_error++;
    }

    if(num_method==0)
    {
        puts("Error: no scoring function has been chosen.");
        num_error++;
    }

    if(num_error!=0)
    {
        printf("\n");
        printf("%d errors have been detected in the parameter file.\n",num_error);
        printf("Please correct them and try again.\n");
        exit(1);
    }

    fclose(fp); return;
}

void SCORE_Input::Show_Contents() const
{
    printf("RECEPTOR_PDB_FILE = %s\n",receptor_file);
    printf("COFACTOR_MOL2_FILE = %s\n",cofactor_file);
    printf("REFERENCE_MOL2_FILE = %s\n",reference_file);
    printf("LIGAND_MOL2_FILE = %s\n",ligand_file);
    printf("OUTPUT_TABLE_FILE = %s\n",output_file);
    printf("OUTPUT_LOG_FILE = %s\n",log_file);

    printf("PARAMETER_DIRECTORY = %s\n",parameter_dir);

    printf("NUMBER_OF_HITS = %d\n",num_hits);
    printf("HITS_DIRECTORY = %s\n",hits_dir);
    printf("SHOW_ATOM_BIND_SCORE = %s\n",show_abs);

    printf("APPLY_HPSCORE = %s\n",apply_hpscore);
    printf("APPLY_HMSCORE = %s\n",apply_hmscore);
    printf("APPLY_HSSCORE = %s\n",apply_hsscore);
    printf("Number of methods = %d\n", num_method);

    return;
}

PRESCREEN_Input::PRESCREEN_Input()
{
    strcpy(function,"PRESCREEN");

    strcpy(apply_chemical_rules,"YES");
    max_weight=600.0; min_weight=200.0;
    max_logp=6.00; min_logp=1.00;
    max_hb_atom=8; min_hb_atom=2;
    max_rotor=10; min_rotor=0;
}

PRESCREEN_Input::~PRESCREEN_Input()
{
}

void PRESCREEN_Input::Read_Inputs(char *filename)
{
    FILE *fp;
    char line[256],head[256];
    int num_error;

    if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

    for(;;)
    {
        if(fgets(line,256,fp)==NULL) break;
        else if(line[0]=='#') continue;
        else if(Blank_Line_Check(line)==TRUE) continue;
        else sscanf(line,"%s",head);

        if(!strcmp(head,"INPUT_MOL2_FILE"))
        {
            sscanf(line,"%*s%s",input_file);
        }
        else if(!strcmp(head,"OUTPUT_MOL2_FILE"))
        {
            sscanf(line,"%*s%s",output_file);
        }
        else if(!strcmp(head,"OUTPUT_LOG_FILE"))
        {
            sscanf(line,"%*s%s",log_file);
        }
        else if(!strcmp(head,"APPLY_CHEMICAL_RULES"))
        {
            sscanf(line,"%*s%s",apply_chemical_rules);
        }
        else if(!strcmp(head,"MAXIMAL_MOLECULAR_WEIGHT"))
        {
            sscanf(line,"%*s%f",&max_weight);
        }
        else if(!strcmp(head,"MINIMAL_MOLECULAR_WEIGHT"))
        {
            sscanf(line,"%*s%f",&min_weight);
        }
        else if(!strcmp(head,"MAXIMAL_LOGP"))
        {
            sscanf(line,"%*s%f",&max_logp);
        }
        else if(!strcmp(head,"MINIMAL_LOGP"))
        {
            sscanf(line,"%*s%f",&min_logp);
        }
        else if(!strcmp(head,"MAXIMAL_HB_ATOM"))
        {
            sscanf(line,"%*s%d",&max_hb_atom);
        }
        else if(!strcmp(head,"MINIMAL_HB_ATOM"))
        {
            sscanf(line,"%*s%d",&min_hb_atom);
        }
        else if(!strcmp(head,"MAXIMAL_ROTOR"))
        {
            sscanf(line,"%*s%d",&max_rotor);
        }
        else if(!strcmp(head,"MINIMAL_ROTOR"))
        {
            sscanf(line,"%*s%d",&min_rotor);
        }
        else continue;
    }

    fclose(fp);

    num_error=0;

    if(!strcmp(input_file,"none"))
    {
        Missing_Parameter_Error("INPUT_MOL2_FILE"); num_error++;
    }

    if(!strcmp(output_file,"none"))
    {
        Missing_Parameter_Error("OUTPUT_MOL2_FILE"); num_error++;
    }

    if(!strcmp(parameter_dir,"none"))
    {
        Missing_Parameter_Error("PARAMETER_DIRECTORY"); num_error++;
    }

    if(!strncasecmp(apply_chemical_rules,"Y",1))
    {
        strcpy(apply_chemical_rules,"YES");
    }
    else if(!strncasecmp(apply_chemical_rules,"N",1))
    {
        strcpy(apply_chemical_rules,"NO");
    }
    else
    {
        Invalid_Parameter_Error("APPLY_CHEMICAL_RULES"); num_error++;
    }

    if(max_weight<min_weight)
    {
        puts("Error: MAXIMAL_MOLECULAR_WEIGHT < MINIMAL_MOLECULAR_WEIGHT");
        num_error++;
    }

    if(max_logp<min_logp)
    {
        puts("Error: MAXIMAL_LOGP < MINIMAL_LOGP");
        num_error++;
    }

    if(max_hb_atom<min_hb_atom)
    {
        puts("Error: MAXIMAL_HB_ATOM < MINIMAL_HB_ATOM");
        num_error++;
    }

    if(max_rotor<min_rotor)
    {
        puts("Error: MAXIMAL_ROTOR < MINIMAL_ROTOR");
        num_error++;
    }

    if(num_error!=0)
    {
        printf("\n");
        printf("%d errors have been detected in the input file.\n",num_error);
        printf("Please correct them and try again.\n");
        exit(1);
    }

    return;
}

void PRESCREEN_Input::Show_Contents() const
{
    printf("INPUT_MOL2_FILE = %s\n",input_file);
    printf("OUTPUT_MOL2_FILE = %s\n",output_file);
    printf("PARAMETER_DIRECTORY = %s\n",parameter_dir);

    printf("APPLY_CHEMICAL_RULES = %s\n",apply_chemical_rules);
    printf("MAXIMAL_MOLECULAR_WEIGHT = %6.1f\n",max_weight);
    printf("MINIMAL_MOLECULAR_WEIGHT = %6.1f\n",min_weight);
    printf("MAXIMAL_LOGP = %-6.2f\n",max_logp);
    printf("MINIMAL_LOGP = %-6.2f\n",min_logp);
    printf("MAXIMAL_HB_ATOM = %d\n",max_hb_atom);
    printf("MINIMAL_HB_ATOM = %d\n",min_hb_atom);

    return;
}

LOGP_Input::LOGP_Input()
{
    strcpy(function,"LOGP");

    strcpy(calculate_logp,"YES");
    strcpy(calculate_mw,"YES");
    strcpy(count_hb_atom,"YES");
    strcpy(count_rotor,"YES");
}

LOGP_Input::~LOGP_Input()
{
}

void LOGP_Input::Read_Inputs(char *filename)
{
    FILE *fp;
    char line[256],head[256];
    int num_error;

    if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

    for(;;)
    {
        if(fgets(line,256,fp)==NULL) break;
        else if(line[0]=='#') continue;
        else if(Blank_Line_Check(line)==TRUE) continue;
        else sscanf(line,"%s",head);

        if(!strcmp(head,"INPUT_MOL2_FILE"))
        {
            sscanf(line,"%*s%s", input_file);
        }
        else if(!strcmp(head,"OUTPUT_LOG_FILE"))
        {
            sscanf(line,"%*s%s", output_file);
        }
        else if(!strcmp(head,"CALCULATE_LOGP"))
        {
            sscanf(line,"%*s%s", calculate_logp);
        }
        else if(!strcmp(head,"CALCULATE_MW"))
        {
            sscanf(line,"%*s%s", calculate_mw);
        }
        else if(!strcmp(head,"COUNT_HB_ATOM"))
        {
            sscanf(line,"%*s%s", count_hb_atom);
        }
        else if(!strcmp(head,"COUNT_ROTOR"))
        {
            sscanf(line,"%*s%s", count_rotor);
        }
        else continue;
    }

    fclose(fp);

    num_error=0;

    if(!strcmp(input_file,"none"))
    {
        Missing_Parameter_Error("INPUT_MOL2_FILE"); num_error++;
    }

    if(!strcmp(parameter_dir,"none"))
    {
        Missing_Parameter_Error("PARAMETER_DIRECTORY"); num_error++;
    }

    if(!strncasecmp(calculate_logp,"N",1)) strcpy(calculate_logp,"NO");
    else strcpy(calculate_logp,"YES");

    if(!strncasecmp(calculate_mw,"N",1)) strcpy(calculate_mw,"NO");
    else strcpy(calculate_mw,"YES");

    if(!strncasecmp(count_hb_atom,"N",1)) strcpy(count_hb_atom,"NO");
    else strcpy(count_hb_atom,"YES");

    if(!strncasecmp(count_rotor,"N",1)) strcpy(count_rotor,"NO");
    else strcpy(count_rotor,"YES");

    if(num_error!=0)
    {
        printf("\n");
        printf("%d errors have been detected in the input file.\n",num_error);
        printf("Please correct them and try again.\n");
        exit(1);
    }

    return;
}

void LOGP_Input::Show_Contents() const
{
    printf("INPUT_MOL2_FILE = %s\n",input_file);
    printf("OUTPUT_LOG_FILE = %s\n",output_file);
    printf("PARAMETER_DIRECTORY = %s\n",parameter_dir);

    return;
}
