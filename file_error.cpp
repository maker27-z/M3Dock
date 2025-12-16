#include <iostream>
#include <cstring>
#include <cstdlib>
#include "xtools.h"
# define TRUE 1
# define FALSE 0
int Blank_Line_Check(const char *line)
{
    int i,mark;
    int len=0;

    len=strlen(line);
    if(len<1) return TRUE;	// blank line!

    mark=0;

    for(i=0;i<len;i++)
    {
        if(line[i]==' ') continue;
        else if(line[i]=='\t') continue;
        else if(line[i]=='\n') break;
        else mark++;
    }

    if(mark==0) return TRUE;	// blank line
    else return FALSE;		// not a blank line
}

void Memory_Allocation_Error(char *position)
{
    printf("\n");
    printf("Memory allocation error!\n");
    printf("This is usually caused by memory overflow.\n");
    printf("Error happened at %s\n", position);
    exit(1);
}

void Open_File_Error(const char *filename)
{
    printf("\n");
    printf("Error: cannot open the file %s\n", filename);
    printf("Please make sure it exists.\n");
    exit(1);
}

void Read_File_Error(const char *filename)
{
    printf("\n");
    printf("Error: something wrong with %s\n", filename);
    printf("It may not have the correct format.\n");
    printf("Please check this file and try again.\n");
    exit(1);
}

void PDB_Format_Error(const char *filename)
{
    printf("\n");
    printf("Error: %s lacks necessary information.\n", filename);
    printf("This file may not be in PDB format.\n");
    printf("Please check it and try again.\n");
    exit(1);
}

void Mol2_Format_Error(const char *filename)
{
    printf("\n");
    printf("Error: %s lacks necessary information.\n", filename);
    printf("This file may not be in Mol2 format.\n");
    printf("Please check it and try again.\n");
    exit(1);
}

void Lig_Format_Error(const char *filename)
{
    printf("\n");
    printf("Error: %s lacks necessary information.\n", filename);
    printf("This file may not be in Lig format.\n");
    printf("Please check it and try again.\n");
    exit(1);
}

int Check_Mol2_File(const char *filename)
// return the number of molecules in the given file.
// if it is not a valid Mol2 file, return 0.
{
    FILE *fp;
    int count;
    char buf[256],head[256];

    if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

    count=0;

    for(;;)
    {
        if(fgets(buf,256,fp)==NULL) break;
        else {strcpy(head,""); sscanf(buf,"%s",head);}

        if(strcmp(head,"@<TRIPOS>MOLECULE")) continue;
        else count++;
    }

    fclose(fp);

    return count;
}


float Distance(const float a[3], const float b[3])
{
    double d,tmpx,tmpy,tmpz;

    tmpx=(a[0]-b[0])*(a[0]-b[0]);
    tmpy=(a[1]-b[1])*(a[1]-b[1]);
    tmpz=(a[2]-b[2])*(a[2]-b[2]);

    d=sqrt(tmpx+tmpy+tmpz);

    return (float)d;
}

char *Chomp(char *string)  // remove the new-line at the end of a string
{
    int i,len;

    len=strlen(string); if(len<1) return NULL;

    for(i=len-1;i>=0;i--)
    {
        if(string[i]!='\n') continue;
        else {string[i]='\0'; break;}
    }

    return string;
}

char *Get_Time()  // get local time in ascii format
{
    extern char *timestring;
    struct tm *ptr;
    time_t lt;

    lt=time(NULL); ptr=localtime(&lt);
    timestring=asctime(ptr);

    Chomp(timestring); return timestring;
}

float Angle_of_Two_Vectors(const float v1[3], const float v2[3])
{
    double angle;
    double l1,l2,tmp1,tmp2;

    l1=sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
    l2=sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);

    tmp1=v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
    tmp2=l1*l2;

    angle=acos(tmp1/tmp2); angle=angle/(PI)*180.0;

    return (float)angle;  // return angle in degree, 0-180
}

float Angle(const float a[3], const float b[3], const float c[3])
{
    int i;
    float angle,v1[3],v2[3];

    for(i=0;i<3;i++) {v1[i]=b[i]-a[i]; v2[i]=b[i]-c[i];}

    angle=Angle_of_Two_Vectors(v1,v2);

    return angle;  // the a-b-c angle in degree 0-180
}
