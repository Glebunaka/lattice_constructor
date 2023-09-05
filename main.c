#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#define Ar 6
#define MAX 4

struct Atom{
    int N;
    int Type;
    int Num1;
    char Num2[4];
    double Pos[3];
    double Mass;
    double Q;
};
struct Dia{
    int a,b,c,d,e,f,in,out;
    char F[MAX][4];
    int F1[MAX];
    double M[MAX],B[3],C[3],D[3];
    double aa;
    int Ntypes[MAX];
};
void Dialog(int T, struct Dia **tm);
void Lattice(struct Atom **pm, struct Dia **tm);
void Printer(struct Atom *pm, struct Dia *tm1);
double Table(char *pst, int *ps2);
int clear_input_buffer(void);
void Reader(struct Atom **pm, struct Dia **tm);
void Charge (struct Atom **pm,struct Dia **tm);
void Charge2 (struct Atom **pm,struct Dia **tm);
void Charge3 (struct Atom **pm,struct Dia **tm);


int main(void){
	struct Atom *pp;
	struct Dia *k;
	k=malloc(sizeof(struct Dia));
	/*char El[4];
	int d;
	scanf("%s",El);
	printf("%lf",Table(El,&d));
	printf("\n%d\n",d);
	puts(El);*/
    Dialog(0,&k);
	Lattice(&pp, &k);
	Charge3(&pp, &k);
	/*Charge2(&pp, &k);*/
	Dialog(2,&k);
	Printer(pp,k);
	free(pp);
	free(k);
	/*Dialog(1,&k);
	Reader(&pp,&k);
	Dialog(2,&k);
	Printer(pp,k);
	free(pp);
	free(k);*/
	return 0;
}
void Charge (struct Atom **pm, struct Dia **tm){
    int n,p,z;
    double k;
    const double L=RAND_MAX;
    p=0;
    for (n=0;n<(**tm).e;n++){
        if((*(*pm+n)).Type==1){
            (*(*pm+n)).Q=-2;
        }
        else if ((*(*pm+n)).Type==3){
            (*(*pm+n)).Q=3;
        }
        else{
            k=rand()/L;
            /*printf("\nk=%lf",k);*/
            if (k<0.5){
                (*(*pm+n)).Q=2;
                p++;
            }
            else{
                (*(*pm+n)).Q=3;
            };
        }
    };
    /*printf("\n1=%d\n2=%d\n3=%d\np=%d",(**tm).Ntypes[0],(**tm).Ntypes[1],(**tm).Ntypes[2],p);*/
    z=((-2)*(**tm).Ntypes[0]+p*2+((**tm).Ntypes[1]-p)*3+(**tm).Ntypes[2]*3);
    /*printf("\nz=%d",z);*/
    while (z!=0){
        int r = rand() % (**tm).e;
        /*printf("\nr=%d",r);*/
        if ((*(*pm+r)).Type==2){
            if (z<0){
                if((*(*pm+r)).Q==2){
                    (*(*pm+r)).Q=3;
                    z++;
                };
            }
            else{
                if((*(*pm+r)).Q==3){
                    (*(*pm+r)).Q=2;
                    z--;
                };
            };
        };
        /*printf("\nz=%d",z);*/
    };
    for (n=0;n<(**tm).e;n++){
        if((*(*pm+n)).Q==3){
            (*(*pm+n)).Type=3;
            (**tm).Ntypes[1]--;
            (**tm).Ntypes[2]++;
        };
    };
    z=0;
    for (n=0;n<(**tm).e;n++){
        z=z+(*(*pm+n)).Q;
        printf("\nz=%d",z);
    };
}

void Charge2 (struct Atom **pm, struct Dia **tm){
    int n;
    for (n=0;n<(**tm).e;n++){
        if((*(*pm+n)).Type==1){
            (*(*pm+n)).Q=-2;
        }
        else if ((*(*pm+n)).Type==3){
            (*(*pm+n)).Q=3;
        }
        else{
            (*(*pm+n)).Q=2.5;
        };
    };
}

void Charge3 (struct Atom **pm, struct Dia **tm){
    int n;
    for (n=0;n<(**tm).e;n++){
        (*(*pm+n)).Q=0;
    };
}

void Reader(struct Atom **pm, struct Dia **tm){
    struct Atom *q;
    int d;
    d=(**tm).in;
    if (d==1){
        int k,m,ch,n=0;
        char T[10];
        FILE *fp;
        fp = fopen("PDB1.pdb","r");
        if ( (fp = fopen("PDB1.pdb","r")) == NULL ) {
            printf("\nFile Error\n");
            exit(1);
        };
        (**tm).B[0]=0;
        (**tm).C[1]=0;
        (**tm).D[2]=0;
        q=malloc(sizeof(struct Atom));
        fscanf(fp,"%s",T);
        while (strcmp(T,"END")){
            if (!strcmp(T,"REMARK")){
                while(!strcmp(T,"REMARK")){
                    fscanf(fp,"%s",T);
                    if (!strcmp(T,"285")){
                        fscanf(fp,"%lf",&(**tm).aa);
                        while (((ch = getc(fp)) != '\n'));
                        fseek(fp,10,SEEK_CUR);
                        fscanf(fp,"%lf %lf %lf",&(**tm).B[0],&(**tm).B[1],&(**tm).B[2]);
                        while (((ch = getc(fp)) != '\n'));
                        fseek(fp,10,SEEK_CUR);
                        fscanf(fp,"%lf %lf %lf",&(**tm).C[0],&(**tm).C[1],&(**tm).C[2]);
                        while (((ch = getc(fp)) != '\n'));
                        fseek(fp,10,SEEK_CUR);
                        fscanf(fp,"%lf %lf %lf",&(**tm).D[0],&(**tm).D[1],&(**tm).D[2]);
                        while (((ch = getc(fp)) != '\n'));
                        fscanf(fp,"%s",T);
                    }
                    else{
                        while (((ch = getc(fp)) != '\n'));
                        fscanf(fp,"%s",T);
                    };
                };
            }
            else if (!strcmp(T,"ATOM")){
                int s=0,u;
                while(!strcmp(T,"ATOM")){
                    int ss=0;
                    q=realloc(q,(n+1)*sizeof(struct Atom));
                    fscanf(fp,"%5d  %s",&k,(*(q+n)).Num2);
                    (*(q+n)).Mass=Table((*(q+n)).Num2,&(*(q+n)).Num1);
                    if (n==0){
                        strcpy((**tm).F[0],(*(q+n)).Num2);
                        (*(q+n)).Type=1;
                        (**tm).Ntypes[0]=1;
                        (**tm).f=1;
                        (**tm).M[0]=(*(q+n)).Mass;
                        (**tm).F1[0]=(*(q+n)).Num1;
                        s=1;
                    };
                    for (u=0;u<s;u++){
                        if (!strcmp((*(q+n)).Num2,(**tm).F[u])){
                            (**tm).Ntypes[u]++;
                            (*(q+n)).Type=u+1;
                            ss++;
                        };
                    };
                    if (ss==0){
                        (**tm).Ntypes[s]=1;
                        strcpy((**tm).F[s],(*(q+n)).Num2);
                        (**tm).M[s]=(*(q+n)).Mass;
                        (**tm).F1[s]=(*(q+n)).Num1;
                        (**tm).f++;
                        s++;
                        (*(q+n)).Type=s;
                    };
                    fseek(fp,10,SEEK_CUR);
                    fscanf(fp,"%d%lf%lf%lf",&m,&(*(q+n)).Pos[0],&(*(q+n)).Pos[1],&(*(q+n)).Pos[2]);
                    if((*(q+n)).Pos[0]>(**tm).B[0])
                        (**tm).B[0]=(*(q+n)).Pos[0];
                    if((*(q+n)).Pos[1]>(**tm).C[1])
                        (**tm).C[1]=(*(q+n)).Pos[1];
                    if((*(q+n)).Pos[2]>(**tm).D[2])
                        (**tm).D[2]=(*(q+n)).Pos[2];
                    (*(q+n)).N=(m-1)*99999+k;
                    /*printf("%s %d %lf %lf %lf\n",(*(q+n)).Num2,(*(q+n)).N,(*(q+n)).Pos[0],(*(q+n)).Pos[1],(*(q+n)).Pos[2]);*/
                    while (((ch = getc(fp)) != '\n'));
                    fscanf(fp,"%s",T);
                    n++;
                }
                (**tm).Ntypes[0]--;
            }
            else{
                while (((ch = getc(fp)) != '\n'));
                fscanf(fp,"%s",T);
            };
        };
        (**tm).e=n;
        fclose(fp);
    }
    else if (d==2){
        int k=0,m,ch,n=0,p=0;
        char CH;
        char T[20];
        FILE *fp;
        fp = fopen("Poscar1","r");
        if ( (fp = fopen("Poscar1","r")) == NULL ) {
            printf("\nFile Error\n");
            exit(1);
        };
        (**tm).B[0]=0;
        (**tm).C[1]=0;
        (**tm).D[2]=0;
        q=malloc(sizeof(struct Atom));
        while (((ch = getc(fp)) != '\n'));
        fscanf(fp,"%lf",&(**tm).aa);
        while (((ch = getc(fp)) != '\n'));
        fscanf(fp,"%lf %lf %lf",&(**tm).B[0],&(**tm).B[1],&(**tm).B[2]);
        while (((ch = getc(fp)) != '\n'));
        fscanf(fp,"%lf %lf %lf",&(**tm).C[0],&(**tm).C[1],&(**tm).C[2]);
        while (((ch = getc(fp)) != '\n'));
        fscanf(fp,"%lf %lf %lf",&(**tm).D[0],&(**tm).D[1],&(**tm).D[2]);
        while (((ch = getc(fp)) != '\n'));
        CH=getc(fp);
        if (isdigit(CH)){
            fseek(fp,-1,SEEK_CUR);
            do{
                fscanf(fp,"%d",&(**tm).Ntypes[k]);
                (**tm).F[k][0]=(k+1);
                (**tm).M[k]=Table((**tm).F[k],&(**tm).F1[k]);
                k++;
            }while ((ch=getc(fp))!='\n');
            (**tm).f=k;
        }
        else{
            fseek(fp,-1,SEEK_CUR);
            do{
                fscanf(fp,"%s",(**tm).F[k]);
                printf("%s",(**tm).F[k]);
                (**tm).M[k]=Table((**tm).F[k],&(**tm).F1[k]);
                k++;
            }while ((ch=getc(fp))!='\n');
            (**tm).f=k;
            for(m=0;m<k;m++){
                fscanf(fp,"%d",&(**tm).Ntypes[m]);
            }
            while (((ch = getc(fp)) != '\n'));
        };
        fscanf(fp,"%s",T);
        while (((ch = getc(fp)) != '\n'));
        for (m=0;m<(**tm).f;m++)
            p=p+(**tm).Ntypes[m];
        (**tm).e=p;
        q=realloc(q,p*sizeof(struct Atom));
        int v=0;
        for (m=0;m<(**tm).f;m++){
            for (n=0;n<(**tm).Ntypes[m];n++){
                fscanf(fp,"%lf %lf %lf\n",&(*(q+v)).Pos[0],&(*(q+v)).Pos[1],&(*(q+v)).Pos[2]);
                (*(q+v)).Mass=(**tm).M[m];
                (*(q+v)).Num1=(**tm).F1[m];
                strcpy((*(q+v)).Num2,(**tm).F[m]);
                (*(q+v)).Type=m+1;
                printf("\n%d",v);
                v++;
            };
            getc(fp);
        };
    }
    else if (d==3){
        int k=0,m,ch,n=0;
        char CH;
        char T[20];
        FILE *fp;
        fp = fopen("LAMMPS1.lmp","r");
        if ( (fp = fopen("LAMMPS1.lmp","r")) == NULL ) {
            printf("\nFile Error\n");
            exit(1);
        };
        (**tm).B[0]=0;
        (**tm).C[1]=0;
        (**tm).D[2]=0;
        q=malloc(sizeof(struct Atom));
        CH=getc(fp);
        if (isdigit(CH)){
            fseek(fp,-1,SEEK_CUR);
            fscanf(fp,"%lf",&(**tm).aa);
            while (((ch = getc(fp)) != '\n'));
        }
        else{
                printf("kk");
            fseek(fp,-1,SEEK_CUR);
            while (((ch = getc(fp)) != '\n'));
        };
        while (((ch = getc(fp)) != '\n'));
        while (k<2){
            fscanf(fp,"%d",&m);
            fscanf(fp,"%s",T);
            if (!strcmp(T,"atoms")){
                (**tm).e=m;
                k++;
                while (((ch = getc(fp)) != '\n'));
            }
            else if (!strcmp(T,"atom")){
                (**tm).f=m;
                k++;
                while (((ch = getc(fp)) != '\n'));
            }
            else
                while (((ch = getc(fp)) != '\n'));
        };
        while (((ch = getc(fp)) != '\n'));
        q=realloc(q,(**tm).e*sizeof(struct Atom));
        double h1,h2;
        CH=getc(fp);
        if (isdigit(CH)){
            fseek(fp,-1,SEEK_CUR);
            fscanf(fp,"%lf",&h1);
            fscanf(fp,"%lf",&h2);
            (**tm).B[0]=h2+fabs(h1);
            while (((ch = getc(fp)) != '\n'));
            fscanf(fp,"%lf",&h1);
            fscanf(fp,"%lf",&h2);
            (**tm).C[1]=h2+fabs(h1);
            while (((ch = getc(fp)) != '\n'));
            fscanf(fp,"%lf",&h1);
            fscanf(fp,"%lf",&h2);
            (**tm).D[2]=h2+fabs(h1);
            while (((ch = getc(fp)) != '\n'));
            (**tm).B[1]=0;
            (**tm).B[2]=0;
            (**tm).C[0]=0;
            (**tm).C[2]=0;
            (**tm).D[1]=0;
            (**tm).D[0]=0;
        }
        else{
            while (((ch = getc(fp)) != '\n'));
        };
        k=0;
        while (k<2){
            fscanf(fp,"%s\n",T);
            if (!strcmp(T,"Atoms")){
                int c;
                for(c=0;c<(**tm).e;c++){
                    int s,s1;
                    fscanf(fp,"%d",&s);
                    fscanf(fp,"%d",&s1);
                    fscanf(fp,"%lf",&(*(q+n)).Pos[0]);
                    fscanf(fp,"%lf",&(*(q+n)).Pos[1]);
                    fscanf(fp,"%lf",&(*(q+n)).Pos[2]);
                    (**tm).Ntypes[s1-1]++;
                    (*(q+n)).Mass=(**tm).M[s1-1];
                    (*(q+n)).Num1=(**tm).F1[s1-1];
                    strcpy((*(q+n)).Num2,(**tm).F[s1-1]);
                    (*(q+n)).Type=s1;
                    n++;
                };
                k++;
            }
            else if (!strcmp(T,"Masses")){
                int c;
                for(c=0;c<(**tm).f;c++){
                    int s;
                    double M;
                    fscanf(fp,"%d",&s);
                    fscanf(fp,"%lf",&M);
                    printf("\n%lf\n",M);
                    for(m=1;m<=103;m++){
                        char V[3],*v;
                        int y;
                        v=V;
                        sprintf(v,"%d",m);
                        if (fabs(M-Table(v,&y))<0.2){
                            (**tm).F1[s-1]=y;
                            strcpy((**tm).F[s-1],v);
                            (**tm).M[s-1]=M;
                        };
                    };
                };
                while (((ch = getc(fp)) != '\n'));
                k++;
            }
            else
                while (((ch = getc(fp)) != '\n'));
        };

    }
    *pm=q;

}
void Dialog(int T, struct Dia **tm){
    if (T==0){
        char str_array [Ar] [10]={"1:FCC","2:BCC","3:NaCl","4:CsCl","5:HCP","6:Spinel"};
        int n,t;
        struct Dia *q;
        q=malloc(sizeof(struct Dia));
        for (n=0;n<Ar;n++){
            puts(str_array [n]);
        };
        scanf("%d",&((*q).a));
        printf("Lattice parameter:\n");
        scanf("%lf",&((*q).aa));
        printf("Amount of cells along x y z:\n");
        scanf("%d %d %d",&((*q).b),&((*q).c),&((*q).d));
        (*q).B[0]=((*q).b)*((*q).aa);
        (*q).B[1]=0.0;
        (*q).B[2]=0.0;
        (*q).C[0]=0.0;
        (*q).C[1]=((*q).c)*((*q).aa);
        (*q).C[2]=0.0;
        (*q).D[0]=0.0;
        (*q).D[1]=0.0;
        (*q).D[2]=((*q).d)*((*q).aa);
        clear_input_buffer();
        for (t=1; t<=MAX; t++){
            printf("\nElement for type %d(Enter if there is no need for more types):\n",t);
            gets((*q).F[t-1]);
            if(!*(*q).F[t-1]) break;
        };
        (*q).f=t-1;
        *tm=q;
    }
    else if (T==1){
        printf("From:\n1:File(.pdb)\n2:File(VASP)\n3:File(Lammps)\n");
        scanf("%d",&(**tm).in);
    }
    else if (T==2){
        printf("To:\n1:File(.pdb)\n2:File(VASP)\n3:File(Lammps)\n4:Screen\n");
        scanf("%d",&(**tm).out);
    }
}

void Lattice(struct Atom **pm, struct Dia **tm){
    struct Atom *q;
    int d,x,y,z,X,Y,Z,N=0,zch,ych;
    double a;
    x=(**tm).b;
    y=(**tm).c;
    z=(**tm).d;
    d=(**tm).a;
    a=(**tm).aa;

    if (d==1){
        int gn;
        double M;
        M=Table((**tm).F[0],&gn);
        (**tm).F1[0]=gn;
        (**tm).M[0]=M;
        q=malloc(4*x*y*z*sizeof(struct Atom));
        for (Z=0;Z<(2*z);Z++){
            zch=Z%2;
            for (Y=0;Y<(2*y);Y++){
                ych=Y%2;
                for (X=0;X<x;X++){
                    if (((zch+ych)%2)==0) {
                        (*(q+N)).Pos[0]=(double)X*a;
                        (*(q+N)).Pos[1]=(double)Y*0.5*a;
                        (*(q+N)).Pos[2]=(double)Z*0.5*a;
                        (*(q+N)).Num1=gn;
                        strcpy((*(q+N)).Num2,(**tm).F[0]);
                        (*(q+N)).Type=1;
                        (*(q+N)).Mass=M;
                        N++;
                    }
                    else{
                        (*(q+N)).Pos[0]=((double)X+0.5)*a;
                        (*(q+N)).Pos[1]=(double)Y*0.5*a;
                        (*(q+N)).Pos[2]=(double)Z*0.5*a;
                        (*(q+N)).Num1=gn;
                        strcpy((*(q+N)).Num2,(**tm).F[0]);
                        (*(q+N)).Type=1;
                        (*(q+N)).Mass=M;
                        N++;
                    };
                };
            };
        };
    (**tm).Ntypes[0]=N;
    }
    else if (d==2){
        int gn;
        double M;
        M=Table((**tm).F[0],&gn);
        (**tm).F1[0]=gn;
        (**tm).M[0]=M;
        q=malloc(2*x*y*z*sizeof(struct Atom));
        for (Z=0;Z<(2*z);Z++){
            if (Z%2){
                for (Y=0;Y<y;Y++){
                    for (X=0;X<x;X++){
                        (*(q+N)).Pos[0]=((double)X+0.5)*a;
                        (*(q+N)).Pos[1]=((double)Y+0.5)*a;
                        (*(q+N)).Pos[2]=(double)Z*0.5*a;
                        (*(q+N)).Num1=gn;
                        strcpy((*(q+N)).Num2,(**tm).F[0]);
                        (*(q+N)).Type=1;
                        (*(q+N)).Mass=M;
                        N++;
                    };
                };
            }
            else{
                for (Y=0;Y<y;Y++){
                    for (X=0;X<x;X++){
                        (*(q+N)).Pos[0]=(double)X*a;
                        (*(q+N)).Pos[1]=(double)Y*a;
                        (*(q+N)).Pos[2]=(double)Z*0.5*a;
                        (*(q+N)).Num1=gn;
                        strcpy((*(q+N)).Num2,(**tm).F[0]);
                        (*(q+N)).Type=1;
                        (*(q+N)).Mass=M;
                        N++;
                    };
                };
            };
        };
    (**tm).Ntypes[0]=N;
    }
    else if (d==3){
        int gn1,gn2;
        double M1,M2;
        M1=Table((**tm).F[0],&gn1);
        M2=Table((**tm).F[1],&gn2);
        (**tm).F1[0]=gn1;
        (**tm).M[0]=M1;
        (**tm).F1[1]=gn2;
        (**tm).M[1]=M2;
        (**tm).Ntypes[0]=0;
        (**tm).Ntypes[1]=0;
        q=malloc(8*x*y*z*sizeof(struct Atom));
        for (Z=0;Z<(2*z);Z++){
            zch=Z%2;
            for (Y=0;Y<(2*y);Y++){
                ych=Y%2;
                for (X=0;X<2*x;X++){
                    if (((zch+ych)%2)==0) {
                        (*(q+N)).Pos[0]=(double)X*a*0.5;
                        (*(q+N)).Pos[1]=(double)Y*0.5*a;
                        (*(q+N)).Pos[2]=(double)Z*0.5*a;
                        if (X%2){
                            (*(q+N)).Num1=gn2;
                            strcpy((*(q+N)).Num2,(**tm).F[1]);
                            (*(q+N)).Type=2;
                            (*(q+N)).Mass=M2;
                            (**tm).Ntypes[1]++;
                        }
                        else{
                            (*(q+N)).Num1=gn1;
                            strcpy((*(q+N)).Num2,(**tm).F[0]);
                            (*(q+N)).Type=1;
                            (*(q+N)).Mass=M1;
                            (**tm).Ntypes[0]++;
                        };
                        N++;
                    }
                    else{
                        (*(q+N)).Pos[0]=(double)X*0.5*a;
                        (*(q+N)).Pos[1]=(double)Y*0.5*a;
                        (*(q+N)).Pos[2]=(double)Z*0.5*a;
                        if (X%2){
                            (*(q+N)).Num1=gn1;
                            strcpy((*(q+N)).Num2,(**tm).F[0]);
                            (*(q+N)).Type=1;
                            (*(q+N)).Mass=M1;
                            (**tm).Ntypes[0]++;
                        }
                        else{
                            (*(q+N)).Num1=gn2;
                            strcpy((*(q+N)).Num2,(**tm).F[1]);
                            (*(q+N)).Type=2;
                            (*(q+N)).Mass=M2;
                            (**tm).Ntypes[1]++;
                        };
                        N++;
                    };
                };
            };
        };
    }
    else if (d==4){
        int gn1,gn2;
        double M1,M2;
        M1=Table((**tm).F[0],&gn1);
        M2=Table((**tm).F[1],&gn2);
        (**tm).F1[0]=gn1;
        (**tm).M[0]=M1;
        (**tm).F1[1]=gn2;
        (**tm).M[1]=M2;
        (**tm).Ntypes[0]=0;
        (**tm).Ntypes[1]=0;
        q=malloc(4*x*y*z*sizeof(struct Atom));
        for (Z=0;Z<(2*z);Z++){
            if (Z%2){
                for (Y=0;Y<y;Y++){
                    for (X=0;X<x;X++){
                        (*(q+N)).Pos[0]=((double)X+0.5)*a;
                        (*(q+N)).Pos[1]=((double)Y+0.5)*a;
                        (*(q+N)).Pos[2]=(double)Z*0.5*a;
                        (*(q+N)).Num1=gn2;
                        strcpy((*(q+N)).Num2,(**tm).F[1]);
                        (*(q+N)).Type=2;
                        (*(q+N)).Mass=M2;
                        (**tm).Ntypes[1]++;
                        N++;
                    };
                };
            }
            else{
                for (Y=0;Y<y;Y++){
                    for (X=0;X<x;X++){
                        (*(q+N)).Pos[0]=(double)X*a;
                        (*(q+N)).Pos[1]=(double)Y*a;
                        (*(q+N)).Pos[2]=(double)Z*0.5*a;
                        (*(q+N)).Num1=gn1;
                        strcpy((*(q+N)).Num2,(**tm).F[0]);
                        (*(q+N)).Type=1;
                        (*(q+N)).Mass=M1;
                        (**tm).Ntypes[0]++;
                        N++;
                    };
                };
            };
        };
    }
    else if(d==5){
        int gn;
        double M;
        M=Table((**tm).F[0],&gn);
        (**tm).F1[0]=gn;
        (**tm).M[0]=M;
        q=malloc(2*x*y*z*sizeof(struct Atom));
        for (Z=0;Z<(2*z);Z++){
            if(!(Z%2)){
                for (Y=0;Y<(y);Y++){
                    if (Y%2){
                        for (X=0;X<x;X++){
                            (*(q+N)).Pos[0]=((double)X+0.5)*a;
                            (*(q+N)).Pos[1]=(double)Y*a*cos(3.14/6);
                            (*(q+N)).Pos[2]=(double)Z*0.5*a*1.63;
                            (*(q+N)).Num1=gn;
                            strcpy((*(q+N)).Num2,(**tm).F[0]);
                            (*(q+N)).Type=1;
                            (*(q+N)).Mass=M;
                            N++;
                        };
                    }
                    else{
                        for (X=0;X<x;X++){
                            (*(q+N)).Pos[0]=((double)X)*a;
                            (*(q+N)).Pos[1]=(double)Y*a*cos(3.14/6);
                            (*(q+N)).Pos[2]=(double)Z*0.5*a*1.63;
                            (*(q+N)).Num1=gn;
                            strcpy((*(q+N)).Num2,(**tm).F[0]);
                            (*(q+N)).Type=1;
                            (*(q+N)).Mass=M;
                            N++;
                        };
                    };
                };
            }
            else{
                for (Y=0;Y<(y);Y++){
                    if (Y%2){
                        for (X=0;X<x;X++){
                            (*(q+N)).Pos[0]=((double)X+0.5)*a;
                            (*(q+N)).Pos[1]=((double)Y+0.66667)*a*cos(3.14/6);
                            (*(q+N)).Pos[2]=(double)Z*0.5*a*1.63;
                            (*(q+N)).Num1=gn;
                            strcpy((*(q+N)).Num2,(**tm).F[0]);
                            (*(q+N)).Type=1;
                            (*(q+N)).Mass=M;
                            N++;
                        };
                    }
                    else{
                        for (X=0;X<x;X++){
                            (*(q+N)).Pos[0]=((double)X)*a;
                            (*(q+N)).Pos[1]=((double)Y+0.66667)*a*cos(3.14/6);
                            (*(q+N)).Pos[2]=(double)Z*0.5*a*1.63;
                            (*(q+N)).Num1=gn;
                            strcpy((*(q+N)).Num2,(**tm).F[0]);
                            (*(q+N)).Type=1;
                            (*(q+N)).Mass=M;
                            N++;
                        };
                    };
                };
            };
        };
    (**tm).Ntypes[0]=N;
    (**tm).C[1]=(y)*a*cos(3.14/6);
    (**tm).D[2]=z*a*1.63;
    }
    else if (d==6){
        int gn1,gn2,gn3;
        double M1,M2,M3;
        M1=Table((**tm).F[0],&gn1);
        M2=Table((**tm).F[1],&gn2);
        M3=Table((**tm).F[2],&gn3);
        (**tm).F1[0]=gn1;
        (**tm).M[0]=M1;
        (**tm).F1[1]=gn2;
        (**tm).M[1]=M2;
        (**tm).F1[2]=gn3;
        (**tm).M[2]=M3;
        (**tm).Ntypes[0]=0;
        (**tm).Ntypes[1]=0;
        (**tm).Ntypes[2]=0;
        q=malloc(56*x*y*z*sizeof(struct Atom));
        if (q==NULL) {
            printf("Memory Error");
        };
        for (Z=0;Z<4;Z++){
            zch=Z%2;
            for (Y=0;Y<4;Y++){
                ych=Y%2;
                for (X=0;X<2;X++){
                    if (((zch+ych)%2)==0) {
                        (*(q+N)).Pos[0]=(double)X*a*0.5;
                        (*(q+N)).Pos[1]=(double)Y*0.25*a;
                        (*(q+N)).Pos[2]=(double)Z*0.25*a;
                        (*(q+N)).Num1=gn1;
                        strcpy((*(q+N)).Num2,(**tm).F[0]);
                        (*(q+N)).Type=1;
                        (*(q+N)).Mass=M1;
                        N++;
                    }
                    else{
                        (*(q+N)).Pos[0]=((double)X+0.5)*0.5*a;
                        (*(q+N)).Pos[1]=(double)Y*0.25*a;
                        (*(q+N)).Pos[2]=(double)Z*0.25*a;
                        (*(q+N)).Num1=gn1;
                        strcpy((*(q+N)).Num2,(**tm).F[0]);
                        (*(q+N)).Type=1;
                        (*(q+N)).Mass=M1;
                        N++;
                    };
                };
            };
        };
        Z=0;
        for (X=0, Y=3;X<4;X++,Y--){
            (*(q+N)).Pos[0]=(double)X*a*0.25;
            (*(q+N)).Pos[1]=(double)Y*0.25*a;
            (*(q+N)).Pos[2]=(double)Z*0.25*a;
            (*(q+N)).Num1=gn2;
            strcpy((*(q+N)).Num2,(**tm).F[1]);
            (*(q+N)).Type=2;
            (*(q+N)).Mass=M2;
            N++;
        };
        Z=1;
        for (X=0,Y=0;X<4;X++,Y++){
            (*(q+N)).Pos[0]=(double)X*a*0.25;
            (*(q+N)).Pos[1]=(double)Y*0.25*a;
            (*(q+N)).Pos[2]=(double)Z*0.25*a;
            (*(q+N)).Num1=gn2;
            strcpy((*(q+N)).Num2,(**tm).F[1]);
            (*(q+N)).Type=2;
            (*(q+N)).Mass=M2;
            N++;
        };
        Z=2;
        for (X=0,Y=0;X<4;X++,Y++){
            if (!(X%2)){
                (*(q+N)).Pos[0]=(double)X*a*0.25;
                (*(q+N)).Pos[1]=(double)(0.25*Y+0.25)*a;
                (*(q+N)).Pos[2]=(double)Z*0.25*a;
                (*(q+N)).Num1=gn2;
                strcpy((*(q+N)).Num2,(**tm).F[1]);
                (*(q+N)).Type=2;
                (*(q+N)).Mass=M2;
                N++;
            }
            else{
                (*(q+N)).Pos[0]=(double)X*a*0.25;
                (*(q+N)).Pos[1]=(double)(0.25*Y-0.25)*a;
                (*(q+N)).Pos[2]=(double)Z*0.25*a;
                (*(q+N)).Num1=gn2;
                strcpy((*(q+N)).Num2,(**tm).F[1]);
                (*(q+N)).Type=2;
                (*(q+N)).Mass=M2;
                N++;
            };
        };
        Z=3;
        for (X=0,Y=3;X<4;X++,Y--){
            if (!(X%2)){
                (*(q+N)).Pos[0]=(double)X*a*0.25;
                (*(q+N)).Pos[1]=(double)(0.25*Y-0.25)*a;
                (*(q+N)).Pos[2]=(double)Z*0.25*a;
                (*(q+N)).Num1=gn2;
                strcpy((*(q+N)).Num2,(**tm).F[1]);
                (*(q+N)).Type=2;
                (*(q+N)).Mass=M2;
                N++;
            }
            else{
                (*(q+N)).Pos[0]=(double)X*a*0.25;
                (*(q+N)).Pos[1]=(double)(0.25*Y+0.25)*a;
                (*(q+N)).Pos[2]=(double)Z*0.25*a;
                (*(q+N)).Num1=gn2;
                strcpy((*(q+N)).Num2,(**tm).F[1]);
                (*(q+N)).Type=2;
                (*(q+N)).Mass=M2;
                N++;
            };
        };
        Z=0;
        for (X=1, Y=2;X<=2;X++,Y--){
            (*(q+N)).Pos[0]=((double)X*0.5-0.125)*a;
            (*(q+N)).Pos[1]=((double)Y*0.5-0.125)*a;
            (*(q+N)).Pos[2]=((double)Z*0.5+0.125)*a;
            (*(q+N)).Num1=gn3;
            strcpy((*(q+N)).Num2,(**tm).F[2]);
            (*(q+N)).Type=3;
            (*(q+N)).Mass=M3;
            N++;
        };
        Z=1;
        for (X=0, Y=1;X<2;X++,Y--){
            (*(q+N)).Pos[0]=((double)X*0.5+0.125)*a;
            (*(q+N)).Pos[1]=((double)Y*0.5+0.125)*a;
            (*(q+N)).Pos[2]=((double)Z*0.5-0.125)*a;
            (*(q+N)).Num1=gn3;
            strcpy((*(q+N)).Num2,(**tm).F[2]);
            (*(q+N)).Type=3;
            (*(q+N)).Mass=M3;
            N++;
        };
        Z=1;
        for (X=1, Y=1;X<=2;X++,Y++){
            (*(q+N)).Pos[0]=((double)X*0.5-0.125)*a;
            (*(q+N)).Pos[1]=((double)Y*0.5-0.125)*a;
            (*(q+N)).Pos[2]=((double)Z*0.5+0.125)*a;
            (*(q+N)).Num1=gn3;
            strcpy((*(q+N)).Num2,(**tm).F[2]);
            (*(q+N)).Type=3;
            (*(q+N)).Mass=M3;
            N++;
        };
        Z=2;
        for (X=0, Y=0;X<2;X++,Y++){
            (*(q+N)).Pos[0]=((double)X*0.5+0.125)*a;
            (*(q+N)).Pos[1]=((double)Y*0.5+0.125)*a;
            (*(q+N)).Pos[2]=((double)Z*0.5-0.125)*a;
            (*(q+N)).Num1=gn3;
            strcpy((*(q+N)).Num2,(**tm).F[2]);
            (*(q+N)).Type=3;
            (*(q+N)).Mass=M3;
            N++;
        };
        int L,h;
        L=N;
        N=0;
        for (Z=0;Z<z;Z++){
            for(Y=0;Y<y;Y++){
                for(X=0;X<x;X++){
                    for(h=0;h<L;h++){
                        (*(q+N)).Pos[0]=(*(q+h)).Pos[0]+X*a;
                        (*(q+N)).Pos[1]=(*(q+h)).Pos[1]+Y*a;
                        (*(q+N)).Pos[2]=(*(q+h)).Pos[2]+Z*a;
                        (*(q+N)).Num1=(*(q+h)).Num1;
                        strcpy((*(q+N)).Num2,(*(q+h)).Num2);
                        (*(q+N)).Type=(*(q+h)).Type;
                        (*(q+N)).Mass=(*(q+h)).Mass;
                        N++;
                        (**tm).Ntypes[(*(q+h)).Type-1]++;
                    };
                };
            };
        };
    };
    *pm=q;
    (**tm).e=N;
}

void Printer(struct Atom *pm, struct Dia *tm1){
    struct Atom *q;
    int d,k,K,m;
    K=(*tm1).e;
    q=pm;
    d=(*tm1).out;
    if (d==1){
        FILE *fp;
        fp = fopen("PDB.pdb","w");
        if ( (fp = fopen("PDB.pdb","w")) == NULL ) {
            printf("\nFile Error\n");
            exit(1);
        };
        int f,p;
        fprintf(fp,"\nREMARK 285 %.3lf\nREMARK 285 %.3lf %.3lf %.3lf\nREMARK 285 %.3lf %.3lf %.3lf\nREMARK 285 %.3lf %.3lf %.3lf\n", (*tm1).aa, (*tm1).B[0], (*tm1).B[1], (*tm1).B[2],(*tm1).C[0], (*tm1).C[1], (*tm1).C[2], (*tm1).D[0], (*tm1).D[1], (*tm1).D[2]);
        printf("%d",K);
        for(f=0;f<=9; f++){
            p=f*99999;
            for(k=0;((p+k<K)&&(k<99999));k++){
                fprintf(fp,"ATOM  %5d  %2s       %4d    %8.3lf%8.3lf%8.3lf\n",k+1,(*(q+k+p)).Num2 ,f+1,(*(q+k+p)).Pos[0],(*(q+k+p)).Pos[1],(*(q+k+p)).Pos[2]);};
            if (p+k==K-1) break;
        }
        fprintf(fp,"END");
    fclose(fp);}
    else if (d==2){
        FILE *fp;
        fp = fopen("POSCAR","w");
        if ( (fp = fopen("POSCAR","w")) == NULL ) {
            printf("\nFile Error\n");
            exit(1);
        };
        fprintf(fp,"Cubic");
        fprintf(fp,"\n%5.3lf\n",(*tm1).aa);
        fprintf(fp,"%.3lf %.3lf %.3lf\n%.3lf %.3lf %.3lf\n%.3lf %.3lf %.3lf\n", (*tm1).B[0], (*tm1).B[1], (*tm1).B[2],(*tm1).C[0], (*tm1).C[1], (*tm1).C[2], (*tm1).D[0], (*tm1).D[1], (*tm1).D[2]);
        for(k=0;k<(*tm1).f;k++)
            fprintf(fp,"%s ",(*tm1).F[k]);
        fprintf(fp,"\n");
        for(k=0;k<(*tm1).f;k++)
            fprintf(fp,"%d ",(*tm1).Ntypes[k]);
        fprintf(fp,"\n");
        fprintf(fp,"Cartesian\n");
        for(m=0;m<(*tm1).f;m++){
            for(k=0;k<K;k++){
                if ((*(q+k)).Type==m+1)
                    fprintf(fp,"%.2lf %.2lf %.2lf\n",(*(q+k)).Pos[0],(*(q+k)).Pos[1],(*(q+k)).Pos[2]);
            };
             fprintf(fp,"\n");
        };
    fclose(fp);
    }
    else if (d==3){
        FILE *fp;
        int r;
        fp = fopen("LAMMPS.lmp","w");
        if ( (fp = fopen("LAMMPS.lmp","w")) == NULL ) {
            printf("\nFile Error\n");
            exit(1);
        };
        fprintf(fp,"\n\n%d atoms\n",K);
        fprintf(fp,"%d atom types\n\n",(*tm1).f);
        fprintf(fp,"0.0 %.6lf xlo xhi\n",(*tm1).B[0]);
        fprintf(fp,"0.0 %.6lf ylo yhi\n",(*tm1).C[1]);
        fprintf(fp,"0.0 %.6lf zlo zhi\n",(*tm1).D[2]);
        fprintf(fp,"0.0 0.0 0.0 xy xz yz\n\n");
        fprintf(fp,"Masses\n");
        for (r=0;r<(*tm1).f;r++)
            fprintf(fp,"\n%d %.2lf",r+1,(*tm1).M[r]);
        fprintf(fp,"\n\nAtoms\n\n");
        for(k=0;k<K;k++)
            fprintf(fp,"%-6d %-2d %5.1lf %11.6lf %11.6lf %11.6lf\n",k+1,(*(q+k)).Type,(*(q+k)).Q,(*(q+k)).Pos[0],(*(q+k)).Pos[1],(*(q+k)).Pos[2]);
        fclose(fp);
    }
    else{
        printf("    N   Type   Num.           X           Y           Z\n");
        for(k=0;k<K;k++)
            printf("%5d|%6d|%6d|%11.5lf|%11.5lf|%11.5lf|\n",k+1,(*(q+k)).Type,(*(q+k)).Num1,(*(q+k)).Pos[0],(*(q+k)).Pos[1],(*(q+k)).Pos[2]);
    }
}

double Table(char *pst, int *ps2){
    if ((strcmp(pst,"H")==0)||(strcmp(pst,"1")==0)){
        strcpy(pst,"H");
        *ps2=1;
        return 1.008;
    }
    else if ((strcmp(pst,"He")==0)||(strcmp(pst,"2")==0)){
        strcpy(pst,"He");
        *ps2=2;
        return 4.003;
    }
    else if ((strcmp(pst,"Li")==0)||(strcmp(pst,"3")==0)){
        strcpy(pst,"Li");
        *ps2=3;
        return 6.941;
    }
    else if ((strcmp(pst,"Be")==0)||(strcmp(pst,"4")==0)){
        strcpy(pst,"Be");
        *ps2=4;
        return 9.012;
    }
    else if ((strcmp(pst,"B")==0)||(strcmp(pst,"5")==0)){
        strcpy(pst,"B");
        *ps2=5;
        return 10.811;
    }
    else if ((strcmp(pst,"C")==0)||(strcmp(pst,"6")==0)){
        strcpy(pst,"C");
        *ps2=6;
        return 12.011;
    }
    else if ((strcmp(pst,"N")==0)||(strcmp(pst,"7")==0)){
        strcpy(pst,"N");
        *ps2=7;
        return 14.007;
    }
    else if ((strcmp(pst,"O")==0)||(strcmp(pst,"8")==0)){
        strcpy(pst,"O");
        *ps2=8;
        return 15.999;
    }
    else if ((strcmp(pst,"F")==0)||(strcmp(pst,"9")==0)){
        strcpy(pst,"F");
        *ps2=9;
        return 18.998;
    }
    else if ((strcmp(pst,"Ne")==0)||(strcmp(pst,"10")==0)){
        strcpy(pst,"Ne");
        *ps2=10;
        return 20.179;
    }
    else if ((strcmp(pst,"Na")==0)||(strcmp(pst,"11")==0)){
        strcpy(pst,"Na");
        *ps2=11;
        return 22.990;
    }
    else if ((strcmp(pst,"Mg")==0)||(strcmp(pst,"12")==0)){
        strcpy(pst,"Mg");
        *ps2=12;
        return 24.305;
    }
    else if ((strcmp(pst,"Al")==0)||(strcmp(pst,"13")==0)){
        strcpy(pst,"Al");
        *ps2=13;
        return 26.980;
    }
    else if ((strcmp(pst,"Si")==0)||(strcmp(pst,"14")==0)){
        strcpy(pst,"Si");
        *ps2=14;
        return 28.086;
    }
    else if ((strcmp(pst,"P")==0)||(strcmp(pst,"15")==0)){
        strcpy(pst,"P");
        *ps2=15;
        return 30.974;
    }
    else if ((strcmp(pst,"S")==0)||(strcmp(pst,"16")==0)){
        strcpy(pst,"S");
        *ps2=16;
        return 32.066;
    }
    else if ((strcmp(pst,"Cl")==0)||(strcmp(pst,"17")==0)){
        strcpy(pst,"Cl");
        *ps2=17;
        return 35.453;
    }
    else if ((strcmp(pst,"Ar")==0)||(strcmp(pst,"18")==0)){
        strcpy(pst,"Ar");
        *ps2=18;
        return 39.948;
    }
    else if ((strcmp(pst,"K")==0)||(strcmp(pst,"19")==0)){
        strcpy(pst,"K");
        *ps2=19;
        return 39.098;
    }
    else if ((strcmp(pst,"Ca")==0)||(strcmp(pst,"20")==0)){
        strcpy(pst,"Ca");
        *ps2=20;
        return 40.078;
    }
    else if ((strcmp(pst,"Sc")==0)||(strcmp(pst,"21")==0)){
        strcpy(pst,"Sc");
        *ps2=21;
        return 44.956;
    }
    else if ((strcmp(pst,"Ti")==0)||(strcmp(pst,"22")==0)){
        strcpy(pst,"Ti");
        *ps2=22;
        return 47.880;
    }
    else if ((strcmp(pst,"V")==0)||(strcmp(pst,"23")==0)){
        strcpy(pst,"V");
        *ps2=23;
        return 50.942;
    }
    else if ((strcmp(pst,"Cr")==0)||(strcmp(pst,"24")==0)){
        strcpy(pst,"Cr");
        *ps2=24;
        return 51.996;
    }
    else if ((strcmp(pst,"Mn")==0)||(strcmp(pst,"25")==0)){
        strcpy(pst,"Mn");
        *ps2=25;
        return 54.938;
    }
    else if ((strcmp(pst,"Fe")==0)||(strcmp(pst,"26")==0)){
        strcpy(pst,"Fe");
        *ps2=26;
        return 55.847;
    }
    else if ((strcmp(pst,"Co")==0)||(strcmp(pst,"27")==0)){
        strcpy(pst,"Co");
        *ps2=27;
        return 58.933;
    }
    else if ((strcmp(pst,"Ni")==0)||(strcmp(pst,"28")==0)){
        strcpy(pst,"Ni");
        *ps2=28;
        return 58.690;
    }
    else if ((strcmp(pst,"Cu")==0)||(strcmp(pst,"29")==0)){
        strcpy(pst,"Cu");
        *ps2=29;
        return 63.546;
    }
    else if ((strcmp(pst,"Zn")==0)||(strcmp(pst,"30")==0)){
        strcpy(pst,"Zn");
        *ps2=30;
        return 65.390;
    }
    else if ((strcmp(pst,"Ga")==0)||(strcmp(pst,"31")==0)){
        strcpy(pst,"Ga");
        *ps2=31;
        return 69.723;
    }
    else if ((strcmp(pst,"Ge")==0)||(strcmp(pst,"32")==0)){
        strcpy(pst,"Ge");
        *ps2=32;
        return 72.590;
    }
    else if ((strcmp(pst,"As")==0)||(strcmp(pst,"33")==0)){
        strcpy(pst,"As");
        *ps2=33;
        return 74.922;
    }
    else if ((strcmp(pst,"Se")==0)||(strcmp(pst,"34")==0)){
        strcpy(pst,"Se");
        *ps2=34;
        return 78.960;
    }
    else if ((strcmp(pst,"Br")==0)||(strcmp(pst,"35")==0)){
        strcpy(pst,"Br");
        *ps2=35;
        return 79.904;
    }
    else if ((strcmp(pst,"Kr")==0)||(strcmp(pst,"36")==0)){
        strcpy(pst,"Kr");
        *ps2=36;
        return 83.800;
    }
    else if ((strcmp(pst,"Rb")==0)||(strcmp(pst,"37")==0)){
        strcpy(pst,"Rb");
        *ps2=37;
        return 85.468;
    }
    else if ((strcmp(pst,"Sr")==0)||(strcmp(pst,"38")==0)){
        strcpy(pst,"Sr");
        *ps2=38;
        return 87.620;
    }
    else if ((strcmp(pst,"Y")==0)||(strcmp(pst,"39")==0)){
        strcpy(pst,"Y");
        *ps2=39;
        return 88.906;
    }
    else if ((strcmp(pst,"Zr")==0)||(strcmp(pst,"40")==0)){
        strcpy(pst,"Zr");
        *ps2=40;
        return 91.224;
    }
    else if ((strcmp(pst,"Nb")==0)||(strcmp(pst,"41")==0)){
        strcpy(pst,"Nb");
        *ps2=41;
        return 92.906;
    }
    else if ((strcmp(pst,"Mo")==0)||(strcmp(pst,"42")==0)){
        strcpy(pst,"Mo");
        *ps2=42;
        return 95.940;
    }
    else if ((strcmp(pst,"Tc")==0)||(strcmp(pst,"43")==0)){
        strcpy(pst,"Tc");
        *ps2=43;
        return 97.907;
    }
    else if ((strcmp(pst,"Ru")==0)||(strcmp(pst,"44")==0)){
        strcpy(pst,"Ru");
        *ps2=44;
        return 101.070;
    }
    else if ((strcmp(pst,"Rh")==0)||(strcmp(pst,"45")==0)){
        strcpy(pst,"Rh");
        *ps2=45;
        return 102.906;
    }
    else if ((strcmp(pst,"Pd")==0)||(strcmp(pst,"46")==0)){
        strcpy(pst,"Pd");
        *ps2=46;
        return 106.420;
    }
    else if ((strcmp(pst,"Ag")==0)||(strcmp(pst,"47")==0)){
        strcpy(pst,"Ag");
        *ps2=47;
        return 107.868;
    }
    else if ((strcmp(pst,"Cd")==0)||(strcmp(pst,"48")==0)){
        strcpy(pst,"Cd");
        *ps2=48;
        return 112.410;
    }
    else if ((strcmp(pst,"In")==0)||(strcmp(pst,"49")==0)){
        strcpy(pst,"In");
        *ps2=49;
        return 114.820;
    }
    else if ((strcmp(pst,"Sn")==0)||(strcmp(pst,"50")==0)){
        strcpy(pst,"Sn");
        *ps2=50;
        return 118.710;
    }
    else if ((strcmp(pst,"Sb")==0)||(strcmp(pst,"51")==0)){
        strcpy(pst,"Sb");
        *ps2=51;
        return 121.750;
    }
    else if ((strcmp(pst,"Te")==0)||(strcmp(pst,"52")==0)){
        strcpy(pst,"Te");
        *ps2=52;
        return 127.600;
    }
    else if ((strcmp(pst,"I")==0)||(strcmp(pst,"53")==0)){
        strcpy(pst,"I");
        *ps2=53;
        return 126.905;
    }
    else if ((strcmp(pst,"Xe")==0)||(strcmp(pst,"54")==0)){
        strcpy(pst,"Xe");
        *ps2=54;
        return 131.290;
    }
    else if ((strcmp(pst,"Cs")==0)||(strcmp(pst,"55")==0)){
        strcpy(pst,"Cs");
        *ps2=55;
        return 132.905;
    }
    else if ((strcmp(pst,"Ba")==0)||(strcmp(pst,"56")==0)){
        strcpy(pst,"Ba");
        *ps2=56;
        return 137.330;
    }
    else if ((strcmp(pst,"La")==0)||(strcmp(pst,"57")==0)){
        strcpy(pst,"La");
        *ps2=57;
        return 138.906;
    }
    else if ((strcmp(pst,"Ce")==0)||(strcmp(pst,"58")==0)){
        strcpy(pst,"Ce");
        *ps2=58;
        return 140.120;
    }
    else if ((strcmp(pst,"Pr")==0)||(strcmp(pst,"59")==0)){
        strcpy(pst,"Pr");
        *ps2=59;
        return 140.908;
    }
    else if ((strcmp(pst,"Nd")==0)||(strcmp(pst,"60")==0)){
        strcpy(pst,"Nd");
        *ps2=60;
        return 144.240;
    }
    else if ((strcmp(pst,"Pm")==0)||(strcmp(pst,"61")==0)){
        strcpy(pst,"Pm");
        *ps2=61;
        return 144.913;
    }
    else if ((strcmp(pst,"Sm")==0)||(strcmp(pst,"62")==0)){
        strcpy(pst,"Sm");
        *ps2=62;
        return 150.360;
    }
    else if ((strcmp(pst,"Eu")==0)||(strcmp(pst,"63")==0)){
        strcpy(pst,"Eu");
        *ps2=63;
        return 151.960;
    }
    else if ((strcmp(pst,"Gd")==0)||(strcmp(pst,"64")==0)){
        strcpy(pst,"Gd");
        *ps2=64;
        return 157.250;
    }
    else if ((strcmp(pst,"Tb")==0)||(strcmp(pst,"65")==0)){
        strcpy(pst,"Tb");
        *ps2=65;
        return 158.925;
    }
    else if ((strcmp(pst,"Dy")==0)||(strcmp(pst,"66")==0)){
        strcpy(pst,"Dy");
        *ps2=66;
        return 162.500;
    }
    else if ((strcmp(pst,"Ho")==0)||(strcmp(pst,"67")==0)){
        strcpy(pst,"Ho");
        *ps2=67;
        return 164.930;
    }
    else if ((strcmp(pst,"Er")==0)||(strcmp(pst,"68")==0)){
        strcpy(pst,"Er");
        *ps2=68;
        return 167.260;
    }
    else if ((strcmp(pst,"Tm")==0)||(strcmp(pst,"69")==0)){
        strcpy(pst,"Tm");
        *ps2=69;
        return 168.934;
    }
    else if ((strcmp(pst,"Yb")==0)||(strcmp(pst,"70")==0)){
        strcpy(pst,"Yb");
        *ps2=70;
        return 173.040;
    }
    else if ((strcmp(pst,"Lu")==0)||(strcmp(pst,"71")==0)){
        strcpy(pst,"Lu");
        *ps2=71;
        return 174.967;
    }
    else if ((strcmp(pst,"Hf")==0)||(strcmp(pst,"72")==0)){
        strcpy(pst,"Hf");
        *ps2=72;
        return 178.490;
    }
    else if ((strcmp(pst,"Ta")==0)||(strcmp(pst,"73")==0)){
        strcpy(pst,"Ta");
        *ps2=73;
        return 180.948;
    }
    else if ((strcmp(pst,"W")==0)||(strcmp(pst,"74")==0)){
        strcpy(pst,"W");
        *ps2=74;
        return 183.850;
    }
    else if ((strcmp(pst,"Re")==0)||(strcmp(pst,"75")==0)){
        strcpy(pst,"Re");
        *ps2=75;
        return 186.207;
    }
    else if ((strcmp(pst,"Os")==0)||(strcmp(pst,"76")==0)){
        strcpy(pst,"Os");
        *ps2=76;
        return 190.200;
    }
    else if ((strcmp(pst,"Ir")==0)||(strcmp(pst,"77")==0)){
        strcpy(pst,"Ir");
        *ps2=77;
        return 192.220;
    }
    else if ((strcmp(pst,"Pt")==0)||(strcmp(pst,"78")==0)){
        strcpy(pst,"Pt");
        *ps2=78;
        return 195.080;
    }
    else if ((strcmp(pst,"Au")==0)||(strcmp(pst,"79")==0)){
        strcpy(pst,"Au");
        *ps2=79;
        return 196.967;
    }
    else if ((strcmp(pst,"Hg")==0)||(strcmp(pst,"80")==0)){
        strcpy(pst,"Hg");
        *ps2=80;
        return 200.590;
    }
    else if ((strcmp(pst,"Tl")==0)||(strcmp(pst,"81")==0)){
        strcpy(pst,"Tl");
        *ps2=81;
        return 204.383;
    }
    else if ((strcmp(pst,"Pb")==0)||(strcmp(pst,"82")==0)){
        strcpy(pst,"Pb");
        *ps2=1;
        return 207.200;
    }
    else if ((strcmp(pst,"Bi")==0)||(strcmp(pst,"83")==0)){
        strcpy(pst,"Bi");
        *ps2=83;
        return 208.980;
    }
    else if ((strcmp(pst,"Po")==0)||(strcmp(pst,"84")==0)){
        strcpy(pst,"Po");
        *ps2=84;
        return 208.982;
    }
    else if ((strcmp(pst,"At")==0)||(strcmp(pst,"85")==0)){
        strcpy(pst,"At");
        *ps2=85;
        return 209.987;
    }
    else if ((strcmp(pst,"Rn")==0)||(strcmp(pst,"86")==0)){
        strcpy(pst,"Rn");
        *ps2=86;
        return 222.018;
    }
    else if ((strcmp(pst,"Fr")==0)||(strcmp(pst,"87")==0)){
        strcpy(pst,"Fr");
        *ps2=87;
        return 223.020;
    }
    else if ((strcmp(pst,"Ra")==0)||(strcmp(pst,"88")==0)){
        strcpy(pst,"Ra");
        *ps2=88;
        return 226.025;
    }
    else if ((strcmp(pst,"Ac")==0)||(strcmp(pst,"89")==0)){
        strcpy(pst,"Ac");
        *ps2=89;
        return 227.028;
    }
    else if ((strcmp(pst,"Th")==0)||(strcmp(pst,"90")==0)){
        strcpy(pst,"Th");
        *ps2=90;
        return 232.038;
    }
    else if ((strcmp(pst,"Pa")==0)||(strcmp(pst,"91")==0)){
        strcpy(pst,"Pa");
        *ps2=91;
        return 231.036;
    }
    else if ((strcmp(pst,"U")==0)||(strcmp(pst,"92")==0)){
        strcpy(pst,"U");
        *ps2=92;
        return 238.029;
    }
    else if ((strcmp(pst,"Ne")==0)||(strcmp(pst,"93")==0)){
        strcpy(pst,"Ne");
        *ps2=93;
        return 237.048;
    }
    else if ((strcmp(pst,"Pu")==0)||(strcmp(pst,"94")==0)){
        strcpy(pst,"Pu");
        *ps2=94;
        return 244.064;
    }
    else if ((strcmp(pst,"Am")==0)||(strcmp(pst,"95")==0)){
        strcpy(pst,"Am");
        *ps2=95;
        return 243.061;
    }
    else if ((strcmp(pst,"Cm")==0)||(strcmp(pst,"96")==0)){
        strcpy(pst,"Cm");
        *ps2=96;
        return 247.070;
    }
    else if ((strcmp(pst,"Bk")==0)||(strcmp(pst,"97")==0)){
        strcpy(pst,"Bk");
        *ps2=97;
        return 247.070;
    }
    else if ((strcmp(pst,"Cf")==0)||(strcmp(pst,"98")==0)){
        strcpy(pst,"Cf");
        *ps2=98;
        return 251.080;
    }
    else if ((strcmp(pst,"Es")==0)||(strcmp(pst,"99")==0)){
        strcpy(pst,"Es");
        *ps2=99;
        return 252.083;
    }
    else if ((strcmp(pst,"Fm")==0)||(strcmp(pst,"100")==0)){
        strcpy(pst,"Fm");
        *ps2=100;
        return 257.095;
    }
    else if ((strcmp(pst,"Md")==0)||(strcmp(pst,"101")==0)){
        strcpy(pst,"Md");
        *ps2=101;
        return 258.099;
    }
    else if ((strcmp(pst,"No")==0)||(strcmp(pst,"102")==0)){
        strcpy(pst,"No");
        *ps2=102;
        return 259.101;
    }
    else if ((strcmp(pst,"Lr")==0)||(strcmp(pst,"103")==0)){
        strcpy(pst,"Lr");
        *ps2=103;
        return 260.105;
    }
    else{
        printf("There is no such element");
        strcpy(pst,"NaN");
        *ps2=0;
        return 0.0;
    };
}

int clear_input_buffer(void) {
    int ch;
    while (((ch = getchar()) != EOF) && (ch != '\n')) /* void */;
    return ch;
}
