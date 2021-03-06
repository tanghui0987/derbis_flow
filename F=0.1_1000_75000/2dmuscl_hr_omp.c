#include<omp.h>
#include<math.h>
#include<stdio.h>
#include <time.h>
#include<stdlib.h>

#define FREE_ARG char*
#define NR_END 1

#define PI 3.141592653589793
#define sqrt2 1.414213562373
#define oneoversqrt2 0.707106781186

/*
INPUTS:
1. FILE 'input' formatted as row vector with the following parameter values:
[nx ny dx F tend a ad h0 mtstar0 rnum rint (# of particle size classes) Jentrain hcfriction depthdependentexponent pi Si Ki gi vegcover]
2. FILES for initial conditions:
    "./topoin" - initial topography
    "./solidin" - mask for solid boundaries inside computational domain
    "./roughnessin" - minimum manning roughness coefficient
    "./cmaskin" - mask for bedrock channels
    "./ksin" - saturated hydraulic conductivity
    "./theta0in" - initial volumetric soil moisture content
    "./thetasin" - volumetric soil moisture content at saturation
    "./hfin" - suction head (Green-Ampt)
    "./vinfin" - initial volume infiltrated (m)
    "./erodibilitymaskin" - mask for spatially variable soil erodibility
    "./omegacin" - critical stream power for entrainment (slope dependent, but based on d50)
    "./rain" - rainfall intensity
    "./raindiam" - raindrop diameter 
    "./rainvel" - raindrop velocity
    "./particlesizein" - vector of representative particle diameter for each size class
    "./particlepercentin" - fraction of particle size distribution in each class (vector)
 */

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

double *vector(nl,nh)
long nh,nl;
/* allocate a double vector with subscript range v[nl..nh] */
{
        double *v;

        v=(double *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(double)));
        return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

void free_vector(double *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

double **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	m += NR_END;
	m -= nrl;

	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	m += NR_END;
	m -= nrl;

	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	return m;
}

double ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;

	/* allocate pointers to pointers to rows */
	t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

double ****f4tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh, long nwl, long nwh)
/* allocate a double 4tensor with range t[nrl..nrh][ncl..nch][ndl..ndh][nwl..nwh] */
{
	long i,j,k, nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1,nwid=nwh-nwl+1 ;
	double ****t;

	/* allocate pointers to pointers to pointers to rows */
	t=(double ****) malloc((size_t)((nrow+NR_END)*sizeof(double***)));
	if (!t) nrerror("allocation failure 1 in f4tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to pointers to rows and set pointers to them */
	t[nrl]=(double ***) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double**)));
	if (!t[nrl]) nrerror("allocation failure 2 in f4tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl][ncl]=(double **) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double*)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f4tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl][ndl]=(double *) malloc((size_t)((nrow*ncol*ndep*nwid+NR_END)*sizeof(double)));
	if (!t[nrl][ncl][ndl]) nrerror("allocation failure 4 in f4tensor()");
	t[nrl][ncl][ndl] += NR_END;
	t[nrl][ncl][ndl] -= nwl;

    for(i=nrl;i<=nrh;i++)
	{
		if (i > nrl)
		{
			t[i] = t[i-1] + ncol ;
		    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		    t[i][ncl][ndl] = t[i-1][ncl][ndl] + ncol*ndep*nwid ;
		}
		for(j=ncl;j<=nch;j++)
		{
			if (j > ncl)
			{
				t[i][j]=t[i][j-1] + ndep ;
				t[i][j][ndl] = t[i][j-1][ndl] + ndep*nwid ;
			}

			for(k=ndl;k<=ndh;k++)
			{
				if (k > ndl) t[i][j][k] = t[i][j][k-1] + nwid ;
			}
		}
	}

	/* return pointer to pointer to array of pointers to rows */
	return t;
}

void free_f4tensor(double ****t, long nrl, long nrh, long ncl, long nch, 
	long ndl, long ndh, long nwl, long nwh)
/* free a double f4tensor allocated by f4tensor() */
{
	free((FREE_ARG) (t[nrl][ncl][ndl]+nwl-NR_END));
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

void initialcondition(double **TOPO,double ***VEL,double ***U,double ***VELOLD,double ***UOLD,double ***UOLD2,double **DEPTHINIT,double **TOPOINIT,double ***SED,double **RHO,double **MANNING,double **ROUGHNESS,int **CMASK,double ***BETA,double **MAXVEL,int **SOLID,int **NEARSOLIDX,int **NEARSOLIDY,double ***M,double **KS,double **THETA0,double **THETAS,double **HF,double **VINF,double **ERODIBILITYMASK,double **UC,double **UC2,double rhow,double rhos,int nx,int ny,double hfilm,double afave,double afrave,int nn)
{   int i,j,k,iup,idown,jup,jdown,iup2,jup2,slopex,slopy,ki,kj,ilow,ihigh,jlow,jhigh,temp,sclass;
    int **SOLID2;
    double tempmax;
    
    FILE *fr0,*fr1,*fr2,*fr3,*fr4,*fr5,*fr6,*fr7,*fr8,*fr9,*fr12;
    
    fr0=fopen("./topoin","r");
    fr1=fopen("./solidin","r");
    fr2=fopen("./roughnessin","r");
    fr3=fopen("./cmaskin","r");
    fr4=fopen("./ksin","r");
    fr5=fopen("./theta0in","r");
    fr6=fopen("./thetasin","r");
    fr7=fopen("./hfin","r");
    fr8=fopen("./vinfin","r");
    fr9=fopen("./erodibilitymaskin","r");
    fr12=fopen("./omegacin","r");

    
    sclass=nn-3;
    SOLID2=imatrix(1,nx,1,ny);
        
    for (i=1;i<=nx;i++)
        for (j=1;j<=ny;j++)
        {
            
            fscanf(fr0,"%lf",&TOPO[i][j]);          // TOPOGRAPHY
            VEL[1][i][j]=0;                         // X-VELOCITY
            VEL[2][i][j]=0;                         // Y-VELOCITY
            U[1][i][j]=1e-4;                        // H
            U[2][i][j]=0;                           // UH
            U[3][i][j]=0;                           // VH
            U[4][i][j]=0;                           // CH
            fscanf(fr1,"%d",&SOLID[i][j]);          // MASK FOR SOLID INTERIOR B.C.
            fscanf(fr2,"%lf",&ROUGHNESS[i][j]);     // ROUGHNESS LENGTH SCALE FOR FRICTION IF USING VPE FRICTION (or MANNING's N otherwised)
            fscanf(fr3,"%d",&CMASK[i][j]);          // MASK FOR CHANNEL LOCATIONS
            fscanf(fr4,"%lf",&KS[i][j]);            // SATURATED HYDRAULIC CONDUCTIVITY
            fscanf(fr5,"%lf",&THETA0[i][j]);
            fscanf(fr6,"%lf",&THETAS[i][j]);
            fscanf(fr7,"%lf",&HF[i][j]);
            fscanf(fr8,"%lf",&VINF[i][j]);
            fscanf(fr9,"%lf",&ERODIBILITYMASK[i][j]);
            fscanf(fr12,"%lf",&UC[i][j]);
                        
            TOPOINIT[i][j]=TOPO[i][j];
            DEPTHINIT[i][j]=U[1][i][j];
            NEARSOLIDX[i][j]=SOLID[i][j];
            NEARSOLIDY[i][j]=SOLID[i][j];
            MANNING[i][j]=ROUGHNESS[i][j];
            UC2[i][j]=UC[i][j];
            
            if (U[1][i][j]>hfilm)
            {
                VEL[1][i][j]=U[2][i][j]/U[1][i][j];
                VEL[2][i][j]=U[3][i][j]/U[1][i][j];
                SED[1][i][j]=U[4][i][j]/U[1][i][j];
            }
            if (U[1][i][j]<=hfilm)
            {
                U[2][i][j]=0;
                U[3][i][j]=0;
                U[4][i][j]=0;
                VEL[1][i][j]=0;
                VEL[2][i][j]=0;
                SED[1][i][j]=0;
            }
            
            for (k=4;k<=nn;k++)
            {
                U[k][i][j]=U[4][i][j];
                SED[k-3][i][j]=SED[1][i][j];
            }
            
            RHO[i][j]=rhow*(1-SED[1][i][j])+rhos*SED[1][i][j];
            BETA[1][i][j]=afave;
            BETA[2][i][j]=afrave;
            
            UOLD[1][i][j]=U[1][i][j];
            UOLD[2][i][j]=U[2][i][j];
            UOLD[3][i][j]=U[3][i][j];
            UOLD[4][i][j]=U[4][i][j];
            UOLD2[1][i][j]=U[1][i][j];
            UOLD2[2][i][j]=U[2][i][j];
            UOLD2[3][i][j]=U[3][i][j];
            UOLD2[4][i][j]=U[4][i][j];
            VELOLD[1][i][j]=VEL[1][i][j];
            VELOLD[2][i][j]=VEL[2][i][j];
            MAXVEL[i][j]=0;
            
            //VINF[i][j]=1e-4;                  // Cumulative Infiltrated Depth [m]
        }

    for (i=1;i<=nx;i++)
        for (j=1;j<=ny;j++)
        {
            iup=fmin(i+1,nx);
            iup2=fmin(i+2,nx);
            idown=fmax(i-1,1);
    
            jup=fmin(j+1,ny);
            jup2=fmin(j+2,ny);
            jdown=fmax(j-1,1);
                        
            if (SOLID[idown][j]+SOLID[i][j]+SOLID[iup][j]+SOLID[iup2][j]>0) NEARSOLIDX[i][j]=1;
            if (SOLID[i][jdown]+SOLID[i][j]+SOLID[i][jup]+SOLID[i][jup2]>0) NEARSOLIDY[i][j]=1;
            
            if (SOLID[i][j]!=0)
            {
                if (SOLID[iup][j]==0) TOPO[i][j]=TOPO[iup][j];
                if (SOLID[idown][j]==0) TOPO[i][j]=TOPO[idown][j];
                if (SOLID[iup][jup]==0) TOPO[i][j]=TOPO[iup][jup];
                if (SOLID[iup][jdown]==0) TOPO[i][j]=TOPO[iup][jdown];
                if (SOLID[idown][jup]==0) TOPO[i][j]=TOPO[idown][jup];
                if (SOLID[idown][jdown]==0) TOPO[i][j]=TOPO[idown][jdown];
                if (SOLID[i][jup]==0) TOPO[i][j]=TOPO[i][jup];
                if (SOLID[i][jdown]==0) TOPO[i][j]=TOPO[i][jdown];
            }
            
            if (SOLID[iup][j]!=0 && SOLID[i][j]==0)
                TOPO[iup][j]=TOPO[i][j];
            if (SOLID[idown][j]!=0 && SOLID[i][j]==0)
                TOPO[idown][j]=TOPO[i][j];
            if (SOLID[i][jup]!=0 && SOLID[i][j]==0)
                TOPO[i][jup]=TOPO[i][j];
            if (SOLID[i][jdown]!=0 && SOLID[i][j]==0)
                TOPO[i][jdown]=TOPO[i][j];
        }
    
    fclose(fr0);
    fclose(fr1);
    fclose(fr2);
    fclose(fr3);
    fclose(fr4);
    fclose(fr5);
    fclose(fr6);
    fclose(fr7);
    fclose(fr8);
    fclose(fr9);
    fclose(fr12);
}

double hllc(double hl,double hr,double uhl,double uhr,double vhl,double vhr,double TFLUX[3],double g,double hfilm)
{   double hstar,ustar,Sr,Sl,Sm,S,ul,ur,vl,vr;
    
    S=0;
    
    // Find Velocity and Sediment Concentration From Conserved Variables
    if (hl>hfilm)
    {
        ul=uhl/hl;
        vl=vhl/hl;
    }
    else
    {
        ul=0;
        vl=0;
    }
    if (hr>hfilm)
    {
        ur=uhr/hr;
        vr=vhr/hr;
    }
    else
    {
        ur=0;
        vr=0;
    }
    
    hstar=pow(ul+2*pow(g*hl,.5)-ur-2*pow(g*hr,.5),2)/(16*g);
    ustar=.5*(ul+ur)+pow(g*hl,.5)-pow(g*hr,.5);
    
    // Wave Speeds
    Sl=ul-pow(g*hl,.5);
    if (ustar-pow(g*hstar,.5)<Sl) Sl=ustar-pow(g*hstar,.5);
    Sm=ustar;
    Sr=ur+pow(g*hr,.5);
    if (ustar+pow(g*hstar,.5)>Sr) Sr=ustar+pow(g*hstar,.5);
    
    // Special Cases For Dry Bed On Left/Right Of Cell
    if (hl<=hfilm)
    {
        Sl=ur-2*pow(g*hr,.5);
        Sm=Sl;
        Sr=ur+pow(g*hr,.5);
    }
    if (hr<=hfilm)
    {
        Sl=ul-pow(g*hl,.5);
        Sr=ul+2*pow(g*hl,.5);
        Sm=Sr;
    }
    
    // Flux
    if (Sl>=0)
    {
        TFLUX[0]=hl*ul;
        TFLUX[1]=hl*ul*ul+.5*g*hl*hl;
    }
    if (Sl<=0 && Sr>=0)
    {
        TFLUX[0]=(Sr*(hl*ul)-Sl*(hr*ur)+Sl*Sr*(hr-hl))/(Sr-Sl);
        TFLUX[1]=(Sr*(hl*ul*ul+.5*g*hl*hl)-Sl*(hr*ur*ur+.5*g*hr*hr)+Sl*Sr*(hr*ur-hl*ul))/(Sr-Sl);
    }
    if (Sr<=0)
    {
        TFLUX[0]=hr*ur;
        TFLUX[1]=hr*ur*ur+.5*g*hr*hr;
    }
    if (Sm>=0)
    {
        TFLUX[2]=TFLUX[0]*vl;
    }
    else
    {
        TFLUX[2]=TFLUX[0]*vr;
    }
    
    // Max Wave Speed For Stability Criteria
    if (fabs(Sl)>S) S=fabs(Sl);
    if (fabs(Sr)>S) S=fabs(Sr);
    if (fabs(Sm)>S) S=fabs(Sm);
    
    return S;
}

// Stream-power-dependent erosion of cohesive sediment
double streampowererosion(double uh,double vh,double h,double beta,double rhow,double g,double sedconc,double sf,double minmorph,double uc)
{
    double E,omega,q;
    
    q=pow(uh*uh+vh*vh,0.5);                 // Unit discharge
    omega=rhow*g*sf*q;                    // Stream power
    
    E=0;
    if (omega>uc && sedconc<0.6 && h>minmorph)
    {
        E=beta*(omega-uc);
        //if (omega>1e-5) printf("%lf\t %lf\t %lf\t %lf\n",E,omega,beta,sf);
    }
    
    return E;
}

// Deposition of non-cohesive sediment
double settlingvelocity(double h,double d,double c,double nu,double rhos,double rhow,double g,double minmorph)
{
    double omega,s;
    
    s=rhos/rhow-1;
    
    omega=pow(pow(13.95*nu/d,2)+1.09*s*g*d,0.5)-13.95*nu/d;           // Settling velocity of particle in tranquil water (Cao 2006)
    //omega=omega*pow(1-fmin(c,0.99),5);
    
    //if (h<minmorph) {omega=0;}
    
    return omega;
}

double minmod(double a,double b)
/* minmod function, input two real numbers, output one real number */
{
    double xx;
    if (a*b>0)
    {
        if (fabs(a)<=fabs(b))
        {
            xx=a;
        }
        else
        {
            xx=b;
        }
    }
    else
    {
        xx=0;
    }
    
    return xx;
}

void MUSCLextrap(double valminus, double val, double valplus, double valplus2, double ULR[8][2],int kk,int order)
{
    // 1st Order Godunov
    ULR[kk][0]=val;
    ULR[kk][1]=valplus;
}

int main()
{   FILE *fpi0,*fp0,*fp1,*fp2,*fp3,*fp4,*fp5,*fp6,*fp7,*fp8,*fp9,*fp10,*fp11,*fp12,*fp13,*fp14,*fp15,*fp16,*fp17,*fp18,*fp19,*fp20,*fp21,*fp22,*fp23,*fp24,*fp25,*fp26,*fp27,*fp28,*fp29,*fp30,*fp31,*fp32;
    FILE *fpm0,*fpm1,*fpm2,*fpm3,*fpm4,*fpm5,*fpm6,*fpm7,*fr10,*fr11,*fr12,*fr13,*fr14,*fr15;
    int i,j,k,nx,ny,cflstat,order,ptcnt,ptcntq,ptcntc,steady,*BNDRY,rnum,rind,nn,stagecnt,outlet,outlettopo,sclass,outlety,outlety2,idum,outlet2;
    double **TOPO,**ITOPO,**SX,***BETA,**BETAX,**INITDEPTH,**DEPTHINIT,*INPUT,**TOPOINIT,***D,****E,**RHO,***UOLD2,**TOTSED,**TOTSEDOLD2,**SF,**ROUGHNESS,**FAIL;
    double ***U,***VEL,***SED,***FLUXX,***FLUXY,**SY,**BETAY,**GAMMAX,**GAMMAY,***UOLD,***VELOLD,**MAXVEL,**MAXDEPTH,***M,**TOTM,**DEP,*DS,*PTEMP,***P,***VF,**TOPOOLD,*RAIN,*RAINDIAM,*RAINVEL,**INFL,**TRACER;
    double Smax,Smaxtemp,thetac,cstable,oneoverdx,oneoverdx2,maxstep,L,Dstar,qchange,sedchange,q,rint,**FAILCNT,pgfx,pgfy,cthreshold,cohesion;
    double **ZF,**THETA0,**THETAS,**HF,**VINF,**KS,**SAVEFLOW,**AMAP,**ADMASK,**ERODIBILITYMASK,**JMAP,**ADMAP,**GX,**GY,**GZ;
    double **STAGE,**DFMASK,ucstar,ucstar2,**UC,**UC2,**SLOPE,**SLOPEX,**SLOPEY,d,**MANNING,dt,dx,t,tend,printinterval,printint,phi,velmax,epsilon;
    double hchange,beta,a,ad,af,afr,mtstar,**H,maxtopochange,lambda,maxsoilthickness;
    double timestep,K,rhos,rhow,rho0,R,g,courant,flowvel,hfilm,minmorph,minstep,nu,chtemp,redetach,vs,h0,b,mtstar0,outletelev,outletdepth;
    double stageint,stageinterval,frictionangle,F,Jentrain,s,depthdependentexponent,hcfriction,magdischarge,pi,gi,Si,Ki,Cv,CS,Di,R1;
    int **SOLID,**NEARSOLIDX,**NEARSOLIDY,**CMASK,*IUP,*IDOWN,*JUP,*JDOWN,*XPOINT,*YPOINT;
    
    idum=-2999;
    
    fp0=fopen("./cflstatus","w");
    
    fpm0=fopen("./topomovie","w");
    fpm1=fopen("./depthmovie","w");
    fpm2=fopen("./velocitymovie","w");
    fpm3=fopen("./cmovie","w");
    fpm4=fopen("./Emovie","w");
    fpm5=fopen("./Dmovie","w");
    fpm6=fopen("./Mmovie","w");
    
    fr10=fopen("./input","r");
    fr11=fopen("./rain","r");
    fr12=fopen("./particlesizein","r");
    fr13=fopen("./particlepercentin","r");
    fr14=fopen("./raindiam","r");
    fr15=fopen("./rainvel","r");
    
    // Input Parameters
    INPUT=vector(1,21);
    for (i=1;i<=21;i++)
        fscanf(fr10,"%lf",&INPUT[i]);
    nx=INPUT[1];
    ny=INPUT[2];                            // [m]
    dx=INPUT[3];                            // [m]
    tend=INPUT[5];                          // Total Simulation Time [s]
    
    // Parameters For HR Sediment Transport
    a=INPUT[6];                             // Detachability of Original Soil (Used to determine rainsplash detachment rate)
    ad=INPUT[7];                            // Detachability of Deposited Sediment (Used to determine rainsplash detachment rate for deposited sediment)
    h0=INPUT[8];                            // Critical depth used in equation to determine rainsplash detachment as funciton of water depth
    mtstar0=INPUT[9];                       // Mass of Deposited Sediment Needed to Shield Original Soil From Erosion (From Heng el al. (2011))
    F=INPUT[4];                             // Effective Fraction of Stream Power In Entrainment
    Jentrain=INPUT[13];                     // Energy Expended in Entraining a Unit Mass of Cohesive Sediment (see Kim et al, 2013) [m^2]/[s^2]
    hcfriction=INPUT[14];                   // Critical Overland Flow Depth In Depth-Dependent Manning Eqn
    depthdependentexponent=INPUT[15];       // Exponent In Depth-Dependent Manning Eqn
    b=1;                                    // Exponent in Rainsplash Model (determines rainsplash detachment as funciton of water depth)
    af=1;
    afr=1;
    
    // Parameters For Rutter Interception Model
    pi=INPUT[16];                           // Free Throughfall Coefficient
    Si=INPUT[17];                           // Canopy capacity (m)
    Ki=INPUT[18];                           // Canopy Drainage Rate Coefficient (m/s)
    gi=INPUT[19];                           // Canopy Drainage Rate Exponent (1/m)
    Cv=INPUT[20];                           // Percent Vegetation Cover [-]
    CS=0;
    Di=0;
    
    // Rainfall Input Data
    rnum=INPUT[10];                         // Total number of entries in input rainfall time series
    rint=INPUT[11];                         // Time [s] between rainfall entries in input rainfall time series
    RAIN=vector(1,rnum);
    RAINDIAM=vector(1,rnum);
    RAINVEL=vector(1,rnum);
    for (i=1;i<=rnum;i++)
    {
        fscanf(fr11,"%lf",&RAIN[i]);
        fscanf(fr14,"%lf",&RAINDIAM[i]);
        fscanf(fr15,"%lf",&RAINVEL[i]);
    }
    
    // Particle Size Input Data
    sclass=INPUT[12];
    nn=3+sclass;
    DS=vector(1,sclass);                                // Particle Diameter for Each Sediment Class
    PTEMP=vector(1,sclass);                             // Fraction of Each Sediment Class in Original Soil
    P=f3tensor(1,sclass,1,nx,1,ny);                     // Fraction of Each Sediment Class in Original Soil
    for (k=1;k<=sclass;k++)
    {
        fscanf(fr12,"%lf",&DS[k]);                     // Particle Diameter for Each Sediment Class
        fscanf(fr13,"%lf",&PTEMP[k]);                      // Fraction of Sediment Class k in Original Soil
    }
    
    // Initial Soil Particle Size Distribution In Each Cell
    for (i=1;i<=nx;i++)
        for (j=1;j<=ny;j++)
        {
            for (k=1;k<=sclass;k++)
            {
                P[k][i][j]=PTEMP[k];
            }
        }
    
    printf("Grid Size (Spacing): %d x %d (%lf)\n",nx,ny,dx);
    printf("Number of conserved variables: %d\n",nn);
    printf("Number of particle size classes: %d\n",sclass);
    
    printf("\n");
    printf("H-R Parameters: \n");
    printf("a0: %lf\n",a);
    printf("ad0: %lf\n",ad);
    printf("F: %lf\n",F);
    printf("J: %lf\n",Jentrain);
    printf("h0: %lf\n",h0);
    printf("mtstar0: %lf\n",mtstar0);
    printf("b: %lf\n",b);
    
    printf("\n");
    printf("Interception Parameters: \n");
    printf("Throughfall coefficient: %lf\n",pi);
    printf("Canopy storage: %lf\n",Si);
    printf("Canopy drainage rate coefficient: %lf\n",Ki);
    printf("Canopy drainage rate exponent: %lf\n",gi);
    printf("Fraction vegetation cover: %lf\n",Cv);
    printf("\n");
    
    oneoverdx=1/dx;
    oneoverdx2=1/(dx*dx);
    
    U=f3tensor(1,nn,-1,nx+2,-1,ny+2);                   // Matrix of Conserved Variables
    UOLD=f3tensor(1,nn,-1,nx+2,-1,ny+2);                // Matrix of Conserved Variables
    TOTSEDOLD2=matrix(1,nx,1,ny);                       // Stores Vol Sed Conc Used to Determine Steady State
    UOLD2=f3tensor(1,nn,-1,nx+2,-1,ny+2);               // Matrix of Conserved Variables
    TOPO=matrix(-1,nx+2,-1,ny+2);                       // Topographic Elevation [m]
    ITOPO=matrix(-1,nx+2,-1,ny+2);
    TOPOOLD=matrix(-1,nx+2,-1,ny+2);                    // Topographic Elevation At Previous Time Step [m]
    VEL=f3tensor(1,2,-1,nx+2,-1,ny+2);                  // Velocity [m]/[s]
    VELOLD=f3tensor(1,2,-1,nx+2,-1,ny+2);               // Velocity [m]/[s]
    MAXVEL=matrix(1,nx,1,ny);                           // Velocity [m]/[s]
    MAXDEPTH=matrix(1,nx,1,ny);                         // Maximum Flow Depth [m]
    GX=matrix(-1,nx+2,-1,ny+2);                         // Component of Gravity in X-Direction
    GY=matrix(-1,nx+2,-1,ny+2);                         // Component of Gravity in Y-Direction
    GZ=matrix(-1,nx+2,-1,ny+2);                         // Component of Gravity in Z-Direction
    SED=f3tensor(1,sclass,-1,nx+2,-1,ny+2);             // Sediment Concentration [kg]/[m^3]
    M=f3tensor(1,sclass,-1,nx+2,-1,ny+2);               // Mass in Sediment Class i In Deposited Layer [kg]
    E=f4tensor(1,4,1,sclass,1,nx,1,ny);                 // Values for Erosion Source Terms
    D=f3tensor(1,sclass,1,nx,1,ny);                     // Values for Deposition Source Terms
    VF=f3tensor(1,sclass,1,nx,1,ny);                    // Settling Velocity for Particle Class k
    TRACER=matrix(-1,nx+2,-1,ny+2);                     // Passize Tracer Used to Track Debris Flows
    RHO=matrix(1,nx,1,ny);                              // Fluid Density
    SLOPE=matrix(1,nx,1,ny);                            // Topographic Slope
    SLOPEX=matrix(1,nx,1,ny);                           // Topographic Slope in X-Direction
    SLOPEY=matrix(1,nx,1,ny);                           // Topographic Slope in Y-Direction
    UC=matrix(1,nx,1,ny);                               // Critical Stream Power For Original Soil
    UC2=matrix(1,nx,1,ny);                              // Critical Stream Power For Detached Soil
    INFL=matrix(1,nx,1,ny);                             // Infiltration
    KS=matrix(1,nx,1,ny);                               // Saturated Hydraulic Conductivy
    THETA0=matrix(1,nx,1,ny);                           // Initial Soil Moisture
    THETAS=matrix(1,nx,1,ny);                           // Saturated Zone Soil Moisture
    HF=matrix(1,nx,1,ny);                               // Wetting Front Capillary Pressure Head
    ZF=matrix(1,nx,1,ny);                               // Depth of Wetting Front [m]
    VINF=matrix(1,nx,1,ny);                             // Cumulative Infiltrated Depth [m]
    BETAX=matrix(1,nx,1,ny);                            // X-direction Momentum Forcing
    SX=matrix(1,nx,1,ny);                               // Topographic Forcing
    BETAY=matrix(1,nx,1,ny);                            // Y-direction Momentum Forcing
    SY=matrix(1,nx,1,ny);                               // Topographic Forcing
    SF=matrix(1,nx,1,ny);                               // Friction Slope
    TOPOINIT=matrix(1,nx,1,ny);                         // Initial Topography [m]
    BETA=f3tensor(1,2,1,nx,1,ny);                       // Substrate Erodibility and Erodibility of Deposited Sediment
    DEPTHINIT=matrix(1,nx,1,ny);                        // Initial Fluid Depth [m]
    TOTSED=matrix(1,nx,1,ny);                           // Total Volumetric Sediment Concentration
    TOTM=matrix(-1,nx+2,-1,ny+2);                       // Total Mass of Sediment Per Unit Area in Deposited Layer
    DEP=matrix(1,nx,1,ny);                              // Total Thickness of Sediment in Deposited Layer [m]
    H=matrix(1,nx,1,ny);                                // Fractional Shielding of Original Soil
    ROUGHNESS=matrix(1,nx,1,ny);                        // Roughness Length Scale Used to Compute Friciton
    CMASK=imatrix(1,nx,1,ny);                           // CMASK[i][j]==1 Indicates a Channel Within the Pixel at (i,j)
    ERODIBILITYMASK=matrix(1,nx,1,ny);                  // Can be Used to Create Spatial Variability in HR model parameters
    AMAP=matrix(1,nx,1,ny);                             // Spatially Variable `a_0' coefficient in HR model
    ADMAP=matrix(1,nx,1,ny);                            // Spatially Variable `ad_0' coefficient in HR model
    JMAP=matrix(1,nx,1,ny);                             // Spatially Variable `J' coefficient in HR model
    DFMASK=matrix(1,nx,1,ny);                           // Mask to Quantify Importance of Debris Flow Resistance Terms
    GAMMAX=matrix(1,nx,1,ny);                           // X-direction Momentum Forcing
    GAMMAY=matrix(1,nx,1,ny);                           // Y-direction Momentum Forcing
    MANNING=matrix(1,nx,1,ny);                          // Manning Roughness Coefficient
    FAIL=matrix(1,nx,1,ny);                             // 0==Sediment is Stable, Nonzero==Thickness of Unstable Sediment
    FAILCNT=matrix(1,nx,1,ny);
    SOLID=imatrix(1,nx,1,ny);                           // SOLID[i][j]==1 Defines an Interior B.C.
    NEARSOLIDX=imatrix(1,nx,1,ny);                      // NEARSOLIDX[i][j]==1 Implies Interior Boundary At Neighboring Cell in X-Direction
    NEARSOLIDY=imatrix(1,nx,1,ny);                      // NEARSOLIDY[i][j]==1 Implies Interior Boundary At Neighboring Cell in Y-Direction
    STAGE=matrix(1,27,1,30000);                         // Stores Flow Information At Basin Outlet
    SAVEFLOW=matrix(1,18,1,30000);                      // Stores Flow Information At Interior Points Defined by (XPOINT,YPOINT)
    XPOINT=ivector(1,6);
    YPOINT=ivector(1,6);
    double ULR[nn][2];                                  // Stores Left and Right Hand Values of Variables For Riemann Problem
    double TFLUX[3];                                    // Temp Stores Flux
    IUP=ivector(1,nx);
    IDOWN=ivector(1,nx);
    JUP=ivector(1,ny);
    JDOWN=ivector(1,ny);
    
    FLUXX=f3tensor(1,nn+1,1,nx+1,1,ny);                   // X-direction Flux
    FLUXY=f3tensor(1,nn+1,1,nx,1,ny+1);                   // Y-direction Flux
    
    // Boudnary Conditions
    BNDRY=ivector(1,4);                         // 0==Transmissive, 1==Solid
    BNDRY[1]=0;                                 // B.C. at i==1
    BNDRY[2]=0;                                 // B.C. at j==1
    BNDRY[3]=0;                                 // B.C. at i==nx
    BNDRY[4]=0;                                 // B.C. at j==ny
    
    // Time Parameters
    t=0;                                        // [s]
    dt=0.0001;                                  // Initial Time Step for Fluid Equations [s]
    printinterval=fmax(1,tend/5);
    printint=fmax(1,tend/5);
    stagecnt=0;
    stageint=5;
    stageinterval=0;
    cflstat=0;                                  // clfstat==1 indicates cfl error
    cstable=0.05;                               // Adaptive time step will adjust to make courant=cstable
    maxstep=0.25;                               // Maximum Time Step Regardless of Courant Number
    minstep=1e-5;                               // Simulation will abort if time step is too small
    order=1;                                    // 1=Godunov
    
    // Additional Model Parameters
    g=9.81;                                 // [m]/[s^2]
    hfilm=1e-5;                             // Water Depth < hfilm is Treated as Dry
    phi=0.4;                                // Bed Sediment Porosity
    rhos=2600;                              // Density of sediment [kg]/[m^3]
    rhow=1000;                              // Density of water [kg]/[m^3]
    s=(rhos-rhow)/rhow;
    rho0=rhow*phi+rhos*(1-phi);             // Density of saturated bed [kg]/[m^3]
    nu=1.2e-6;                              // Kinematic viscosity of fluid [Pa s]
    frictionangle=0.6745;                   // Static Angle of Friction
    cohesion=1e10;                          // Cohesion of Sediment in Deposited Layer
    lambda=0.65;                            // Ratio of Basal Pore Pressure to Total Basal Normal Stress in Debris Flow
    cthreshold=0.4;                         // Threshold volumetric sediment concentration for debris flow
    maxsoilthickness=0.5;                   // Maximum Soil Thickness (Prevents Further Erosion in a Pixel if this Limit is Reached in that Pixel)
    
    // Set Initial Conditions
    initialcondition(TOPO,VEL,U,VELOLD,UOLD,UOLD2,DEPTHINIT,TOPOINIT,SED,RHO,MANNING,ROUGHNESS,CMASK,BETA,MAXVEL,SOLID,NEARSOLIDX,NEARSOLIDY,M,KS,THETA0,THETAS,HF,VINF,ERODIBILITYMASK,UC,UC2,rhow,rhos,nx,ny,hfilm,af,afr,nn);
    
    // Print Initial Topography
    fpi0=fopen("./initialtopo","w");
    for (j=1;j<=ny;j++)
        for (i=1;i<=nx;i++)
        {
            if (i!=nx) fprintf(fpi0,"%10.10lf\t",TOPO[i][j]);
            if (i==nx) fprintf(fpi0,"%10.10lf\n",TOPO[i][j]);
            if (i==1 && j==1) printf("%lf\n",TOPO[i][j]);
            MAXDEPTH[i][j]=0;
            ITOPO[i][j]=TOPO[i][j];
            FAILCNT[i][j]=0;
            AMAP[i][j]=a;
            //AMAP[i][j]=a*ERODIBILITYMASK[i][j];
            ADMAP[i][j]=ad;
            //ADMAP[i][j]=ad*ERODIBILITYMASK[i][j];
            //JMAP[i][j]=Jentrain;
            //JMAP[i][j]=0.5*1000*5.5*5.5/AMAP[i][j];         // Energy Expended in Entraining a Unit Mass of Cohesive Sediment (see Kim et al, 2013) [m^2]/[s^2]
        }
    
    for (j=1;j<=ny;j++)
    {
        JUP[j]=j+1;
        JDOWN[j]=j-1;
    }
    JDOWN[1]=1;
    JUP[ny]=ny;
    
    for (i=1;i<=nx;i++)
    {
        IUP[i]=i+1;
        IDOWN[i]=i-1;
    }
    IUP[nx]=nx;
    IDOWN[1]=1;
    
    // Find Basin Outlet
    outlet=1; outlety=1; outletelev=99999;
    for (i=1;i<=nx;i++)
    {
        if (SOLID[i][outlety]==0 && TOPO[i][outlety]<outletelev)
        {
            outletelev=TOPO[i][outlety];
            outlet=i;                                       // Save outlet location
        }
    }
    // Find Basin Outlet 2
    outlet2=1; outlety2=10; outletelev=99999;
    for (i=1;i<=nx;i++)
    {
        if (SOLID[i][outlety2]==0 && TOPO[i][outlety2]<outletelev)
        {
            outletelev=TOPO[i][outlety2];
            outlet2=i;                                       // Save outlet location
        }
    }
    printf("Basin Outlet Pixel: %d\n",outlet);
    printf("Basin Outlet Elevation (m): %lf\n",outletelev);
    
    // Define Interior Points Where Flow Information Will Be Saved in the 'SAVEFLOW' array
    XPOINT[1]=1; YPOINT[1]=1;
    XPOINT[2]=1; YPOINT[2]=2;
    XPOINT[3]=2; YPOINT[3]=2;
    XPOINT[4]=2; YPOINT[4]=1;
    XPOINT[5]=1; YPOINT[5]=3;
    XPOINT[6]=3; YPOINT[6]=1;
    
    // Define GX, GY, GZ
    for (i=1;i<=nx;i++)
        for (j=1;j<=ny;j++)
        {
            // Topographic Forcing in X-Direction
            if (i!=1 && i!=nx)
            {
                SLOPEX[i][j]=0.5*oneoverdx*(TOPO[i+1][j]-TOPO[i-1][j]);
                GX[i][j]=g*SLOPEX[i][j]/pow(SLOPEX[i][j]*SLOPEX[i][j]+SLOPEY[i][j]*SLOPEY[i][j]+1,0.5);
                SX[i][j]=-GX[i][j]*U[1][i][j];
            }
            else if (i==nx)
            {
                SLOPEX[i][j]=oneoverdx*(TOPO[i][j]-TOPO[i-1][j]);
                GX[i][j]=g*SLOPEX[i][j]/pow(SLOPEX[i][j]*SLOPEX[i][j]+SLOPEY[i][j]*SLOPEY[i][j]+1,0.5);
                SX[i][j]=-GX[i][j]*U[1][i][j];
            }
            else
            {
                SLOPEX[i][j]=oneoverdx*(TOPO[i+1][j]-TOPO[i][j]);
                GX[i][j]=g*SLOPEX[i][j]/pow(SLOPEX[i][j]*SLOPEX[i][j]+SLOPEY[i][j]*SLOPEY[i][j]+1,0.5);
                SX[i][j]=-GX[i][j]*U[1][i][j];
            }
            
            // Topographic Forcing in Y-Direction
            if (j!=1 && j!=ny)
            {
                SLOPEY[i][j]=0.5*oneoverdx*(TOPO[i][j+1]-TOPO[i][j-1]);
                GY[i][j]=g*SLOPEY[i][j]/pow(SLOPEX[i][j]*SLOPEX[i][j]+SLOPEY[i][j]*SLOPEY[i][j]+1,0.5);
                SY[i][j]=-GY[i][j]*U[1][i][j];
            }
            else if (j==ny)
            {
                SLOPEY[i][j]=oneoverdx*(TOPO[i][j]-TOPO[i][j-1]);
                GY[i][j]=g*SLOPEY[i][j]/pow(SLOPEX[i][j]*SLOPEX[i][j]+SLOPEY[i][j]*SLOPEY[i][j]+1,0.5);
                SY[i][j]=-GY[i][j]*U[1][i][j];
            }
            else
            {
                SLOPEY[i][j]=oneoverdx*(TOPO[i][j+1]-TOPO[i][j]);
                GY[i][j]=g*SLOPEY[i][j]/pow(SLOPEX[i][j]*SLOPEX[i][j]+SLOPEY[i][j]*SLOPEY[i][j]+1,0.5);
                SY[i][j]=-GY[i][j]*U[1][i][j];
            }
            
            GZ[i][j]=g/pow(SLOPEX[i][j]*SLOPEX[i][j]+SLOPEY[i][j]*SLOPEY[i][j]+1,0.5);
            SLOPE[i][j]=pow(SLOPEX[i][j]*SLOPEX[i][j]+SLOPEY[i][j]*SLOPEY[i][j],0.5);
        }
    
    for (j=-1;j<=ny+2;j++)
    {
        GZ[-1][j]=GZ[1][j];
        GZ[0][j]=GZ[1][j];
        GZ[nx+1][j]=GZ[nx][j];
        GZ[nx+2][j]=GZ[nx][j];
        
        GX[-1][j]=GX[1][j];
        GX[0][j]=GX[1][j];
        GX[nx+1][j]=GX[nx][j];
        GX[nx+2][j]=GX[nx][j];
        
        GY[-1][j]=GY[1][j];
        GY[0][j]=GY[1][j];
        GY[nx+1][j]=GY[nx][j];
        GY[nx+2][j]=GY[nx][j];
    }
    
    for (i=-1;i<=nx+2;i++)
    {
        GZ[i][0]=GZ[i][1];
        GZ[i][-1]=GZ[i][1];
        GZ[i][ny+1]=GZ[i][ny];
        GZ[i][ny+2]=GZ[i][ny];
        
        GX[i][0]=GX[i][1];
        GX[i][-1]=GX[i][1];
        GX[i][ny+1]=GX[i][ny];
        GX[i][ny+2]=GX[i][ny];
        
        GY[i][0]=GY[i][1];
        GY[i][-1]=GY[i][1];
        GY[i][ny+1]=GY[i][ny];
        GY[i][ny+2]=GY[i][ny];
    }
    
    omp_set_num_threads(28);
    
    time_t t0 = time(NULL);
    while (t<tend)
    {
        // Solve SWE For Steady State
        Smax=0; Smaxtemp=0; velmax=0;
        
        // Determine Rainfall Rate From Input Data
        rind=fmin(ceil(t/rint+0.01),rnum);
        R1=RAIN[rind];
        h0=0.33*RAINDIAM[rind];
        for (j=1;j<=ny;j++)
            {
                for (i=1;i<=nx;i++)
                {
                    if (t > 9180)
                    {
                        ROUGHNESS[i][j]=0.07;
                    }
                    JMAP[i][j]=0.5*1000*RAINVEL[rind]*RAINVEL[rind]/AMAP[i][j];
                }
            }

        // Rutter Interception Model
        R=(1-Cv)*R1+Cv*(pi*R1+Di);                      // Effective Rainfall Rate [m/s]
        
        // Update Boundary Conditions
        if (BNDRY[1]==0)
        {
            // Transmissive B.C. at i=1
            for (j=-1;j<=ny+2;j++)
            {
                VEL[1][0][j]=VEL[1][1][j];
                VEL[1][-1][j]=VEL[1][1][j];
                VEL[2][0][j]=VEL[2][1][j];
                VEL[2][-1][j]=VEL[2][1][j];
                TRACER[0][j]=TRACER[1][j];
                TRACER[-1][j]=TRACER[1][j];
                for (k=1;k<=nn;k++)
                {
                    U[k][0][j]=U[k][1][j];
                    U[k][-1][j]=U[k][1][j];
                }
                for (k=1;k<=sclass;k++)
                {
                    SED[k][0][j]=SED[k][1][j];
                    SED[k][-1][j]=SED[k][1][j];
                }
            }
        }
        else if (BNDRY[1]==1)
        {
            // Solid B.C. at i=1
            for (j=-1;j<=ny+2;j++)
            {
                VEL[1][0][j]=-VEL[1][1][j];
                VEL[1][-1][j]=-VEL[1][1][j];
                VEL[2][0][j]=VEL[2][1][j];
                VEL[2][-1][j]=VEL[2][1][j];
                TRACER[0][j]=TRACER[1][j];
                TRACER[-1][j]=TRACER[1][j];
                for (k=1;k<=nn;k++)
                {
                    U[k][0][j]=U[k][1][j];
                    U[k][-1][j]=U[k][1][j];
                }
                U[2][0][j]=-U[2][1][j];
                U[2][-1][j]=-U[2][1][j];
                for (k=1;k<=sclass;k++)
                {
                    SED[k][0][j]=SED[k][1][j];
                    SED[k][-1][j]=SED[k][1][j];
                }
            }
        }
        if (BNDRY[2]==0)
        {
            // Transmissive B.C. at j=1
            for (i=-1;i<=nx+2;i++)
            {
                VEL[1][i][0]=VEL[1][i][1];
                VEL[1][i][-1]=VEL[1][i][1];
                VEL[2][i][0]=VEL[2][i][1];
                VEL[2][i][-1]=VEL[2][i][1];
                TRACER[i][0]=TRACER[i][1];
                TRACER[i][-1]=TRACER[i][1];
                for (k=1;k<=nn;k++)
                {
                    U[k][i][0]=U[k][i][1];
                    U[k][i][-1]=U[k][i][1];
                }
                for (k=1;k<=sclass;k++)
                {
                    SED[k][i][0]=SED[k][i][1];
                    SED[k][i][-1]=SED[k][i][1];
                    //M[k][i][1]=0;
                    //TOTM[i][1]=0;
                }
            }
        }
        if (BNDRY[2]==1)
        {
            // Solid B.C. at j=1
            for (i=-1;i<=nx+2;i++)
            {
                VEL[1][i][0]=VEL[1][i][1];
                VEL[1][i][-1]=VEL[1][i][1];
                VEL[2][i][0]=-VEL[2][i][1];
                VEL[2][i][-1]=-VEL[2][i][1];
                TRACER[i][0]=TRACER[i][1];
                TRACER[i][-1]=TRACER[i][1];
                for (k=1;k<=nn;k++)
                {
                    U[k][i][0]=U[k][i][1];
                    U[k][i][-1]=U[k][i][1];
                }
                U[3][i][0]=-U[3][i][1];
                U[3][i][-1]=-U[3][i][1];
                for (k=1;k<=sclass;k++)
                {
                    SED[k][i][0]=SED[k][i][1];
                    SED[k][i][-1]=SED[k][i][1];
                }
            }
        }
        if (BNDRY[3]==0)
        {
            // Transmissive B.C. at i=nx
            for (j=-1;j<=ny+2;j++)
            {
                VEL[1][nx+1][j]=VEL[1][nx][j];
                VEL[1][nx+2][j]=VEL[1][nx][j];
                VEL[2][nx+1][j]=VEL[2][nx][j];
                VEL[2][nx+2][j]=VEL[2][nx][j];
                TRACER[nx+1][j]=TRACER[nx][j];
                TRACER[nx+2][j]=TRACER[nx][j];
                for (k=1;k<=nn;k++)
                {
                    U[k][nx+1][j]=U[k][nx][j];
                    U[k][nx+2][j]=U[k][nx][j];
                }
                for (k=1;k<=sclass;k++)
                {
                    SED[k][nx+1][j]=SED[k][nx][j];
                    SED[k][nx+2][j]=SED[k][nx][j];
                }
            }
        }
        else if (BNDRY[3]==1)
        {
            // Solid B.C. at i=nx
            for (j=-1;j<=ny+2;j++)
            {
                VEL[1][nx+1][j]=-VEL[1][nx][j];
                VEL[1][nx+2][j]=-VEL[1][nx][j];
                VEL[2][nx+1][j]=VEL[2][nx][j];
                VEL[2][nx+2][j]=VEL[2][nx][j];
                TRACER[nx+1][j]=TRACER[nx][j];
                TRACER[nx+2][j]=TRACER[nx][j];
                for (k=1;k<=nn;k++)
                {
                    U[k][nx+1][j]=U[k][nx][j];
                    U[k][nx+2][j]=U[k][nx][j];
                }
                U[2][nx+1][j]=-U[2][nx][j];
                U[2][nx+2][j]=-U[2][nx][j];
                for (k=1;k<=sclass;k++)
                {
                    SED[k][nx+1][j]=SED[k][nx][j];
                    SED[k][nx+2][j]=SED[k][nx][j];
                }
            }
        }
        if (BNDRY[4]==0)
        {
            // Transmissive B.C. at j=ny
            for (i=-1;i<=nx+2;i++)
            {
                VEL[1][i][ny+1]=VEL[1][i][ny];
                VEL[1][i][ny+2]=VEL[1][i][ny];
                VEL[2][i][ny+1]=VEL[2][i][ny];
                VEL[2][i][ny+2]=VEL[2][i][ny];
                TRACER[i][ny+1]=TRACER[i][ny];
                TRACER[i][ny+2]=TRACER[i][ny];
                for (k=1;k<=nn;k++)
                {
                    U[k][i][ny+1]=U[k][i][ny];
                    U[k][i][ny+2]=U[k][i][ny];
                }
                for (k=1;k<=sclass;k++)
                {
                    SED[k][i][ny+1]=SED[k][i][ny];
                    SED[k][i][ny+2]=SED[k][i][ny];
                }
            }
        }
        else if (BNDRY[4]==1)
        {
            // Solid B.C. at j=ny
            for (i=-1;i<=nx+2;i++)
            {
                VEL[1][i][ny+1]=VEL[1][i][ny];
                VEL[1][i][ny+2]=VEL[1][i][ny];
                VEL[2][i][ny+1]=-VEL[2][i][ny];
                VEL[2][i][ny+2]=-VEL[2][i][ny];
                TRACER[i][ny+1]=TRACER[i][ny];
                TRACER[i][ny+2]=TRACER[i][ny];
                for (k=1;k<=nn;k++)
                {
                    U[k][i][ny+1]=U[k][i][ny];
                    U[k][i][ny+2]=U[k][i][ny];
                }
                U[3][i][ny+1]=-U[3][i][ny];
                U[3][i][ny+2]=-U[3][i][ny];
                for (k=1;k<=sclass;k++)
                {
                    SED[k][i][ny+1]=SED[k][i][ny];
                    SED[k][i][ny+2]=SED[k][i][ny];
                }
            }
        }
        
        // Topography B.C. at i=1 and i=nx
        for (j=-1;j<=ny+2;j++)
        {
            TOPO[-1][j]=TOPO[1][j];
            TOPO[0][j]=TOPO[1][j];
            TOPO[nx+1][j]=TOPO[nx][j];
            TOPO[nx+2][j]=TOPO[nx][j];
        }
        // Topography B.C. at j=1 and j=ny
        for (i=-1;i<=nx+2;i++)
        {
            TOPO[i][-1]=TOPO[i][1];
            TOPO[i][0]=TOPO[i][1];
            TOPO[i][ny+1]=TOPO[i][ny];
            TOPO[i][ny+2]=TOPO[i][ny];
        }
        
        for (j=-1;j<=ny+2;j++)
        {
            GZ[-1][j]=GZ[1][j];
            GZ[0][j]=GZ[1][j];
            GZ[nx+1][j]=GZ[nx][j];
            GZ[nx+2][j]=GZ[nx][j];
            
            GX[-1][j]=GX[1][j];
            GX[0][j]=GX[1][j];
            GX[nx+1][j]=GX[nx][j];
            GX[nx+2][j]=GX[nx][j];
            
            GY[-1][j]=GY[1][j];
            GY[0][j]=GY[1][j];
            GY[nx+1][j]=GY[nx][j];
            GY[nx+2][j]=GY[nx][j];
        }
        
        for (i=-1;i<=nx+2;i++)
        {
            GZ[i][0]=GZ[i][1];
            GZ[i][-1]=GZ[i][1];
            GZ[i][ny+1]=GZ[i][ny];
            GZ[i][ny+2]=GZ[i][ny];
            
            GX[i][0]=GX[i][1];
            GX[i][-1]=GX[i][1];
            GX[i][ny+1]=GX[i][ny];
            GX[i][ny+2]=GX[i][ny];
            
            GY[i][0]=GY[i][1];
            GY[i][-1]=GY[i][1];
            GY[i][ny+1]=GY[i][ny];
            GY[i][ny+2]=GY[i][ny];
        }
        
        // Compute Fluxes
        #pragma omp parallel for private(j,k,Smaxtemp,ULR,TFLUX,q,beta,mtstar,pgfx,pgfy)
        for (i=0;i<=nx;i++)
            for (j=0;j<=ny;j++)
            {
                Smaxtemp=0;
                
                // Flux In X-Direction
                // Flux Across Interface at (i+1/2)
                if (j!=0)
                {
                    for (k=1;k<=nn;k++)
                        FLUXX[k][i+1][j]=0;
                    
                    // Conserved Variables at Left/Right Sides of (i+1/2)
                    for (k=1;k<=3;k++)
                    {
                        if (i!=0 && NEARSOLIDX[i][j]!=0)
                        {
                            // Reduce to 1st Order Flux if Near Interior Boundary
                            MUSCLextrap(U[k][i-1][j],U[k][i][j],U[k][i+1][j],U[k][i+2][j],ULR,k-1,1);
                        }
                        else
                        {
                            // Not Near Interior Boundary
                            MUSCLextrap(U[k][i-1][j],U[k][i][j],U[k][i+1][j],U[k][i+2][j],ULR,k-1,order);
                        }
                    }
                    
                    // Sediment Concentration at Left/Right Sides of (i+1/2)
                    for (k=4;k<=nn;k++)
                    {
                        if (i!=0 && NEARSOLIDX[i][j]!=0)
                        {
                            // Reduce to 1st Order Flux if Near Interior Boundary
                            MUSCLextrap(SED[k-3][i-1][j],SED[k-3][i][j],SED[k-3][i+1][j],SED[k-3][i+2][j],ULR,k-1,1);
                        }
                        else
                        {
                            // Not Near Interior Boundary
                            MUSCLextrap(SED[k-3][i-1][j],SED[k-3][i][j],SED[k-3][i+1][j],SED[k-3][i+2][j],ULR,k-1,order);
                        }
                    }
                    
                    // HLLC Flux
                    Smaxtemp=hllc(ULR[0][0],ULR[0][1],ULR[1][0],ULR[1][1],ULR[2][0],ULR[2][1],TFLUX,GZ[i][j],hfilm);
                    
                    // Transparent Interior B.C. On the Left: Reduce to 1st Order Flux
                    if (i!=0 && SOLID[i][j]!=0 && i!=nx && SOLID[i+1][j]==0)
                    {
                        if ( SOLID[i][j]==1)
                        {
                            Smaxtemp=hllc(U[1][i+1][j],U[1][i+1][j],U[2][i+1][j],U[2][i+1][j],U[3][i+1][j],U[3][i+1][j],TFLUX,GZ[i][j],hfilm);
                        }
                        if ( SOLID[i][j]==2)
                        {
                            Smaxtemp=hllc(U[1][i+1][j],U[1][i+1][j],-U[2][i+1][j],U[2][i+1][j],U[3][i+1][j],U[3][i+1][j],TFLUX,GZ[i][j],hfilm);
                        }
                    }
                    // Transparent Interior B.C. On the Right: Reduce to 1st Order Flux
                    else if (i!=0 && SOLID[i][j]==0 && i!=nx && SOLID[i+1][j]!=0)
                    {
                        if (SOLID[i+1][j]==1)
                        {
                            Smaxtemp=hllc(U[1][i][j],U[1][i][j],U[2][i][j],U[2][i][j],U[3][i][j],U[3][i][j],TFLUX,GZ[i][j],hfilm);
                        }
                        if (SOLID[i+1][j]==2)
                        {
                            Smaxtemp=hllc(U[1][i][j],U[1][i][j],U[2][i][j],-U[2][i][j],U[3][i][j],U[3][i][j],TFLUX,GZ[i][j],hfilm);
                        }
                    }
                    
                    // Store Flux
                    FLUXX[1][i+1][j]=TFLUX[0];
                    FLUXX[2][i+1][j]=TFLUX[1];
                    FLUXX[3][i+1][j]=TFLUX[2];
                    for (k=4;k<=nn;k++)
                    {
                        // Upwind Flux for Each Particle Size Class
                        if (TFLUX[0]>0) FLUXX[k][i+1][j]=TFLUX[0]*ULR[k-1][0];
                        else FLUXX[k][i+1][j]=TFLUX[0]*ULR[k-1][1];
                    }
                    if (TFLUX[0]>0) FLUXX[nn+1][i+1][j]=TFLUX[0]*TRACER[i][j];
                    else FLUXX[nn+1][i+1][j]=TFLUX[0]*TRACER[i+1][j];
                    
                    if (Smaxtemp>Smax) Smax=Smaxtemp;
                    
                }
                
                // Flux In Y-Direction
                // Flux Across Interface at (j+1/2)
                if (i!=0)
                {
                    for (k=1;k<=nn;k++)
                        FLUXY[k][i][j+1]=0;
                    
                    // Conserved Variables at Left/Right Sides of (j+1/2)
                    for (k=1;k<=3;k++)
                    {
                        if (j!=0 && NEARSOLIDY[i][j]!=0)
                        {
                            // Reduce to 1st Order Flux if Near Interior Boundary
                            MUSCLextrap(U[k][i][j-1],U[k][i][j],U[k][i][j+1],U[k][i][j+2],ULR,k-1,1);
                        }
                        else
                        {
                            // Not Near Interior Boundary
                            MUSCLextrap(U[k][i][j-1],U[k][i][j],U[k][i][j+1],U[k][i][j+2],ULR,k-1,order);
                        }
                    }
                    
                    // Sediment Concentration at Left/Right Sides of (j+1/2)
                    for (k=4;k<=nn;k++)
                    {
                        if (j!=0 && NEARSOLIDY[i][j]!=0)
                        {
                            // Reduce to 1st Order Flux if Near Interior Boundary
                            MUSCLextrap(SED[k-3][i][j-1],SED[k-3][i][j],SED[k-3][i][j+1],SED[k-3][i][j+2],ULR,k-1,1);
                        }
                        else
                        {
                            // Not Near Interior Boundary
                            MUSCLextrap(SED[k-3][i][j-1],SED[k-3][i][j],SED[k-3][i][j+1],SED[k-3][i][j+2],ULR,k-1,order);
                        }
                    }
                    
                    // HLLC Flux
                    Smaxtemp=hllc(ULR[0][0],ULR[0][1],ULR[2][0],ULR[2][1],ULR[1][0],ULR[1][1],TFLUX,GZ[i][j],hfilm);
                    
                    // Transparent Interior B.C. On the Left: Reduce to 1st Order Flux
                    if (j>0 && SOLID[i][j]!=0 && SOLID[i][j+1]==0)
                    {
                        if (SOLID[i][j]==1)
                        {
                            Smaxtemp=hllc(U[1][i][j+1],U[1][i][j+1],U[3][i][j+1],U[3][i][j+1],U[2][i][j+1],U[2][i][j+1],TFLUX,GZ[i][j],hfilm);
                        }
                        if (SOLID[i][j]==2)
                        {
                            Smaxtemp=hllc(U[1][i][j+1],U[1][i][j+1],-U[3][i][j+1],U[3][i][j+1],U[2][i][j+1],U[2][i][j+1],TFLUX,GZ[i][j],hfilm);
                        }
                    }
                    // Transparent Interior B.C. On the Right: Reduce to 1st Order Flux
                    if (j>0 && SOLID[i][j]==0 && SOLID[i][j+1]!=0)
                    {
                        if (SOLID[i][j+1]==1)
                        {
                            Smaxtemp=hllc(U[1][i][j],U[1][i][j],U[3][i][j],U[3][i][j],U[2][i][j],U[2][i][j],TFLUX,GZ[i][j],hfilm);
                        }
                         if (SOLID[i][j+1]==2)
                         {
                             Smaxtemp=hllc(U[1][i][j],U[1][i][j],U[3][i][j],-U[3][i][j],U[2][i][j],U[2][i][j],TFLUX,GZ[i][j],hfilm);
                         }
                    }
                    
                    FLUXY[1][i][j+1]=TFLUX[0];
                    FLUXY[3][i][j+1]=TFLUX[1];
                    FLUXY[2][i][j+1]=TFLUX[2];
                    for (k=4;k<=nn;k++)
                    {
                        // Upwind Flux for Each Particle Size Class
                        if (TFLUX[0]>0) FLUXY[k][i][j+1]=TFLUX[0]*ULR[k-1][0];
                        else FLUXY[k][i][j+1]=TFLUX[0]*ULR[k-1][1];
                    }
                    if (TFLUX[0]>0) FLUXY[k][i][j+1]=TFLUX[0]*TRACER[i][j];
                    else FLUXY[k][i][j+1]=TFLUX[0]*TRACER[i][j+1];
                    
                    if (Smaxtemp>Smax) Smax=Smaxtemp;
                    
                }
                
                // Source Terms
                if (i!=0 && j!=0 && SOLID[i][j]==0)
                {
                    BETAX[i][j]=0;
                    BETAY[i][j]=0;
                    FAIL[i][j]=0;
                    
                    // Topographic Forcing in X-Direction
                    if (i!=1 && i!=nx)
                    {
                        SLOPEX[i][j]=0.5*oneoverdx*(TOPO[i+1][j]-TOPO[i-1][j]);
                        GX[i][j]=g*SLOPEX[i][j]/pow(SLOPEX[i][j]*SLOPEX[i][j]+SLOPEY[i][j]*SLOPEY[i][j]+1,0.5);
                        SX[i][j]=-GX[i][j]*U[1][i][j];
                    }
                    else if (i==nx)
                    {
                        SLOPEX[i][j]=oneoverdx*(TOPO[i][j]-TOPO[i-1][j]);
                        GX[i][j]=g*SLOPEX[i][j]/pow(SLOPEX[i][j]*SLOPEX[i][j]+SLOPEY[i][j]*SLOPEY[i][j]+1,0.5);
                        SX[i][j]=-GX[i][j]*U[1][i][j];
                    }
                    else
                    {
                        SLOPEX[i][j]=oneoverdx*(TOPO[i+1][j]-TOPO[i][j]);
                        GX[i][j]=g*SLOPEX[i][j]/pow(SLOPEX[i][j]*SLOPEX[i][j]+SLOPEY[i][j]*SLOPEY[i][j]+1,0.5);
                        SX[i][j]=-GX[i][j]*U[1][i][j];
                    }
                    
                    // Topographic Forcing in Y-Direction
                    if (j!=1 && j!=ny)
                    {
                        SLOPEY[i][j]=0.5*oneoverdx*(TOPO[i][j+1]-TOPO[i][j-1]);
                        GY[i][j]=g*SLOPEY[i][j]/pow(SLOPEX[i][j]*SLOPEX[i][j]+SLOPEY[i][j]*SLOPEY[i][j]+1,0.5);
                        SY[i][j]=-GY[i][j]*U[1][i][j];
                    }
                    else if (j==ny)
                    {
                        SLOPEY[i][j]=oneoverdx*(TOPO[i][j]-TOPO[i][j-1]);
                        GY[i][j]=g*SLOPEY[i][j]/pow(SLOPEX[i][j]*SLOPEX[i][j]+SLOPEY[i][j]*SLOPEY[i][j]+1,0.5);
                        SY[i][j]=-GY[i][j]*U[1][i][j];
                    }
                    else
                    {
                        SLOPEY[i][j]=oneoverdx*(TOPO[i][j+1]-TOPO[i][j]);
                        GY[i][j]=g*SLOPEY[i][j]/pow(SLOPEX[i][j]*SLOPEX[i][j]+SLOPEY[i][j]*SLOPEY[i][j]+1,0.5);
                        SY[i][j]=-GY[i][j]*U[1][i][j];
                    }
                    
                    GZ[i][j]=g/pow(SLOPEX[i][j]*SLOPEX[i][j]+SLOPEY[i][j]*SLOPEY[i][j]+1,0.5);
                    SLOPE[i][j]=pow(SLOPEX[i][j]*SLOPEX[i][j]+SLOPEY[i][j]*SLOPEY[i][j],0.5);
                    
                    // Momentum Forcing Terms From Variable Density
                    GAMMAX[i][j]=0;
                    
                    if (i>1 && i<nx) GAMMAX[i][j]-=(rhos-RHO[i][j])*GZ[i][j]*U[1][i][j]*U[1][i][j]/(2*RHO[i][j])*.5*oneoverdx*(TOTSED[i+1][j]-TOTSED[i-1][j]);
                    
                    // Momentum Forcing Terms From Variable Density
                    GAMMAY[i][j]=0;
                    
                    if (j>1 && j<ny) GAMMAY[i][j]-=(rhos-RHO[i][j])*GZ[i][j]*U[1][i][j]*U[1][i][j]/(2*RHO[i][j])*.5*oneoverdx*(TOTSED[i][j+1]-TOTSED[i][j-1]);
                    
                    // Determine if Deposited Sediment is Stable in X-Direction
                    pgfx=DEP[i][j]*rho0*GZ[i][j]*0.5*oneoverdx*(DEP[IUP[i]][j]*cos(atan(fabs(SLOPEX[IUP[i]][j])))-DEP[IDOWN[i]][j]*cos(atan(fabs(SLOPEX[IDOWN[i]][j]))));
                    // PGF and Slope have same sign then PGF acts as a driving force
                    if (pgfx*SLOPEX[i][j]>0) pgfx=-1*fabs(pgfx);
                    // PGF and Slope have opposite sign then PGF acts as a resisting force
                    if (pgfx*SLOPEX[i][j]<0) pgfx=fabs(pgfx);
                    if (FAIL[i][j]==0 && TOTM[i][j]>26 && ((rho0*GZ[i][j]*DEP[i][j]+RHO[i][j]*GZ[i][j]*U[1][i][j])*sin(atan(fabs(SLOPEY[i][j])))>(cohesion+((rho0*GZ[i][j]*DEP[i][j]+RHO[i][j]*GZ[i][j]*U[1][i][j])*cos(atan(fabs(SLOPEY[i][j])))-rhow*GZ[i][j]*(DEP[i][j]+U[1][i][j]))*frictionangle+pgfx)))
                    {
                        //printf("!! Slope failure at (i=%d,j=%d): %lf\n",i,j,DEP[i][j]);
                        FAIL[i][j]=1;
                        FAILCNT[i][j]+=DEP[i][j];
                    }
                    
                    // Determine if Deposited Sediment is Stable in Y-Direction
                    pgfy=DEP[i][j]*rho0*GZ[i][j]*0.5*oneoverdx*(DEP[i][JUP[j]]*cos(atan(fabs(SLOPEY[i][JUP[j]])))-DEP[i][JDOWN[j]]*cos(atan(fabs(SLOPEY[i][JDOWN[j]]))));
                    // PGF and Slope have same sign then PGF acts as a driving force
                    if (pgfy*SLOPEY[i][j]>0) pgfy=-1*fabs(pgfy);
                    // PGF and Slope have opposite sign then PGF acts as a resisting force
                    if (pgfy*SLOPEY[i][j]<0) pgfy=fabs(pgfy);
                    if (FAIL[i][j]==0 && TOTM[i][j]>26 && ((rho0*GZ[i][j]*DEP[i][j]+RHO[i][j]*GZ[i][j]*U[1][i][j])*sin(atan(fabs(SLOPEY[i][j])))>(cohesion+((rho0*GZ[i][j]*DEP[i][j]+RHO[i][j]*GZ[i][j]*U[1][i][j])*cos(atan(fabs(SLOPEY[i][j])))-rhow*GZ[i][j]*(DEP[i][j]+U[1][i][j]))*frictionangle+pgfy)))
                    {
                        //printf("!! Slope failure at (i=%d,j=%d): %lf\n",i,j,DEP[i][j]);
                        FAIL[i][j]=1;
                        FAILCNT[i][j]+=DEP[i][j];
                    }
                    
                    // Debris Flow Forcing Terms
                    DFMASK[i][j]=fmin(1,(TOTSED[i][j]-0.2)/0.2);
                    if (U[1][i][j]<9*h0 || TOTSED[i][j]<0.2) DFMASK[i][j]=0;
                    BETAX[i][j]=-(1-lambda)*GZ[i][j]*U[1][i][j]*frictionangle;
                    BETAY[i][j]=-(1-lambda)*GZ[i][j]*U[1][i][j]*frictionangle;
                    if (U[2][i][j]<0) BETAX[i][j]=-1*BETAX[i][j];
                    if (U[3][i][j]<0) BETAY[i][j]=-1*BETAY[i][j];
                    
                    for (k=1;k<=sclass;k++)
                    {
                        E[1][k][i][j]=0;
                        E[2][k][i][j]=0;
                        E[3][k][i][j]=0;
                        E[4][k][i][j]=0;
                        
                        // Fractional Shielding of Original Soil by Deposited Sediment
                        mtstar=mtstar0;
                        if (U[1][i][j]/cos(atan(fabs(SLOPE[i][j])))>h0) mtstar=mtstar0*pow(h0/(U[1][i][j]/cos(atan(fabs(SLOPE[i][j])))),b);
                        H[i][j]=fmin(1,TOTM[i][j]/mtstar);
                        
                        // Rainsplash (Note: No Rainasplash When Sediment Concentrations Are High)
                        // NOTE: No Rainsplash if grain diameter exceeds critical grain diameter based on rainfall intensity
                        if (U[1][i][j]>10*hfilm && TOTSED[i][j]<0.6  && DS[k]<0.004)
                        {
                            // Rainsplash Detachment
                            E[1][k][i][j]=(1-H[i][j])*P[k][i][j]*AMAP[i][j]*(Cv*pi+(1-Cv))*R1;
                            
                            // Rainsplash Re-Detachment
                            if (TOTM[i][j]>1e-8 && M[k][i][j]>1e-8)
                            {
                                E[2][k][i][j]=H[i][j]*M[k][i][j]/TOTM[i][j]*ADMAP[i][j]*(Cv*pi+(1-Cv))*R1;
                            }
                            
                            // Modification of Rainsplash Detachment/Re-Detachment To Account For Flow Depth
                            // Need to Convert Flow Depth From Slope Normal to Vertical
                            if ((U[1][i][j]/cos(atan(fabs(SLOPE[i][j]))))>h0)
                            {
                                E[1][k][i][j]=E[1][k][i][j]*pow(h0/(U[1][i][j]/cos(atan(fabs(SLOPE[i][j])))),b);
                                E[2][k][i][j]=E[2][k][i][j]*pow(h0/(U[1][i][j]/cos(atan(fabs(SLOPE[i][j])))),b);
                            }
                        }
                        
                        // Friciton Slope
                        if (U[1][i][j]<=hcfriction) MANNING[i][j]=ROUGHNESS[i][j]*pow(U[1][i][j]/hcfriction,depthdependentexponent);
                        else MANNING[i][j]=ROUGHNESS[i][j];
                        
                        q=pow(U[2][i][j]*U[2][i][j]+U[3][i][j]*U[3][i][j],0.5);
                        SF[i][j]=MANNING[i][j]*MANNING[i][j]*q*q/pow(U[1][i][j],3.333);                      // Manning friciton slope
                        
                        if (U[1][i][j]>10*hfilm && TOTSED[i][j]<0.6)
                        {
                            // Entrainment Due to Runoff
                            E[3][k][i][j]=streampowererosion(U[2][i][j],U[3][i][j],U[1][i][j],BETA[1][i][j],RHO[i][j],GZ[i][j],TOTSED[i][j],SF[i][j],DS[k],UC[i][j]);
                            E[3][k][i][j]=(1-H[i][j])*P[k][i][j]*F/JMAP[i][j]*E[3][k][i][j];
                            
                            // Re-Entrainment Due to Runoff
                            if (TOTM[i][j]>1e-6 && M[k][i][j]>0)
                            {
                                E[4][k][i][j]=streampowererosion(U[2][i][j],U[3][i][j],U[1][i][j],BETA[2][i][j],RHO[i][j],GZ[i][j],TOTSED[i][j],SF[i][j],DS[k],UC2[i][j]);
                                E[4][k][i][j]=H[i][j]*M[k][i][j]/TOTM[i][j]*F/((rhos-RHO[i][j])/rhos*GZ[i][j]*U[1][i][j])*E[4][k][i][j];
                            }
                            else
                            {
                                E[4][k][i][j]=0;
                            }
                        }
                        
                        // Enforce Limited Soil Depth If Desired
                        if (TOPO[i][j]-ITOPO[i][j]<-1*maxsoilthickness)
                        {
                            E[1][k][i][j]=0;
                            E[3][k][i][j]=0;
                        }
                        
                        // Deposition Rate
                        VF[k][i][j]=settlingvelocity(U[1][i][j],DS[k],TOTSED[i][j],nu,rhos,rhow,GZ[i][j],DS[k]);
                        VF[k][i][j]=VF[k][i][j]*pow(1-fmin(TOTSED[i][j],0.99),4);
                        
                        // Modify settling velocity in shallow interill areas (e.g. McGuire et al., 2016)
                        if (DS[k]>0.0001 && U[1][i][j]<9*h0) VF[k][i][j]=2e-3+fmax(U[1][i][j]-3*h0,0)*(VF[k][i][j]-2e-3)/(6*h0);
                        if (DS[k]>0.002 && U[1][i][j]<9*h0) VF[k][i][j]=2e-3+fmax(U[1][i][j]-3*h0,0)*(VF[k][i][j]-2e-3)/(6*h0);
                        
                        D[k][i][j]=VF[k][i][j]*SED[k][i][j];
                        
                        // No Entrainment or Deposition In Debris Flows
                        if (U[1][i][j]>9*h0 && TOTSED[i][j]>cthreshold)
                        {
                            E[1][k][i][j]=0;
                            E[2][k][i][j]=0;
                            E[3][k][i][j]=0;
                            E[4][k][i][j]=0;
                            VF[k][i][j]=0;
                            D[k][i][j]=0;
                        }
                        
                    }
                    if (fabs(VEL[1][i][j])>velmax) velmax=fabs(VEL[1][i][j]);
                    if (fabs(VEL[2][i][j])>velmax) velmax=fabs(VEL[2][i][j]);
                }
                
            }
        
        // Adaptive Time Step
        dt=cstable*dx/(Smax);
        if (dt>maxstep) dt=maxstep;
        if (t+dt>tend) dt=tend-t;
        courant=(Smax)*dt*oneoverdx;
        if (courant>1)  {printf("!!! CFL Error at t=%lf\n",t); t=2*tend;}
        
        // Rutter Interception Model
        Di=Ki*exp(gi*(CS-Si));
        CS=CS+dt*(1-pi)*R1-dt*Di;        // Canopy Storage
        
        // Time Update
        maxtopochange=0;
        #pragma omp parallel for private(j,k,chtemp,redetach,beta)
        for (i=1;i<=nx;i++)
            for (j=1;j<=ny;j++)
            {
                chtemp=0; redetach=0; SF[i][j]=0; beta=0;
                if (SOLID[i][j]==0)
                {
                    // Time Update
                    U[1][i][j] = U[1][i][j]
                            - dt*oneoverdx*(FLUXX[1][i+1][j]-FLUXX[1][i][j])
                            - dt*oneoverdx*(FLUXY[1][i][j+1]-FLUXY[1][i][j])
                            + dt*R;
                    
                    U[2][i][j] = U[2][i][j]
                            - dt*oneoverdx*(FLUXX[2][i+1][j]-FLUXX[2][i][j])
                            - dt*oneoverdx*(FLUXY[2][i][j+1]-FLUXY[2][i][j])
                            + dt*SX[i][j]
                            + dt*GAMMAX[i][j]
                            + dt*DFMASK[i][j]*BETAX[i][j];
                    
                    U[3][i][j] = U[3][i][j]
                            - dt*oneoverdx*(FLUXX[3][i+1][j]-FLUXX[3][i][j])
                            - dt*oneoverdx*(FLUXY[3][i][j+1]-FLUXY[3][i][j])
                            + dt*SY[i][j]
                            + dt*GAMMAY[i][j]
                            + dt*DFMASK[i][j]*BETAY[i][j];
                    
                    TRACER[i][j] = TRACER[i][j]
                            - dt*oneoverdx*(FLUXX[nn+1][i+1][j]-FLUXX[nn+1][i][j])
                            - dt*oneoverdx*(FLUXY[nn+1][i][j+1]-FLUXY[nn+1][i][j])
                            + FAIL[i][j]*U[1][i][j];
                    
                    // Infiltration
                    ZF[i][j]=VINF[i][j]/(THETAS[i][j]-THETA0[i][j]);                                        // Depth of wetting front
                    INFL[i][j]=fmin(U[1][i][j]-hfilm,dt*KS[i][j]*(ZF[i][j]+HF[i][j]+U[1][i][j])/ZF[i][j]);  // Infiltration Capacity
                    VINF[i][j]+=INFL[i][j];                                                                 // Infiltrated depth
                    U[1][i][j]-=INFL[i][j];
                    
                    for (k=4;k<=nn;k++)
                    {
                        U[k][i][j] = U[k][i][j]
                                - dt*oneoverdx*(FLUXX[k][i+1][j]-FLUXX[k][i][j])
                                - dt*oneoverdx*(FLUXY[k][i][j+1]-FLUXY[k][i][j]);
                    }
                    
                    TOTM[i][j]=0;
                    DEP[i][j]=0;
                    
                    for (k=4;k<=nn;k++)
                    {
                        if (U[1][i][j]<10*hfilm)
                        {
                            // Flow Too Shallow To Transport Sediment For Particle Size Class k-3
                            M[k-3][i][j]+=dt*E[1][k-3][i][j]+U[k][i][j];
                            U[k][i][j]=0;
                        }
                        else
                        {
                            chtemp=U[k][i][j];                                                      // Store current value of ch
                            redetach=fmin(M[k-3][i][j],dt*(E[2][k-3][i][j]+E[4][k-3][i][j]));       // Total Mass of Re-detached Sediment <= Total Mass of Deposited Sediment
                            
                            U[k][i][j] = exp(-VF[k-3][i][j]/U[1][i][j]*(0.5*dt))*U[k][i][j];
                            
                            U[k][i][j] = U[k][i][j]
                                    + dt*(E[1][k-3][i][j]+E[3][k-3][i][j])
                                    + redetach;
                            
                            U[k][i][j] = exp(-VF[k-3][i][j]/U[1][i][j]*(0.5*dt))*U[k][i][j];
                            
                            TOPO[i][j] = TOPO[i][j]
                                    - 1/(rhos*(1-phi))*(U[k][i][j]-chtemp);
                            
                            M[k-3][i][j] = M[k-3][i][j]
                                    + dt*(E[1][k-3][i][j]+E[3][k-3][i][j])
                                    - (U[k][i][j]-chtemp);
                            
                            if (M[k-3][i][j]<-1e-6) printf("!! M_k < 0 %d\t %d\t %lf\t %lf\n",i,j,M[k-3][i][j],E[4][k-3][i][j]);
                            
                            D[k-3][i][j]=(chtemp-U[k][i][j])/dt;
                            
                            // If Sediment is Unstable, Add it to Water Column and Remove From Deposited Layer
                            if (FAIL[i][j]==1)
                            {
                                TOPO[i][j] = TOPO[i][j]-0.99*M[k-3][i][j]/(rhos*(1-phi));
                                U[k][i][j] = U[k][i][j]+0.99*M[k-3][i][j];
                                TOTM[i][j] = 0;
                                M[k-3][i][j] = 0.01*M[k-3][i][j];
                                DEP[i][j] = 0;
                            }
                        }
                        
                        U[1][i][j] = U[1][i][j]+(U[k][i][j]-chtemp)/(rhos*(1-phi));
                        
                        TOTM[i][j]+=M[k-3][i][j];                                          // Total Thickness of Deposited Sediment
                        DEP[i][j]+=M[k-3][i][j]/(rhos);                                    // Total Thickness of Deposited Sediment
                    }
                    
                    TOPOOLD[i][j]=TOPO[i][j];
                    TOTSED[i][j]=0;
                    
                    if (U[1][i][j]>hfilm)
                    {
                        VEL[1][i][j]=U[2][i][j]/U[1][i][j];
                        VEL[2][i][j]=U[3][i][j]/U[1][i][j];
                        for (k=4;k<=nn;k++)
                        {
                            SED[k-3][i][j]=U[k][i][j]/U[1][i][j];
                            TOTSED[i][j]+=SED[k-3][i][j]/rhos;                                  // Total Volumetric Sediment Concentration
                            if (SED[k-3][i][j]<-1e-5)
                            {
                                SED[k-3][i][j]=0;
                                printf("!! Sed < 0 \n");
                                printf("%lf\t %lf\t %lf\t %lf\t %lf\n",SED[k-3][i][j],E[1][k-3][i][j],E[2][k-3][i][j],E[3][k-3][i][j],E[4][k-3][i][j]);
                            }
                        }
                    }
                    else
                    {
                        U[1][i][j]=hfilm;
                        U[2][i][j]=0;
                        U[3][i][j]=0;
                        VEL[1][i][j]=0;
                        VEL[2][i][j]=0;
                        for (k=4;k<=nn;k++)
                        {
                            // Add Sediment From Water Column To Bed
                            TOPO[i][j] = TOPO[i][j]
                                    + 1/(rhos*(1-phi))*U[k][i][j];
                            
                            M[k-3][i][j] = M[k-3][i][j]
                                    + U[k][i][j];
                            
                            U[k][i][j]=0;
                            SED[k-3][i][j]=0;
                        }
                    }
                    
                    // Implicit Friction Update (LeVeque, 2011)
                    if (U[1][i][j]>hfilm)
                    {
                        if (U[1][i][j]<=hcfriction) MANNING[i][j]=ROUGHNESS[i][j]*pow(U[1][i][j]/hcfriction,depthdependentexponent);
                        else MANNING[i][j]=ROUGHNESS[i][j];
                        
                        SF[i][j]=MANNING[i][j]*MANNING[i][j]*GZ[i][j]*pow(U[2][i][j]*U[2][i][j]+U[3][i][j]*U[3][i][j],0.5)*1/pow(U[1][i][j],2.333);
                        U[2][i][j]=U[2][i][j]/(1+dt*SF[i][j]);
                        U[3][i][j]=U[3][i][j]/(1+dt*SF[i][j]);
                        
                        VEL[1][i][j]=U[2][i][j]/U[1][i][j];
                        VEL[2][i][j]=U[3][i][j]/U[1][i][j];
                    }
                }
                else
                {
                    // Pixel Is Not In Computational Domain
                    VEL[1][i][j]=0;
                    VEL[2][i][j]=0;
                    for (k=1;k<=nn;k++)
                        U[k][i][j]=0;
                    for (k=1;k<=sclass;k++)
                        SED[k][i][j]=0;
                }
                
                RHO[i][j]=rhow*(1-TOTSED[i][j])+rhos*TOTSED[i][j];
                if (RHO[i][j]>0.8*rhos) RHO[i][j]=0.8*rhos;
                
                for (k=1;k<=nn;k++)
                    UOLD[k][i][j]=U[k][i][j];
                
                VELOLD[1][i][j]=VEL[1][i][j];
                VELOLD[2][i][j]=VEL[2][i][j];
                
                if (pow(VEL[1][i][j]*VEL[1][i][j]+VEL[2][i][j]*VEL[2][i][j],0.5)>MAXVEL[i][j]) MAXVEL[i][j]=pow(VEL[1][i][j]*VEL[1][i][j]+VEL[2][i][j]*VEL[2][i][j],0.5);
                if (U[1][i][j]>MAXDEPTH[i][j]) MAXDEPTH[i][j]=U[1][i][j];
            }

        t+=dt;
        
        if (t>stageinterval)
        {
            // Find Basin Outlet
            outlet=1; outlettopo=1; outlety=1; outletelev=99999; outletdepth=0;
            for (i=1;i<=nx;i++)
            {
                if (SOLID[i][outlety]==0 && U[1][i][outlety]>outletdepth)
                {
                    outletdepth=U[1][i][outlety];
                    outlet=i;                                       // Save outlet location
                }
                if (SOLID[i][outlety]==0 && TOPO[i][outlety]<outletelev)
                {
                    outletelev=TOPO[i][outlety];
                    outlettopo=i;                                       // Save outlet location
                }
            }
            
            stageinterval+=stageint;
            stagecnt++;
            STAGE[1][stagecnt]=t;
            STAGE[2][stagecnt]=U[1][outlet][outlety];
            STAGE[3][stagecnt]=pow(VEL[1][outlet][outlety]*VEL[1][outlet][outlety]+VEL[2][outlet][outlety]*VEL[2][outlet][outlety],0.5);
            STAGE[4][stagecnt]=TOTSED[outlet][outlety];
            STAGE[5][stagecnt]=U[1][outlettopo][outlety];
            STAGE[6][stagecnt]=pow(VEL[1][outlettopo][outlety]*VEL[1][outlettopo][outlety]+VEL[2][outlettopo][outlety]*VEL[2][outlettopo][outlety],0.5);
            STAGE[7][stagecnt]=TOTSED[outlettopo][outlety];
            
            if (nn>=4) STAGE[8][stagecnt]=M[1][outlet][outlety];
            else STAGE[8][stagecnt]=0;
            if (nn>=5) STAGE[9][stagecnt]=M[2][outlet][outlety];
            else STAGE[9][stagecnt]=0;
            if (nn>=6) STAGE[10][stagecnt]=M[3][outlet][outlety];
            else STAGE[10][stagecnt]=0;
            if (nn>=7) STAGE[11][stagecnt]=M[4][outlet][outlety];
            else STAGE[11][stagecnt]=0;
            if (nn>=8) STAGE[12][stagecnt]=M[5][outlet][outlety];
            else STAGE[12][stagecnt]=0;
            if (nn>=9) STAGE[13][stagecnt]=M[6][outlet][outlety];
            else STAGE[13][stagecnt]=0;
            if (nn>=10) STAGE[14][stagecnt]=M[7][outlet][outlety];
            else STAGE[14][stagecnt]=0;
            if (nn>=11) STAGE[15][stagecnt]=M[8][outlet][outlety];
            else STAGE[15][stagecnt]=0;
            if (nn>=12) STAGE[16][stagecnt]=M[9][outlet][outlety];
            else STAGE[16][stagecnt]=0;
            if (nn>=13) STAGE[17][stagecnt]=M[10][outlet][outlety];
            else STAGE[17][stagecnt]=0;
            if (nn>=4) STAGE[18][stagecnt]=U[4][outlet][outlety];
            else STAGE[18][stagecnt]=0;
            if (nn>=5) STAGE[19][stagecnt]=U[5][outlet][outlety];
            else STAGE[19][stagecnt]=0;
            if (nn>=6) STAGE[20][stagecnt]=U[6][outlet][outlety];
            else STAGE[20][stagecnt]=0;
            if (nn>=7) STAGE[21][stagecnt]=U[7][outlet][outlety];
            else STAGE[21][stagecnt]=0;
            if (nn>=8) STAGE[22][stagecnt]=U[8][outlet][outlety];
            else STAGE[22][stagecnt]=0;
            if (nn>=9) STAGE[23][stagecnt]=U[9][outlet][outlety];
            else STAGE[23][stagecnt]=0;
            if (nn>=10) STAGE[24][stagecnt]=U[10][outlet][outlety];
            else STAGE[24][stagecnt]=0;
            if (nn>=11) STAGE[25][stagecnt]=U[11][outlet][outlety];
            else STAGE[25][stagecnt]=0;
            if (nn>=12) STAGE[26][stagecnt]=U[12][outlet][outlety];
            else STAGE[26][stagecnt]=0;
            if (nn>=13) STAGE[27][stagecnt]=U[13][outlet][outlety];
            else STAGE[27][stagecnt]=0;
            
            SAVEFLOW[1][stagecnt]=U[1][XPOINT[1]][YPOINT[1]];
            SAVEFLOW[2][stagecnt]=pow(VEL[1][XPOINT[1]][YPOINT[1]]*VEL[1][XPOINT[1]][YPOINT[1]]+VEL[2][XPOINT[1]][YPOINT[1]]*VEL[2][XPOINT[1]][YPOINT[1]],0.5);
            SAVEFLOW[3][stagecnt]=TOTSED[XPOINT[1]][YPOINT[1]];
            
            SAVEFLOW[4][stagecnt]=U[1][XPOINT[2]][YPOINT[2]];
            SAVEFLOW[5][stagecnt]=pow(VEL[1][XPOINT[2]][YPOINT[2]]*VEL[1][XPOINT[2]][YPOINT[2]]+VEL[2][XPOINT[2]][YPOINT[2]]*VEL[2][XPOINT[2]][YPOINT[2]],0.5);
            SAVEFLOW[6][stagecnt]=TOTSED[XPOINT[2]][YPOINT[2]];
            
            SAVEFLOW[7][stagecnt]=U[1][XPOINT[3]][YPOINT[3]];
            SAVEFLOW[8][stagecnt]=pow(VEL[1][XPOINT[3]][YPOINT[3]]*VEL[1][XPOINT[3]][YPOINT[3]]+VEL[2][XPOINT[3]][YPOINT[3]]*VEL[2][XPOINT[3]][YPOINT[3]],0.5);
            SAVEFLOW[9][stagecnt]=TOTSED[XPOINT[3]][YPOINT[3]];
            
            SAVEFLOW[10][stagecnt]=U[1][XPOINT[4]][YPOINT[4]];
            SAVEFLOW[11][stagecnt]=pow(VEL[1][XPOINT[4]][YPOINT[4]]*VEL[1][XPOINT[4]][YPOINT[4]]+VEL[2][XPOINT[4]][YPOINT[4]]*VEL[2][XPOINT[4]][YPOINT[4]],0.5);
            SAVEFLOW[12][stagecnt]=TOTSED[XPOINT[4]][YPOINT[4]];
            
            SAVEFLOW[13][stagecnt]=U[1][XPOINT[5]][YPOINT[5]];
            SAVEFLOW[14][stagecnt]=pow(VEL[1][XPOINT[5]][YPOINT[5]]*VEL[1][XPOINT[5]][YPOINT[5]]+VEL[2][XPOINT[5]][YPOINT[5]]*VEL[2][XPOINT[5]][YPOINT[5]],0.5);
            SAVEFLOW[15][stagecnt]=TOTSED[XPOINT[5]][YPOINT[5]];
            
            SAVEFLOW[16][stagecnt]=U[1][XPOINT[6]][YPOINT[6]];
            SAVEFLOW[17][stagecnt]=pow(VEL[1][XPOINT[6]][YPOINT[6]]*VEL[1][XPOINT[6]][YPOINT[6]]+VEL[2][XPOINT[6]][YPOINT[6]]*VEL[2][XPOINT[6]][YPOINT[6]],0.5);
            SAVEFLOW[18][stagecnt]=TOTSED[XPOINT[6]][YPOINT[6]];
        }

        // Intermediate Results
        if (t>=printinterval)
        {
            printf("%lf\t %lf\t %lf\t %lf\n",t,courant,dt,R1*3600*1000);
            printinterval+=printint;
            for (j=1;j<=ny;j++)
                for (i=1;i<=nx;i++)
                {
                    if (i!=nx) fprintf(fpm0,"%10.10lf\t",TOPO[i][j]);
                    if (i==nx) fprintf(fpm0,"%10.10lf\n",TOPO[i][j]);
                    if (i!=nx) fprintf(fpm1,"%10.10lf\t",U[1][i][j]);
                    if (i==nx) fprintf(fpm1,"%10.10lf\n",U[1][i][j]);
                    if (i!=nx) fprintf(fpm2,"%10.10lf\t",pow(VEL[1][i][j]*VEL[1][i][j]+VEL[2][i][j]*VEL[2][i][j],0.5));
                    if (i==nx) fprintf(fpm2,"%10.10lf\n",pow(VEL[1][i][j]*VEL[1][i][j]+VEL[2][i][j]*VEL[2][i][j],0.5));
                    if (i!=nx) fprintf(fpm3,"%10.10lf\t",TOTSED[i][j]);
                    if (i==nx) fprintf(fpm3,"%10.10lf\n",TOTSED[i][j]);
                    if (i!=nx) fprintf(fpm4,"%10.10lf\t",E[3][1][i][j]);
                    if (i==nx) fprintf(fpm4,"%10.10lf\n",E[3][1][i][j]);
                    if (i!=nx) fprintf(fpm5,"%10.10lf\t",FAILCNT[i][j]);
                    if (i==nx) fprintf(fpm5,"%10.10lf\n",FAILCNT[i][j]);
                    if (i!=nx) fprintf(fpm6,"%10.10lf\t",TOTM[i][j]);
                    if (i==nx) fprintf(fpm6,"%10.10lf\n",TOTM[i][j]);
                }
            
        }

    }
    
    time_t t1 = time(NULL);
    printf ("Elapsed wall clock time: %ld\n", (long) (t1 - t0));
    
    // Print Final Results
    fprintf(fp0,"%d\n",cflstat);
    
    fp1=fopen("./topo","w");
    fp2=fopen("./depth","w");
    fp3=fopen("./uh","w");
    fp4=fopen("./vh","w");
    fp5=fopen("./ch","w");
    fp6=fopen("./m","w");
    fp7=fopen("./c","w");
    fp8=fopen("./vel","w");
    fp9=fopen("./E1","w");
    fp10=fopen("./E2","w");
    fp11=fopen("./E3","w");
    fp12=fopen("./E4","w");
    fp13=fopen("./D","w");
    fp14=fopen("./c1","w");
    fp15=fopen("./c2","w");
    fp16=fopen("./c3","w");
    fp17=fopen("./m1","w");
    fp18=fopen("./m2","w");
    fp19=fopen("./m3","w");
    fp20=fopen("./solid","w");
    fp21=fopen("./stage","w");
    fp22=fopen("./maxvel","w");
    fp23=fopen("./m4","w");
    fp24=fopen("./m5","w");
    fp25=fopen("./m6","w");
    fp26=fopen("./m7","w");
    fp27=fopen("./m8","w");
    fp28=fopen("./m9","w");
    fp29=fopen("./m10","w");
    fp30=fopen("./infl","w");
    fp31=fopen("./saveflow","w");
    fp32=fopen("./maxdepth","w");
    
    for (j=1;j<=ny;j++)
        for (i=1;i<=nx;i++)
        {
            if (i!=nx)
            {
                fprintf(fp1,"%10.10lf\t",TOPO[i][j]);
                fprintf(fp2,"%10.10lf\t",U[1][i][j]);
                fprintf(fp3,"%10.10lf\t",U[2][i][j]);
                fprintf(fp4,"%10.10lf\t",U[3][i][j]);
                fprintf(fp5,"%10.10lf\t",U[4][i][j]);
                fprintf(fp6,"%10.10lf\t",TOTM[i][j]);
                fprintf(fp7,"%10.10lf\t",TOTSED[i][j]);
                fprintf(fp8,"%10.10lf\t",pow(VEL[1][i][j]*VEL[1][i][j]+VEL[2][i][j]*VEL[2][i][j],0.5));
                fprintf(fp9,"%10.10lf\t",GX[i][j]);
                fprintf(fp10,"%10.10lf\t",GY[i][j]);
                fprintf(fp11,"%10.10lf\t",GZ[i][j]);
                fprintf(fp12,"%10.10lf\t",E[4][1][i][j]);
                fprintf(fp13,"%10.10lf\t",FAILCNT[i][j]);
                fprintf(fp22,"%10.10lf\t",MAXVEL[i][j]);
                fprintf(fp32,"%10.10lf\t",MAXDEPTH[i][j]);
                fprintf(fp30,"%10.10lf\t",INFL[i][j]);
                if (nn>=4)
                {
                    fprintf(fp14,"%10.10lf\t",U[4][i][j]);
                    fprintf(fp17,"%10.10lf\t",M[1][i][j]);
                }
                else
                {
                    fprintf(fp14,"%10.10lf\t",0.0);
                    fprintf(fp17,"%10.10lf\t",0.0);
                }
                if (nn>=5)
                {
                    fprintf(fp15,"%10.10lf\t",U[5][i][j]);
                    fprintf(fp18,"%10.10lf\t",M[2][i][j]);
                }
                else
                {
                    fprintf(fp15,"%10.10lf\t",0.0);
                    fprintf(fp18,"%10.10lf\t",0.0);
                }
                if (nn>=6)
                {
                    fprintf(fp16,"%10.10lf\t",U[6][i][j]);
                    fprintf(fp19,"%10.10lf\t",M[3][i][j]);
                }
                else
                {
                    fprintf(fp16,"%10.10lf\t",0.0);
                    fprintf(fp19,"%10.10lf\t",0.0);
                }
                if (nn>=7)
                {
                    fprintf(fp23,"%10.10lf\t",M[4][i][j]);
                }
                else
                {
                    fprintf(fp23,"%10.10lf\t",0.0);
                }
                if (nn>=8)
                {
                    fprintf(fp24,"%10.10lf\t",M[5][i][j]);
                }
                else
                {
                    fprintf(fp24,"%10.10lf\t",0.0);
                }
                if (nn>=9)
                {
                    fprintf(fp25,"%10.10lf\t",M[6][i][j]);
                }
                else
                {
                    fprintf(fp25,"%10.10lf\t",0.0);
                }
                if (nn>=10)
                {
                    fprintf(fp26,"%10.10lf\t",M[7][i][j]);
                }
                else
                {
                    fprintf(fp26,"%10.10lf\t",0.0);
                }
                if (nn>=11)
                {
                    fprintf(fp27,"%10.10lf\t",M[8][i][j]);
                }
                else
                {
                    fprintf(fp27,"%10.10lf\t",0.0);
                }
                if (nn>=12)
                {
                    fprintf(fp28,"%10.10lf\t",M[9][i][j]);
                }
                else
                {
                    fprintf(fp28,"%10.10lf\t",0.0);
                }
                if (nn>=13)
                {
                    fprintf(fp29,"%10.10lf\t",M[10][i][j]);
                }
                else
                {
                    fprintf(fp29,"%10.10lf\t",0.0);
                }
                fprintf(fp20,"%d\t",SOLID[i][j]);
            }
            if (i==nx)
            {
                fprintf(fp1,"%10.10lf\n",TOPO[i][j]);
                fprintf(fp2,"%10.10lf\n",U[1][i][j]);
                fprintf(fp3,"%10.10lf\n",U[2][i][j]);
                fprintf(fp4,"%10.10lf\n",U[3][i][j]);
                fprintf(fp5,"%10.10lf\n",U[4][i][j]);
                fprintf(fp6,"%10.10lf\n",TOTM[i][j]);
                fprintf(fp7,"%10.10lf\n",TOTSED[i][j]);
                fprintf(fp8,"%10.10lf\n",pow(VEL[1][i][j]*VEL[1][i][j]+VEL[2][i][j]*VEL[2][i][j],0.5));
                fprintf(fp9,"%10.10lf\n",GX[i][j]);
                fprintf(fp10,"%10.10lf\n",GY[i][j]);
                fprintf(fp11,"%10.10lf\n",GZ[i][j]);
                fprintf(fp12,"%10.10lf\n",E[4][1][i][j]);
                fprintf(fp13,"%10.10lf\n",FAILCNT[i][j]);
                fprintf(fp22,"%10.10lf\n",MAXVEL[i][j]);
                fprintf(fp32,"%10.10lf\n",MAXDEPTH[i][j]);
                fprintf(fp30,"%10.10lf\n",INFL[i][j]);
                if (nn>=4)
                {
                    fprintf(fp14,"%10.10lf\n",U[4][i][j]);
                    fprintf(fp17,"%10.10lf\n",M[1][i][j]);
                }
                else
                {
                    fprintf(fp14,"%10.10lf\n",0.0);
                    fprintf(fp17,"%10.10lf\n",0.0);
                }
                if (nn>=5)
                {
                    fprintf(fp15,"%10.10lf\n",U[5][i][j]);
                    fprintf(fp18,"%10.10lf\n",M[2][i][j]);
                }
                else
                {
                    fprintf(fp15,"%10.10lf\n",0.0);
                    fprintf(fp18,"%10.10lf\n",0.0);
                }
                if (nn>=6)
                {
                    fprintf(fp16,"%10.10lf\n",U[6][i][j]);
                    fprintf(fp19,"%10.10lf\n",M[3][i][j]);
                }
                else
                {
                    fprintf(fp16,"%10.10lf\n",0.0);
                    fprintf(fp19,"%10.10lf\n",0.0);
                }
                if (nn>=7)
                {
                    fprintf(fp23,"%10.10lf\n",M[4][i][j]);
                }
                else
                {
                    fprintf(fp23,"%10.10lf\n",0.0);
                }
                if (nn>=8)
                {
                    fprintf(fp24,"%10.10lf\n",M[5][i][j]);
                }
                else
                {
                    fprintf(fp24,"%10.10lf\n",0.0);
                }
                if (nn>=9)
                {
                    fprintf(fp25,"%10.10lf\n",M[6][i][j]);
                }
                else
                {
                    fprintf(fp25,"%10.10lf\n",0.0);
                }
                if (nn>=10)
                {
                    fprintf(fp26,"%10.10lf\n",M[7][i][j]);
                }
                else
                {
                    fprintf(fp26,"%10.10lf\n",0.0);
                }
                if (nn>=11)
                {
                    fprintf(fp27,"%10.10lf\n",M[8][i][j]);
                }
                else
                {
                    fprintf(fp27,"%10.10lf\n",0.0);
                }
                if (nn>=12)
                {
                    fprintf(fp28,"%10.10lf\n",M[9][i][j]);
                }
                else
                {
                    fprintf(fp28,"%10.10lf\n",0.0);
                }
                if (nn>=13)
                {
                    fprintf(fp29,"%10.10lf\n",M[10][i][j]);
                }
                else
                {
                    fprintf(fp29,"%10.10lf\n",0.0);
                }
                fprintf(fp20,"%d\n",SOLID[i][j]);
            }
        }
    
    for (j=1;j<=stagecnt;j++)
        for (i=1;i<=27;i++)
        {
            if (i!=27) fprintf(fp21,"%lf\t",STAGE[i][j]);
            if (i==27) fprintf(fp21,"%lf\n",STAGE[i][j]);
        }
    
    for (j=1;j<=stagecnt;j++)
        for (i=1;i<=18;i++)
        {
            if (i!=18) fprintf(fp31,"%lf\t",SAVEFLOW[i][j]);
            if (i==18) fprintf(fp31,"%lf\n",SAVEFLOW[i][j]);
        }
    
    
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);
    fclose(fp6);
    fclose(fp7);
    fclose(fp8);
    fclose(fp9);
    fclose(fp10);
    fclose(fp11);
    fclose(fp12);
    fclose(fp13);
    fclose(fp14);
    fclose(fp15);
    fclose(fp16);
    fclose(fp17);
    fclose(fp18);
    fclose(fp19);
    fclose(fp20);
    fclose(fp21);
    fclose(fp22);
    fclose(fp30);
    fclose(fp31);
    fclose(fp32);
    
    fclose(fpm0);
    fclose(fpm1);
    fclose(fpm2);
    fclose(fpm3);
    fclose(fpm4);
    fclose(fpm5);
    fclose(fpm6);
    
    fclose(fr10);
    
}
