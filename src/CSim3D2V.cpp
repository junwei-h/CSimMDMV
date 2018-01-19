/*This is program is to simulate Multi-dimensional, multi-variable and multi-cale heterogeneous reservoirs according to:
	Huang, J.-W., G. Bellefleur, and B. Milkereit, 2009. Seismic modeling of multidimensional heterogeneity scales of Mallik gas hydrate reservoirs, Northwest Territories of Canada, J. Geophys. Res., 114, B07306, doi:10.1029/2008JB006172.
	Huang, J-W, G. Bellefleur, and B. Milkereit, 2010. Stochastic Characterization of Multi-dimensional, Multi-variant, and Multi-scale Distribution of Heterogeneous Reservoir Rock Properties, Geophyscial Prospecting, submitted.
  
PLEASE DO NOT DISTRIBUTE. PLEASE REFER OTHER PEOPLE TO THE AUTHOR. 

Author	
Dr. Jun-Wei Huang, University of Toronto, Physics Department
422-60 St. George Street, Toronto, Ontario,Canada
ph: +1 416 978 2676,
fax: +1 416 978 7606, mailto:jhuang@physics.utoronto.ca
Homepage: http://www.physics.utoronto.ca/~jhuang/

If you want to publish results calculated with this program please
give a reference to the following papers:
	Huang, J.-W., G. Bellefleur, and B. Milkereit, 2009. Seismic modeling of multidimensional heterogeneity scales of Mallik gas hydrate reservoirs, Northwest Territories of Canada, J. Geophys. Res., 114, B07306, doi:10.1029/2008JB006172.
	Huang, J.-W., G. Bellefleur, and B. Milkereit, 2010. Stochastic Characterization of Multi-dimensional, Multi-variant, and Multi-scale Distribution of Heterogeneous Reservoir Rock Properties, Geophyscial Prospecting, submitted. 

----------------------------------------------------------------

History of important modifications:

15.01.2009	Version 1.0	original implementation
				submitted to JGR	
				Jun-Wei Huang
16.06.2009  Version 1.1 conditional and unconditional Monte Carlo simulation was added 	
				by Jun-Wei Huang
26.06.2010	Version 1.1 was cleaned for publication on Computers and Geosciences
				by Jun-Wei Huang
*/

#include "csim.h"

/*--------maximum array size: MaxA----------*/
# ifndef MaxA
# define MaxA	10000
# endif

/*--------FFT Accuracy Control Switch----------*/
#define FFT_NOFLOAT
#undef REAL
#ifdef TEST_FLOAT
# define REAL	float
# define fftn	fftnf
#else
# define REAL	double
#endif

/*--------FFT or IFFT Switch----------*/
#if 1
#define FORWARD_SCALE	-1.0
#define INVERSE_SCALE	0.0
#endif	 

int main(int argc, char *argv[])
{	
	clock_t Tstart;//start time recording
	Tstart = clock();
	int lenk[3],M=0,Nax,Nay,Nvar,condflag,my_rank,p,MaxMC;
	double HOZx,HOZy,simds,VET,kmax,kxs,kys,kzs,rNB,rNB2,rNB3,sum=0.0,sum2=0.0,facDes=0.0,Telapse,Mstart,temptime,MCtime;
	double shiftR,stdR,mymean,mystd;
	char *inpf,*outf, line[200],inp[20][80],logname[10], tempar[8][80];
	char rbinfile[100],ibinfile[100],cmd[100],outfile[100];
	double par[8],*cdfpar;
	double Global_facDes=0.0,Global_sum=0.0,Global_sum2=0.0; 
	double boreX0,boreY0,boreZ0,GHCon=0.0, Poro=0.0,GH[2],GHVol;
	
	float temp;
	
	double *PP=new double[5];
	double *SS=new double[5],*SP=new double[5];
	double S1;
	
	clock_t start;
	FILE *stream,*istream,*rstream;

	//===========start multi-processes========
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	
	MPI_Status status;
	if (my_rank==0)	{Telapse=MPI_Wtime();info();}
	
	if(argc<1){
		if(my_rank==0) {usage();}
		return 1;
	}
	
	inpf=argv[1];
	
    if( (stream  = fopen( inpf, "r" )) == NULL ) 
		if(my_rank==0){printf( "The file %10s was not opened\n",inpf );}
   else
		if(my_rank==0){printf( "The file %10s was opened\n",inpf );}

	while(fgets(line,100,stream)){
		if(line[0]!='#'){
			strcpy(inp[M++],line);
		}
	}
	fclose(stream);
	//==========Read input parameters==============	
	sscanf(inp[0],"%s%s%s%s%s",tempar,tempar+1,tempar+2,tempar+3,tempar+4);
	HOZx=atof(tempar[0]);HOZy=atof(tempar[1]);
	GH[0]=atof(tempar[2]);GH[1]=atof(tempar[3]);simds=atof(tempar[4]);
	VET=GH[1]-GH[0]+simds;
	lenk[0]=(int)round(pow(simds,-1)*HOZx);lenk[1]=(int)round(pow(simds,-1)*HOZy);lenk[2]=(int)round(pow(simds,-1)*VET);
	double s=lenk[0]*lenk[1]*lenk[2]*sizeof(double)/1024.0/1024.0;

	if(my_rank==0){
		if (lenk[2] % p != 0){
			printf("Vertical Sampling(%d) and Number of Processor (%d) does not match\n",lenk[2],p);
			return 1;}
		else if (fabs(round(pow(simds,-1)*VET)-VET/simds) >1e-9) {
			printf("Depth Interval (%f) and Sampling interval (%f) does not match\n",VET,simds);
			return 1;}
	}			
 
	sscanf(inp[1],"%s",logname);
	
	sscanf(inp[2],"%s",tempar); Nvar=atoi(tempar[0]);
	cdfpar=new (nothrow) double [Nvar*(Nvar+1)/2*3];
	
	sscanf(inp[3],"%s",tempar);
	condflag=atoi(tempar[0]);
	if(my_rank==0){
		if (condflag) 
			printf("Conditional Simulation\n");
		else
			printf("Nonconditional Simulation\n");
		
	}
	
	sscanf(inp[4],"%s%s%s",tempar,tempar+1,tempar+2);boreX0=atof(tempar[0]);boreY0=atof(tempar[1]);boreZ0=atof(tempar[2]);//origin point of UTM coordinate
	
	sscanf(inp[5],"%s%s%s",tempar,tempar+1,tempar+2);
	Nax=(int)round((atof(tempar[2])-atof(tempar[0]))/atof(tempar[1]))+1;
	if (my_rank==0)printf("ax=%f:%f:%f, total: %d\n",atof(tempar[0]),atof(tempar[1]),atof(tempar[2]),Nax);
	double *ax=new double[Nax];
	ax[0]=atof(tempar[0]);
	for (int i=1;i<Nax;i++)
		ax[i]=ax[0]+i*atof(tempar[1]);
		
	sscanf(inp[6],"%s%s%s",tempar,tempar+1,tempar+2);
	Nay=(int)round((atof(tempar[2])-atof(tempar[0]))/atof(tempar[1]))+1;
	if (my_rank==0) printf("ay=%f:%f:%f, total: %d\n",atof(tempar[0]),atof(tempar[1]),atof(tempar[2]),Nay);
	double *ay=new double[Nay];
	ay[0]=atof(tempar[0]);
	for (int i=1;i<Nay;i++)
		ay[i]=ay[0]+i*atof(tempar[1]);	
	if (my_rank==0) printf("total iteration: %d\n",Nax*Nay);
	
	sscanf(inp[7],"%s",tempar);MaxMC=atoi(tempar[0]);
	double Global_GHVol[MaxMC];
	
	for (int j=0;j<Nvar*(Nvar+1)/2;j++){
		sscanf(inp[8+j],"%s%s%s",tempar,tempar+1,tempar+2);
		for(int i=0;i<3;i++)
			cdfpar[j*3+i]=atof(tempar[i]);
		if(my_rank==0){printf("%d: az=%f,v=%f, n0=%f\n",j, cdfpar[j*3+0],cdfpar[j*3+1],cdfpar[j*3+2]);}
	}
	
	sscanf(inp[8+Nvar*(Nvar+1)/2],"%s",outfile);
	if(my_rank==0){
		if (condflag) 
			printf("Output file name: %s3D*_cond_ax*.bin\n",outfile);
		else
			printf("Output file name: %s3D*_uncond_ax*.bin\n",outfile);}
	
	//=============Read logs===========
	if ((stream=fopen(logname,"r"))==NULL)
		if(my_rank==0){printf("The file %s cannot be opened\n",logname);}
	else
		if(my_rank==0){printf("The file %s was opened\n",logname);}	
	int count=0;
	while(fgets(line,200,stream)!=NULL){
		count++;
	}//find out the total lines of logs
	fclose(stream);
	
	int Lenz=count;	
	
	/*---Define X,Y,Z coordinates and Nvar Variables---*/
	vector<double> boreX(Lenz),boreY(Lenz),boreZ(Lenz);
	vector<double> P(Lenz),lP(Lenz),S(Lenz),lS(Lenz);

	if ((stream=fopen(logname,"r"))==NULL)
		if(my_rank==0){printf("The file %s cannot be opened\n",logname);}
	else
		if(my_rank==0){printf("The file %s was opened\n",logname);}	
	count=0;
	while(fgets(line,200,stream)!=NULL){
		sscanf(line,"%s%s%s%s%s%s%s",tempar,tempar+1,tempar+2,tempar+3,tempar+4,tempar+5,tempar+6);
		boreX[count]=atof(tempar[0])-boreX0;boreY[count]=atof(tempar[1])-boreY0;boreZ[count]=atof(tempar[2])-boreZ0;		
		P[count]=atof(tempar[3]);lP[count]=atof(tempar[4]);
		S[count]=atof(tempar[5]);lS[count]=atof(tempar[6]);
		count++;
	}
	fclose(stream);

	if(my_rank==0){
		printf("Depth=[%f %f]\n",GH[0],GH[1]);printf("Hoz X=%f m; Hoz Y=%f m; Dep Z=%f m\n",HOZx,HOZy,VET);
		printf("simds=%f m\n",simds);printf("input file: ' %10s '\n",inpf);
	}

	//double logds=boreZ[1]-boreZ[0];//logs sampling interval;
	vector<double> fp(Lenz),xp(Lenz),fs(Lenz),xs(Lenz);
	vector<double> vp(Lenz),gvp(Lenz),vs(Lenz),gvs(Lenz);
	vector<int> pindex(Lenz),sindex(Lenz);
	//============calculating ecdf =====
	
	for (int i=0;i<Lenz;i++){
		vp[i]=(P[i]-lP[i]);pindex[i]=i;
		vs[i]=(S[i]-lS[i]);sindex[i]=i;
	}//detrend
	ecdf(vp,fp,xp,pindex);
	ecdf(vs,fs,xs,sindex);
	//====pdf trasform=====
	for (int i=0;i<Lenz;i++){
		gvp[pindex[i]]=inormcdf(fp[i],0,1);
		gvs[sindex[i]]=inormcdf(fs[i],0,1);
	}
	/*
	if( (stream  = fopen( "fx.txt", "w" )) == NULL ) 
		printf( "The file fx.txt was not opened\n");
	for (int i=0;i<Lenz;i++){
		fprintf(stream,"%6.5f	%6.5f\n",xp[i],fp[i]);
		}
	fclose(stream);
		if( (stream  = fopen( "gp.txt", "w" )) == NULL ) 
		printf( "The file g.txt was not opened\n");
	for (int i=0;i<Lenz;i++){
		fprintf(stream,"%6.5f\n",gvp[i]);
		}
	fclose(stream);
	cin>>temp;*/
	//====pdf trasform=====
	//Poly_interp myfunR(fr,xr,4),myfunP(fp,xp,4),myfunS(fs,xs,4);
	//Fitab myfitR(boreZ,lR), myfitP(boreZ,lP), myfitS(boreZ,lS);//a linear fit, lR=a+b*boreZ;
	double myfitP_b=(lP[1]-lP[0])/(boreZ[1]-boreZ[0]),  myfitS_b=(lS[1]-lS[0])/(boreZ[1]-boreZ[0]);
	double myfitP_a=-myfitP_b*boreZ[0]+lP[0],myfitS_a=-myfitS_b*boreZ[0]+lS[0];
	/*if (my_rank==0){
		printf("R,(b,a)=(%f,%f)\n",myfitR_b,myfitR_a);
		printf("P,(b,a)=(%f,%f)\n",myfitP_b,myfitP_a);
		printf("S,(b,a)=(%f,%f)\n",myfitS_b,myfitS_a);
	}
	cin>>temp;*/
	//====end of ecdf=====	
	kmax=M_PI/simds;kxs=2*M_PI/HOZx;kys=2*M_PI/HOZy;kzs=2*M_PI/VET;
	
	if(my_rank==0)
	printf("Size of Model: %f Mb (Code may demand: >%5.2f Mb) \n",s,6*s);
		
	int dims2[]={lenk[0],lenk[1]},dims1[]={lenk[2]},ret;
	int Loc_k=lenk[2]/p;
	
	double *kx=new double[lenk[0]],*ky=new double[lenk[1]],*kz=new double[lenk[2]];
		
	for(int i=0;i<lenk[0];i++)
		kx[i]=-kmax+kxs*i;
	for(int i=0;i<lenk[1];i++)	
		ky[i]=-kmax+kys*i;
	for(int i=0;i<lenk[2];i++)	
		kz[i]=-kmax+kzs*i;
		
	for(int i=0;i<lenk[0]/2;i++)
		swap(kx[i],kx[lenk[0]/2+i]);
	for(int i=0;i<lenk[1]/2;i++)	
		swap(ky[i],ky[lenk[1]/2+i]);
	for(int i=0;i<lenk[2]/2;i++)	
		swap(kz[i],kz[lenk[2]/2+i]);		
		
	stream=fopen("kx3.dat","w");
	for(int i=0;i<lenk[0];i++)
		fprintf(stream,"%f\n",kx[i]);
	fclose(stream);
	stream=fopen("ky3.dat","w");
	for(int i=0;i<lenk[1];i++)
		fprintf(stream,"%f\n",ky[i]);
	fclose(stream);
	stream=fopen("kz3.dat","w");
	for(int i=0;i<lenk[2];i++)
		fprintf(stream,"%f\n",kz[i]);
	fclose(stream);
	
	if(my_rank==0){
		printf("Total kx: %d\n",lenk[0]);
		printf("Total ky: %d\n",lenk[1]);
		printf("Total kz: %d\n",lenk[2]);
	}

	double *Zmat_P=new(std::nothrow)double [lenk[0]*lenk[1]*Loc_k];
	double *Zmat_S=new(std::nothrow)double [lenk[0]*lenk[1]*Loc_k];//3D raw data kriging subvolume
	
	for (int axi=0;axi<Nax;axi++){
		for(int ayi=0;ayi<Nay;ayi++){	
			Ran ranP(my_rank+10), ranP2(my_rank+10), ranP3(my_rank+10);//random number seed;
			Ran ranS(my_rank+20), ranS3(my_rank+20);	
	
	int MCNUM=0;
	while (MCNUM<MaxMC){
	MCtime=MPI_Wtime();	
		
// ==========generate P-wave in parallel============
	double *rP=new(std::nothrow)double [lenk[0]*lenk[1]*Loc_k];//3D subvolume
	double *iP=new(std::nothrow)double [lenk[0]*lenk[1]*Loc_k];
	
	if (rP && iP){
		if(my_rank==0) {printf("\n%d Realization of Vp, ax=%6.5f, ay=%6.5f...\n",MCNUM+1,ax[axi],ay[ayi]);Mstart=MPI_Wtime();}
//  0=az vertical,1=ax horizontal 1, 2=ay horizontal 2,3=Hurst, 4=nugget 
	PP[0]=cdfpar[0*3+0];PP[1]=ax[axi];PP[2]=ay[ayi];	
	PP[3]=cdfpar[0*3+1];PP[4]=cdfpar[0*3+2];

	facDes=0.0;
	for(int k=0;k<Loc_k;k++){
		for(int i=0;i<lenk[1];i++){
			for(int j=0;j<lenk[0];j++){
				rNB=ranP.doub();
				rP[k*lenk[0]*lenk[1]+i*lenk[0]+j]=lenk[0]*lenk[1]*lenk[2]*cos(2*M_PI*rNB)*sqrt(pskarman3d(PP,kx[j],ky[i],kz[k+my_rank*Loc_k])/HOZx/HOZy/VET);
				iP[k*lenk[0]*lenk[1]+i*lenk[0]+j]=lenk[0]*lenk[1]*lenk[2]*sin(2*M_PI*rNB)*sqrt(pskarman3d(PP,kx[j],ky[i],kz[k+my_rank*Loc_k])/HOZx/HOZy/VET);
				facDes=facDes+pskarman3d(PP,kx[j],ky[i],kz[k+my_rank*Loc_k])/(1-PP[4]);
			}
		}
		ret = fftn(2, dims2, rP+k*lenk[0]*lenk[1], iP+k*lenk[0]*lenk[1],  1, FORWARD_SCALE);
		if (ret) return 1;
	}
	
	for(int i=0;i<lenk[1];i++){
		if (my_rank ==p-1) temptime=MPI_Wtime();
		for(int j=0;j<lenk[0];j++){
				double *rVec=new(std::nothrow)double [Loc_k];
				double *iVec=new(std::nothrow)double [Loc_k];
				for(int k=0;k<Loc_k;k++){
					rVec[k]=rP[k*lenk[0]*lenk[1]+i*lenk[0]+j];			
					iVec[k]=iP[k*lenk[0]*lenk[1]+i*lenk[0]+j];
				}							
				if (my_rank !=p-1){
					MPI_Send(rVec,Loc_k,MPI_DOUBLE,p-1,0,MPI_COMM_WORLD);
					MPI_Send(iVec,Loc_k,MPI_DOUBLE,p-1,2,MPI_COMM_WORLD);
					MPI_Recv(rVec,Loc_k,MPI_DOUBLE,p-1,1,MPI_COMM_WORLD,&status);					
					}
				else{
					double *Global_rVec=new(std::nothrow)double [lenk[2]];
					double *Global_iVec=new(std::nothrow)double [lenk[2]];
					for (int r=0;r<p-1;r++)
						MPI_Recv(Global_rVec+r*Loc_k,Loc_k,MPI_DOUBLE,r,0,MPI_COMM_WORLD,&status);
					for (int r=0;r<p-1;r++)
						MPI_Recv(Global_iVec+r*Loc_k,Loc_k,MPI_DOUBLE,r,2,MPI_COMM_WORLD,&status);
						
					for (int k=(p-1)*Loc_k;k<lenk[2];k++){
						Global_rVec[k]=rVec[k-(p-1)*Loc_k];
						Global_iVec[k]=iVec[k-(p-1)*Loc_k];
					}
					
					ret = fftn(1, dims1, Global_rVec, Global_iVec,  1, FORWARD_SCALE);
					if (ret) return 1;
					for (int r=0;r<p-1;r++)
						MPI_Send(Global_rVec+r*Loc_k,Loc_k,MPI_DOUBLE,r,1,MPI_COMM_WORLD);

					for (int k=(p-1)*Loc_k;k<lenk[2];k++)
						rVec[k-(p-1)*Loc_k]=Global_rVec[k];
						
					free(Global_rVec);free(Global_iVec);
					}

				for(int k=0;k<Loc_k;k++)
					rP[k*lenk[0]*lenk[1]+i*lenk[0]+j]=rVec[k];
				delete [] rVec;delete [] iVec;					
			}
			//if (my_rank ==p-1) printf("P: fft for traces at %d row is completed, using %5.3f sec\n",i,MPI_Wtime()-temptime);
	}
	delete [] iP;
	sum=0.0;sum2=0.0;
	for(int k=0;k<Loc_k;k++){ 
		for(int i=0;i<lenk[1];i++){
			for(int j=0;j<lenk[0];j++){
				sum=sum+rP[k*lenk[0]*lenk[1]+i*lenk[0]+j];
				sum2=sum2+rP[k*lenk[0]*lenk[1]+i*lenk[0]+j]*rP[k*lenk[0]*lenk[1]+i*lenk[0]+j];
			}
		}	
	} 
	
   	MPI_Allreduce( &sum, &Global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &sum2, &Global_sum2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &facDes, &Global_facDes, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	mymean=Global_sum/lenk[0]/lenk[1]/lenk[2];mystd=sqrt((Global_sum2-Global_sum*Global_sum/lenk[0]/lenk[1]/lenk[2])/(lenk[0]*lenk[1]*lenk[2]-1));
	Global_facDes=sqrt(Global_facDes/HOZx/HOZy/VET);

	for(int k=0;k<Loc_k;k++){
		for(int i=0;i<lenk[1];i++)
			for(int j=0;j<lenk[0];j++)
				rP[k*lenk[0]*lenk[1]+i*lenk[0]+j]=Global_facDes*(rP[k*lenk[0]*lenk[1]+i*lenk[0]+j]-mymean)/mystd;
	}
			
	if(my_rank==0) 
		printf("Time for Nonconditional Simulation of Vp: %6.4f sec\n",(MPI_Wtime()-Mstart));
	}		
	else {printf("Not enough processors!\n");return 1;}
	
	// ==========generate S-wave in parallel============
	double SS1,SS2;
	double *rS=new(std::nothrow)double [lenk[0]*lenk[1]*Loc_k];//3D subvolume
	double *iS=new(std::nothrow)double [lenk[0]*lenk[1]*Loc_k];
	
	if (rS && iS) {
		if(my_rank==0) {printf("\n%d Realization of Vs, ax=%6.5f, ay=%6.5f...\n",MCNUM+1,ax[axi],ay[ayi]);Mstart=MPI_Wtime();}
//  0=az vertical,1=ax horizontal 1, 2=ay horizontal 2,3=Hurst, 4=nugget 
	PP[0]=cdfpar[0*3+0];PP[1]=ax[axi];PP[2]=ay[ayi];	
	PP[3]=cdfpar[0*3+1];PP[4]=cdfpar[0*3+2];//PP
	
	SP[0]=cdfpar[1*3+0];SP[1]=ax[axi];SP[2]=ay[ayi];
	SP[3]=cdfpar[1*3+1];SP[4]=cdfpar[1*3+2];//SP
	
	SS[0]=cdfpar[2*3+0];SS[1]=ax[axi];SS[2]=ay[ayi];
	SS[3]=cdfpar[2*3+1];SS[4]=cdfpar[2*3+2];//SS

	facDes=0.0;
	for(int k=0;k<Loc_k;k++){
		for(int i=0;i<lenk[1];i++){
			for(int j=0;j<lenk[0];j++){						
				rNB2=ranP2.doub();rNB=ranS.doub();
				SS1=pskarman3d(SP,kx[j],ky[i],kz[k+my_rank*Loc_k])*pskarman3d(SP,kx[j],ky[i],kz[k+my_rank*Loc_k])/pskarman3d(PP,kx[j],ky[i],kz[k+my_rank*Loc_k])/HOZx/HOZy/VET;
				SS2=(pskarman3d(SS,kx[j],ky[i],kz[k+my_rank*Loc_k])-pskarman3d(SP,kx[j],ky[i],kz[k+my_rank*Loc_k])*pskarman3d(SP,kx[j],ky[i],kz[k+my_rank*Loc_k])/pskarman3d(PP,kx[j],ky[i],kz[k+my_rank*Loc_k]))/HOZx/HOZy/VET;				
				if(SS2>0){
					rS[k*lenk[0]*lenk[1]+i*lenk[0]+j]=lenk[0]*lenk[1]*lenk[2]*(cos(2*M_PI*rNB2)*sqrt(SS1)+cos(2*M_PI*rNB)*sqrt(SS2));
					iS[k*lenk[0]*lenk[1]+i*lenk[0]+j]=lenk[0]*lenk[1]*lenk[2]*(sin(2*M_PI*rNB2)*sqrt(SS1)+sin(2*M_PI*rNB)*sqrt(SS2));
				}
				else{
					rS[k*lenk[0]*lenk[1]+i*lenk[0]+j]=lenk[0]*lenk[1]*lenk[2]*(cos(2*M_PI*rNB2)*sqrt(SS1)+cos(2*M_PI*rNB+M_PI/2)*sqrt(-SS2));
					iS[k*lenk[0]*lenk[1]+i*lenk[0]+j]=lenk[0]*lenk[1]*lenk[2]*(sin(2*M_PI*rNB2)*sqrt(SS1)+sin(2*M_PI*rNB+M_PI/2)*sqrt(-SS2));
				}
				facDes=facDes+pskarman3d(SS,kx[j],ky[i],kz[k+my_rank*Loc_k])/(1-SS[4]);
			}
		}
		ret = fftn(2, dims2, rS+k*lenk[0]*lenk[1], iS+k*lenk[0]*lenk[1],  1, FORWARD_SCALE);
		if (ret) return 1;
	}
	
	for(int i=0;i<lenk[1];i++){
		if (my_rank ==p-1) temptime=MPI_Wtime();
		for(int j=0;j<lenk[0];j++){
				double *rVec=new(std::nothrow)double [Loc_k];
				double *iVec=new(std::nothrow)double [Loc_k];
				for(int k=0;k<Loc_k;k++){
					rVec[k]=rS[k*lenk[0]*lenk[1]+i*lenk[0]+j];			
					iVec[k]=iS[k*lenk[0]*lenk[1]+i*lenk[0]+j];
				}							
				if (my_rank !=p-1){
					MPI_Send(rVec,Loc_k,MPI_DOUBLE,p-1,0,MPI_COMM_WORLD);
					MPI_Send(iVec,Loc_k,MPI_DOUBLE,p-1,2,MPI_COMM_WORLD);
					MPI_Recv(rVec,Loc_k,MPI_DOUBLE,p-1,1,MPI_COMM_WORLD,&status);					
					}
				else{
					double *Global_rVec=new(std::nothrow)double [lenk[2]];
					double *Global_iVec=new(std::nothrow)double [lenk[2]];
					for (int r=0;r<p-1;r++)
						MPI_Recv(Global_rVec+r*Loc_k,Loc_k,MPI_DOUBLE,r,0,MPI_COMM_WORLD,&status);
					for (int r=0;r<p-1;r++)
						MPI_Recv(Global_iVec+r*Loc_k,Loc_k,MPI_DOUBLE,r,2,MPI_COMM_WORLD,&status);
						
					for (int k=(p-1)*Loc_k;k<lenk[2];k++){
						Global_rVec[k]=rVec[k-(p-1)*Loc_k];
						Global_iVec[k]=iVec[k-(p-1)*Loc_k];
					}
					
					ret = fftn(1, dims1, Global_rVec, Global_iVec,  1, FORWARD_SCALE);
					if (ret) return 1;
					for (int r=0;r<p-1;r++)
						MPI_Send(Global_rVec+r*Loc_k,Loc_k,MPI_DOUBLE,r,1,MPI_COMM_WORLD);

					for (int k=(p-1)*Loc_k;k<lenk[2];k++)
						rVec[k-(p-1)*Loc_k]=Global_rVec[k];
						
					free(Global_rVec);free(Global_iVec);
					}

				for(int k=0;k<Loc_k;k++)
					rS[k*lenk[0]*lenk[1]+i*lenk[0]+j]=rVec[k];
				delete [] rVec;delete [] iVec;					
			}
			//if (my_rank ==p-1) printf("R: fft for traces at %d row is completed, using %5.3f sec\n",i,MPI_Wtime()-temptime);
	}
	delete [] iS;
	sum=0.0;sum2=0.0;
	for(int k=0;k<Loc_k;k++){ 
		for(int i=0;i<lenk[1];i++){
			for(int j=0;j<lenk[0];j++){
				sum=sum+rS[k*lenk[0]*lenk[1]+i*lenk[0]+j];
				sum2=sum2+rS[k*lenk[0]*lenk[1]+i*lenk[0]+j]*rS[k*lenk[0]*lenk[1]+i*lenk[0]+j];
			}
		}	
	} 
	
   	MPI_Allreduce( &sum, &Global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &sum2, &Global_sum2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &facDes, &Global_facDes, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	mymean=Global_sum/lenk[0]/lenk[1]/lenk[2];mystd=sqrt((Global_sum2-Global_sum*Global_sum/lenk[0]/lenk[1]/lenk[2])/(lenk[0]*lenk[1]*lenk[2]-1));
	Global_facDes=sqrt(Global_facDes/HOZx/HOZy/VET);

	for(int k=0;k<Loc_k;k++){
		for(int i=0;i<lenk[1];i++)
			for(int j=0;j<lenk[0];j++)
				rS[k*lenk[0]*lenk[1]+i*lenk[0]+j]=Global_facDes*(rS[k*lenk[0]*lenk[1]+i*lenk[0]+j]-mymean)/mystd;
	}
	if(my_rank==0) printf("Time for Nonconditional Simulation of Vs: %6.4f sec\n",(MPI_Wtime()-Mstart));
	}		
	else {printf("Not enough processors!\n");return 1;}
	//cout<<mymean<<' '<<mystd<<' '<<Global_facDes<<endl;
	/*sprintf(rbinfile,"gR.bin.%d",my_rank);
	if( (rstream  = fopen( rbinfile, "wb" )) == NULL ) 
		printf( "The file %s was not opened\n", rbinfile );
	fwrite(rS,sizeof(double),lenk[0]*lenk[1]*Loc_k,rstream);fclose(rstream);	
	MPI_Barrier( MPI_COMM_WORLD );
	if(my_rank==0) {
	for (int i=0;i<p;i++){
		sprintf(cmd,"cat gR.bin.%d >> gNCR_ax%d_ay%d.bin ",i,int(ax[axi]),int(ay[ayi]));system(cmd);	
		sprintf(cmd,"rm gR.bin.%d",i);system(cmd);	
	}
	}
cout<<"stop"<<endl;cin>>temp;*/		
	
	//======Start kriging part for three variables====
				
	if (condflag){
		if(MCNUM==0){			
			vector<double> ZCmat_P(boreX.size()*boreX.size()),ZCmat_S(boreX.size()*boreX.size());
			vector<double> ZLam_P(boreX.size()),ZLam_S(boreX.size());
			if (my_rank==0) {printf("Data Kriging...\n");printf("Size of Data: %d X 1 Column\n",boreX.size());}
			temptime=MPI_Wtime();
			for (int i=0;i<boreX.size();i++){
				for (int j=0;j<boreX.size();j++){
					ZCmat_P[i*boreX.size()+j]=karman3d(PP,boreX[i]-boreX[j],boreY[i]-boreY[j],boreZ[i]-boreZ[j]);
					ZCmat_S[i*boreX.size()+j]=karman3d(SS,boreX[i]-boreX[j],boreY[i]-boreY[j],boreZ[i]-boreZ[j]);
				}
			}

			LinearSolver(ZCmat_P,gvp,ZLam_P);LinearSolver(ZCmat_S,gvs,ZLam_S);
			//printf("PE%d: size of gvs=%d, size of Lam= %d\n",my_rank,gvs.size(),Lam.size());
			for(int k=0;k<Loc_k;k++){					
				for(int i=0;i<lenk[1];i++){				
					for(int j=0;j<lenk[0];j++){	
					
						S1=0.0;
						for (int ikar=0;ikar<boreX.size();ikar++){
							S1=S1+ZLam_P[ikar]*karman3d(PP,j*simds-boreX[ikar],i*simds-boreY[ikar],(k+Loc_k*my_rank)*simds+GH[0]-boreZ[ikar]);
						}
						Zmat_P[k*lenk[0]*lenk[1]+i*lenk[0]+j]=S1;
					
						S1=0.0;
						for (int ikar=0;ikar<boreX.size();ikar++){
							S1=S1+ZLam_S[ikar]*karman3d(SS,j*simds-boreX[ikar],i*simds-boreY[ikar],(k+Loc_k*my_rank)*simds+GH[0]-boreZ[ikar]);
						}
						Zmat_S[k*lenk[0]*lenk[1]+i*lenk[0]+j]=S1;					
					}				
				}			
			}
			if (my_rank==0) printf("Data Kriging completed; takes time %f sec\n",MPI_Wtime()-temptime);	
		}
		else if (my_rank==0)
			printf("Not the first realization, recycle the Data kriging model\n");
	sprintf(rbinfile,"Zmat_P.bin.%d",my_rank);
	if( (rstream  = fopen( rbinfile, "wb" )) == NULL ) 
		printf( "The file %s was not opened\n", rbinfile );
	fwrite(Zmat_P,sizeof(double),lenk[0]*lenk[1]*Loc_k,rstream);fclose(rstream);	
	MPI_Barrier( MPI_COMM_WORLD );
	if(my_rank==0) {
	for (int i=0;i<p;i++){
		sprintf(cmd,"cat Zmat_P.bin.%d >> Zmat_P_ax%d_ay%d.bin ",i,int(ax[axi]),int(ay[ayi]));system(cmd);	
		sprintf(cmd,"rm Zmat_P.bin.%d",i);system(cmd);	
	}
	}	
	//MPI_Barrier( MPI_COMM_WORLD );
		
		if (my_rank==0) printf("Model Kriging...\n");		
		temptime=MPI_Wtime();
		double *Smat_P=new(std::nothrow)double [lenk[0]*lenk[1]*Loc_k];
		double *Smat_S=new(std::nothrow)double [lenk[0]*lenk[1]*Loc_k];//3D model data kriging subvolume
		vector<double> modX,modY,modZ;
		bore2mod(boreX,boreY,boreZ,modX,modY,modZ,simds);//assign likely irregularly sampled raw data to regular grid of the model
		if (my_rank==0) printf("Size of Model: %d X 1 Column\n",modX.size());
		
		vector<double> SCmat_P(modX.size()*modX.size()),SCmat_S(modX.size()*modX.size());
		vector<double>	SLam_P(modX.size()), Globalgvp(modX.size());
		vector<double> SLam_S(modX.size()), Globalgvs(modX.size());
		
		//printf("PE%d: size of Mat=%d, size of Lam= %d, size of Globalgvs=%d\n",my_rank,Mat.nrows(),Lam.size(),Globalgvs.size());
		int k,rankID;
		for (int i=0;i<modX.size();i++){
			k=int((modZ[i]-GH[0])/simds) % Loc_k;
			rankID=(int((modZ[i]-GH[0])/simds)-k)/Loc_k;
			//printf("k=%d,rankID=%d\n",k,rankID);
			if (my_rank==rankID) {
				Globalgvp[i]=rP[k*lenk[0]*lenk[1]+int(modY[i]/simds)*lenk[0]+int(modX[i]/simds)];
				Globalgvs[i]=rS[k*lenk[0]*lenk[1]+int(modY[i]/simds)*lenk[0]+int(modX[i]/simds)];
			}
				MPI_Bcast(&Globalgvp[i], 1, MPI_DOUBLE, rankID, MPI_COMM_WORLD);
				MPI_Bcast(&Globalgvs[i], 1, MPI_DOUBLE, rankID, MPI_COMM_WORLD);
			
		}//all process share the same model data
		
		if (my_rank==0) {
			if( (stream  = fopen( "modXYZ.txt", "w" )) == NULL ) 
				printf( "The file fx.txt was not opened\n");
			fprintf(stream,"X	Y	Z	grho	gvp		gvs\n");
			for (int i=0;i<modX.size();i++){
				fprintf(stream,"%6.5f	%6.5f	%6.5f	%6.5f	%6.5f\n",modX[i],modY[i],modZ[i],Globalgvp[i],Globalgvs[i]);
			}
			fclose(stream);
		}
		
		for (int i=0;i<modX.size();i++){
			for (int j=0;j<modX.size();j++){
				SCmat_P[i*modX.size()+j]=karman3d(PP,modX[i]-modX[j],modY[i]-modY[j],modZ[i]-modZ[j]);
				SCmat_S[i*modX.size()+j]=karman3d(SS,modX[i]-modX[j],modY[i]-modY[j],modZ[i]-modZ[j]);
			}
		}

		LinearSolver(SCmat_P,Globalgvp,SLam_P);LinearSolver(SCmat_S,Globalgvs,SLam_S);
		for(int k=0;k<Loc_k;k++){		
			for(int i=0;i<lenk[1];i++){
				for(int j=0;j<lenk[0];j++){				
					
					S1=0.0;
					for (int ikar=0;ikar<modX.size();ikar++){
						S1=S1+SLam_P[ikar]*karman3d(PP,j*simds-modX[ikar],i*simds-modY[ikar],(k+Loc_k*my_rank)*simds+GH[0]-modZ[ikar]);
					}
					Smat_P[k*lenk[0]*lenk[1]+i*lenk[0]+j]=S1;
					
					S1=0.0;
					for (int ikar=0;ikar<modX.size();ikar++){
						S1=S1+SLam_S[ikar]*karman3d(SS,j*simds-modX[ikar],i*simds-modY[ikar],(k+Loc_k*my_rank)*simds+GH[0]-modZ[ikar]);
					}
					Smat_S[k*lenk[0]*lenk[1]+i*lenk[0]+j]=S1;						
				}
			}			
		}
		if (my_rank==0) printf("Model Kriging completed; takes time %f sec\n",MPI_Wtime()-temptime);
	sprintf(rbinfile,"Smat_P.bin.%d",my_rank);
	if( (rstream  = fopen( rbinfile, "wb" )) == NULL ) 
		printf( "The file %s was not opened\n", rbinfile );
	fwrite(Smat_P,sizeof(double),lenk[0]*lenk[1]*Loc_k,rstream);fclose(rstream);	
	MPI_Barrier( MPI_COMM_WORLD );
	if(my_rank==0) {
	for (int i=0;i<p;i++){
		sprintf(cmd,"cat Smat_P.bin.%d >> Smat_P_ax%d_ay%d.bin ",i,int(ax[axi]),int(ay[ayi]));system(cmd);	
		sprintf(cmd,"rm Smat_P.bin.%d",i);system(cmd);	
	}
	}
		for(int k=0;k<Loc_k;k++){
			for(int i=0;i<lenk[1];i++)
				for(int j=0;j<lenk[0];j++){	
					rP[k*lenk[0]*lenk[1]+i*lenk[0]+j]=myInterp(fp,xp,normcdf(Zmat_P[k*lenk[0]*lenk[1]+i*lenk[0]+j]+rP[k*lenk[0]*lenk[1]+i*lenk[0]+j]-Smat_P[k*lenk[0]*lenk[1]+i*lenk[0]+j],0,1))+\
						myfitP_a+myfitP_b*((k+Loc_k*my_rank)*simds+GH[0]);
					rS[k*lenk[0]*lenk[1]+i*lenk[0]+j]=myInterp(fs,xs,normcdf(Zmat_S[k*lenk[0]*lenk[1]+i*lenk[0]+j]+rS[k*lenk[0]*lenk[1]+i*lenk[0]+j]-Smat_S[k*lenk[0]*lenk[1]+i*lenk[0]+j],0,1))+\
						myfitS_a+myfitS_b*((k+Loc_k*my_rank)*simds+GH[0]);
				}
		}
		delete [] Smat_P; delete [] Smat_S; 
	}
	else{
	//Nonconditional
		for(int k=0;k<Loc_k;k++){ 
		temptime=MPI_Wtime();		
		for(int i=0;i<lenk[1];i++)
			for(int j=0;j<lenk[0];j++){
				rP[k*lenk[0]*lenk[1]+i*lenk[0]+j]=myInterp(fp,xp,normcdf(rP[k*lenk[0]*lenk[1]+i*lenk[0]+j],0,1))+\
					myfitP_a+myfitP_b*((k+Loc_k*my_rank)*simds+GH[0]);		
				rS[k*lenk[0]*lenk[1]+i*lenk[0]+j]=myInterp(fs,xs,normcdf(rS[k*lenk[0]*lenk[1]+i*lenk[0]+j],0,1))+\
					myfitS_a+myfitS_b*((k+Loc_k*my_rank)*simds+GH[0]);
				
			}
		}
	}
	// =====kriging end====
	
	if (MaxMC==1 || ax[axi]==ay[ayi]){

		//===== Output Vp============	
		sprintf(rbinfile,"rP.bin.%d",my_rank);
		if( (rstream  = fopen( rbinfile, "wb" )) == NULL ) 
			printf( "The file %s was not opened\n", rbinfile );
		fwrite(rP,sizeof(double),lenk[0]*lenk[1]*Loc_k,rstream);fclose(rstream);			
		
		//===== Output Vs============	
		sprintf(rbinfile,"rS.bin.%d",my_rank);
		if( (rstream  = fopen( rbinfile, "wb" )) == NULL ) 
			printf( "The file %s was not opened\n", rbinfile );
		fwrite(rS,sizeof(double),lenk[0]*lenk[1]*Loc_k,rstream);fclose(rstream);	
		MPI_Barrier( MPI_COMM_WORLD );
		
		//=====Merge all files ==========
		if(my_rank==0) {
			if (condflag){
				for (int i=0;i<p;i++){
					sprintf(cmd,"cat rP.bin.%d >> %s3DVp_cond_ax%d_ay%d.bin ",i,outfile,int(ax[axi]),int(ay[ayi]));system(cmd);	
					sprintf(cmd,"rm rP.bin.%d",i);system(cmd);	
					sprintf(cmd,"cat rS.bin.%d >> %s3DVs_cond_ax%d_ay%d.bin ",i,outfile,int(ax[axi]),int(ay[ayi]));system(cmd);	
					sprintf(cmd,"rm rS.bin.%d",i);system(cmd);
				}			
			}
			else {
				for (int i=0;i<p;i++){
					sprintf(cmd,"cat rP.bin.%d >> %s3DVp_uncond_ax%d_ay%d.bin ",i,outfile,int(ax[axi]),int(ay[ayi]));system(cmd);	
					sprintf(cmd,"rm rP.bin.%d",i);system(cmd);	
					sprintf(cmd,"cat rS.bin.%d >> %s3DVs_uncond_ax%d_ay%d.bin ",i,outfile,int(ax[axi]),int(ay[ayi]));system(cmd);	
					sprintf(cmd,"rm rS.bin.%d",i);system(cmd);	
				}
			}
			sprintf(cmd,"mv *.bin ../models/");system(cmd);
		}
				
	}
	
	MCNUM=MCNUM+1;
	if(my_rank==0) printf("Monto Carlo Realization %d, using %5.3f sec \n",MCNUM,MPI_Wtime()-MCtime);
	
	}//MCNUM
	 
	MPI_Barrier( MPI_COMM_WORLD );	
	} //ay
	} //ax
	if(my_rank==0) {printf("Total Elapse Time: %6.4f hours \n", (MPI_Wtime()-Telapse)/3600);}
	MPI_Finalize();
}	
