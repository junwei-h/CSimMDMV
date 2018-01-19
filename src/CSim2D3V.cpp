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
give a reference to the following paper:
	Huang, J.W., Gilles Bellefleur and Bernd Milkereit, Seismic Modeling of Multi-dimensional Heterogeneity Scales of Mallik Gas Hydrates Reservoirs, Northwest Territories of Canada, Journal of Geophysical Research-Solid Earth, (submitted). 

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
	int lenk[3],M=0,Nax,Nay,Nvar,condflag,RockPhyTag=0, my_rank,p,MaxMC;
	double HOZx,HOZy,simds,VET,kmax,kxs,kys,kzs,rNB,rNB2,rNB3,sum=0.0,sum2=0.0,facDes=0.0,Telapse,Mstart,temptime,MCtime;
	double shiftR,stdR,mymean,mystd;
	char *inpf,*outf, line[200],inp[90][80],logname[10], tempar[8][80];
	char rbinfile[100],ibinfile[100],cmd[100],outfile[100];
	//double *boreX=new double[MaxA],*boreY=new double[MaxA],*boreZ=new double[MaxA], *R=new double[MaxA],*lR=new double[MaxA];
	double par[8],*cdfpar;
	double Global_facDes=0.0,Global_sum=0.0,Global_sum2=0.0; 
	double boreX0,boreY0,boreZ0,GHCon=0.0, Poro=0.0,GH[2],GHVol;
	
	float temp;
	
	double *PP=new double[5];
	double *SS=new double[5],*SP=new double[5];
	double *RP=new double[5],*RS=new double[5],*RR=new double[5];
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
			strcpy(inp[M],line);
			/*if (my_rank==0){
				printf("inp[%d]=%s",M,inp[M]);
				printf("line=%s\n",line);
				}*/
			M=M+1;
		}
	}
	fclose(stream);
	
	//==========Read input parameters==============	
	sscanf(inp[0],"%s%s%s%s%s",tempar,tempar+1,tempar+2,tempar+3,tempar+4);
	HOZx=atof(tempar[0]);
	HOZy=atof(tempar[1]);//HOZy and lenk[1] are ignored in this 2D program!
	GH[0]=atof(tempar[2]);GH[1]=atof(tempar[3]);simds=atof(tempar[4]);
	VET=GH[1]-GH[0]+simds;
	lenk[0]=(int)round(pow(simds,-1)*HOZx);lenk[1]=(int)round(pow(simds,-1)*HOZy);lenk[2]=(int)round(pow(simds,-1)*VET);
	double s=lenk[0]*lenk[2]*sizeof(double)/1024.0/1024.0;

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
		
	/*sscanf(inp[5],"%s%s%s",tempar,tempar+1,tempar+2);
	Nay=int((atof(tempar[2])-atof(tempar[0]))/atof(tempar[1]))+1;
	if (my_rank==0) printf("ay=%f:%f:%f, total: %d\n",atof(tempar[0]),atof(tempar[1]),atof(tempar[2]),Nay);
	double *ay=new double[Nay];
	ay[0]=atof(tempar[0]);
	for (int i=1;i<Nay;i++)
		ay[i]=ay[0]+i*atof(tempar[1]);	*/
	if (my_rank==0) printf("total iteration: %d\n",Nax*1);
	
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
			printf("Output file name: %s2D*_cond_ax*.bin\n",outfile);
		else
			printf("Output file name: %s2D*_uncond_ax*.bin\n",outfile);
	}
	sscanf(inp[9+Nvar*(Nvar+1)/2],"%s",tempar);RockPhyTag=atoi(tempar[0]);
	if (my_rank==0 && RockPhyTag) {
		printf("RockPhyTag=%d, tempar=%s\n",RockPhyTag,tempar[0]);
		printf("Use Modified Biot-Gassmann Rock Physics Theory! Check if you need it....\n");
	}
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
	
	/*---Define X,Y,Z coordinates and 3 Variables---*/
	vector<double> boreX(Lenz),boreY(Lenz),boreZ(Lenz),boreR(Lenz);
	vector<double> R(Lenz),lR(Lenz),P(Lenz),lP(Lenz),S(Lenz),lS(Lenz);
	
	if ((stream=fopen(logname,"r"))==NULL)
		if(my_rank==0){printf("The file %s cannot be opened\n",logname);}
	else
		if(my_rank==0){printf("The file %s was opened\n",logname);}	
	count=0;
	while(fgets(line,200,stream)!=NULL){
		sscanf(line,"%s%s%s%s%s%s%s%s%s",tempar,tempar+1,tempar+2,tempar+3,tempar+4,tempar+5,tempar+6,tempar+7,tempar+8);
		boreX[count]=atof(tempar[0])-boreX0;boreY[count]=atof(tempar[1])-boreY0;boreZ[count]=atof(tempar[2])-boreZ0;		
		P[count]=atof(tempar[3]);lP[count]=atof(tempar[4]);
		S[count]=atof(tempar[5]);lS[count]=atof(tempar[6]);
		R[count]=atof(tempar[7]);lR[count]=atof(tempar[8]);
		count++;
	}
	fclose(stream);

	if(my_rank==0){
		printf("Depth=[%f %f]\n",GH[0],GH[1]);printf("Hoz X=%f m; Dep Z=%f m\n",HOZx,VET);
		printf("simds=%f m\n",simds);printf("input file: ' %10s '\n",inpf);
	}

	//double logds=boreZ[1]-boreZ[0];//logs sampling interval;
	vector<double> fr(Lenz),xr(Lenz),fp(Lenz),xp(Lenz),fs(Lenz),xs(Lenz);
	vector<double> rho(Lenz),grho(Lenz),vp(Lenz),gvp(Lenz),vs(Lenz),gvs(Lenz);
	vector<int> rindex(Lenz),pindex(Lenz),sindex(Lenz);
	//============calculating ecdf =====
	
	for (int i=0;i<Lenz;i++){
		rho[i]=(R[i]-lR[i]);rindex[i]=i;
		vp[i]=(P[i]-lP[i]);pindex[i]=i;
		vs[i]=(S[i]-lS[i]);sindex[i]=i;
	}//detrend
	ecdf(rho,fr,xr,rindex);
	ecdf(vp,fp,xp,pindex);
	ecdf(vs,fs,xs,sindex);
	//====pdf trasform=====
	for (int i=0;i<Lenz;i++){
		grho[rindex[i]]=inormcdf(fr[i],0,1);
		gvp[pindex[i]]=inormcdf(fp[i],0,1);
		gvs[sindex[i]]=inormcdf(fs[i],0,1);
	}
	/*if( (stream  = fopen( "fx.txt", "w" )) == NULL ) 
		printf( "The file fx.txt was not opened\n");
	for (int i=0;i<Lenz;i++){
		fprintf(stream,"%6.5f	%6.5f\n",xs[i],fs[i]);
		}
	fclose(stream);*/
		if( (stream  = fopen( "g.txt", "w" )) == NULL ) 
		printf( "The file g.txt was not opened\n");
	for (int i=0;i<Lenz;i++){
		fprintf(stream,"%6.5f\n",grho[i]);
		}
	fclose(stream);
	//cin>>temp;
	//====pdf trasform=====
	//a linear fit, lR=a+b*boreZ;
	double myfitR_b=(lR[1]-lR[0])/(boreZ[1]-boreZ[0]),  myfitP_b=(lP[1]-lP[0])/(boreZ[1]-boreZ[0]),  myfitS_b=(lS[1]-lS[0])/(boreZ[1]-boreZ[0]);
	double myfitR_a=-myfitR_b*boreZ[0]+lR[0], myfitP_a=-myfitP_b*boreZ[0]+lP[0],myfitS_a=-myfitS_b*boreZ[0]+lS[0];
	//====end of ecdf=====	
	kmax=M_PI/simds;kxs=2*M_PI/HOZx;kys=2*M_PI/HOZy;kzs=2*M_PI/VET;
	
	if(my_rank==0)
	printf("Size of Model: %f Mb (Code may demand: >%5.2f Mb) \n",s,12*s);
		
	int dims2[]={lenk[0]},dims1[]={lenk[2]},ret;
	int Loc_k=lenk[2]/p;
	
	double *kx=new double[lenk[0]],*kz=new double[lenk[2]];
		
	for(int i=0;i<lenk[0];i++)
		kx[i]=-kmax+kxs*i;
	for(int i=0;i<lenk[2];i++)	
		kz[i]=-kmax+kzs*i;
		
	for(int i=0;i<lenk[0]/2;i++)
		swap(kx[i],kx[lenk[0]/2+i]);
	for(int i=0;i<lenk[2]/2;i++)	
		swap(kz[i],kz[lenk[2]/2+i]);		
		
	stream=fopen("kx2.dat","w");
	for(int i=0;i<lenk[0];i++)
		fprintf(stream,"%f\n",kx[i]);
	fclose(stream);
	stream=fopen("kz2.dat","w");
	for(int i=0;i<lenk[2];i++)
		fprintf(stream,"%f\n",kz[i]);
	fclose(stream);
	
	if(my_rank==0){
		printf("Total kx: %d\n",lenk[0]);
		printf("Total kz: %d\n",lenk[2]);
	}

	double *Zmat_R=new(std::nothrow)double [lenk[0]*Loc_k];
	double *Zmat_P=new(std::nothrow)double [lenk[0]*Loc_k];
	double *Zmat_S=new(std::nothrow)double [lenk[0]*Loc_k];//2D raw data kriging subvolume
	
	for (int axi=0;axi<Nax;axi++){
			Ran ranP(my_rank+10), ranP2(my_rank+10), ranP3(my_rank+10);//random number seed;
			Ran ranS(my_rank+20), ranS3(my_rank+20);	
			Ran ranR(my_rank+1);
	
	int MCNUM=0;
	while (MCNUM<MaxMC){
	MCtime=MPI_Wtime();	
	
// ==========generate P-wave in parallel============
	double *rP=new(std::nothrow)double [lenk[0]*Loc_k];//2D subvolume
	double *iP=new(std::nothrow)double [lenk[0]*Loc_k];
	
	if (rP && iP){
		if(my_rank==0) {printf("\n%d Realization of Vp, ax=%6.5f,...\n",MCNUM+1,ax[axi]);Mstart=MPI_Wtime();}
//  0=az vertical,1=ax horizontal 1, 2=ay horizontal 2,3=Hurst, 4=nugget 
	PP[0]=cdfpar[0*3+0];PP[1]=ax[axi];PP[2]=PP[1];	
	PP[3]=cdfpar[0*3+1];PP[4]=cdfpar[0*3+2];

	facDes=0.0;
	for(int k=0;k<Loc_k;k++){
		for(int j=0;j<lenk[0];j++){
			rNB=ranP.doub();
			rP[k*lenk[0]+j]=lenk[0]*lenk[2]*cos(2*M_PI*rNB)*sqrt(pskarman2d(PP,kx[j],kz[k+my_rank*Loc_k])/HOZx/VET);
			iP[k*lenk[0]+j]=lenk[0]*lenk[2]*sin(2*M_PI*rNB)*sqrt(pskarman2d(PP,kx[j],kz[k+my_rank*Loc_k])/HOZx/VET);
			facDes=facDes+pskarman2d(PP,kx[j],kz[k+my_rank*Loc_k])/(1-PP[4]);
		}		
		ret = fftn(1, dims2, rP+k*lenk[0], iP+k*lenk[0],  1, FORWARD_SCALE);
		if (ret) return 1;
	}
	
		for(int j=0;j<lenk[0];j++){
			if (my_rank ==p-1) temptime=MPI_Wtime();
			double *rVec=new(std::nothrow)double [Loc_k];
			double *iVec=new(std::nothrow)double [Loc_k];
			for(int k=0;k<Loc_k;k++){
				rVec[k]=rP[k*lenk[0]+j];			
				iVec[k]=iP[k*lenk[0]+j];
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
				rP[k*lenk[0]+j]=rVec[k];
			delete [] rVec;delete [] iVec;					
		}
			//if (my_rank ==p-1) printf("P: fft for traces at %d row is completed, using %5.3f sec\n",i,MPI_Wtime()-temptime);
	delete [] iP;
	sum=0.0;sum2=0.0;
	for(int k=0;k<Loc_k;k++){ 
		for(int j=0;j<lenk[0];j++){
			sum=sum+rP[k*lenk[0]+j];
			sum2=sum2+rP[k*lenk[0]+j]*rP[k*lenk[0]+j];
		}
	} 
	
   	MPI_Allreduce( &sum, &Global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &sum2, &Global_sum2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &facDes, &Global_facDes, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	mymean=Global_sum/lenk[0]/lenk[2];mystd=sqrt((Global_sum2-Global_sum*Global_sum/lenk[0]/lenk[2])/(lenk[0]*lenk[2]-1));
	Global_facDes=sqrt(Global_facDes/HOZx/VET);

	for(int k=0;k<Loc_k;k++){
			for(int j=0;j<lenk[0];j++)
				rP[k*lenk[0]+j]=Global_facDes*(rP[k*lenk[0]+j]-mymean)/mystd;
	}
			
	if(my_rank==0) 
		printf("Time for Nonconditional Simulation of Vp: %6.4f sec\n",(MPI_Wtime()-Mstart));
	}		
	else {printf("Not enough processors!\n");return 1;}
	
	// ==========generate S-wave in parallel============
	double rNB2,SS1,SS2;
	double *rS=new(std::nothrow)double [lenk[0]*Loc_k];//2D subvolume
	double *iS=new(std::nothrow)double [lenk[0]*Loc_k];
	
	if (rS && iS) {
		if(my_rank==0) {printf("\n%d Realization of Vs, ax=%6.5f,...\n",MCNUM+1,ax[axi]);Mstart=MPI_Wtime();}
//  0=az vertical,1=ax horizontal 1, 2=ay horizontal 2,3=Hurst, 4=nugget 
	PP[0]=cdfpar[0*3+0];PP[1]=ax[axi];PP[2]=PP[1];	
	PP[3]=cdfpar[0*3+1];PP[4]=cdfpar[0*3+2];//PP
	
	SP[0]=cdfpar[1*3+0];SP[1]=ax[axi];SP[2]=SP[1];
	SP[3]=cdfpar[1*3+1];SP[4]=cdfpar[1*3+2];//SP
	
	SS[0]=cdfpar[2*3+0];SS[1]=ax[axi];SS[2]=SS[1];
	SS[3]=cdfpar[2*3+1];SS[4]=cdfpar[2*3+2];//SS

	facDes=0.0;
	for(int k=0;k<Loc_k;k++){
		for(int j=0;j<lenk[0];j++){						
			rNB2=ranP2.doub();rNB=ranS.doub();
			SS1=pskarman2d(SP,kx[j],kz[k+my_rank*Loc_k])*pskarman2d(SP,kx[j],kz[k+my_rank*Loc_k])/pskarman2d(PP,kx[j],kz[k+my_rank*Loc_k])/HOZx/VET;
			SS2=(pskarman2d(SS,kx[j],kz[k+my_rank*Loc_k])-pskarman2d(SP,kx[j],kz[k+my_rank*Loc_k])*pskarman2d(SP,kx[j],kz[k+my_rank*Loc_k])/pskarman2d(PP,kx[j],kz[k+my_rank*Loc_k]))/HOZx/VET;				
			if(SS2>0){
				rS[k*lenk[0]+j]=lenk[0]*lenk[2]*(cos(2*M_PI*rNB2)*sqrt(SS1)+cos(2*M_PI*rNB)*sqrt(SS2));
				iS[k*lenk[0]+j]=lenk[0]*lenk[2]*(sin(2*M_PI*rNB2)*sqrt(SS1)+sin(2*M_PI*rNB)*sqrt(SS2));
			}
			else{
				rS[k*lenk[0]+j]=lenk[0]*lenk[2]*(cos(2*M_PI*rNB2)*sqrt(SS1)+cos(2*M_PI*rNB+M_PI/2)*sqrt(-SS2));
				iS[k*lenk[0]+j]=lenk[0]*lenk[2]*(sin(2*M_PI*rNB2)*sqrt(SS1)+sin(2*M_PI*rNB+M_PI/2)*sqrt(-SS2));
			}
			facDes=facDes+pskarman2d(SS,kx[j],kz[k+my_rank*Loc_k])/(1-SS[4]);
		}
		ret = fftn(1, dims2, rS+k*lenk[0], iS+k*lenk[0],  1, FORWARD_SCALE);
		if (ret) return 1;
	}
	
	for(int j=0;j<lenk[0];j++){
		if (my_rank ==p-1) temptime=MPI_Wtime();
		double *rVec=new(std::nothrow)double [Loc_k];
		double *iVec=new(std::nothrow)double [Loc_k];
		for(int k=0;k<Loc_k;k++){
			rVec[k]=rS[k*lenk[0]+j];			
			iVec[k]=iS[k*lenk[0]+j];
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
				rS[k*lenk[0]+j]=rVec[k];
			delete [] rVec;delete [] iVec;					
	}
			//if (my_rank ==p-1) printf("R: fft for traces at %d row is completed, using %5.3f sec\n",i,MPI_Wtime()-temptime);

	delete [] iS;
	sum=0.0;sum2=0.0;
	for(int k=0;k<Loc_k;k++){ 
		for(int j=0;j<lenk[0];j++){
			sum=sum+rS[k*lenk[0]+j];
			sum2=sum2+rS[k*lenk[0]+j]*rS[k*lenk[0]+j];
		}	
	} 
	
   	MPI_Allreduce( &sum, &Global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &sum2, &Global_sum2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &facDes, &Global_facDes, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	mymean=Global_sum/lenk[0]/lenk[2];mystd=sqrt((Global_sum2-Global_sum*Global_sum/lenk[0]/lenk[2])/(lenk[0]*lenk[2]-1));
	Global_facDes=sqrt(Global_facDes/HOZx/VET);

	for(int k=0;k<Loc_k;k++){
			for(int j=0;j<lenk[0];j++)
				rS[k*lenk[0]+j]=Global_facDes*(rS[k*lenk[0]+j]-mymean)/mystd;
	}
	if(my_rank==0) printf("Time for Nonconditional Simulation of Vs: %6.4f sec\n",(MPI_Wtime()-Mstart));
	}		
	else {printf("Not enough processors!\n");return 1;}
	//cout<<mymean<<' '<<mystd<<' '<<Global_facDes<<endl;
	/*sprintf(rbinfile,"gR.bin.%d",my_rank);
	if( (rstream  = fopen( rbinfile, "wb" )) == NULL ) 
		printf( "The file %s was not opened\n", rbinfile );
	fwrite(rS,sizeof(double),lenk[0]*Loc_k,rstream);fclose(rstream);	
	MPI_Barrier( MPI_COMM_WORLD );
	if(my_rank==0) {
	for (int i=0;i<p;i++){
		sprintf(cmd,"cat gR.bin.%d >> gNCR_ax%d.bin ",i,int(ax[axi]));system(cmd);	
		sprintf(cmd,"rm gR.bin.%d",i);system(cmd);	
	}
	}
cout<<"stop"<<endl;cin>>temp;*/
	
	// ==========generate density in parallel============
	double RR1,RR2,RR3,t1=0.,t2=0.,t3=0.;
	double *rRho=new(std::nothrow)double [lenk[0]*Loc_k];//2D subvolume
	double *iRho=new(std::nothrow)double [lenk[0]*Loc_k];
	
	if (rRho && iRho) {
		if(my_rank==0) {printf("\n%d Realization of Density, ax=%6.5f,...\n",MCNUM+1,ax[axi]);Mstart=MPI_Wtime();}
//  0=az vertical,1=ax horizontal 1, 2=ay horizontal 2,3=Hurst, 4=nugget 
	PP[0]=cdfpar[0*3+0];PP[1]=ax[axi];PP[2]=PP[1];	
	PP[3]=cdfpar[0*3+1];PP[4]=cdfpar[0*3+2];//PP
	
	SP[0]=cdfpar[1*3+0];SP[1]=ax[axi];SP[2]=SP[1];
	SP[3]=cdfpar[1*3+1];SP[4]=cdfpar[1*3+2];//SP
	
	SS[0]=cdfpar[2*3+0];SS[1]=ax[axi];SS[2]=SS[1];
	SS[3]=cdfpar[2*3+1];SS[4]=cdfpar[2*3+2];//SS
	
	RP[0]=cdfpar[3*3+0];RP[1]=ax[axi];RP[2]=RP[1];
	RP[3]=cdfpar[3*3+1];RP[4]=cdfpar[3*3+2];//RP
	
	RS[0]=cdfpar[4*3+0];RS[1]=ax[axi];RS[2]=RS[1];
	RS[3]=cdfpar[4*3+1];RS[4]=cdfpar[4*3+2];//RS
	
	RR[0]=cdfpar[5*3+0];RR[1]=ax[axi];RR[2]=RR[1];
	RR[3]=cdfpar[5*3+1];RR[4]=cdfpar[5*3+2];//RR	

	facDes=0.0;
	for(int k=0;k<Loc_k;k++){
		for(int j=0;j<lenk[0];j++){
			rNB=ranR.doub();rNB2=ranP3.doub();rNB3=ranS3.doub();
			RR1=pskarman2d(RP,kx[j],kz[k+my_rank*Loc_k])*pskarman2d(RP,kx[j],kz[k+my_rank*Loc_k])/pskarman2d(PP,kx[j],kz[k+my_rank*Loc_k])/HOZx/VET;
			RR2=pskarman2d(RS,kx[j],kz[k+my_rank*Loc_k])-pskarman2d(SP,kx[j],kz[k+my_rank*Loc_k])*pskarman2d(RP,kx[j],kz[k+my_rank*Loc_k])/pskarman2d(PP,kx[j],kz[k+my_rank*Loc_k]);
			RR2=RR2*RR2/(pskarman2d(SS,kx[j],kz[k+my_rank*Loc_k])-pskarman2d(SP,kx[j],kz[k+my_rank*Loc_k])*pskarman2d(SP,kx[j],kz[k+my_rank*Loc_k])/pskarman2d(PP,kx[j],kz[k+my_rank*Loc_k]))/HOZx/VET;
			RR3=pskarman2d(RR,kx[j],kz[k+my_rank*Loc_k])/HOZx/VET-RR1-RR2;
			//if (my_rank==0) {printf("(RR1,RR2,RR3)=(%f,%f,%f)\n",RR1,RR2,RR3);}
			//cin>>temp;
			if (RR1<0) {t1=M_PI/2; RR1=-RR1;} if (RR2<0) {t2=M_PI/2; RR2=-RR2;} if (RR3<0) {t3=M_PI/2; RR3=-RR3;}
			rRho[k*lenk[0]+j]=lenk[0]*lenk[2]*(cos(2*M_PI*rNB2+t1)*sqrt(RR1)+cos(2*M_PI*rNB3+t2)*sqrt(RR2)+cos(2*M_PI*rNB+t3)*sqrt(RR3));
			iRho[k*lenk[0]+j]=lenk[0]*lenk[2]*(sin(2*M_PI*rNB2+t1)*sqrt(RR1)+sin(2*M_PI*rNB3+t2)*sqrt(RR2)+sin(2*M_PI*rNB+t3)*sqrt(RR3));
			facDes=facDes+pskarman2d(RR,kx[j],kz[k+my_rank*Loc_k])/(1-RR[4]);
		}
		ret = fftn(1, dims2, rRho+k*lenk[0], iRho+k*lenk[0],  1, FORWARD_SCALE);
		if (ret) return 1;
	}

	for(int j=0;j<lenk[0];j++){
		if (my_rank ==p-1) temptime=MPI_Wtime();
			double *rVec=new(std::nothrow)double [Loc_k];
			double *iVec=new(std::nothrow)double [Loc_k];
			for(int k=0;k<Loc_k;k++){
				rVec[k]=rRho[k*lenk[0]+j];			
				iVec[k]=iRho[k*lenk[0]+j];
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
				rRho[k*lenk[0]+j]=rVec[k];
			delete [] rVec;delete [] iVec;					
	}
	//if (my_rank ==p-1) printf("R: fft for traces at %d row is completed, using %5.3f sec\n",i,MPI_Wtime()-temptime);

	delete [] iRho;
	sum=0.0;sum2=0.0;
	for(int k=0;k<Loc_k;k++){ 
		for(int j=0;j<lenk[0];j++){
			sum=sum+rRho[k*lenk[0]+j];
			sum2=sum2+rRho[k*lenk[0]+j]*rRho[k*lenk[0]+j];
		}	
	} 
	
   	MPI_Allreduce( &sum, &Global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &sum2, &Global_sum2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &facDes, &Global_facDes, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	mymean=Global_sum/lenk[0]/lenk[2];mystd=sqrt((Global_sum2-Global_sum*Global_sum/lenk[0]/lenk[2])/(lenk[0]*lenk[2]-1));
	Global_facDes=sqrt(Global_facDes/HOZx/VET);

	for(int k=0;k<Loc_k;k++){
			for(int j=0;j<lenk[0];j++)
				rRho[k*lenk[0]+j]=Global_facDes*(rRho[k*lenk[0]+j]-mymean)/mystd;
	}
	//cout<<mymean<<' '<<mystd<<' '<<Global_facDes<<endl;
	/*sprintf(rbinfile,"gR.bin.%d",my_rank);
	if( (rstream  = fopen( rbinfile, "wb" )) == NULL ) 
		printf( "The file %s was not opened\n", rbinfile );
	fwrite(rRho,sizeof(double),lenk[0]*Loc_k,rstream);fclose(rstream);	
	MPI_Barrier( MPI_COMM_WORLD );
	if(my_rank==0) {
	for (int i=0;i<p;i++){
		sprintf(cmd,"cat gR.bin.%d >> gNCR_ax%d.bin ",i,int(ax[axi]));system(cmd);	
		sprintf(cmd,"rm gR.bin.%d",i);system(cmd);	
	}
	}
cout<<"stop"<<endl;cin>>temp;*/
	
	if(my_rank==0) 
		printf("Time for Nonconditional Simulation of Density: %6.4f sec\n",(MPI_Wtime()-Mstart));
	}		
	else {printf("Not enough processors!\n");return 1;}
	
	//======Start kriging part for three variables====

	if (condflag){
		if(MCNUM==0){
			xy2r(boreX,boreY,boreR);//to flat the connecting line of boreholes along x axis
			/*if (my_rank==5) {
					stream=fopen("boreRZ.txt","w");
					for(int i=0;i<boreX.size();i++)
						fprintf(stream,"%f %f\n",boreR[i],boreZ[i]);
					fclose(stream);
			}*/
				
			vector<double> ZCmat_R(boreX.size()*boreX.size()),ZCmat_P(boreX.size()*boreX.size()),ZCmat_S(boreX.size()*boreX.size());
			vector<double> ZLam_R(boreX.size()),ZLam_P(boreX.size()),ZLam_S(boreX.size());
			if (my_rank==0) {printf("\nData Kriging...\n");printf("Size of Data: %d X 1 Column\n",boreX.size());}
			temptime=MPI_Wtime();
			for (int i=0;i<boreX.size();i++){
				for (int j=0;j<boreX.size();j++){
					ZCmat_R[i*boreX.size()+j]=karman2d(RR,boreR[i]-boreR[j],boreZ[i]-boreZ[j]);
					ZCmat_P[i*boreX.size()+j]=karman2d(PP,boreR[i]-boreR[j],boreZ[i]-boreZ[j]);
					ZCmat_S[i*boreX.size()+j]=karman2d(SS,boreR[i]-boreR[j],boreZ[i]-boreZ[j]);
					//if (my_rank==5 && i==j) printf("i=%d,j=%d, karman2d=%f\n",i,j,karman2d(RR,boreR[i]-boreR[j],boreZ[i]-boreZ[j]));
				}
			}			
			
			//printf("PE%d: size of gvs=%d, size of Lam= %d\n",my_rank,gvs.size(),Lam.size());
			if (my_rank==0){
				if( (rstream  = fopen( "ZCmat_R.bin", "wb" )) == NULL ) 
					printf( "The file was not opened\n" );
				fwrite(&ZCmat_R[0],sizeof(double),boreX.size()*boreX.size(),rstream);fclose(rstream);

				if( (rstream  = fopen( "ZCmat_P.bin", "wb" )) == NULL ) 
					printf( "The file was not opened\n" );
				fwrite(&ZCmat_P[0],sizeof(double),boreX.size()*boreX.size(),rstream);fclose(rstream);

			}
			LUSolver(ZCmat_R,grho,ZLam_R);LUSolver(ZCmat_P,gvp,ZLam_P);LUSolver(ZCmat_S,gvs,ZLam_S);
			for(int k=0;k<Loc_k;k++){								
				for(int j=0;j<lenk[0];j++){
					S1=0.0;
					for (int ikar=0;ikar<boreX.size();ikar++){
						S1=S1+ZLam_R[ikar]*karman2d(RR,j*simds-boreR[ikar],(k+Loc_k*my_rank)*simds+GH[0]-boreZ[ikar]);
					}
					Zmat_R[k*lenk[0]+j]=S1;	
					
					S1=0.0;
					for (int ikar=0;ikar<boreX.size();ikar++){
						S1=S1+ZLam_P[ikar]*karman2d(PP,j*simds-boreR[ikar],(k+Loc_k*my_rank)*simds+GH[0]-boreZ[ikar]);
					}
					Zmat_P[k*lenk[0]+j]=S1;
					
					S1=0.0;
					for (int ikar=0;ikar<boreX.size();ikar++){
						S1=S1+ZLam_S[ikar]*karman2d(SS,j*simds-boreR[ikar],(k+Loc_k*my_rank)*simds+GH[0]-boreZ[ikar]);
					}
					Zmat_S[k*lenk[0]+j]=S1;					
				}							
			}
			if (my_rank==0) printf("Data Kriging completed; takes time %f sec\n",MPI_Wtime()-temptime);	
		}
		else if (my_rank==0)
			printf("Not the first realization, recycle the Data kriging model\n");
		sprintf(rbinfile,"Zmat_P.bin.%d",my_rank);
		if( (rstream  = fopen( rbinfile, "wb" )) == NULL ) 
			printf( "The file %s was not opened\n", rbinfile );
		fwrite(Zmat_R,sizeof(double),lenk[0]*Loc_k,rstream);fclose(rstream);	
		MPI_Barrier( MPI_COMM_WORLD );
		if(my_rank==0) {
			for (int i=0;i<p;i++){
				sprintf(cmd,"cat Zmat_P.bin.%d >> Zmat_P_ax%d.bin ",i,int(ax[axi]));system(cmd);	
				sprintf(cmd,"rm Zmat_P.bin.%d",i);system(cmd);	
			}
		}
		//MPI_Barrier( MPI_COMM_WORLD );
		
		if (my_rank==0) printf("Model Kriging...\n");		
		temptime=MPI_Wtime();
		double *Smat_R=new(std::nothrow)double [lenk[0]*Loc_k];
		double *Smat_P=new(std::nothrow)double [lenk[0]*Loc_k];
		double *Smat_S=new(std::nothrow)double [lenk[0]*Loc_k];//2D model data kriging subvolume
		vector<double> modY,modZ,modR;
		bore2mod(boreR,boreY,boreZ,modR,modY,modZ,simds);//assign likely irregularly sampled raw data to regular grid of the model
														// in 2D, boreY and modY are ignored
		
		vector<double> SCmat_R(modR.size()*modR.size()),SCmat_P(modR.size()*modR.size()),SCmat_S(modR.size()*modR.size());
		vector<double> SLam_R(modR.size()), Globalgrho(modR.size());
		vector<double> SLam_P(modR.size()), Globalgvp(modR.size());
		vector<double> SLam_S(modR.size()), Globalgvs(modR.size());
		
		//printf("PE%d: size of Mat=%d, size of Lam= %d, size of Globalgvs=%d\n",my_rank,Mat.nrows(),Lam.size(),Globalgvs.size());
		int k,rankID;
		for (int i=0;i<modR.size();i++){
			k=int((modZ[i]-GH[0])/simds) % Loc_k;
			rankID=(int((modZ[i]-GH[0])/simds)-k)/Loc_k;
			//printf("k=%d,rankID=%d\n",k,rankID);
			if (my_rank==rankID){
				Globalgrho[i]=rRho[k*lenk[0]+int(modR[i]/simds)];
				Globalgvp[i]=rP[k*lenk[0]+int(modR[i]/simds)];
				Globalgvs[i]=rS[k*lenk[0]+int(modR[i]/simds)];
			}
				MPI_Bcast(&Globalgrho[i], 1, MPI_DOUBLE, rankID, MPI_COMM_WORLD);
				MPI_Bcast(&Globalgvp[i], 1, MPI_DOUBLE, rankID, MPI_COMM_WORLD);
				MPI_Bcast(&Globalgvs[i], 1, MPI_DOUBLE, rankID, MPI_COMM_WORLD);
		}//all process share the same model data
		
		if (my_rank==0) {
			if( (stream  = fopen( "modRZ.txt", "w" )) == NULL ) 
				printf( "The file modRZ.txt was not opened\n");
			fprintf(stream,"R	Z	grho	gvp		gvs\n");
			for (int i=0;i<modR.size();i++){
				fprintf(stream,"%6.5f	%6.5f	%6.5f	%6.5f	%6.5f\n",modR[i],modZ[i],Globalgrho[i],Globalgvp[i],Globalgvs[i]);
			}
			fclose(stream);
			/*if( (stream  = fopen( "modRYR.txt", "w" )) == NULL ) 
				printf( "The file fx.txt was not opened\n");
			for (int i=0;i<modR.size();i++){
				fprintf(stream,"%6.5f	%6.5f	%6.5f\n",modR[i],modY[i],modR[i]);
			}
			fclose(stream);*/
		}
		
		for (int i=0;i<modR.size();i++){
			for (int j=0;j<modR.size();j++){
				SCmat_R[i*modR.size()+j]=karman2d(RR,modR[i]-modR[j],modZ[i]-modZ[j]);
				SCmat_P[i*modR.size()+j]=karman2d(PP,modR[i]-modR[j],modZ[i]-modZ[j]);
				SCmat_S[i*modR.size()+j]=karman2d(SS,modR[i]-modR[j],modZ[i]-modZ[j]);
			}
		}
		
			if (my_rank==0){
				if( (rstream  = fopen( "SCmat_R.bin", "wb" )) == NULL ) 
					printf( "The file was not opened\n" );
				fwrite(&SCmat_R[0],sizeof(double),modR.size()*modR.size(),rstream);fclose(rstream);

				if( (rstream  = fopen( "SCmat_P.bin", "wb" )) == NULL ) 
					printf( "The file was not opened\n" );
				fwrite(&SCmat_P[0],sizeof(double),modR.size()*modR.size(),rstream);fclose(rstream);

			}
		LUSolver(SCmat_R,Globalgrho,SLam_R);LUSolver(SCmat_P,Globalgvp,SLam_P);LUSolver(SCmat_S,Globalgvs,SLam_S);
		for(int k=0;k<Loc_k;k++){		
			for(int j=0;j<lenk[0];j++){				
				S1=0.0;
				for (int ikar=0;ikar<modR.size();ikar++){
					S1=S1+SLam_R[ikar]*karman2d(RR,j*simds-modR[ikar],(k+Loc_k*my_rank)*simds+GH[0]-modZ[ikar]);
				}
				Smat_R[k*lenk[0]+j]=S1;	
					
				S1=0.0;
				for (int ikar=0;ikar<modR.size();ikar++){
					S1=S1+SLam_P[ikar]*karman2d(PP,j*simds-modR[ikar],(k+Loc_k*my_rank)*simds+GH[0]-modZ[ikar]);
				}
				Smat_P[k*lenk[0]+j]=S1;
					
				S1=0.0;
				for (int ikar=0;ikar<modR.size();ikar++){
					S1=S1+SLam_S[ikar]*karman2d(SS,j*simds-modR[ikar],(k+Loc_k*my_rank)*simds+GH[0]-modZ[ikar]);
				}
				Smat_S[k*lenk[0]+j]=S1;						
			}		
		}
		if (my_rank==0) printf("Model Kriging completed; takes time %f sec\n",MPI_Wtime()-temptime);
		sprintf(rbinfile,"Smat_P.bin.%d",my_rank);
		if( (rstream  = fopen( rbinfile, "wb" )) == NULL ) 
			printf( "The file %s was not opened\n", rbinfile );
		fwrite(Smat_R,sizeof(double),lenk[0]*Loc_k,rstream);fclose(rstream);	
		MPI_Barrier( MPI_COMM_WORLD );
		if(my_rank==0) {
			for (int i=0;i<p;i++){
				sprintf(cmd,"cat Smat_P.bin.%d >> Smat_P_ax%d.bin ",i,int(ax[axi]));system(cmd);	
				sprintf(cmd,"rm Smat_P.bin.%d",i);system(cmd);	
			}
		}
		for(int k=0;k<Loc_k;k++){
			for(int j=0;j<lenk[0];j++){	
				rRho[k*lenk[0]+j]=myInterp(fr,xr,normcdf(Zmat_R[k*lenk[0]+j]+rRho[k*lenk[0]+j]-Smat_R[k*lenk[0]+j],0,1))+\
					myfitR_a+myfitR_b*((k+Loc_k*my_rank)*simds+GH[0]);
				rP[k*lenk[0]+j]=myInterp(fp,xp,normcdf(Zmat_P[k*lenk[0]+j]+rP[k*lenk[0]+j]-Smat_P[k*lenk[0]+j],0,1))+\
					myfitP_a+myfitP_b*((k+Loc_k*my_rank)*simds+GH[0]);
				rS[k*lenk[0]+j]=myInterp(fs,xs,normcdf(Zmat_S[k*lenk[0]+j]+rS[k*lenk[0]+j]-Smat_S[k*lenk[0]+j],0,1))+\
					myfitS_a+myfitS_b*((k+Loc_k*my_rank)*simds+GH[0]);
					//trendR.interp((k+Loc_k*my_rank)*simds);
			}
		//	printf("kriging condioning R: slice  at depth %d is completed, using %5.3f sec\n",k+Loc_k*my_rank,MPI_Wtime()-temptime);
		}
		delete [] Smat_R; delete [] Smat_P; delete [] Smat_S; 
	}
	else{
	//unconditional
		for(int k=0;k<Loc_k;k++){ 
			temptime=MPI_Wtime();		
			for(int j=0;j<lenk[0];j++){
				rRho[k*lenk[0]+j]=myInterp(fr,xr,normcdf(rRho[k*lenk[0]+j],0,1))+\
					myfitR_a+myfitR_b*((k+Loc_k*my_rank)*simds+GH[0]);
				rP[k*lenk[0]+j]=myInterp(fp,xp,normcdf(rP[k*lenk[0]+j],0,1))+\
					myfitP_a+myfitP_b*((k+Loc_k*my_rank)*simds+GH[0]);		
				rS[k*lenk[0]+j]=myInterp(fs,xs,normcdf(rS[k*lenk[0]+j],0,1))+\
					myfitS_a+myfitS_b*((k+Loc_k*my_rank)*simds+GH[0]);
			}
		}
	}
	// =====kriging end====
	
	if (MaxMC==1){	
		//===== Output Density Rho============	
		sprintf(rbinfile,"rR.bin.%d",my_rank);
		if( (rstream  = fopen( rbinfile, "wb" )) == NULL ) 
			printf( "The file %s was not opened\n", rbinfile );
		fwrite(rRho,sizeof(double),lenk[0]*Loc_k,rstream);fclose(rstream);	
		//===== Output Vp============	
		sprintf(rbinfile,"rP.bin.%d",my_rank);
		if( (rstream  = fopen( rbinfile, "wb" )) == NULL ) 
			printf( "The file %s was not opened\n", rbinfile );
		fwrite(rP,sizeof(double),lenk[0]*Loc_k,rstream);fclose(rstream);		
		//===== Output Vs============	
		sprintf(rbinfile,"rS.bin.%d",my_rank);
		if( (rstream  = fopen( rbinfile, "wb" )) == NULL ) 
			printf( "The file %s was not opened\n", rbinfile );
		fwrite(rS,sizeof(double),lenk[0]*Loc_k,rstream);fclose(rstream);	
		MPI_Barrier( MPI_COMM_WORLD );		
		
		if(my_rank==0) {
			if (condflag){
				for (int i=0;i<p;i++){
					sprintf(cmd,"cat rR.bin.%d >> %s2DRho_cond_ax%d.bin ",i,outfile,int(ax[axi]));system(cmd);	
					sprintf(cmd,"rm rR.bin.%d",i);system(cmd);	
					sprintf(cmd,"cat rP.bin.%d >> %s2DVp_cond_ax%d.bin ",i,outfile,int(ax[axi]));system(cmd);	
					sprintf(cmd,"rm rP.bin.%d",i);system(cmd);
					sprintf(cmd,"cat rS.bin.%d >> %s2DVs_cond_ax%d.bin ",i,outfile,int(ax[axi]));system(cmd);	
					sprintf(cmd,"rm rS.bin.%d",i);system(cmd);
				}			
			}
			else {
				for (int i=0;i<p;i++){
					sprintf(cmd,"cat rR.bin.%d >> %s2DRho_uncond_ax%d.bin ",i,outfile,int(ax[axi]));system(cmd);	
					sprintf(cmd,"rm rR.bin.%d",i);system(cmd);	
					sprintf(cmd,"cat rP.bin.%d >> %s2DVp_uncond_ax%d.bin ",i,outfile,int(ax[axi]));system(cmd);	
					sprintf(cmd,"rm rP.bin.%d",i);system(cmd);	
					sprintf(cmd,"cat rS.bin.%d >> %s2DVs_uncond_ax%d.bin ",i,outfile,int(ax[axi]));system(cmd);	
					sprintf(cmd,"rm rS.bin.%d",i);system(cmd);
				}
			}
			sprintf(cmd,"mv *.bin ../models/");system(cmd);
		}
	}
	
	if (RockPhyTag && Nvar==3) {
		//====Rock Physics relationship: Modified Biot-Gassmann=========	
		int IntNUM=10;
		vector<double> spratio(IntNUM+1),cc(IntNUM+1);
		for (int i=0;i<=IntNUM;i++)
			cc[i]=double(i)/double(IntNUM);
		
		GHVol=0.0;
		for(int k=0;k<Loc_k;k++){
			temptime=MPI_Wtime();		
			for(int j=0;j<lenk[0];j++){
				for(int pp=0;pp<=IntNUM;pp++)
					spratio[pp]=fsp(cc[pp],rRho[k*lenk[0]+j]/1000);
				if (rS[k*lenk[0]+j]/rP[k*lenk[0]+j]<spratio[0])
					GHCon=0.0;
				else if (rS[k*lenk[0]+j]/rP[k*lenk[0]+j]>spratio[IntNUM])
					GHCon=1.0;
				else
					GHCon=myInterp(spratio,cc,rS[k*lenk[0]+j]/rP[k*lenk[0]+j]);
				if (GHCon<0)
					GHCon=0;
				Poro=fporo(GHCon,rRho[k*lenk[0]+j]/1000);
				//printf("GHCon=%f, Poro=%f\n",GHCon,Poro);
				GHVol=GHVol+GHCon*Poro;
			}
			//if (MCNUM==0){
			//	printf("The 1st realization MBGL with ax=%f ay=%f: slice  at depth %d is completed, using %5.3f sec\n",ax[axi],ay[ayi],k+Loc_k*my_rank,MPI_Wtime()-temptime);							
			//}
		}
		delete [] rRho;	delete [] rP;delete [] rS;
		GHVol=GHVol/lenk[0]/lenk[2];
		MPI_Allreduce( &GHVol, Global_GHVol+MCNUM, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
		//====The End of Rock Physics====
	}
	MCNUM=MCNUM+1;
	if(my_rank==0) printf("Monto Carlo Realization %d, using %5.3f sec \n",MCNUM,MPI_Wtime()-MCtime);
	
	}//MCNUM
	 
//=======Write disk==========	
	if (RockPhyTag && Nvar==3) {
		if (my_rank==0){
			if (condflag)
				sprintf(rbinfile,"%s_GHVol_cond_np%d_x%d_ax%d_ds%3.2f.txt",outfile,p,int(HOZx),int(ax[axi]),simds);
			else
				sprintf(rbinfile,"%s_GHVol_uncond_np%d_x%d_ax%d_ds%3.2f.txt",outfile,p,int(HOZx),int(ax[axi]),simds);
			if( (rstream  = fopen( rbinfile, "w" )) == NULL ) 
				printf( "The file %s was not opened\n", rbinfile );
			for (int i=0;i<MaxMC;i++)
				fprintf(rstream,"%10.9f\n",Global_GHVol[i]);
		fclose(rstream);
		}
	}
	MPI_Barrier( MPI_COMM_WORLD );	
	} //ax
	if(my_rank==0) {printf("Total Elapse Time: %6.4f hours \n", (MPI_Wtime()-Telapse)/3600);}
	MPI_Finalize();
}	
