clear all
close all
inpf='../par/CSimMD3V.inp';
fid=fopen(inpf,'r');

i=1;tline=fgetl(fid);
while tline > 0
    TF=strncmp(tline,'#',1);
    if ~TF
        C(i)=textscan(tline,'%s');i=i+1;
    end
    tline=fgetl(fid);
end
fclose(fid);

if str2num(char(C{4}))==0
    filename='_uncond_'; %Conditional 'cond', or unconditional 'uncond'
else
    filename='_cond_';
    Skrig=['../models/Smat_P_ax',char(C{6}(1)),'_ay',char(C{7}(1)),'.bin'];
    Zkrig=['../models/Zmat_P_ax',char(C{6}(1)),'_ay',char(C{7}(1)),'.bin'];
end

pfile=['../models/',char(C{end-1}),'3DVp',filename,'ax',char(C{6}(1)),'_ay',char(C{7}(1)),'.bin'];
sfile=['../models/',char(C{end-1}),'3DVs',filename,'ax',char(C{6}(1)),'_ay',char(C{7}(1)),'.bin'];
rfile=['../models/',char(C{end-1}),'3DRho',filename,'ax',char(C{6}(1)),'_ay',char(C{7}(1)),'.bin'];

Hozx=str2num(char(C{1}(1)));Hozy=str2num(char(C{1}(2)));
z0=str2num(char(C{1}(3)));zN=str2num(char(C{1}(4)));simds=str2num(char(C{1}(5)));
nx=Hozx/simds;ny=Hozy/simds;nz=(zN-z0+simds)/simds;

savetag=0;

vp3D=zeros(nx,ny,nz);
fid=fopen(pfile,'rb');
vp3D(:)=fread(fid,inf,'double');
fclose(fid);
vp3D=permute(vp3D,[2,1,3]);%matlab is colume major and C++ is row major

vs3D=zeros(nx,ny,nz);
fid=fopen(sfile,'rb');
vs3D(:)=fread(fid,inf,'double');
fclose(fid);
vs3D=permute(vs3D,[2,1,3]);

r3D=zeros(nx,ny,nz);
fid=fopen(rfile,'rb');
r3D(:)=fread(fid,inf,'double');
fclose(fid);
r3D=permute(r3D,[2,1,3]);

if str2num(char(C{4}))~=0
    Z3D=zeros(nx,ny,nz);
    fid=fopen(Zkrig,'rb');
    Z3D(:)=fread(fid,inf,'double');
    fclose(fid);
    Z3D=permute(Z3D,[2,1,3]);

    S3D=zeros(nx,ny,nz);
    fid=fopen(Skrig,'rb');
    S3D(:)=fread(fid,inf,'double');
    fclose(fid);
    S3D=permute(S3D,[2,1,3]);
    
    figure(2)
    subplot(121);
    slice(Z3D,nx/2,ny/2,nz/2);shading interp;set(gca,'zdir','reverse');
    box on;axis tight;colorbar;
    xlabel('x');ylabel('y');title('Data Kriging Model');
    subplot(122);
    slice(S3D,nx/2,ny/2,nz/2);shading interp;set(gca,'zdir','reverse');
    box on;axis tight;colorbar;
    xlabel('x');ylabel('y');title('Model Kriging Model');
end

figure(1);subplot(3,3,[1,4]);slice(vp3D,1,1,1);shading interp;set(gca,'zdir','reverse');
box on;axis tight;colorbar;caxis([1700 4000]);
xlabel('x');ylabel('y');
if str2num(char(C{4}))==0
    title('v_p (m/s): Nonconditional');
else
    title('v_p (m/s): Conditional');
end
subplot(337);hist(vp3D(:),1000);grid on;axis tight;

subplot(3,3,[2,5]);slice(vs3D,1,1,1);shading interp;set(gca,'zdir','reverse');
box on;axis tight;colorbar;caxis([600,2000]);
xlabel('x');ylabel('y');title('v_s (m/s)');
subplot(338);hist(vs3D(:),1000);grid on;axis tight;

subplot(3,3,[3,6]);slice(r3D,1,1,1);shading interp;set(gca,'zdir','reverse');
box on;axis tight;colorbar;xlabel('x');ylabel('y');title('\rho (kg/m^3)');
subplot(339);hist(r3D(:),1000);grid on;axis tight;

if savetag
fid=fopen('./filename.vp','wb');
a=size(vp3D);
 for i=1:a(1)
     for j=1:a(2)
        for k=1:a(3)         
            fwrite(fid,vp3D(i,j,k),'float');
        end
    end
 end
fclose(fid);

fid=fopen('./filename.vs','wb');
 for i=1:a(1)
     for j=1:a(2)
        for k=1:a(3)         
            fwrite(fid,vs3D(i,j,k),'float');
        end
    end
 end
fclose(fid);

fid=fopen('./filename.rho','wb');
 for i=1:a(1)
     for j=1:a(2)
        for k=1:a(3)         
            fwrite(fid,rho3D(i,j,k),'float');
        end
    end
 end
fclose(fid);
end
