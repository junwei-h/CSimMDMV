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
    Skrig=['../models/Smat_P_ax',char(C{6}(1)),'.bin'];
    Zkrig=['../models/Zmat_P_ax',char(C{6}(1)),'.bin'];
end

pfile=['../models/',char(C{end-1}),'2DVp',filename,'ax',char(C{6}(1)),'.bin'];
sfile=['../models/',char(C{end-1}),'2DVs',filename,'ax',char(C{6}(1)),'.bin'];
rfile=['../models/',char(C{end-1}),'2DRho',filename,'ax',char(C{6}(1)),'.bin'];

Hozx=str2num(char(C{1}(1)));
z0=str2num(char(C{1}(3)));zN=str2num(char(C{1}(4)));simds=str2num(char(C{1}(5)));
nx=Hozx/simds;ny=(zN-z0+simds)/simds;

savetag=0;

vp2d=zeros(nx,ny);
fid=fopen(pfile,'r');
vp2d(:)=fread(fid,inf,'double');
vp2d=permute(vp2d,[2,1]);%matlab is colume major and C++ is row major
fclose(fid);

vs2d=zeros(nx,ny);
fid=fopen(sfile,'r');
vs2d(:)=fread(fid,inf,'double');
vs2d=permute(vs2d,[2,1]);
fclose(fid);

r2d=zeros(nx,ny);
fid=fopen(rfile,'r');
r2d(:)=fread(fid,inf,'double');
r2d=permute(r2d,[2,1]);
fclose(fid);

xx=simds:simds:simds*nx;yy=z0:simds:zN;
figure;
subplot(3,3,[1 4]);
imagesc(xx,yy,vp2d);colorbar;xlabel('Horizon (m)');ylabel('Depth (m)');title('v_p (m/s)');caxis([1700 4400]);
subplot(337);hist(vp2d(:),100);grid on;
subplot(3,3,[2 5]);
imagesc(xx,yy,vs2d);colorbar;xlabel('Horizon (m)');
if str2num(char(C{4}))==0
    title('v_s (m/s): Nonconditional');
else
    title('v_s (m/s): Conditional');
end
caxis([800 2200]);
subplot(338);hist(vs2d(:),100);grid on;
subplot(3,3,[3 6]);imagesc(xx,yy,r2d);colorbar;xlabel('Horizon (m)');title('\rho (kg/m^3)');
subplot(339);hist(r2d(:),100);grid on;

if str2num(char(C{4}))~=0
    Z2D=zeros(nx,ny);
    fid=fopen(Zkrig,'rb');
    Z2D(:)=fread(fid,inf,'double');
    fclose(fid);
    Z2D=permute(Z2D,[2,1]);

    S2D=zeros(nx,ny);
    fid=fopen(Skrig,'rb');
    S2D(:)=fread(fid,inf,'double');
    fclose(fid);
    S2D=permute(S2D,[2,1]);
    
    figure(2)
    subplot(121);
    imagesc(xx,yy,Z2D);set(gca,'zdir','reverse');
    box on;axis tight;colorbar;
    xlabel('x');ylabel('y');title('Data Kriging Model');
    subplot(122);
    imagesc(xx,yy,S2D);shading interp;set(gca,'zdir','reverse');
    box on;axis tight;colorbar;
    xlabel('x');ylabel('y');title('Model Kriging Model');
end

if savetag
    a=size(vp2d);
    fid=fopen('./yourfilename.vp','wb');
    for j=1:a(2)
        for i=1:a(1)
            fwrite(fid,vp2d(i,j),'float');
        end
    end
    fclose(fid);
        fid=fopen('./yourfilename.vs','wb');
    for j=1:a(2)
        for i=1:a(1)
            fwrite(fid,vs2d(i,j),'float');
        end
    end
    fclose(fid);
        fid=fopen('./yourfilename.rho','wb');
    for j=1:a(2)
        for i=1:a(1)
            fwrite(fid,r2d(i,j),'float');
        end
    end
    fclose(fid);
end
