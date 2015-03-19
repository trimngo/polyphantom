% function valdiation_fullsampling_UT
%example usage of polyhedralFTscript
%% construct an ellipoid
javaaddpath external/STBBFolder/JarFile/STBB.jar
import STBB.*
%delta_x, delta_y, delta_z,        a,       b,       c,            phi,  theta,  psi,     rho 
eparm =[0,       0,       0,     0.50,    0.75,     0.30,              0,      0,    0,      1];
sl3d = SheppLogan3D(eparm); 

Nx=64;
Ny=64;
Nz=64;
Lx=2;
Ly=2;
Lz=2;

%k-space grid
factor=(Nx*Ny*Nz)/(Lx*Ly*Lz);
k=gengrid(Nx,Ny,Nz,Lx,Ly,Lz);
klin = reshape(k,Nx*Ny*Nz,3);

%generate image domain grid
dx=Lx/Nx;
dy=Ly/Ny;
dz=Lz/Nz;
xcoor=-Lx/2:dx:Lx/2-dx;
ycoor=-Ly/2:dy:Ly/2-dy;
zcoor=-Lz/2:dz:Lz/2-dz;

[X,Y,Z]=ndgrid(xcoor,ycoor,zcoor);
plin=[X(:) Y(:) Z(:)];

idsignal3d = sl3d.ImageDomainSignal(plin); 
idata=reshape(idsignal3d, [Nx Ny Nz]);

% get the Fourier domain signal at (0.25,0.15,0.05).
% 
fdsignal3d = sl3d.FourierDomainSignal(klin); 
data=reshape(complex(fdsignal3d(:,1),fdsignal3d(:,2)),[Nx Ny Nz]);

%% optionally save the gold standard
% fn=['optional/ellipsoidFT_64_GS.mat'];
% save(fn,'data');

%% reconstruct koay
img=factor*fftshift(ifftn(fftshift(data)));
figure; imagescn(abs(img),[],[],[],3);

%% show image space
figure; imagescn(abs(idata),[],[],[],3);
%% loop
% nf=4:2:60;
% nf=60;
nf=10:10:200;
diff={};
ndiff={};
kimgsave={};
kdiff={};
nkdiff={};
diffl2=[];
diffmax=[];
ndiffmax=[];
nfeff=[];

for a=1:length(nf)
    a
    %% make ellipsoid mesh
    [x, y, z] = ellipsoid(0,0,0,eparm(4),eparm(5),eparm(6),nf(a));
%     figure;surfl(x, y, z)
%     colormap copper
%     axis equal

    k=convhull(x,y,z);
    nfeff(a)=size(k,1);
%     figure; trisurf(k,x,y,z, 'Facecolor','cyan'); axis equal;
   %print('-painter','-deps','-r600',['images/ellipsoid_' num2str(nfeff(a)) '.eps']);
    %% compute polyft
    cverts=[x(:) y(:) z(:)]';
    cfacescell=k';
    singleflag=false;
    debug=true;
%     [kdata]=polyhedralFT({cverts},{cfacescell},1,klin',singleflag,debug);
    %[kdata]=polyhedralFT({cverts},{cfacescell},1,klin',singleflag,false);
    
    %% reconstruct ours
    %kdata2=reshape(kdata,[Nx Ny Nz]);
        fn=['optional/ellipsoidFT_64_' num2str(nfeff(a)) '.mat'];
%     fn=['optional/ellipsoidFT_64_GS.mat'];
    %save(fn,'kdata2');
    kdata2=importdata(fn);
    kimg=factor*fftshift(ifftn(fftshift(kdata2)));
    kimgsave{a}=kimg;
    %     figure; imagescn(abs(kimg),[],[],[],3);
    
%     kdata2_fast=reshape(kdata_fast,[Nx Ny Nz]);
    
    %% calculate difference image
    diff{a}=img-kimg;
    
    %% compute normalized difference image
    ndiff{a}=abs(img-kimg)./abs(img);
    
    %% calc error
    diffl2(a)=norm(diff{a}(:))/norm(img(:)); %normalized
    diffmax(a)=max(abs(diff{a}(:)));
    
    %% calc max normalized diff
    ndiffmax(a)=max(abs(ndiff{a}(:)));
    
    %% display
%     dp=cat(4,img,kimg,diff,ndiff);
%     figure; imagescn(abs(dp),[],[],[],3);

    %% calculate kspace diff img
    kdiff{a}=kdata2-data;
    kdiffl2(a)=norm(kdiff{a}(:))/norm(data(:));
    nkdiff{a}=abs(kdata2-data)./abs(data);
    
%     kdiff_fast=kdata2_fast-data;
%     nkdiff_fast=(kdata2_fast-data)./data;
    
%     kdp=cat(4,kdata2,data,kdiff,nkdiff);
    
%     kdp=cat(4,kdata2,data,kdata2_fast,kdiff,nkdiff,kdiff_fast,nkdiff_fast);
%     figure; imagescn(abs(kdp),[],[],[],3);
    
    %% compare the debug and non debug methods
%     mdiff=(kdata2-kdata2_fast);
%     nmdiff=mdiff./kdata2;
%     figure; imagescn(abs(mdiff));
end

%% generate ellipsoid mesh
for a=1:4%length(nf)
    a
    %% make ellipsoid mesh
    [x, y, z] = ellipsoid(0,0,0,eparm(4),eparm(5),eparm(6),nf(a));

    k=convhull(x,y,z);
    nfeff(a)=size(k,1);
    figure; trisurf(k,x,y,z, 'Facecolor','white'); axis equal;
   %print('-painter','-deps','-r600',['images/ellipsoid_' num2str(nfeff(a)) '.eps']);
end

%% plot result
% figure; plot(nfeff,diffl2); xlabel('faces'); ylabel('|error|/|gs|');
figure; plot(nfeff,log10(kdiffl2),'-o'); 
xlabel('faces','FontSize',20);
% ylabel('log_{10}(|error|/|gs|)','FontSize',16);
ylabel('$$\|s_{gs}-s_{3D}\|/\|s_{gs}\| $$','FontSize',20,'interpreter','latex');
% title('log of normalized l_2 norm of error versus number of faces','FontSize',16);
set(gca,'FontSize',20)
axis square;
print('-painter','-deps','-r600',['images/ellipsoid_normalizedl2error.eps']);
%figure; plot(nfeff,diffmax);
%figure; plot(nfeff,ndiffmax);

%% show imspace reconstruction
kzs=(Nz/2)+1;
dp=cat(3,kimgsave{1}(:,:,kzs),...
    kimgsave{2}(:,:,kzs),...
    kimgsave{3}(:,:,kzs),...
    kimgsave{4}(:,:,kzs));
figure; imagescn(abs(dp),[],[1 4]);

%% line profiles
kzs=(Nz/2)+1;
kxs=(Nx/2)+1;
dist=-Ly/2:Ly/Ny:Ly/2-Ly/Ny;
gsp=abs(img(kxs,:,kzs));
p{1}=abs(kimgsave{1}(kxs,:,kzs));
p{2}=abs(kimgsave{2}(kxs,:,kzs));
p{3}=abs(kimgsave{3}(kxs,:,kzs));
p{4}=abs(kimgsave{20}(kxs,:,kzs));
figure; plot(dist,gsp,'-o'); hold all;
plot(dist,p{1});
plot(dist,p{2});
plot(dist,p{3});
plot(dist,p{4});
axis equal;
legend('GS','180','760', '1,740', '79,600');
ylabel('intensity','FontSize',16);
xlabel('position','FontSize',16);
set(gca,'FontSize',16);

%% show where error located in kspace
cr=linspace(0,1,256);
cb=[linspace(0,1,128) linspace(1,0,128)];
cg=linspace(1,0,256);
cm=[cr(:) cb(:) cg(:)];

kzs=(Nz/2)+1;
dp=cat(3,kdiff{1}(:,:,kzs),...
    kdiff{2}(:,:,kzs),...
    kdiff{3}(:,:,kzs),...
    kdiff{20}(:,:,kzs));
figure; imagescn(log(abs(dp)),[],[1 4]);
colormap(cm);

kzs=(Nz/2)+1;
dp=cat(3,nkdiff{1}(:,:,kzs),...
    nkdiff{2}(:,:,kzs),...
    nkdiff{3}(:,:,kzs),...
    nkdiff{20}(:,:,kzs));
figure; imagescn(log(abs(dp)),[],[1 4]);
colormap(cm);

%% show where error is in imagespace
kzs=(Nz/2)+1;
dp=cat(3,diff{1}(:,:,kzs),...
    diff{2}(:,:,kzs),...
    diff{3}(:,:,kzs),...
    diff{20}(:,:,kzs));
figure; imagescn(abs(dp),[],[1 4]);

kzs=(Nz/2)+1;
dp=cat(3,ndiff{1}(:,:,kzs),...
    ndiff{2}(:,:,kzs),...
    ndiff{3}(:,:,kzs),...
    ndiff{20}(:,:,kzs));
figure; imagescn(log(abs(dp)),[],[1 4]);

% print('-painter','-deps','-r600',['images/ellipsoid_larger_normabsl2error_bar.eps']);
%figure; imagescn(abs(diff),[],[],[],3); %abs error
%figure; imagescn(abs(ndiff),[],[],[],3); %relative error

%figure; imagescn(abs(kdiff),[],[],[],3); %abs error
%figure; imagescn(abs(nkdiff),[],[],[],3); %relative error

% %% show ellipsoid mesh with different numbers of faces
% for a=[1 ceil(length(nf)/2) length(nf)]
%     [x, y, z] = ellipsoid(0,0,0,0.69,0.92,0.9,nf(a));
% %     figure;surfl(x, y, z)
% %     colormap copper
% %     axis equal
% 
%     k=convhull(x,y,z);
%     nfeff(a)=size(k,1);
%     figure; trisurf(k,x,y,z, 'Facecolor','white'); axis equal;
% end