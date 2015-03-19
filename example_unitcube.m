%% make cube
nfaces=12;
   
Wx=1; %total width of cube
Wy=1;
Wz=1;
[cverts,~]=makeCube(Wx/2,Wy/2,Wz/2);

%compute the convex hull
[cfaces]=convhull(cverts(1,:),cverts(2,:),cverts(3,:));

%% create kspace sampling matrix
aNx=[64];
aNy=[64];
aNz=[64];
kdata3out={};
sout={};
kcoors={};

Nx=aNx(1);
Ny=aNy(1);
Nz=aNz(1);

Lx=2*Wx;
Ly=2*Wy;
Lz=2*Wz;

if mod(Nx,2)==0
    tNx=Nx+1;
    kxcoors=(-floor(tNx/2):floor(tNx/2))*(1/Lx);
    kxcoors=kxcoors(1:end-1);
else
    kxcoors=(-floor(Nx/2):floor(Nx/2))*(1/Lx);
end
if mod(Ny,2)==0
    tNy=Ny+1;
    kycoors=(-floor(tNy/2):floor(tNy/2))*(1/Ly);
    kycoors=kycoors(1:end-1);
else
    kycoors=(-floor(Ny/2):floor(Ny/2))*(1/Ly);
end
if mod(Nz,2)==0
    tNz=Nz+1;
    kzcoors=(-floor(tNz/2):floor(tNz/2))*(1/Lz);
    kzcoors=kzcoors(1:end-1);
else
    kzcoors=(-floor(Nz/2):floor(Nz/2))*(1/Lz);
end

clear k;
[k(:,:,:,1),k(:,:,:,2),k(:,:,:,3)]=ndgrid(kxcoors,kycoors,kzcoors);
    
%% compute ground truth fourier transform of cube
s=Wx*Wy*Wz*sinc(Wx*k(:,:,:,1)).*...
    sinc(Wy*k(:,:,:,2)).*...
    sinc(Wz*k(:,:,:,3));

%% compute via polyhedral transform
klin=reshape(k,[size(k,1)*size(k,2)*size(k,3),3])';
cfacescell=mat2cell(cfaces,ones(1,size(cfaces,1)),[3]);
singleflag=0;
[kdata]=polyhedralFT({cverts},{cfacescell},1,klin,singleflag);
kdata3=reshape(kdata,Nx,Ny,Ny);

%% display results
figure; imagesc(abs(kdata3(:,:,32)).^(1/3));
rec=fftshift(ifftn(fftshift(kdata3)));
figure; imagesc(abs(rec(:,:,32)));
figure; imagesc(angle(rec(:,:,32)));

