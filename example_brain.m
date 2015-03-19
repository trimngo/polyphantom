%% Load in full brain mesh
prefix='braindataoff/';
fname{1}='centered_OAS1_0001_WMSurface.off';
fname{2}='centered_OAS1_0001_GMSurface.off';
fname{3}='centered_OAS1_0001_centralSurface.off';

%% read in all meshes
for i=1:length(fname)
    [cverts{i}, cfaces{i}]=tri_readoff([prefix fname{i}]);
end

%% set the intensities of each surface
colors(1)=38; %white matter
colors(2)=74;   %grey matter
colors(3)=10;   %in between

%% compute the extent for the entire dataset
extents=extent([cverts{:}]);
a=range(extents,2)*1.1;
Lx=a(1);
Ly=a(2);
Lz=a(3);

%% create kspace sampling matrix
%number of samples in each dimension (must be positive)
Nx=[128];
Ny=[128];
Nz=[128];

%there is one extra negative freq. sample
tNx=Nx+1;
kxcoors=(-floor(tNx/2):floor(tNx/2))*(1/Lx);
kxcoors=kxcoors(1:end-1);

tNy=Ny+1;
kycoors=(-floor(tNy/2):floor(tNy/2))*(1/Ly);
kycoors=kycoors(1:end-1);

tNz=Nz+1;
kzcoors=(-floor(tNz/2):floor(tNz/2))*(1/Lz);
kzcoors=kzcoors(1:end-1);

[k(:,:,:,1),k(:,:,:,2),k(:,:,:,3)]=ndgrid(kxcoors,kycoors,kzcoors);

%% compute via polyhedral transform
klin=reshape(k,[size(k,1)*size(k,2)*size(k,3),3])';
cfacescell=mat2cell(cfaces,ones(1,size(cfaces,1)),[3]);
[kdata]=polyhedralFT(cverts,cfaces,colors,klin);
kdata3=reshape(kdata,Nx,Ny,Ny);

%% reconstruct
rec=fftshift(ifftn(ifftshift(kdata3)));

%% show
figure; imagescn(abs(rec),[],[],[],3);
figure; imagescn(abs(recomb),[],[],[],3);