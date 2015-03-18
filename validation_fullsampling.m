function [sout,kdata3out,kcoors]=validation_fullsampling(Wx,Wy,Wz,shift,cverts,cfaces,singleflag)
%% create kspace sampling matrix
aNx=[64];
aNy=[64];
aNz=[64];
kdata3out={};
sout={};
kcoors={};
for t=1:length(aNx)
    %TODO: replace the below with gengrid function
    t
    Nx=aNx(t);
    Ny=aNy(t);
    Nz=aNz(t);
    
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
    
    % compute ground truth fourier transform of cube
    %s=(Nx*Ny*Nz)*(1/(Lx*Ly*Lz))*Wx*Wy*Wz*sinc(Wx*k(:,:,:,1)).*sinc(Wy*k(:,:,:,2)).*sinc(Wz*k(:,:,:,3));
    s=Wx*Wy*Wz*sinc(Wx*k(:,:,:,1)).*...
        sinc(Wy*k(:,:,:,2)).*...
        sinc(Wz*k(:,:,:,3)).*...
        exp(-1i*2*pi*k(:,:,:,1)*shift(1)).*...
        exp(-1i*2*pi*k(:,:,:,2)*shift(2)).*...
        exp(-1i*2*pi*k(:,:,:,3)*shift(3));
    % volume_browser(abs(s));

    
    % compute via polyhedral transform
    klin=reshape(k,[size(k,1)*size(k,2)*size(k,3),3])';
    cfacescell=mat2cell(cfaces,ones(1,size(cfaces,1)),[3]);
    
    [kdata]=polyhedralFT({cverts},{cfacescell},1,klin,singleflag);
    
%     [kdata]=polyhedralFTscriptsingle(cverts,cfacescell,klin);
    kdata3=reshape(kdata,Nx,Ny,Ny);
    % volume_browser(real(double(p_pt)));
    % volume_browser(imag(double(p_pt)));
    
    sout{t}=s;
    kdata3out{t}=kdata3;
    kcoors{t}=k;
end
