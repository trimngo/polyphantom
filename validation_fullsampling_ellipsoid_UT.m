%% ellipsoid params
%delta_x, delta_y, delta_z,        a,       b,       c,            phi,  theta,  psi,     rho 
eparm =[0,       0,       0,     0.50,    0.75,     0.30,              0,      0,    0,      1];

%% construct kspace sampling 
Nx=64;
Ny=64;
Nz=64;
Lx=2;
Ly=2;
Lz=2;

factor=(Nx*Ny*Nz)/(Lx*Ly*Lz);
k=gengrid(Nx,Ny,Nz,Lx,Ly,Lz);
klin = reshape(k,Nx*Ny*Nz,3);

%% setup loop
% nf=4:2:60;
% nf=60;
nf=10:10:200;
diffl2=[];
diffmax=[];
ndiffmax=[];
nfeff=[];

%% data generation loop
for a=1:length(nf)
    a
    %% make ellipsoid mesh
    [x, y, z] = ellipsoid(0,0,0,eparm(4),eparm(5),eparm(6),nf(a));
%     figure;surfl(x, y, z)
%     colormap copper
%     axis equal

    k=convhull(x,y,z);
    nfeff(a)=size(k,1);
    figure; trisurf(k,x,y,z, 'Facecolor','cyan'); axis equal;

    %% write the mesh
    writeoff(['ellipsoid_' num2str(nfeff(a)) '.off'],verts,faces);
    
    %% compute polyft
    cverts=[x(:) y(:) z(:)]';
    cfacescell=k';
    singleflag=false;
    debug=true;
%     [kdata]=polyhedralFT({cverts},{cfacescell},1,klin',singleflag,debug);
    [kdata]=polyhedralFT({cverts},{cfacescell},1,klin',singleflag,false);
    
    %% reconstruct ours
    kdata2=reshape(kdata,[Nx Ny Nz]);
    fn=['optional/ellipsoidFT_64_' num2str(nfeff(a)) '.mat'];
    save(fn,'kdata2');
end

%% plot faces verus nf
figure; plot(nf,nfeff);

%% write out the first three and last mesh
for a=[1:3 length(nf)]
    a
    %% make ellipsoid mesh
    [x, y, z] = ellipsoid(0,0,0,eparm(4),eparm(5),eparm(6),nf(a));
%     figure;surfl(x, y, z)
%     colormap copper
%     axis equal

    k=convhull(x,y,z);
    nfeff(a)=size(k,1);
%     figure; trisurf(k,x,y,z, 'Facecolor','cyan'); axis equal;

    %% write the mesh
    cverts=[x(:) y(:) z(:)]';
    cfaces=k';
    writeoff(['ellipsoid_' num2str(nfeff(a)) '.off'],cverts,cfaces);

end

