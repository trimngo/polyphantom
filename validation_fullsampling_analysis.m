%% setup parameters
nfaces=12*2.^[0:13];
% nfaces=[12 3072 6144 12288 98304];
shift=[0,0,0;...
    0.25,0.25,0.25;...
    -0.2,-0.25,-0.3;...
    0.1,-0.2,0.3;...
    0.1,0.01,0.4;...
    0.4,0.3,0.2;...
    0.4,0.3,0.2;...
    0.4,0.3,0.2;...
    0.4,0.3,0.2;...
    0.4,0.3,0.2;...
    0.4,0.3,0.2]';
shift(:,6)=shift(:,6)*-1;
shift(:,7)=shift(:,7)*-0.6;
shift(:,8)=shift(:,8)*-0.2;
shift(:,9)=shift(:,9)*0.2;
shift(:,10)=shift(:,10)*-0.6;
shift(:,11)=shift(:,11);
%% setup vars
l2diff=[];
norml2diff=[];

%% analyze
for a=1:size(shift,2)
    %% generate the gold standard
    Wx=1; %total width of cube
    Wy=1;
    Wz=1;
    [cverts,cfaces]=makeCube(Wx/2,Wy/2,Wz/2);
    cverts=cverts+repmat(shift(:,a),[1 size(cverts,2)]); %not really neccesary
    %compute the convex hull
    [cfaces]=convhull(cverts(1,:),cverts(2,:),cverts(3,:));
    % scface{1}=cfaces;
    % scverts{1}=cverts;
    
    [s,~,~]=validation_fullsampling(Wx,Wy,Wz,shift(:,a),cverts,cfaces,0);
    snorm(a)=norm(abs(s{1}(:)));
    %% analysis loop
    data={};
    kdata3={};
    kdiff={};
    for j=1:length(nfaces)
        %% load in kspace data
        dfn=sprintf('optional/validationdata_64_%g_%g_%g_%ifaces_double.mat',...
            shift(1,a),shift(2,a),shift(3,a),nfaces(j));
        data{j}=importdata(dfn);
        
        %% analyze results
        kdata3{j}=data{j}{1};
        kdiff{j}=s{1}-kdata3{j};
        l2diff(a,j)=norm(abs(kdiff{j}(:)));
        norml2diff(a,j)=l2diff(a,j)/snorm(a);
    end
end

%% graph  l2 error
figure;
plot(nfaces(1:end),l2diff','-o');
xlabel('faces','FontSize',16);
ylabel('|error|','FontSize',16);
title('l^2 norm of error versus number of faces','FontSize',16);
set(gca,'FontSize',16)
legend('0,0,0',...
    '0.25,0.25,0.25',...
    '-0.2,-0.25,-0.3',...
    '0.1,-0.2,0.3',...
    '0.1,0.01,0.4');
%% graph normalized l2 error (edges)
figure;
plot(nfaces(1:end)*3/2,norml2diff','-o','Color',[0.8 0.8 0.8]);hold all;
errorbar(nfaces*3/2, mean(norml2diff,1),std(norml2diff,0,1),'Color',[ 0 0 0],...
    'LineWidth',2);
xlabel('edges','FontSize',14);
% ylabel('s_gs-s_3D|/|s_gs|','FontSize',14,'interpreter','latex');
ylabel('$$\|s_{gs}-s_{3D}\|/\|s_{gs}\| $$','FontSize',14,'interpreter','latex');
% title('normalized l^2 norm of error versus edges','FontSize',14);
set(gca,'FontSize',16);
xlim([0 (nfaces(end)*3/2)*1.1]);
axis square;
print('-painter','-depsc2','-r600',['images/cube_normalizedl2error_edges.eps']);

%% graph normalized l2 error (faces)
figure;
plot(nfaces(1:end),norml2diff','-o','Color',[0.8 0.8 0.8]);hold all;
errorbar(nfaces, mean(norml2diff,1),std(norml2diff,0,1),'Color',[ 0 0 0],...
    'LineWidth',2);
xlabel('faces','FontSize',14);
% ylabel('s_gs-s_3D|/|s_gs|','FontSize',14,'interpreter','latex');
ylabel('$$\|s_{gs}-s_{3D}\|/\|s_{gs}\| $$','FontSize',14,'interpreter','latex');
% title('normalized l^2 norm of error versus edges','FontSize',14);
set(gca,'FontSize',16);
xlim([0 nfaces(end)*1.1]);
axis square;
print('-painter','-depsc2','-r600',['images/cube_normalizedl2error_faces.eps']);

%% plot norm of gold standard
figure;
plot(1:size(shift,2),snorm);