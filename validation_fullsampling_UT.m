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
%% compute the data
for a=1%:size(shift,2)
    for j=1%:length(nfaces)
%         nfaces(j)
%         dfn=sprintf('optional/validationdata_64_%g_%g_%g_%ifaces_double.mat',...
%             shift(1,a),shift(2,a),shift(3,a),nfaces(j));
%         if exist(dfn, 'file')
%             display('skipping');
%             continue;
%         end
        
        %% construct a cube
        Wx=1; %total width of cube
        Wy=1;
        Wz=1;
        [cverts,cfaces]=makeCube(Wx/2,Wy/2,Wz/2);
        cverts=cverts+repmat(shift(:,a),[1 size(cverts,2)]);
        %compute the convex hull
        [cfaces]=convhull(cverts(1,:),cverts(2,:),cverts(3,:));
        % scface{1}=cfaces;
        % scverts{1}=cverts;
        
        %% show untesselated cube
%         figure; grid on; drawMesh(cverts',cfaces,'FaceColor','white','FaceAlpha',0.7,'EdgeAlpha',0.25); axis equal; camproj('perspective');
%         campos([0.5,-2,1.5]*3);
%         camtarget([0 0 0]);
%         xlim([-1,1]);ylim([-1,1]);zlim([-1,1]);
%         camproj('perspective');
%         set(gca,'FontUnits','inches','FontSize',0.15);
        % % xlabel('x');
        % % ylabel('y');
        % % zlabel('z');
        
        %% print image
        % http://www.mathworks.com/matlabcentral/newsreader/view_thread/157417
        % print -painter -deps -r600 images/untesscube.eps
        
        %% sample
        % [kdatags,kdata_12_double,k]=validation_fullsampling(Wx,Wy,Wz,cverts,cfaces,0);
        % [kdatags,kdata_12_single,k]=validation_fullsampling(Wx,Wy,Wz,cverts,cfaces,1);
        
        %% tesselate
%         while size(cfaces,1)<nfaces(j)
%             [cverts,cfaces]=tesselate(cverts,cfaces);
%             size(cfaces,1)
%         end
        
        %% show tesselated cube
%         figure; grid on; drawMesh(cverts',cfaces,'FaceColor','white','FaceAlpha',0.7,'EdgeAlpha',0.25); axis equal; camproj('perspective');
%         campos([0.5,-2,1.5]*3);
%         camtarget([0 0 0]);
%         xlim([-1,1]);ylim([-1,1]);zlim([-1,1]);
%         camproj('perspective');
%         set(gca,'FontUnits','inches','FontSize',0.15);
        % xlabel('x');
        % ylabel('y');
        % zlabel('z');
        
        %% print image
        % http://www.mathworks.com/matlabcentral/newsreader/view_thread/157417
        % print -painter -deps -r600 images/tesscube_600.eps
        
        %% sample
        [kdatags,kdata,k]=validation_fullsampling(Wx,Wy,Wz,shift(:,a),cverts,cfaces,0);
        
        %% save results
        save(dfn,'kdata');
        
    end
end