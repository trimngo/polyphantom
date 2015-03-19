function kdata=polyhedralFTscript(r_v,face_node,kvec)
%r_v: an 3xP (where P is number of points) matrix describing vertices of
%surface
%face_node: a cell array where each element is a vector containing the
%vertex indices describing each face. The cell array should be a vector
%kvec:a 3xN matrix containing the kspace points to be sampled

%replicate the first vertex
face_node_X=cellfun(@(x) [x(:);x(1)],face_node,'UniformOutput',false);
%k is in linear form

%% Parallel computing
% tmp  = ver;
% nCores = any(cellfun(@any,  strfind({tmp.Name}, 'Parallel')));
% if nCores
%     nCores = findResource('scheduler', 'type', 'local');
%     nCores = nCores.ClusterSize;   
% 
% 	if ~matlabpool('size') > 0
%         fprintf(1,'  ');
%         matlabpool('open', nCores)
% 		disp('  Matlab pool opening...');
%     else
%         disp('  Matlab pool already open...');
%     end
% else
%     disp('  Not using parallel computing...');
% end

F=size(face_node_X,2);
%% allocate output
kdata=zeros(1,size(kvec,2));

%% calculate associated parameters
%need to adjust this to work with polyhedron, make them into cell arrays
%

% L=cell(1,F);
% t=cell(1,F);
% %N_f=cell(1,F);
% n=cell(1,F);
% r_c=cell(1,F);

for f=1:F
    E=length(face_node_X{f})-1;
    L=zeros(1,E);
    t=zeros(3,E);
    n=zeros(3,E);
    r_c=zeros(3,E);
    
    for e=1:E
        L_hat=r_v(:,face_node_X{f}(e+1))-r_v(:,face_node_X{f}(e));
        L(e)=norm(L_hat);
        t(:,e)=L_hat/L(e);
        r_c(:,e)=r_v(:,face_node_X{f}(e))+t(:,e)*L(e)/2;
    end
    N_f=cross(t(:,1),t(:,2))/norm(cross(t(:,1),t(:,2)));
    for e=1:E
        n(:,e)=cross(t(:,e),N_f);
    end
    for i=1:size(kvec,2)
        %calculate appropriate face contribution
        k=kvec(:,i);
        if isequal(k,[0 0 0]')
            %volume of polyhedron
            a=[0 0 0]';
            for e=1:E
                a=a+cross(r_v(:,face_node_X{f}(e)),r_v(:,face_node_X{f}(e+1)));
            end
            kdata(i)=kdata(i)+dot(r_v(:,face_node_X{f}(1)),N_f)*norm(dot(N_f,a));
        else
            %check to see if vectors are almost perpenicular, notice that
            %k and N_f are unit length so the length of crossprod is
            %sin of angle between them
            if norm(cross(k,N_f))<1e-6
                a=0;
                %k perpedicular to plane of face
                for e=1:length(face_node_X{f})-1
                    a=a+cross(r_v(:,face_node_X{f}(e)),r_v(:,face_node_X{f}(e+1)));
                end
                %note below that factor of two is gone because expression 
                %for P_f contains (1/2) which cancels out factor of 2
                kdata(i)=kdata(i)-1i*pi*dot(k,N_f)*norm(dot(N_f,a))*exp(-pi*2i*dot(k,r_v(:,face_node_X{f}(1))));
            else
                %the usual contribution
                b=0;
                for e=1:length(face_node_X{f})-1
                    b=b+L(e)*dot(k,n(:,e))*sinc(dot(k,t(:,e))*L(e))*exp(-pi*2i*dot(k,r_c(:,e)));
                end
                kdata(i)=kdata(i)+(dot(k,N_f)/(norm(k)^2-dot(k,N_f)^2))*b;
            end
        end
    end
end

%% multiply with appropriate coefficient
for i=1:size(kvec,2)
    k=kvec(:,i);
    if ~isequal(k,[0 0 0]')
        kdata(i)=-kdata(i)/(2*pi*norm(k))^2;
    else
        kdata(i)=abs(kdata(i))/6;
    end
end
