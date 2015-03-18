function [v,f]=makeCube(xw,yw,zw,triflag)
if ~exist('triflag','var')
    triflag=false;
end

top=[-1, -1, 1;...
    1,-1,1;...
    1,1,1;...
    -1,1,1];
bottom=top+repmat([0,0,-2],[4,1]);
v=[top;bottom]';
f={[1 2 3 4];[8 7 6 5];[1 5 6 2];[2 6 7 3];[3 7 8 4]; [4 8 5 1]};
v(1,:)=v(1,:)*xw;
v(2,:)=v(2,:)*yw;
v(3,:)=v(3,:)*zw;

%% triangulate
if triflag
    f=convhull(v(1,:),v(2,:),v(3,:));
end
