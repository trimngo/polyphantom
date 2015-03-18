function [kdata,removed]=polyhedralFT(vertices,faces,values,kvec,singleflag,debug)
% PURPOSE
% Computes the continuous Fourier Transform of a set of polyhedral surfaces
%
% INPUTS
% [vertices] a cell array where each element are the vertices of a single
% surface
% [faces] a cell array where each element are the faces of a single surface.
% Each element of the cell array can either be a cell array or a 3xE matrix
% In the cell array case, each element should be a column vector describing
% the connectivity between each vertex in the face. Each face can have a 
% different number of vertices. If each surface is described as a 3xE
% matrix, where E is the number of vertices in the face. The matrix is 
% converted to an equivalent cell array.
% [values] an array containing the intensity value for each surface
% [kvec] a 3xN matrix containing the kspace points to be sampled. N is the
% number of k-space samples 
% [singleflag TRUE: compute using single or double precision. default
% default=FALSE;
% [debug] TRUE: use the slow matlab version of the script. default=FALSE
%
% OUTPUTS
% kdata
% removed: (optional) The [vertices] parameter after removing surfaces with 
% no vertices

kdata=[];
%check input variables
if ~iscell(vertices)
    error('Error: vertices must be a cell array')
end
if ~iscell(faces)
    error('Error: faces must be a cell array')
end
if length(faces)~=length(vertices)
    error('Error:faces and vertices must have the same number of components')
end
if length(faces)~=length(values)
    error('Error: values vector must have the same length as number of components in faces and vertices')
end
if ~exist('singleflag','var')
    singleflag=false;
    display('No single flag specified, defaulting to doubles');
end
if ~exist('debug','var')
    debug=false;
end
%should check to see that faces array are unsigned integers, otherwise
%inputing the vertex array where the faces array should be will cause
%program to crash

%should also optimize so that the single thread at the end isn't running by
%itself, figure out what that single thread is doing.

%remove all components whose vertex array is empty.  Remove corresponding
%value faces
%known issue:
kept=[];
for c=1:length(vertices)
    if ~isempty(vertices{c})
        kept=[kept c];
    end
end
removed=1:length(vertices);
removed(kept)=[]

if isempty(kept)
    kdata=zeros(1,length(kvec));
    return
end

kept=1:length(vertices);

if ~iscell(faces{1})
    display('we have a triangular mesh->convert to cell array');
    for k=1:length(faces)
        cfaces{k}=mat2cell(faces{k},[3],ones(1,size(faces{k},2)));
    end
    faces=cfaces;
end

if debug
    display('debug mode on!');
    kdata=zeros(1,length(kvec));
    for i=1:length(kept)
        kdata=kdata+values(kept(i))*polyhedralFTscript(vertices{kept(i)},faces{kept(i)},kvec);
    end
else
    if singleflag
        if iscell(faces{1})
            %we have a polyhedron
            display('Using singles!');
            kdata=fastmtpolyhedralFTcellpoly_single(vertices(kept),faces(kept),values(kept),kvec);
        else
            display('we have a triangular mesh->convert to cell array');
            kdata=fastmtpolyhedralFTcell(vertices(kept),faces(kept),values(kept),kvec);
        end
    else
        if iscell(faces{1})
            %we have a polyhedron
            tic
            kdata=fastmtpolyhedralFTcellpoly(vertices(kept),faces(kept),values(kept),kvec);
            toc
        else
            display('we have a triangular mesh->convert to cell array');
            tic
            kdata=fastmtpolyhedralFTcell(vertices(kept),faces(kept),values(kept),kvec);
            toc
        end
    end
    
end