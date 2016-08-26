function M=GetCellOptions(further,opt)
% get the parameter to the right of the options
ind=strcmpi(further,opt);
if max(ind)==0 || length(ind(ind==1))>1
    M=Nan;
    return
end

[~,loc]=max(ind);
if loc==length(further)
    M=Nan;
    return
end

M=further{loc+1};

