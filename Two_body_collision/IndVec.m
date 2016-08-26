
function v=IndVec(A)
iUnderConsideration = 1;

[q,r] = qr(A);
rInd = [];
for j = 1:size(r,2),
    if(r(iUnderConsideration,j) ~= 0)
        rInd = [rInd r(:,j)];
        iUnderConsideration = iUnderConsideration + 1;
    end
    if(iUnderConsideration > size(r,1))
        break;
    end
end


end