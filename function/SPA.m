function Ind = SPA(X,R )
%SPA 此处显示有关此函数的摘要
%   此处显示详细说明
[M, N] = size(X);
maxX = max(max(X));
threshold  = maxX/1e6;
SumX = sum(X.^2);
Z = zeros(M,N);
for nn = 1:N
    if sum(X(:,nn))<threshold
        Z(:,nn) = X(:,nn);
    else
        Z(:,nn) = X(:,nn)/sum(X(:,nn));
    end
end
% SPA
Rmat = Z;
Ind =[];
for rr = 1:R
    RmatL2 = sum(Rmat.^2);
    maxR = max(RmatL2);
    index = find(RmatL2==maxR);
    
    
    selectMat = SumX(index);
    indexmax = find( selectMat ==max(selectMat));
    indexrr = index(indexmax(1));
    Ind = [Ind;indexrr];
    Rmat = (eye(M)-Rmat(:,indexrr)*Rmat(:,indexrr)'/norm(Rmat(:,indexrr))^2)*Rmat;
end

end

