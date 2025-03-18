function [P,S,Congruence,Factor2_calibrate] = Unifying_Ambiguities(Factor1, Factor2)
%Unifying permutation and scaling uncertainties given a pair of matrices

    if size(Factor1,1) ~= size(Factor2,1)
        error('The row dimensions of factor one and two are not identical.');
    end
    
    R1 = size(Factor1,2);
    R2 = size(Factor2,2);




    Factor1_norm = bsxfun(@rdivide,Factor1,sqrt(dot(Factor1,Factor1,1)));
    Factor2_norm = bsxfun(@rdivide,Factor2,sqrt(dot(Factor2,Factor2,1)));

    C = abs(Factor2_norm'*Factor1_norm);
    
    Congruence = C;
   

    % permutation uncertainty
    P = zeros(R2,R1);

    

    for r = 1:R1
        [Cr,i] = max(C,[],1);
        [~,j] = max(Cr);
        P(i(j),j) = 1;   
        C(i(j),:) = 0;
        C(:,j) = 0;
    end

    


     
    Factor2_permute = Factor2 * P;


    S = diag(conj(dot(Factor1,Factor2_permute,1)./dot(Factor2_permute,Factor2_permute,1)));
    S(~isfinite(S)) = 1;
    Factor2_calibrate = Factor2_permute * S;
    
end