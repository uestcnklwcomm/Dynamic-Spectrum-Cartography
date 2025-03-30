function [Si,Ci] = SMF_LL1(I,J,Yob3,L,Rest,SampleIndex)
    M = size(Yob3,1);
    P = zeros(M,I*J);
    for mm = 1:M
        row_selection = SampleIndex(mm);
        P(mm,row_selection) = 1;
    end
    lambda = 1e-8;
%     svalue = 0;
    tic;
    iter = 50;
%     [S_init,C] = NMF_HALS(Yob3,Rest,50);%Initialized with HALS;
    SelectInd = SPA(Yob3,Rest);
    S_init = Yob3(:,SelectInd);
    C_init = (S_init\Yob3)';

    
    S = zeros(I*J,Rest);
    S(SampleIndex,:) = S_init;
    C = C_init;
    %%Projection Gradient
    for ii = 1:iter    
        PS = P*S;
        LS = 1/norm(PS'*PS);
        C = C - rand*LS * (lambda*C + (C*PS' - Yob3') * PS);
        LC = 1/norm(C'*C);
        S = S - rand*LC * (P'*(PS*C'- Yob3)*C); 
        for rr = 1:Rest
            sr = reshape(S(:,rr),[],1);
            SrMat = reshape(sr,I,J);
            [Us,Ss,Vs] = svds(SrMat,L);
            SrMat = Us*Ss*Vs';
            sr = reshape(SrMat,[],1);
            S(:,rr) = sr;
        end
    end
    Ci = C;
    Si = S(SampleIndex,:);
end