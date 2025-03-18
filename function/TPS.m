function z = TPS( X,Y,Z,x,y,lambda,dB )
%	Input data: (X,Y)-->Z; vector-form or matrix-form.
%   We use thin-plate splines to interpolate (X,Y)-->Z with smoothing
%   parameter lambda.
%   We want to predict z date given x and y.
%% demo  %%
% % x_grid = 0:2:10;
% % y_grid = 0:2:10;
% % [X,Y] = meshgrid(x_grid,y_grid);
% % Z = 1*sin((X-Y).^2) + 1*randn(size(X));
% % xi_grid = 0:0.2:10;
% % yi_grid = 0:0.2:10;
% % [x,y] = meshgrid(xi_grid,yi_grid);
% % lambda = 1e-2;
% % z = TPS( X,Y,Z,x,y,lambda );
% % mesh(x,y,z);
% % hold on;plot3(X,Y,Z,'o')

%% function begin 
% kernel function

if nargin<7
    dB=0;
end

kf=@(r) r.^2.*log(abs(r)+eps);

% vectorize
Xvec=X(:);
Yvec=Y(:);
Zvec=Z(:);
if dB
    Zvec=log10(Z(:));
end
xvec = x(:);
yvec = y(:);
len=length(Xvec);
% K=zeros(len,len);

dist = sqrt((bsxfun(@minus,Xvec,Xvec')).^2+(bsxfun(@minus,Yvec,Yvec')).^2);
K = kf(dist);
% for ii=1:len
%     for jj=1:len
%         K(ii,jj)=kf(norm([Xvec(ii)-Xvec(jj),Yvec(ii)-Yvec(jj)]));
%     end
% end
P=[ones(size(Xvec)),Xvec,Yvec];
Amat=[K+lambda*eye(size(K,1)),P;P',zeros(3,3)];
bmat=[Zvec;zeros(3,1)];
% xmat=pinv(Amat)*bmat;
xmat = Amat\bmat;
w=xmat(1:len);
a=xmat(end-2:end);
GK = kf( sqrt(  (bsxfun(@minus,xvec,Xvec')).^2+(bsxfun(@minus,yvec,Yvec')).^2  ) );
N = [ones(size(xvec)),xvec,yvec];
zvec = GK*w + N*a;
z = reshape(zvec,size(x));
if dB
    z = 10.^(z);
end

% f=@(x,y) [1 x y]*a+kf(sqrt((x-Xvec).^2+(y-Yvec).^2))'*w;
% for ii=1:size(x,1)
%     for jj=1:size(x,2)
%         z(ii,jj)=f(x(ii,jj),y(ii,jj));
%     end
% end
%% plot
% h=pcolor(XX,YY,f(XX,YY));
% set(h, 'LineStyle','none');
% colorbar;
end

