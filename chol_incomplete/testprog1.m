n=500; d=5; maxd=20;
hyppars=[0.25; 0.5; 0.1; 0.2; 0.3; 1; 0];

data=randn(n,d);
dvec=sum(data.*muldiag(data,hyppars(1:d)),2);
[lfact,pind]=chol_incomplete(n,1,maxd,0,data,dvec,hyppars);

covmat=sqexpcf(data,data,hyppars(d+1),hyppars(1:d)*d);
lf2=chol(covmat(pind,pind))';
lf3=lf2(:,1:maxd);
temp=abs(lfact-lf3);
max(max(temp))
figure(1);
imagesc(temp(1:20,:));
temp=abs(lfact(1:20,:)*lfact(1:20,:)'-covmat(pind(1:20),...
					     pind(1:20)));
figure(2);
imagesc(temp);
