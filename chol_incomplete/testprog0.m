n=500; d=5; maxd=20;
hyppars=[0.25; 1; 0];

data=randn(n,d);
dvec=sum(data.*data,2);
[lfact,pind]=chol_incomplete(n,0,maxd,0,data,dvec,hyppars);

covmat=radialcf(data,data,hyppars(2),hyppars(1)*d);
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
