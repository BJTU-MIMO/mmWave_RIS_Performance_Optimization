clc ; 
clear all;
Size = 1e6;
KK =3;
mm1=5;
mm2=5;
hh=0.5;
uu=log10(10)*10;
%ee为yigama^2
f122=FTR(KK,mm1,uu,hh,Size);
f111=FTR(KK,mm2,uu,hh,Size);
fpro=f122.*f111;
xi = 0.1:2:10.1;
for ii=1:length(xi)
   CDF(ii)=sum(fpro<xi(ii))/Size;
end
plot(xi,CDF,'r*');hold on;
% 下面是mathematic的计算结果：
mathematica=[0.000313764, 0.0729952, 0.217457, 0.385434, 0.543984, 0.675935]
plot(xi,mathematica,'b');hold on;

% xi = linspace(0,5,501);
% f = ksdensity(f1,xi,'function','pdf');
% plot(xi,f);hold on;
% grid on;