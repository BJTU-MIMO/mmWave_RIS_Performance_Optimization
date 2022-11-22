clc ; 
clear all;

Size = 1e6;
KK = 1;
%KK=15;
mm=5;
uu=12 ;
hh=0.5;
ee=uu/(2*(KK+1));

syms V1 V2;
S=solve((V1^2+V2^2)/(2*ee)==KK,(2*V1*V2)/(V1^2+V2^2)==hh);
S=[S.V1 S.V2];
V1=S(3,1);
V2=S(3,2);
V1=double(V1);
V2=double(V2);

d1=unifrnd(0,2*pi,[1,Size]);
d2=unifrnd(0,2*pi,[1,Size]);

X=normrnd(0,sqrt(ee),[1,Size]);
Y=normrnd(0,sqrt(ee),[1,Size]);
aa=gamrnd(mm,1/mm,[1,Size]);
V=V1*(aa.^0.5).*exp(i*d1)+V2*(aa.^0.5).*exp(i*d2)+X+i*Y;
f1=(abs(V)).^2;%squared FTR RVs
f2=abs(V);%FTR RVs

xi = 0:1:10;
for i=1:length(xi)
   CDF(i)=sum(f2<xi(i))/Size;
end
plot(xi,CDF,'r*');hold on;
% 下面是mathematic的计算结果：
mathematica=[0., 0.0690291, 0.256668, 0.503385, 0.729581, 0.882408, 0.95971, 0.989214, 0.997755, 0.999638, 0.999955]
plot(xi,mathematica,'b');hold on;

% xi = linspace(0,5,501);
% f = ksdensity(f1,xi,'function','pdf');
% plot(xi,f);hold on;
% grid on;