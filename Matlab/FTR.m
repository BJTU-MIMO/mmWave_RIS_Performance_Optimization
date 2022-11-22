function  Y =  FTRS(K,m,SNR_mean_dB,h,Num)
Size = Num;
KK = K;
mm=m;
SNR_mean=10^(SNR_mean_dB/10);
uu=SNR_mean ;
hh=h;
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
Y=(abs(V));%squared FTR RVs
end