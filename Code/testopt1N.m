function[absoluteErr,relativeErr,solntime]=testopt1(N,tf,coefs)
dx=2*pi/N;
dy=2*pi/N;
x=dx*(0:N-1);
y=dy*(0:N-1);
[x2,y2]=meshgrid(x,y);
x2=reshape(x2,numel(x2),1);
y2=reshape(y2,numel(y2),1);
eN2=ones(N^2,1);
kb=1.0;         
I0=1.0;   
w1=[0:N-1]/2;
w2=[0:N-1]/2;
ci=1.0;    
if coefs==1
    kon=10^(-0.5);    %e-0.5 (1./sqrt(10))1/s   kon*re^2 ~ D1 %correct reaction dominant
    koff=1.e-1;   % rn=0.5 \mu m
    D1 = 30;
    D2 = 1.e-4;
    rn=0.5;
elseif coefs==2
    kon=10^(3.5);    % 1/s kon*re^2 ~ D1  %correct effective diffusion
    koff=1;  % rn=0.5 \mu m
    D1 = 30;
    D2 = 1.e-4;
    rn=0.5;
elseif coefs==3
    kon=255.0;    % 1/s   kon*re^2 ~ D1  %correct-Irn diffusion dominant
    koff=31.0;  % rn=0.6 \mu m
    D1 = 45;
    D2 = 2.5;
    rn=0.6;
elseif coefs==4
    kon=1.e-2;    % 1/s   kon*re^2 ~ D1  %correct pure diffusion dominant
    koff=10;     % rn=0.5 \mu m
    D1 = 30;
    D2 = 1.e-4;
    rn=0.5;
else
    kon=10^2;    % 1/s   kon*re^2 ~ D1  %correct full model
    koff=10^(-1);
    D1 = 30;
    D2 = 10^(-1);
    rn=0.5;
end
Irn=(2.0*I0./pi.*rn.^2).*exp(-2.0*((x2-pi).^2+(y2-pi).^2)./rn.^2);
B12=(ci.*koff)./(koff+kon)*eN2;
B22=(ci.*kon)./(koff+kon)*eN2;
%L12=koff*speye(N^2);
%L21=kon*speye(N^2);
coef=-kb*Irn;
%L11=-kb*Irn-kon;
%L22=-kb*Irn-koff;

tic
[w1m,w2m]=meshgrid(w1,w2);
w1m=reshape(w1m,numel(w1m),1);
w2m=reshape(w2m,numel(w2m),1);
u11=1.0.*koff*eN2;
u12=1.0.*koff*eN2;
w12m2=(w1m.^2+w2m.^2);
u21=0.5.*(w12m2.*(D1-D2)+(kon-koff))+0.5.*sqrt(w12m2.^2.*(D1-D2).^2+...
    w12m2.*(2.0.*(D1+D2).*(kon+koff)-4.0.*((D1.*koff)+(D2.*kon)))+(kon+koff).^2);
u22=0.5.*(w12m2.*(D1-D2)+(kon-koff))-0.5.*sqrt(w12m2.^2.*(D1-D2).^2+...
    w12m2.*(2.0.*(D1+D2).*(kon+koff)-4.0.*((D1.*koff)+(D2.*kon)))+(kon+koff).^2);
cons=koff*(u22-u21);
v11=(u22./cons);   
v21=(-koff./cons);
v12=(-u21./cons);
v22=(koff./cons);
dt=tf;

u2=B12;
b2=B22;
Irn0=mean(coef);
RM11=Irn0+(((koff.*conj(v11).*v21)+(kon.*conj(v21).*v11)-(((D1.*v11.^2)+(D2.*v21.^2)).*w12m2)-(koff.*v21.^2)-(kon.*v11.^2))./(v11.^2+v21.^2));
M11=Irn0+(((koff.*conj(v12).*v22)+(kon.*conj(v22).*v12)-(((D1.*v12.^2)+(D2.*v22.^2)).*w12m2)-(koff.*v22.^2)-(kon.*v12.^2))./(v12.^2+v22.^2));
Rlam1=RM11;
lam1=M11;

p=coef.*u2-kon*u2+koff*b2;
q=kon*u2-koff*b2+coef.*b2;  %   Q=[p;q];  %Q=L*B
%p2=reshape(p,N,N);
%q2=reshape(q,N,N);
TIrn=reshape(mydct2(reshape(coef,N,N)),N^2,1);
%Tpm=mydct2(p2);
%Tqm=mydct2(q2);   %   Tpq=[Tp;Tq];
%Tpm=reshape(Tpm,N^2,1);
%Tqm=reshape(Tqm,N^2,1);
Tqm=b2(1)*TIrn;
Tqm(1)=Tqm(1)+N*(kon*u2(1)-koff*b2(1));
Tpm=u2(1)*TIrn;
Tpm(1)=Tpm(1)+N*(koff*b2(1)-kon*u2(1));
Rlam2=(u2'*p+b2'*q)/(u2'*u2+b2'*b2);
lam2=Rlam1;
Rc1=((exp(Rlam2*dt)-exp(Rlam1*dt))./(Rlam2-Rlam1));
Rc0=exp(Rlam2*dt)-Rc1.*Rlam2;
c1=((exp(lam2*dt)-exp(lam1*dt))./(lam2-lam1));
c0=exp(lam2*dt)-c1.*lam2;
Tv1=koff*(Rc1.*((u22).*Tpm-koff*Tqm))+...
    koff.*(c1.*((-u21).*Tpm+koff*Tqm));
Tv2=u21.*(Rc1.*((u22).*Tpm-koff*Tqm))+...
    u22.*(c1.*((-u21).*Tpm+koff*Tqm));
Tum1=sum(B12)/N;
Tbm1=sum(B22)/N;
Tv1(1)=Tv1(1)+u11(1)*Rc0(1)*u22(1)*Tum1-u11(1)*Rc0(1)*u12(1)*Tbm1-...
    u12(1)*c0(1)*u21(1)*Tum1+u12(1)*c0(1)*u11(1)*Tbm1;
Tv2(1)=Tv2(1)+u21(1)*Rc0(1)*u22(1)*Tum1-u21(1)*Rc0(1)*u12(1)*Tbm1-...
    u22(1)*c0(1)*u21(1)*Tum1+u22(1)*c0(1)*u11(1)*Tbm1;
Tv1=reshape(Tv1./cons,N,N); 
Tv2=reshape(Tv2./cons,N,N);

v1=myidct2(Tv1);
v2=myidct2(Tv2);
v1=reshape(v1,N^2,1);
v2=reshape(v2,N^2,1);
Sf=[v1;v2];
solntime=toc;

[exact1,exact2]=refsolnN(coefs,64,tf);
exact1=reshape(exact1,64,64);
exact2=reshape(exact2,64,64);

r=64/N;
Sf=reshape(Sf,numel(Sf),1);

Sf1=Sf(1:end/2);
Sf2=Sf(end/2+1:end);
Sf1=reshape(Sf1,N,N);
Sf2=reshape(Sf2,N,N);

norm1=norm(exact1(1:r:64,1:r:64),'inf');
norm2=norm(exact2(1:r:64,1:r:64),'inf');
absoluteErr1=norm(Sf1-exact1(1:r:64,1:r:64),'inf');
absoluteErr2=norm(Sf2-exact2(1:r:64,1:r:64),'inf');
absoluteErr=max([ absoluteErr1 absoluteErr2 ])
relativeErr=absoluteErr/max([ norm1 norm2 ])
