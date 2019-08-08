function Sf=kssmain(w1,w2,kon,koff,D1,D2,eN2,tf,nsteps,B12,B22,N,kb,Irn,L11,L12,L21,L22)
[w1m,w2m]=meshgrid(w1,w2);
w1m=reshape(w1m,numel(w1m),1);
w2m=reshape(w2m,numel(w2m),1);
u11=1.0.*koff*eN2;
u12=1.0.*koff*eN2;
u21=0.5.*((w1m.^2+w2m.^2).*(D1-D2)+(kon-koff))+0.5.*sqrt((w1m.^2+w2m.^2).^2.*(D1-D2).^2+...
    (w1m.^2+w2m.^2).*(2.0.*(D1+D2).*(kon+koff)-4.0.*((D1.*koff)+(D2.*kon)))+(kon+koff).^2);
u22=0.5.*((w1m.^2+w2m.^2).*(D1-D2)+(kon-koff))-0.5.*sqrt((w1m.^2+w2m.^2).^2.*(D1-D2).^2+...
    (w1m.^2+w2m.^2).*(2.0.*(D1+D2).*(kon+koff)-4.0.*((D1.*koff)+(D2.*kon)))+(kon+koff).^2);
cons=u11.*u22-u12.*u21;
v11=(u22./cons);   
v21=(-u12./cons);
v12=(-u21./cons);
v22=(u11./cons);
dt=tf/nsteps;
nstep=1;

for t=dt:dt:tf
    u2=reshape(B12,N,N);
    b2=reshape(B22,N,N);
    Tum=fft2(u2);  %Tu=\hat{u}
    Tbm=fft2(b2);   %  Tub=[Tu;Tb];
    u2=reshape(u2,N^2,1);    
    b2=reshape(b2,N^2,1);
    Tum=reshape(Tum,N^2,1);
    Tbm=reshape(Tbm,N^2,1);
    RM11=-(kb.*mean(Irn))+(((koff.*conj(v11).*v21)+(kon.*conj(v21).*v11)-(((D1.*v11.^2)+(D2.*v21.^2)).*(w1m.^2+w2m.^2))-(koff.*v21.^2)-(kon.*v11.^2))./(v11.^2+v21.^2));
    M11=-(kb.*mean(Irn))+(((koff.*conj(v12).*v22)+(kon.*conj(v22).*v12)-(((D1.*v12.^2)+(D2.*v22.^2)).*(w1m.^2+w2m.^2))-(koff.*v22.^2)-(kon.*v12.^2))./(v12.^2+v22.^2));
    Rlam1=RM11;
    lam1=M11;
    p=L11*u2+L12*b2;
    q=L21*u2+L22*b2;  %   Q=[p;q];  %Q=L*B
    p2=reshape(p,N,N);
    q2=reshape(q,N,N);
    Tpm=fft2(p2);
    Tqm=fft2(q2);   %   Tpq=[Tp;Tq];
    Tpm=reshape(Tpm,N^2,1);
    Tqm=reshape(Tqm,N^2,1);
    Rlam2=(u2'*p+b2'*q)/(u2'*u2+b2'*b2);
    lam2=Rlam1;
    Rc1=((exp(Rlam2*dt)-exp(Rlam1*dt))./(Rlam2-Rlam1));
    Rc0=exp(Rlam2*dt)-Rc1.*Rlam2;
    %RC=[Rc0;Rc1];
    c1=((exp(lam2*dt)-exp(lam1*dt))./(lam2-lam1));
    c0=exp(lam2*dt)-c1.*lam2;
    %C=[c0;c1];
    Tv1=u11.*(Rc0.*((u22).*Tum+(-u12).*Tbm))+...
        u11.*(Rc1.*((u22).*Tpm+(-u12).*Tqm))+...
        u12.*(c0.*((-u21).*Tum+(u11).*Tbm))+...
        u12.*(c1.*((-u21).*Tpm+(u11).*Tqm));
    Tv2=u21.*(Rc0.*((u22).*Tum+(-u12).*Tbm))+...
        u21.*(Rc1.*((u22).*Tpm+(-u12).*Tqm))+...
        u22.*(c0.*((-u21).*Tum+(u11).*Tbm))+...
        u22.*(c1.*((-u21).*Tpm+(u11).*Tqm));
    Tv1=reshape(Tv1./cons,N,N); 
    Tv2=reshape(Tv2./cons,N,N);
    
    v1=ifft2(Tv1);
    v2=ifft2(Tv2);
    v1=reshape(v1,N^2,1);
    v2=reshape(v2,N^2,1);
    vv=[v1;v2];
    B12=v1;
    B22=v2;
    Sf=vv;
    nstep=nstep+1;
end
