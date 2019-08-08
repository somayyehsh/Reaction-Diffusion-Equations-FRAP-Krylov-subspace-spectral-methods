% This can be used for updating the excution times fro theprevious results,
% where tiem-stepping is done
function[absoluteErr,relativeErr,solntime]=testopt2(N,nsteps,tf,coefs)
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
w1=[0:N/2 -N/2+1:-1];
w2=[0:N/2 -N/2+1:-1];
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
D=-2*eye(N)+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);
D(1,N)=1;
D(N,1)=1;
D=D/dx^2;
D=sparse(D);
Lapl=kron(D,eye(N))+kron(eye(N),D);
B12=(ci.*koff)./(koff+kon)*eN2;
B22=(ci.*kon)./(koff+kon)*eN2;
L12=koff*speye(N^2);
L21=kon*speye(N^2);
L11=-kb*spdiags(Irn,0,N^2,N^2)+D1*Lapl-L21;
L22=-kb*spdiags(Irn,0,N^2,N^2)+D2*Lapl-L12;

tic
Sf=kssmain(w1,w2,kon,koff,D1,D2,eN2,tf,nsteps,B12,B22,N,kb,Irn,L11,L12,L21,L22);
solntime=toc;

[exact1,exact2]=refsoln(coefs,64,tf);
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
absoluteErr=max([ absoluteErr1 absoluteErr2 ]);
relativeErr=absoluteErr/max([ norm1 norm2 ]);
