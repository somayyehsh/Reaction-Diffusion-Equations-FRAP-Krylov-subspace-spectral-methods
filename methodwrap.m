function [solntime,absoluteErr,relativeErr]=methodwrap(methfun,N,nsteps,tf,coefs)
% N=64;
% nsteps=1000;
% tf=1;
% coefs=1;
%S=[];
%Sf=0;
dx=2*pi/N;
dy=2*pi/N;
x=dx*(0:N-1);
y=dy*(0:N-1);
[x2,y2]=meshgrid(x,y);
x2=reshape(x2,numel(x2),1);
y2=reshape(y2,numel(y2),1);
%kon=255.0*ones(size(x2));    % 1/s   kon*re^2 ~ D1
%%% SS: don't use vectors for constant coefficients, use scalars instead
%%% also, define ones(N^2,1) once and re-use
eN2=ones(N^2,1);
kb=1.0;         % 1/s, doubt on its value
%kb=reshape(kb,N,N);
I0=1.0;   % micro A, doubt this to have micro or not
%rn=0.5*ones(N^2,1);    % micro m
%x=1.0*1.e-6*ones(size(x2));
%y=1.0*1.e-6*ones(size(y2));
%w1=[0:N/2 -N/2+1:-1];
%w2=[0:N/2 -N/2+1:-1];
%w1=[-N/2+1:N/2];
%I=eye(N^2);
%TIrn=fft2(reshape(Irn,N,N));
% D1 = 45*ones(N^2,1);
% D2 = 2.5*ones(N^2,1);
ci=1.0;    %doubt on unit and value
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
%     kon=10.0*ones(N^2,1);    % 1/s   kon*re^2 ~ D1 %nc reaction dominant
%     koff=1.e-2*ones(N^2,1);
%     D1 = 30*1.e-6*ones(N^2,1); %removed 1.e-6 kang
%     D2 = 1.e-6*ones(N^2,1); %removed 1.e-6  kang
D=-2*eye(N)+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);
D(1,N)=1;
D(N,1)=1;
D=D/dx^2;
D=sparse(D);
Lapl=kron(D,eye(N))+kron(eye(N),D);
%D3=diag(ones(N-1,1),1);
%D3(N,1)=1;
%D3=D3-D3';
%D3=D3/(2*dx);
%D3=sparse(D3);
%Dx=kron(D3,eye(N));
%Dy=kron(eye(N),D3);
%%% SS: because ci, koff, kon now scalars, need vector of 1's here
B12=(ci.*koff)./(koff+kon)*eN2;
B22=(ci.*kon)./(koff+kon)*eN2;
%B0=[B12;B22];
%%% SS: made matrices completely sparse, don't use eye or diag,
%%% use speye or spdiags instead
L12=koff*speye(N^2);
L21=kon*speye(N^2);
L11=-kb*spdiags(Irn,0,N^2,N^2)+D1*Lapl-L21;
L22=-kb*spdiags(Irn,0,N^2,N^2)+D2*Lapl-L12;
A=[ L11 L12; L21 L22 ];
y0=[ B12; B22 ];
tic
y=methfun(A,y0,tf,nsteps);
solntime=toc;


[exact1,exact2]=refsoln(coefs,64,tf);
%exact=exact1+exact2;
exact1=reshape(exact1,64,64);
exact2=reshape(exact2,64,64);
%[exact1,exact2,~]=test64(64,1000,1,coefs);
r=64/N;
%Sf=reshape(Sf,numel(Sf),1);
Sf=y;
%norm(Sf,'inf')
% exact=expm(L*tf)*B0;
% exact1=exact(1:end/2);
% exact2=exact(end/2+1:end);
% exact1=reshape(exact1,N,N);
% exact2=reshape(exact2,N,N);
Sf1=Sf(1:end/2);
Sf2=Sf(end/2+1:end);
Sf1=reshape(Sf1,N,N);
Sf2=reshape(Sf2,N,N);
% exact=exact1+exact2;
% exact=reshape(exact,numel(exact),1);
% absoluteErr=norm(Sf-exact,'inf')
% relativeErr=norm(Sf-exact,'inf')./norm(exact,'inf')
absoluteErr1=norm(Sf1-exact1(1:r:64,1:r:64),'inf');  % max is infinity norm in the whole solution
relativeErr1=norm(Sf1-exact1(1:r:64,1:r:64),'inf')./norm(exact1(1:r:64,1:r:64),'inf');
absoluteErr2=norm(Sf2-exact2(1:r:64,1:r:64),'inf');
relativeErr2=norm(Sf2-exact2(1:r:64,1:r:64),'inf')./norm(exact2(1:r:64,1:r:64),'inf');
absoluteErr=absoluteErr1+absoluteErr2;
relativeErr=relativeErr1+relativeErr2;
