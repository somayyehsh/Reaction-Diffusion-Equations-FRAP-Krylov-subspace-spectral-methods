function [exact1,exact2]=refsolnN(coefs,N,tf)
% set final time
%T=tf(coefs);
% construct filename to save exact solution in
a1=[ 'refsoln1' num2str(coefs) '_' num2str(N) 'N.mat' ];
a2=[ 'refsoln2' num2str(coefs) '_' num2str(N) 'N.mat' ];
% check if filename already exists
s=dir;
found=false;
for j=1:length(s)
    if strcmp(s(j).name,a1)
        % if file is found, load exact solution from file and exit
        ty1=load(a1);
        exact1=ty1.exact1;
        found=true;
    end
    if strcmp(s(j).name,a2)
        % if file is found, load exact solution from file and exit
        ty2=load(a2);
        exact2=ty2.exact2;
        found=true;
    end
end
if found
    return;
end
N=64;
dx=2*pi/N;
dy=2*pi/N;
x=dx*(0:N-1);
y=dy*(0:N-1);
[x2,y2]=meshgrid(x,y);
x2=reshape(x2,numel(x2),1);
y2=reshape(y2,numel(y2),1);
ci=1.0*ones(N^2,1); 
kb=1.0*ones(N^2,1);
I0=1.0*ones(N^2,1);
if coefs==1
    kon=10^(-0.5)*ones(N^2,1);    %e-0.5 (1./sqrt(10))1/s   kon*re^2 ~ D1 %correct reaction dominant
    koff=1.e-1*ones(N^2,1);   % rn=0.5 \mu m
    D1 = 30*ones(N^2,1);
    D2 = 1.e-4*ones(N^2,1);
    rn=0.5*ones(N^2,1);
elseif coefs==2
    kon=10^(3.5)*ones(N^2,1);    % 1/s kon*re^2 ~ D1  %correct effective diffusion
    koff=1*ones(N^2,1);  % rn=0.5 \mu m
    D1 = 30*ones(N^2,1);
    D2 = 1.e-4*ones(N^2,1);
    rn=0.5*ones(N^2,1);
elseif coefs==3
    kon=255.0*ones(N^2,1);    % 1/s   kon*re^2 ~ D1  %correct-Irn diffusion dominant
    koff=31.0*ones(N^2,1);  % rn=0.6 \mu m
    D1 = 45*ones(N^2,1);
    D2 = 2.5*ones(N^2,1);
    rn=0.6*ones(N^2,1);
elseif coefs==4
    kon=1.e-2*ones(N^2,1);    % 1/s   kon*re^2 ~ D1  %correct pure diffusion dominant
    koff=10*ones(N^2,1);     % rn=0.5 \mu m
    D1 = 30*ones(N^2,1);
    D2 = 1.e-4*ones(N^2,1);
    rn=0.5*ones(N^2,1);
else
    kon=10^2*ones(N^2,1);    % 1/s   kon*re^2 ~ D1  %correct full model
    koff=10^(-1)*ones(N^2,1);
    D1 = 30*ones(N^2,1);
    D2 = 10^(-1)*ones(N^2,1);
    rn=0.5*ones(N^2,1);
end
Irn=(2.0*I0./pi.*rn.^2).*exp(-2.0*((x2-pi).^2+(y2-pi).^2)./rn.^2);
D=-2*eye(N)+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);
%D(1,N)=1;
%D(N,1)=1;
D(1,1)=-1;
D(N,N)=-1;
D=D/dx^2;
D=sparse(D);
Lapl=kron(D,eye(N))+kron(eye(N),D);
D3=diag(ones(N-1,1),1);
D3(N,1)=1;
D3=D3-D3';
D3=D3/(2*dx);
D3=sparse(D3);
Dx=kron(D3,eye(N));
Dy=kron(eye(N),D3);
B12=(ci.*koff)./(koff+kon);
B22=(ci.*kon)./(koff+kon);
B0=[B12;B22];
L11=-diag(kb)*diag(Irn)+diag(D1)*Lapl-diag(kon);
L12=diag(koff)*eye(N^2);
L21=diag(kon)*eye(N^2);
L22=-diag(kb)*diag(Irn)+diag(D2)*Lapl-diag(koff);
L=[L11 L12;L21 L22];
L=sparse(L);
exact64=expm(L*tf)*B0;
exact1=exact64(1:numel(exact64)/2);
exact2=exact64(numel(exact64)/2+1:numel(exact64));
% exact1=reshape(exact1,64,64);
% exact2=reshape(exact2,64,64);
%fname=[ 'refsoln' num2str(coefs) '_' num2str(N) '.mat' ];

%ex=reshape(y,numel(y),1);
        save(a1,'exact1');
        save(a2,'exact2');
% file is not found, so compute exact solution and save it
% (actualy saving is done in myout.m)

% exact1=reshape(exact1,64,64);
% exact2=reshape(exact2,64,64);
% % % A=makemat(coefs,N);
% % % options=odeset('AbsTol',1e-14,'RelTol',1e-14,'Jacobian',@(t,y)fjac(t,y,coefs,A),'OutputFcn',@(t,y,flag)myout(t,y,flag,coefs,N));
% % % y0=initdata(coefs,N);
% % % oda15s(@(t,y)jacobian(t,y,coefs,A),[ 0 T ],y0,options);
% at this point, solution has been computed and saved, so load from
% file and return it
%ty=load(fname);
ty1=load(a1);
ty2=load(a2);
exact1=ty1.exact1;
exact2=ty2.exact2;
