%[solntime,absoluteErr,relativeErr]=methodwrap(methfun,N,nsteps,tf,coefs)
%[exact1,exact2]=refsoln(coefs,N,tf);
%[absoluteErr,relativeErr,solntime]=testopt1(N,tf,coefs)


%cases
% coefs=1;
% coefs=2;
% coefs=3;
% coefs=4;
 coefs=5;

%[absoluteErr,relativeErr,solntime]=testopt1(64,1,coefs)
[solntime,absoluteErr,relativeErr]=methodwrap(@crank,64,1,1,coefs)
[solntime,absoluteErr,relativeErr]=methodwrap(@rk4,64,10000,1,coefs)
[exact1,exact2]=refsoln(1,64,1);
[solntime,absoluteErr,relativeErr]=methodwrap(@eulerfwd,64,100000,1,coefs)

%[absoluteErr,relativeErr,solntime]=testopt1N(64,1,coefs)
[solntime,absoluteErr,relativeErr]=methodwrapN(@crank,64,1,1,coefs)
[solntime,absoluteErr,relativeErr]=methodwrapN(@rk4,64,10000,1,coefs)
[exact1,exact2]=refsolnN(1,64,1);
[solntime,absoluteErr,relativeErr]=methodwrapN(@eulerfwd,64,100000,1,coefs)

%[absoluteErr,relativeErr,solntime]=testopt1(32,1,coefs)
[solntime,absoluteErr,relativeErr]=methodwrap(@crank,32,1,1,coefs)
[solntime,absoluteErr,relativeErr]=methodwrap(@rk4,32,10000,1,coefs)
[exact1,exact2]=refsoln(1,32,1);
[solntime,absoluteErr,relativeErr]=methodwrap(@eulerfwd,32,100000,1,coefs)

%[absoluteErr,relativeErr,solntime]=testopt1N(32,1,coefs)
[solntime,absoluteErr,relativeErr]=methodwrapN(@crank,32,1,1,coefs)
[solntime,absoluteErr,relativeErr]=methodwrapN(@rk4,32,10000,1,coefs)
[exact1,exact2]=refsolnN(1,32,1);
[solntime,absoluteErr,relativeErr]=methodwrapN(@eulerfwd,32,100000,1,coefs)