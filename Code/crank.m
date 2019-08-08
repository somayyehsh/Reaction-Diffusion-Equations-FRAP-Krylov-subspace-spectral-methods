function y=crank(A,y0,tf,nstep)
dt=tf/nstep;
y=y0;
I=speye(size(A));
A1=I+dt/2*A;
A2=I-dt/2*A;
for i=1:nstep
    b=A1*y;
    y=A2\b;
end
