function y=eulerfwd(A,y0,tf,nstep)
dt=tf/nstep;
y=y0;
for i=1:nstep
    k1=A*y;
    y=y+dt*k1;
end
