function y=rk4(A,y0,tf,nstep)
dt=tf/nstep;
y=y0;
for i=1:nstep
    k1=A*y;
    k2=k1+dt/2*A*k1;
    k3=k1+dt/2*A*k2;
    k4=k1+dt*A*k3;
    y=y+dt/6*(k1+2*(k2+k3)+k4);
end
