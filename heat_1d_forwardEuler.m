function heat_1d_forwardEuler()
%%-- get initial wall time:
time0=clock();
%--- Simulation cell parameters:
L=1;
Nx=128;
Nx1=Nx+1;
Nx2=Nx+2;
dx = L/Nx1;
D=1; % Diffusion Coefficient
%--- Time integration parameters:
t_end=1.0e3;

nprint = 200;
dtime  = 1.0e-5;

%-- Initialize temperature field & grid:
u0(1:Nx2)=0;
x(1:Nx2)=(0:Nx1)*dx;


%%
u1=u0;
t=0;
istep=0;

%%
while t < t_end
    
    istep=istep+1;
    t=t+dtime;
    u1(2:Nx1) = u0(2:Nx1) + D*dtime/dx^2 * ( u0(3:Nx2)-2*u0(2:Nx1)+u0(1:Nx)) + dtime * rhs(x(2:Nx1),t,D);
    u1(1)  = u_real(x(1),t);
    u1(Nx2)= u_real(x(Nx2),t);
    
    Error = u1 - u_real(x,t);
    u0 = u1;
    
    if((mod(istep,nprint) == 0) ||(istep == 1))
        subplot(1,3,1);
        plot(x,u1);
        time=sprintf('%f',t);
        title(['time' time]);
        axis([0 L -1 1]);
        xlabel('x');
        ylabel('Temperature');
        subplot(1,3,2);
        plot(x,Error);
        time=sprintf('%f',t);
        title(['time' time]);
        xlabel('x');
        ylabel('Temperature Error');
        subplot(1,3,3);
        plot(x,u1,x,u_real(x,t));
        time=sprintf('%f',t);
        title(['time' time]);
        axis([0 L -1 1]);
        xlabel('x');
        ylabel('Temperature');
        getframe(gcf);
    end
end

%%
compute_time=etime(clock(),time0);
fprintf('Compute Time: %5d\n', compute_time);

%%
    function rhs_value=rhs(x,t,D)
        rhs_value = cos(t).*cos(2*pi.*x)+4*pi^2*D*sin(t).*cos(2*pi.*x);
    end
%%
    function RealU=u_real(x,t)
        RealU = sin(t)*cos(2*pi*x);
    end
end
