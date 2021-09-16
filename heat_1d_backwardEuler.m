function heat_1d_backwardEuler()
%%
%-- get initial wall time:
time0=clock();

L=1;
Nx=256;
Nx1=Nx+1;
Nx2=Nx+2;
dx = L/Nx1;
D=1; % Diffusion Coefficient
t_end=1.0e3;

nprint = 200;
dtime  = 1.0e-3;

%% -- Initialize temperature field & grid:
u0(1:Nx2)=0;
x(1:Nx2)=(0:Nx1)*dx;

x_num = Nx;
cfl = D*dtime/dx^2;
A=sparse_A(x_num,cfl);

%%
u1=u0;
t=0;
istep=0;

while t < t_end
    
    istep=istep+1;
    t=t+dtime;
    my_rhs = u0(2:Nx1) + dtime * rhs(x(2:Nx1),t,D);
    my_rhs(1) = my_rhs(1) + cfl * u_real(x(1),t);
    my_rhs(Nx)= my_rhs(Nx)+ cfl * u_real(x(Nx2),t);
    %size(A),size(my_rhs)
    u1(2:Nx1) = A\my_rhs';
    u1(1)  = u_real(x(1),t);
    u1(Nx2)= u_real(x(Nx2),t);
    
    Error = u1 - u_real(x,t);
    u0 = u1;
    
    if((mod(istep,nprint) == 0) ||(istep == 1))
        subplot(1,2,1);
        plot(x,u1);
        time=sprintf('%f',t);
        title(['time' time]);
        axis([0 L -1 1]);
        xlabel('x');
        ylabel('Temperature');
        subplot(1,2,2);
        plot(x,Error);
        time=sprintf('%f',t);
        title(['time' time]);
        xlabel('x');
        ylabel('Temperature Error');
        getframe(gcf);
    end
end

%%
compute_time=etime(clock(),time0);
fprintf('Compute Time: %5d\n', compute_time);

%%
    function A=sparse_A(x_num, cfl)
        A = sparse ( [], [], [], x_num, x_num );
        A(1,1) = 1.0 + 2.*cfl;
        A(1,2) =  - cfl;
        
        for i = 2 : x_num - 1
            A(i,i-1) =           - cfl;
            A(i,i  ) = 1.0 + 2.0 * cfl;
            A(i,i+1) =           - cfl;
        end
        A(x_num,x_num-1)=-cfl;
        A(x_num,x_num) = 1.0 + 2.*cfl;
    end
%%
    function rhs_value=rhs(x,t,D)
        rhs_value = cos(t).*cos(2*pi.*x)+4*pi^2*D*sin(t).*cos(2*pi.*x);
    end
%%
    function RealU=u_real(x,t)
        RealU = sin(t)*cos(2*pi*x);
    end

end
