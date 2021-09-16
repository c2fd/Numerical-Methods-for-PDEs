function heat_2d_forwardEuler()
%-- get initial wall time:
time0=clock();
%--- Simulation cell parameters:
L=1;
Lx=L;Ly=L;
assert(Lx==Ly);
Nx=128;Ny=128;
Nx1=Nx+1;Ny1=Ny+1;
Nx2=Nx+2;Ny2=Ny+2;
dx = Lx/Nx1; dy=Ly/Ny1; 
D=1; % Diffusion Coefficient
%--- Time integration parameters:
t_end=1.0e3;

nprint = 200;
dtime  = 1.0e-6;

%-- Initialize temperature field & grid:
u0(1:Nx2,1:Ny2)=0;
x(1:Nx2)=(0:Nx1)*dx;
y(1:Ny2)=(0:Ny1)*dy;

%%ncount=0;
u1=u0;
t=0;
istep=0;

while t < t_end

	istep=istep+1;
	t=t+dtime;
	u1(2:Nx1,2:Ny1) = u0(2:Nx1,2:Ny1) + D*dtime *  ...
		( + ( u0(3:Nx2,2:Ny1)-2*u0(2:Nx1,2:Ny1)+u0(1:Nx,2:Ny1))/dx^2  ...
		  + ( u0(2:Nx1,3:Ny2)-2*u0(2:Nx1,2:Ny1)+u0(2:Nx1,1:Ny))/dy^2  ...
		)	+ dtime * rhs(x(2:Nx1),y(2:Ny1),t,D);
	u1(1,1:Ny2)  = u_real(x(1)  ,y(1:Ny2),t);
	u1(Nx2,1:Ny2)= u_real(x(Nx2),y(1:Ny2),t);
	u1(1:Nx2,1)  = u_real(x(1:Nx2),y(1),t);
	u1(1:Nx2,Ny2)= u_real(x(1:Nx2),y(Ny2),t);

	Error = abs(u1 - u_real(x,y,t));
	u0 = u1;

    if((mod(istep,nprint) == 0) ||(istep == 1))
        norm(Error,2)
        subplot(1,2,1);
        %pcolor(u1),colorbar, shading interp,axis('off'), axis('equal'), title('Solution');
        imagesc(u1),colorbar,title('Solution');
        time=sprintf('%f',t);
        title(['time' time]);
        xlabel('x');
        ylabel('Temperature');
        subplot(1,2,2);
        %pcolor(Error),colorbar, shading interp,axis('off'), axis('equal'), title('Error');
        imagesc(Error),colorbar,title('Error');
        time=sprintf('%f',t);
        title(['time' time]);
        xlabel('x');
        ylabel('Temperature Error');
        getframe(gcf);
    end
end

compute_time=etime(clock(),time0);
fprintf('Compute Time: %5d\n', compute_time);


function rhs_value=rhs(x,y,t,D)
	%x_ones = ones(size(x));
	%y_ones = ones(size(y));
	%rhs_value = cos(t)*cos(2*pi.*x)'*cos(2*pi.*y)+4*pi^2*D*sin(t)*(cos(2*pi*x)'*y_ones+x_ones'*cos(2*pi*y));
	 rhs_value = cos(t)*cos(2*pi.*x)'*cos(2*pi.*y)+8*pi^2*D*sin(t)*(cos(2*pi*x)'*cos(2*pi*y));
end

function RealU=u_real(x,y,t)
	RealU = sin(t)*cos(2*pi*x)'*cos(2*pi*y);
end
%movie(gcf,mov,10,3)
end
