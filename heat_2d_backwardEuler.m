function heat_2d_backwardEuler()
%%
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
dtime  = 1.0e-2;

%-- Initialize temperature field & grid:
u0(1:Nx2,1:Ny2)=0;
x(1:Nx2)=(0:Nx1)*dx;
y(1:Ny2)=(0:Ny1)*dy;

cflx = D*dtime/dx^2;
cfly = D*dtime/dy^2;

% A = full(sparse_A(5,5,1,1))
% A = full(rhs_pad(5,5,1,2,1))

A = sparse_A(Nx,Ny,cflx,cfly);

%% ncount=0;
u1=u0;
%Error(1:Nx2,1:Ny2)=0;
t=0;
istep=0;

%%
while t < t_end
    
    istep=istep+1;
    t=t+dtime;
    my_rhs = u0(2:Nx1,2:Ny1) + dtime * rhs(x(2:Nx1),y(2:Ny1),t,D) + rhs_pad(Nx,Ny,cflx,cfly,t);
    my_rhs = reshape(my_rhs',1,Nx*Ny);
    
    vec_u1 = A\my_rhs';
    u1(2:Nx1,2:Ny1) = reshape(vec_u1,Ny,Nx)';
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
        %imagesc(u1),colorbar,title('Solution');
        surf(x,y,u1)
        time=sprintf('%f',t);
        title(['time' time]);
        xlabel('x');
        ylabel('Temperature');
        subplot(1,2,2);
        %pcolor(Error),colorbar, shading interp,axis('off'), axis('equal'), title('Error');
        %imagesc(Error),colorbar,title('Error');
        surf(x,y,Error)
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
    function rhs_value=rhs(x,y,t,D)
        %x_ones = ones(size(x));
        %y_ones = ones(size(y));
        %rhs_value = cos(t)*cos(2*pi.*x)'*cos(2*pi.*y)+4*pi^2*D*sin(t)*(cos(2*pi*x)'*y_ones+x_ones'*cos(2*pi*y));
        rhs_value = cos(t)*cos(2*pi.*x)'*cos(2*pi.*y)+8*pi^2*D*sin(t)*(cos(2*pi*x)'*cos(2*pi*y));
    end
%%
    function RealU=u_real(x,y,t)
        RealU = sin(t)*cos(2*pi*x)'*cos(2*pi*y);
    end

% function id = id_2_to_1(idx,idy,Nx,Ny)
% 	assert(idy>=1);
% 	id = (idy-1)*Nx+idx;
% end
%%
    function A=sparse_A(Nx,Ny,cflx,cfly)
        x_num = Nx * Ny;
        A = sparse ( [], [], [], x_num, x_num );
        
        A(1,1) = 1.0 + 2.*(cflx+cfly);
        A(1,1+1) =-cflx;
        A(1,1+Nx)=-cfly;
        
        A(Nx,Nx) = 1.0 + 2.*(cflx+cfly);
        A(Nx,Nx-1)=-cflx;
        A(Nx,Nx+Nx)=-cfly;
        
        A(Nx*(Ny-1)+1,Nx*(Ny-1)+1) = 1.0 + 2.*(cflx+cfly);
        A(Nx*(Ny-1)+1,Nx*(Ny-1)+2) =-cflx;
        A(Nx*(Ny-1)+1,Nx*(Ny-2)+1) =-cfly;
        
        A(Nx*Ny,Nx*Ny) = 1.0 + 2.*(cflx+cfly);
        A(Nx*Ny,Nx*Ny-1) = -cflx;
        A(Nx*Ny,Nx*(Ny-1))= -cfly;
        
        for j = 2 : Ny - 1
            A(j,j-1) = - cflx;
            A(j,j  ) = 1.0 + 2.0 * (cflx + cfly);
            A(j,j+1) = - cflx;
            A(j,j+Nx)= - cfly;
            
            A(Nx*(Ny-1)+j,Nx*(Ny-1)+j-1) = -cflx;
            A(Nx*(Ny-1)+j,Nx*(Ny-1)+j)   = 1.0 + 2.0*(cflx+cfly);
            A(Nx*(Ny-1)+j,Nx*(Ny-1)+j+1) = -cflx;
            A(Nx*(Ny-1)+j,Nx*(Ny-2)+j)   = -cfly;
        end
        
        for i= 2: Nx - 1
            A(Ny*(i-1)+1,Ny*(i-1)+2)    = -cflx;
            A(Ny*(i-1)+1,Ny*(i-1)+1)    = 1.0 + 2.0 * (cflx+cfly);
            A(Ny*(i-1)+1,Ny*(i-0)+1)    = -cfly;
            A(Ny*(i-1)+1,Ny*(i-2)+1)    = -cfly;
            
            A(Ny*i,Ny*i-1)  = -cflx;
            A(Ny*i,Ny*i)    = 1.0 + 2.0*(cflx+cfly);
            A(Ny*i,Ny*(i+1))= -cfly;
            A(Ny*i,Ny*(i-1))= -cfly;
        end
        
        for i = 2:Nx -1
            for j= 2:Ny-1
                A((i-1)*Ny+j,(i-1)*Ny+j)    = 1. + 2.*(cflx+cfly);
                A((i-1)*Ny+j,(i-2)*Ny+j)    = -cfly;
                A((i-1)*Ny+j,(i-0)*Ny+j)    = -cfly;
                A((i-1)*Ny+j,(i-1)*Ny+j-1)  = -cflx;
                A((i-1)*Ny+j,(i-1)*Ny+j+1)  = -cflx;
            end
        end
        
    end
%%
    function rhs_add_value=rhs_pad(Nx,Ny,cflx,cfly,t)
        rhs_add_value = zeros(Nx,Ny);
        
        rhs_add_value(1,1) 	= cflx * u_real(x(2),y(1),t)   		+ cfly * u_real(x(1),y(2),t);
        rhs_add_value(Nx,1)	= cflx * u_real(x(Nx+1),y(1),t)		+ cfly * u_real(x(Nx+2),y(2),t);
        rhs_add_value(1,Ny) = cflx * u_real(x(2),y(Ny+2),t)     + cfly * u_real(x(1),y(Ny+1),t);
        rhs_add_value(Nx,Ny)= cflx * u_real(x(Nx+1),y(Ny+2),t)	+ cfly * u_real(x(Nx+2),y(Ny+1),t);
        
        rhs_add_value(1,2:Ny-1) = cfly * u_real(x(1),y(3:Ny),t)';
        rhs_add_value(Nx,2:Ny-1)= cfly * u_real(x(Nx+2),y(3:Ny),t)';
        
        rhs_add_value(2:Nx-1,1) = cflx * u_real(x(3:Nx),y(1),t);
        rhs_add_value(2:Nx-1,Ny)= cflx * u_real(x(3:Nx),y(Ny+2),t);
        
        % for i=2:Ny-1
        %     rhs_add_value(1,i) = cfly * u_real(x(1),y(i+1),t);
        %     rhs_add_value(Nx,i)= cfly * u_real(x(Nx+2),y(i+1),t);
        % end
        
        % for j=2:Nx-1
        %     rhs_add_value(j,1) = cflx * u_real(x(j+1),y(1),t);
        %     rhs_add_value(j,Ny)= cflx * u_real(x(j+1),y(Ny+2),t);
        % end
    end

end
