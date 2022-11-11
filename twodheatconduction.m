%Alexander Bazzi
%Center for Advanced Turbomachinery and Energy Research
%Propulsion and Energy Research Laboratory
%10/6/2016

%% 2D STEADY-STATE HEAT DIFFUSION WITH SOURCES %%

clear all;

%% variable declarations
nx = 30;
ny = 30;
dy = 0.1;
dx = dy;
tol = 10^-3;
maxIt = 5000;
guess = 300;
temp = 10^5;

k = 1000; %thermal conductivity of material
t = 0.01; %thickness of plate
areaX = t*dx;
areaY = t*dy;
q = 5*10^7; %heat flux (heat generation) [W/m^2]

count = 1;

%% boundary conditions
N_TEMP = 100;
E_TEMP = 100;
W_TEMP = 100;
S_TEMP = 100;

N_FLUX = 100;
E_FLUX = 0;
W_FLUX = 0;
S_FLUX = 0;


%% matrix initializations
T(nx,ny) = guess;
oldT(nx,ny) = guess - 100;

T(1,:) = N_TEMP;
T(nx,:) = S_TEMP;
T(:,1) = W_TEMP;
T(:,ny) = E_TEMP;
disp(T);

W_COEFF = k*areaX/dx;
E_COEFF = W_COEFF;
N_COEFF = W_COEFF;
S_COEFF = W_COEFF;

%% iterating subroutine
for i = 1:maxIt
    
    %% west boundary
    for i = 1:ny %west boundary
        if W_FLUX == 0
            break;
        end
        w = 0;
        e = E_COEFF;
        s = S_COEFF;
        n = N_COEFF;
        source = W_FLUX*areaX;
        if i == 1 && source ~= 0 %north/west corner
            n = 0;
            T(i,1) = temperature2d(i, 1, T, 'nw', w, e, s, n, source);
        elseif i == ny && source ~= 0 %south/west corner
            s = 0;
            T(i,1) = temperature2d(i, 1, T, 'sw', w, e, s, n, source);
        elseif i ~= 1 && i ~= ny %along west boundary
            T(i,1) = temperature2d(i, 1, T, 'w', w, e, s, n, source);
        end
    end
    
    %% east boundary
    for i = 1:ny
        if E_FLUX == 0
            break;
        end
        w = W_COEFF;
        e = 0;
        s = S_COEFF;
        n = N_COEFF;
        source = E_FLUX*areaX;
        if i == 1 %north/east corner
            n = 0;
            T(i,nx) = temperature2d(i, nx, T, 'ne', w, e, s, n, source);
        elseif i == ny %south/east corner
            s = 0;
            T(i,nx) = temperature2d(i, nx, T, 'se', w, e, s, n, source);
        elseif i ~= 1 && i ~= ny %along east boundary
            T(i,nx) = temperature2d(i, nx, T, 'e', w, e, s, n, source);
        end
    end
    
    %% north boundary
    for j = 1:nx
        if N_FLUX == 0
            break;
        end
        w = W_COEFF;
        e = E_COEFF;
        s = S_COEFF;
        n = 0;
        source = N_FLUX*areaY;
        if j == 1 %north/west corner
            w = 0;
            T(1,j) = temperature2d(1, j, T, 'nw', w, e, s, n, source);
        elseif j == nx %north/east corner
            e = 0;
            T(1,j) = temperature2d(1, j, T, 'ne', w, e, s, n, source);
        elseif j ~= 1 && j ~= ny %along north boundary
            T(1,j) = temperature2d(1, j, T, 'n', w, e, s, n, source);
        end
    end
    
    %% south boundary
    for j = 1:nx
        if S_FLUX == 0
            break;
        end
        w = W_COEFF;
        e = E_COEFF;
        s = 0;
        n = N_COEFF;
        source = S_FLUX*areaY;
        if j == 1 %south/west corner with generation
            w = 0;
            T(ny,j) = temperature2d(ny, j, T, 'sw', w, e, s, n, source);
        elseif j == nx %south/east corner with generation
            e = 0;
            T(ny,j) = temperature2d(ny, j, T, 'se', w, e, s, n, source);
        elseif j ~= 1 && j ~= ny %along south boundary
            T(ny,j) = temperature2d(ny, j, T, 's', w, e, s, n, source);
        end
    end
    
    %% interior nodes
    for i = 2:ny - 1
        for j = 2:nx - 1
            w = W_COEFF;
            e = E_COEFF;
            s = S_COEFF;
            n = N_COEFF;
            source = 0;
           %{
            residual(count, resCount) = abs(temperature2d(i, j, T, 'interior', w, e, s, n, source) - T(i,j));
            resCount = resCount + 1;
            T(i,j) = temperature2d(i, j, T, 'interior', w, e, s, n, source);
             itConv(count + 1) = abs(T(i,j) - oldT(i,j));
            
            if itConv(count + 1) < temp
                itConv(count + 1) = temp;
            end
            
            temp = itConv(count + 1);
            %}
            
            residual(i-1,j-1) = abs(temperature2d(i, j, T, 'interior', w, e, s, n, source) - T(i,j));
            T(i,j) = temperature2d(i, j, T, 'interior', w, e, s, n, source);
            
            itConv(i-1,j-1) = abs(T(i,j) - oldT(i,j));
            oldT(i,j) = T(i,j);

        end
    end
    
    res(count) = norm(residual, 2);
    conv(count) = norm(itConv, 2);
    
    if mod(count,80) == 0
        fprintf('residual:%5.2e\n', res(count));
        fprintf('convergence:%5.2e\n', conv(count));
    end
    
    if res(count) < tol && conv(count) < tol  
        break;
    end
    
    count = count + 1;
end

%% visualizing the results
contourplot(T, 'jet', min(min(T)), max(max(T)));
meshplot(T, nx, ny);