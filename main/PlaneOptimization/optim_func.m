function mi = optim_func(x0,Data3D, hist, alpha_max)

hist = double(hist);
[X_size Y_size Z_size] = size(Data3D);
[x_mesh, y_mesh]=meshgrid(1:Y_size,1:X_size);
norm_al = alpha_max; norm_phi = 2*pi; norm_d = Z_size;
alpha = x0(1)*norm_al;
phi = x0(2)*norm_phi;
d = x0(3)*norm_d;


len = 1.000;
B_ = len * sin(alpha)*sin(phi);
C = len * cos(alpha);
A = len * sin(alpha)*cos(phi);
D = -d * C;
Z_fit=-(y_mesh.*B_+ x_mesh.*A + D)/C;
    
try
   Slice= interp3(single(Data3D), x_mesh, y_mesh,Z_fit);
catch
   Slice= ones(X_size, Y_size);
end

    
    
    

if (isnan(Slice(:))==0)
    %%%%%%%%%% MI
    mi = - similarity(hist,Slice,'NMI');
else
    for i = 1:X_size
        for j = 1:Y_size
            if (isnan(Slice(i,j)) ==1)
                Slice(i,j) = 1;
            end
        end
    end
    mi = - similarity(hist,Slice,'NMI');
end
   
    
    
    
