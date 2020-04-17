function [ravg,rstd] = MOToptim_ForceCalculation_2(det)
    PLOTRESULTS=0;
    if(nargin==0)
        det=-50;
        PLOTRESULTS=1;
    end

    %---------------------------------------------------------------------%
    %---------------------------------------------------------------------%
    
    %--------------------   6 BEAMS   ------------------------------------%
%     BEAMS = [
%         -1	0   0
%         +1  0   0
%         0   -1  0
%         0   +1  0
%         0   0   -1
%         0   0   +1
%     ];
%     LHCP_Beam = [1 1 1 1 0 0];
%     RHCP_Beam = [0 0 0 0 1 1];
%     so = 160*[1 1 1 1 1 1];
%     GB = 0.046*[0.5,0.5,-1.0];

    %--------------------   6 BEAMS - FOUNTAIN GEOMETRY   ----------------%
%     BEAMS = [
%         -1      0           +sqrt(2)/2
%         +1/2    -sqrt(3)/2  +sqrt(2)/2
%         +1/2	+sqrt(3)/2  +sqrt(2)/2
%         +1      0           -sqrt(2)/2
%         -1/2    +sqrt(3)/2  -sqrt(2)/2
%         -1/2	-sqrt(3)/2  -sqrt(2)/2
%     ];
%     LHCP_Beam = [1 1 1 1 1 1];
%     RHCP_Beam = [0 0 0 0 0 0];
%     so = 160*[1 1 1 1 1 1];
%     GB = 0.046*[-1,0.5,0.5];
    
    %------------   6 BEAMS - FOUNTAIN GEOMETRY (30 degrees)   -----------%
%     BEAMS = [
%         -0.866025,-0.500000,+0.707107
%         +0.500000,-0.000000,+0.707107
%         +1.000000,+0.866025,+0.707107
%         +0.866025,+0.500000,-0.707107
%         -0.500000,+0.000000,-0.707107
%         -1.000000,-0.866025,-0.707107
%     ];
%     LHCP_Beam = [1 1 1 1 1 1];
%     RHCP_Beam = [0 0 0 0 0 0];
%     so = 160*[1 1 1 1 1 1];
%     GB = 0.046*[-1,0.5,0.5];
    
    %--------------------   5 BEAMS   ------------------------------------%
%     BEAMS = [
%         -1	0   0
%         +1  0   0
%         0   -1  0
%         0   +1  0
%         0   0   +1
%     ];
%     LHCP_Beam = [1 1 1 1 0];
%     RHCP_Beam = [0 0 0 0 1];
%     so = 160*[1 1 1 1 1];
%     GB = 0.046*[0.5,0.5,-1];

    %--------------------   4 BEAMS   ------------------------------------%
%     BEAMS = [
%         -1	0  +1
%         +1  0  +1
%         0   -1  0
%         0   +1  0
%     ];
%     LHCP_Beam = [1 1 0 0];
%     RHCP_Beam = [0 0 1 1];
%     so = 160*[1 1 1 1];
%     GB = 0.046*[0.5,-1,0.5];
    %--------------------   3 BEAMS   ------------------------------------%
%     BEAMS = [
%         -1      0           sqrt(2)/2
%         +1/2    -sqrt(3)/2  sqrt(2)/2
%         +1/2	+sqrt(3)/2  sqrt(2)/2
%     ];
%     LHCP_Beam = [1 1 1];
%     RHCP_Beam = [0 0 0];
%     so = [160 160 160];   
%     GB = 0.046*[-1,0.5,0.5];
    
    %--------------------   3 BEAMS   ------------------------------------%
    BEAMS = [
        -1      0           0.5*sqrt(2)/2
        +1/2    -sqrt(3)/2  0.5*sqrt(2)/2
        +1/2	+sqrt(3)/2  0.5*sqrt(2)/2
    ];
    LHCP_Beam = [1 1 1];
    RHCP_Beam = [0 0 0];
    so = [160 160 160];   
    GB = 0.046*[0.5,0.5,-1];
    %---------------------------------------------------------------------%
    %---------------------------------------------------------------------%
    
    g = [0,0,-9.8];
    
    h = 6.62607004e-34;
    mub = 9.274009994e-24;
    e = 1.60217662e-19;                 %#ok
    kb = 1.38064852e-23;                %#ok
    c = 2.99792458e8;                   %#ok

    m = (1.660539040e-27)*163.9291748;
    Gamma = 136e3;
    lambda = 626e-9;
    Jgnd = 8;
    Jexc = 9;                           %#ok
    gjg = 1.24;
    gje = 1.29;
    
    detuning = det*Gamma;
    
    mjg = -Jgnd;
    
    %---------------------------------------------------------------------%
    %---------------------------------------------------------------------%
    
    
    LHCP_Beam = repmat(LHCP_Beam',[1,1,1,3]);
    RHCP_Beam = repmat(RHCP_Beam',[1,1,1,3]);
    g = reshape(g,[1 1 1 3]);


    BEAMS = normalize(permute(BEAMS,[1 3 4 2]));
    rvb = (rand(size(BEAMS,1),1,1,3)-0.5)*2;
    %x,y for beams
    bX = rvb - BEAMS.*repmat(sum(BEAMS.*rvb,4),[1,1,1,3]); %Gram-Schmidt
    bX = normalize(bX);
    bY = CrossProduct(BEAMS,bX);

    % LHCP
    PoL = (1/sqrt(2))*(bX + 1i*bY);
    % RHCP
    PoR = (1/sqrt(2))*(bX - 1i*bY);
    BeamPol = normalize(LHCP_Beam.*PoL + RHCP_Beam.*PoR);

    xyzrange = 30e-3;
    [x,y,z]=meshgrid(linspace(-xyzrange,xyzrange,100));
    
    F = CalculateForce(x,y,z);
    F = F/m;
    divF = divergence(x,y,z,F(:,:,:,1),F(:,:,:,2),F(:,:,:,3));
    absF = sqrt(sum(F.^2,4)).*double(divF<0) + 1000*double(divF>=0);
    
    
    if(PLOTRESULTS)
        figure;
        ISOSURF = isosurface(x*1e3,y*1e3,z*1e3,absF,9.8*0.3);
        ISOSURFX = ISOSURF;
        ISOSURFY = ISOSURF;
        ISOSURFZ = ISOSURF;
        ISOSURFX.vertices(:,1) = -xyzrange*1e3;
        ISOSURFY.vertices(:,2) = +xyzrange*1e3;
        ISOSURFZ.vertices(:,3) = -xyzrange*1e3;
        ph = patch(ISOSURF);
        isonormals(x*1e3,y*1e3,z*1e3,absF,ph);
        ph.FaceColor = [1 0 0];
        ph.EdgeColor = 'none';
        hold on;
        phx = patch(ISOSURFX);
        phy = patch(ISOSURFY);
        phz = patch(ISOSURFZ);
        hold off;
        phx.FaceColor = 0.9*[1 1 1];
        phy.FaceColor = 0.9*[1 1 1];
        phz.FaceColor = 0.9*[1 1 1];
        phx.EdgeColor = 'none';
        phy.EdgeColor = 'none';
        phz.EdgeColor = 'none';
        camlight;
        lighting gouraud;
        daspect([1 1 1]);
        axis image;
        xlim(1e3*xyzrange*[-1 1]);
        ylim(1e3*xyzrange*[-1 1]);
        zlim(1e3*xyzrange*[-1 1]);
        xlabel('x (mm)');
        ylabel('y (mm)');
        zlabel('z (mm)');
        view([45 30]);
        grid on;
        title(sprintf('$\\delta = %+.2f\\ \\Gamma$',detuning/Gamma),'Interpreter','Latex','fontsize',12);    


        figure;
        slch = slice(x*1e3,y*1e3,z*1e3,absF,0,0,0);
        title('$\left| \vec{F} \right| /m$','Interpreter','Latex','fontsize',12);
        set(slch,'Linestyle','none');
        set(gca,'Clim',2*[0 9.8]);
        colorbar;
        colormap(jet);
        axis equal;
        xlabel('x (mm)');
        ylabel('y (mm)');
        zlabel('z (mm)');
        view(3);
        rotate3d on;

        figure;
        subplot(1,3,1);
            slch = slice(x*1e3,y*1e3,z*1e3,F(:,:,:,1),0,0,0);
            set(slch,'Linestyle','none');
            colorbar;
            colormap(jet);
            axis equal;
            xlabel('x (mm)');
            ylabel('y (mm)');
            zlabel('z (mm)');
            view([0,0]);
            rotate3d on;
            title('$\vec{F}\cdot\hat{x}/m\ (m/s^2)$','Interpreter','Latex','Fontsize',12);
        subplot(1,3,2);
            slch = slice(x*1e3,y*1e3,z*1e3,F(:,:,:,2),0,0,0);
            set(slch,'Linestyle','none');
            colorbar;
            colormap(jet);
            axis equal;
            xlabel('x (mm)');
            ylabel('y (mm)');
            zlabel('z (mm)');
            view([90,0]);
            rotate3d on;
            title('$\vec{F}\cdot\hat{y}/m\ (m/s^2)$','Interpreter','Latex','Fontsize',12);
        subplot(1,3,3);
            slch = slice(x*1e3,y*1e3,z*1e3,F(:,:,:,3),0,0,0);
            set(slch,'Linestyle','none');
            colorbar;
            colormap(jet);
            axis equal;
            xlabel('x (mm)');
            ylabel('y (mm)');
            zlabel('z (mm)');
            view([0,0]);
            rotate3d on;
            title('$\vec{F}\cdot\hat{z}/m\ (m/s^2)$','Interpreter','Latex','Fontsize',12);
    end


    BOOL = double(absF<9.8*0.7);
    BOOL = BOOL/sum(BOOL(:));
    
    xavg = sum(sum(sum(x.*BOOL)));
    yavg = sum(sum(sum(y.*BOOL)));
    zavg = sum(sum(sum(z.*BOOL)));
	
    x2avg = sum(sum(sum(x.^2.*BOOL)));
    y2avg = sum(sum(sum(y.^2.*BOOL)));
    z2avg = sum(sum(sum(z.^2.*BOOL)));
    
    xstd = sqrt(x2avg - xavg^2);
    ystd = sqrt(y2avg - yavg^2);
    zstd = sqrt(z2avg - zavg^2);
    
    ravg = [xavg,yavg,zavg];
    rstd = [xstd,ystd,zstd];
    
    1;   
    
    function B = MagneticField(x,y,z,GB)
        B = zeros(size(x,1),size(x,2),size(x,3),3);
        B(:,:,:,1) = GB(1)*x;
        B(:,:,:,2) = GB(2)*y;
        B(:,:,:,3) = GB(3)*z;
    end

    function nV = normalize(V)
        nV = V./repmat(sqrt(sum(abs(V).^2,4)),[1,1,1,3]);
    end

    function V = CrossProduct(A,B)
        V(:,:,:,1) = A(:,:,:,2).*B(:,:,:,3) - A(:,:,:,3).*B(:,:,:,2);
        V(:,:,:,2) = A(:,:,:,3).*B(:,:,:,1) - A(:,:,:,1).*B(:,:,:,3);
        V(:,:,:,3) = A(:,:,:,1).*B(:,:,:,2) - A(:,:,:,2).*B(:,:,:,1);
    end

    function F = CalculateForce(X,Y,Z)
        B = MagneticField(X,Y,Z,GB);

        nB = normalize(B);
        nB(isnan(nB)) = 0;

        rvc = (rand(size(X,1),size(X,2),size(X,3),3)-0.5)*2;
        pX = rvc - nB.*repmat(sum(nB.*rvc,4),[1,1,1,3]); %Gram-Schmidt
        pX = normalize(pX);
        pY = CrossProduct(nB,pX);

        T_sp = (1/sqrt(2))*(pX + 1i*pY);
        T_sm = (1/sqrt(2))*(pX - 1i*pY);
        T_pi = nB;

        Z_sp = mub*sqrt(sum(B.^2,4))*(gjg*mjg - gje*(mjg+1))/h;
        Z_pi = mub*sqrt(sum(B.^2,4))*(gjg*mjg - gje*(mjg+0))/h;
        Z_sm = mub*sqrt(sum(B.^2,4))*(gjg*mjg - gje*(mjg-1))/h;


        Fscatt = zeros(size(X,1),size(X,2),size(X,3),3);
        for bi=1:size(BEAMS,1)
            s_sm = so(bi)*(abs(sum(bsxfun(@times,BeamPol(bi,1,1,:),conj(T_sm)),4)).^2);
            s_pi = so(bi)*(abs(sum(bsxfun(@times,BeamPol(bi,1,1,:),conj(T_pi)),4)).^2);
            s_sp = so(bi)*(abs(sum(bsxfun(@times,BeamPol(bi,1,1,:),conj(T_sp)),4)).^2);
            Fscatt = Fscatt + bsxfun(@times,BEAMS(bi,1,1,:),(h/lambda)*(Gamma/2)*s_sm./(1 + s_sm + 4*((detuning + Z_sm)/Gamma).^2));
            Fscatt = Fscatt + bsxfun(@times,BEAMS(bi,1,1,:),(h/lambda)*(Gamma/2)*s_pi./(1 + s_pi + 4*((detuning + Z_pi)/Gamma).^2));
            Fscatt = Fscatt + bsxfun(@times,BEAMS(bi,1,1,:),(h/lambda)*(Gamma/2)*s_sp./(1 + s_sp + 4*((detuning + Z_sp)/Gamma).^2));
        end

        Fg = m*repmat(g,[size(X,1),size(X,2),size(X,3),1]);
        F = Fscatt + Fg;
    end
end