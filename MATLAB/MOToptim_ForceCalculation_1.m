function MOToptim_ForceCalculation_1()
    %---------------------------------------------------------------------%
    %---------------------------------------------------------------------%
    
    %--------------------   6 BEAMS   ------------------------------------%
    BEAMS = [
        -1	0   0
        +1  0   0
        0   -1  0
        0   +1  0
        0   0   -1
        0   0   +1
    ];
    LHCP_Beam = [1 1 1 1 0 0];
    RHCP_Beam = [0 0 0 0 1 1];
    so = [160 160 160 160 160 160];   
    GB = [0.046/2,0.046/2,-0.046];
    
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
%     so = [160 160 160 160 160];   
%     GB = [0.046/2,0.046/2,-0.046];
    %--------------------   4 BEAMS   ------------------------------------%
%     BEAMS = [
%         -1	0  +1
%         +1  0  +1
%         0   -1  0
%         0   +1  0
%     ];
%     LHCP_Beam = [1 1 0 0];
%     RHCP_Beam = [0 0 1 1];
%     so = [160 160 160 160];   
%     GB = [0.046/2,-0.046,0.046/2];
    %--------------------   3 BEAMS   ------------------------------------%
%     BEAMS = [
%         -1      0           sqrt(2)/2
%         +1/2    -sqrt(3)/2  sqrt(2)/2
%         +1/2	+sqrt(3)/2  sqrt(2)/2
%     ];
%     LHCP_Beam = [0 0 0];
%     RHCP_Beam = [1 1 1];
%     so = [160 160 160];   
%     GB = [-0.046,0.046/2,0.046/2];
    %---------------------------------------------------------------------%
    %---------------------------------------------------------------------%
    
    
    
    g = [0,0,-9.8];
    RADIUS = 5e-3;
    
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
    
    detuning = -70*Gamma;
    
    mjg = -Jgnd;
    
    xc = 0;
    yc = 0;
    zc = 0;
    %---------------------------------------------------------------------%
    %---------------------------------------------------------------------%
    
    
    LHCP_Beam = repmat(LHCP_Beam',[1,1,3]);
    RHCP_Beam = repmat(RHCP_Beam',[1,1,3]);
    g = reshape(g,[1 1 3]);


    BEAMS = normalize(permute(BEAMS,[1 3 2]));
    rvb = (rand(size(BEAMS,1),1,3)-0.5)*2;
    %x,y for beams
    bX = rvb - BEAMS.*repmat(sum(BEAMS.*rvb,3),[1,1,3]); %Gram-Schmidt
    bX = normalize(bX);
    bY = CrossProduct(BEAMS,bX);

    % LHCP
    PoL = (1/sqrt(2))*(bX + 1i*bY);
    % RHCP
    PoR = (1/sqrt(2))*(bX - 1i*bY);   
    BeamPol = normalize(LHCP_Beam.*PoL + RHCP_Beam.*PoR);
    
    SliderRange = 60e-3;
    fh = figure('units','normalized','outerposition',[0 0.1 1 0.7]);
    sldXC = uicontrol('style','slider','max',SliderRange,'min',-SliderRange,'value',0,'sliderstep',[0.0005 0.002]/(2*SliderRange),'units','normalized','position',[0.020 0.2 0.01 0.6],'callback',@recalculate);
    sldYC = uicontrol('style','slider','max',SliderRange,'min',-SliderRange,'value',0,'sliderstep',[0.0005 0.002]/(2*SliderRange),'units','normalized','position',[0.035 0.2 0.01 0.6],'callback',@recalculate);
    sldZC = uicontrol('style','slider','max',SliderRange,'min',-SliderRange,'value',0,'sliderstep',[0.0005 0.002]/(2*SliderRange),'units','normalized','position',[0.050 0.2 0.01 0.6],'callback',@recalculate);
    txtXC = uicontrol('style','text','string',sprintf('xc = %+f mm',xc*1e3),'HorizontalAlignment','Left','units','normalized','position',[0.020 0.90 0.08 0.02]);
    txtYC = uicontrol('style','text','string',sprintf('yc = %+f mm',yc*1e3),'HorizontalAlignment','Left','units','normalized','position',[0.020 0.88 0.08 0.02]);
    txtZC = uicontrol('style','text','string',sprintf('zc = %+f mm',zc*1e3),'HorizontalAlignment','Left','units','normalized','position',[0.020 0.86 0.08 0.02]);
    uicontrol('style','pushbutton','string','Reset','units','normalized','position',[0.020 0.15 0.040 0.04],'callback',@resetxcyczc);
%     uicontrol('style','pushbutton','string','Auto xyz','units','normalized','position',[0.020 0.10 0.040 0.04],'callback',@AutomaticCentering);
%     uicontrol('style','pushbutton','string','x','units','normalized','position',[0.020 0.05 0.010 0.03],'callback',@AutomaticCenteringX);
%     uicontrol('style','pushbutton','string','y','units','normalized','position',[0.035 0.05 0.010 0.03],'callback',@AutomaticCenteringY);
%     uicontrol('style','pushbutton','string','z','units','normalized','position',[0.050 0.05 0.010 0.03],'callback',@AutomaticCenteringZ);
    uicontrol('style','text','string','x','HorizontalAlignment','center','units','normalized','position',[0.0205 0.8 0.01 0.02]);
    uicontrol('style','text','string','y','HorizontalAlignment','center','units','normalized','position',[0.0355 0.8 0.01 0.02]);
    uicontrol('style','text','string','z','HorizontalAlignment','center','units','normalized','position',[0.0505 0.8 0.01 0.02]);
    
    global X Y Z Fproj;
    
    CalculateForce();
    PlotResults();
    
    
    
    
    
    
	1;

    function B = MagneticField(x,y,z,GB)
        B = zeros(size(x,1),size(x,2),3);
        B(:,:,1) = GB(1)*x;
        B(:,:,2) = GB(2)*y;
        B(:,:,3) = GB(3)*z;
    end

    function nV = normalize(V)
        nV = V./repmat(sqrt(sum(abs(V).^2,3)),[1,1,3]);
    end

    function V = CrossProduct(A,B)
        V(:,:,1) = A(:,:,2).*B(:,:,3) - A(:,:,3).*B(:,:,2);
        V(:,:,2) = A(:,:,3).*B(:,:,1) - A(:,:,1).*B(:,:,3);
        V(:,:,3) = A(:,:,1).*B(:,:,2) - A(:,:,2).*B(:,:,1);
    end

    function CalculateForce()
        [nX,nY,nZ] = sphere(200);
        nR(:,:,1) = nX;
        nR(:,:,2) = nY;
        nR(:,:,3) = nZ;
        X = RADIUS*nX;
        Y = RADIUS*nY;
        Z = RADIUS*nZ;

        B = MagneticField(X-xc,Y-yc,Z-zc,GB);

        nB = normalize(B);
        nB(isnan(nB)) = 0;

        rvc = (rand(size(X,1),size(X,2),3)-0.5)*2;
        pX = rvc - nB.*repmat(sum(nB.*rvc,3),[1,1,3]); %Gram-Schmidt
        pX = normalize(pX);
        pY = CrossProduct(nB,pX);

        T_sp = (1/sqrt(2))*(pX + 1i*pY);
        T_sm = (1/sqrt(2))*(pX - 1i*pY);
        T_pi = nB;

        Z_sp = mub*sqrt(sum(B.^2,3))*(gjg*mjg - gje*(mjg+1))/h;
        Z_pi = mub*sqrt(sum(B.^2,3))*(gjg*mjg - gje*(mjg+0))/h;
        Z_sm = mub*sqrt(sum(B.^2,3))*(gjg*mjg - gje*(mjg-1))/h;      


        Fscatt = zeros(size(X,1),size(X,2),3);
        for bi=1:size(BEAMS,1)
            Fscatt = Fscatt + repmat((abs(sum(bsxfun(@times,BeamPol(bi,1,:),conj(T_sm)),3)).^2),[1,1,3]).*bsxfun(@times,BEAMS(bi,1,:),(h/lambda)*(Gamma/2)*so(bi)./(1 + so(bi) + 4*((detuning + Z_sm)/Gamma).^2));
            Fscatt = Fscatt + repmat((abs(sum(bsxfun(@times,BeamPol(bi,1,:),conj(T_pi)),3)).^2),[1,1,3]).*bsxfun(@times,BEAMS(bi,1,:),(h/lambda)*(Gamma/2)*so(bi)./(1 + so(bi) + 4*((detuning + Z_pi)/Gamma).^2));
            Fscatt = Fscatt + repmat((abs(sum(bsxfun(@times,BeamPol(bi,1,:),conj(T_sp)),3)).^2),[1,1,3]).*bsxfun(@times,BEAMS(bi,1,:),(h/lambda)*(Gamma/2)*so(bi)./(1 + so(bi) + 4*((detuning + Z_sp)/Gamma).^2));
        end

        Fg = m*repmat(g,[size(X,1),size(X,2),1]);
        F = Fscatt + Fg;
        Fproj = sum(F.*nR,3);
    end

%     function MaxNormalForce = CalculateMaxNormalForce(r0)
%         [nX,nY,nZ] = sphere(200);
%         nR(:,:,1) = nX;
%         nR(:,:,2) = nY;
%         nR(:,:,3) = nZ;
%         X = RADIUS*nX;
%         Y = RADIUS*nY;
%         Z = RADIUS*nZ;
% 
%         B = MagneticField(X-r0(1),Y-r0(2),Z-r0(3),GB);
% 
%         nB = normalize(B);
%         nB(isnan(nB)) = 0;
% 
%         rvc = (rand(size(X,1),size(X,2),3)-0.5)*2;
%         pX = rvc - nB.*repmat(sum(nB.*rvc,3),[1,1,3]); %Gram-Schmidt
%         pX = normalize(pX);
%         pY = CrossProduct(nB,pX);
% 
%         T_sp = (1/sqrt(2))*(pX - 1i*pY);
%         T_sm = (1/sqrt(2))*(pX + 1i*pY);
%         T_pi = nB;
% 
%         Z_sp = mub*sqrt(sum(B.^2,3))*(gjg*mjg - gje*(mjg+1))/h;
%         Z_pi = mub*sqrt(sum(B.^2,3))*(gjg*mjg - gje*(mjg+0))/h;
%         Z_sm = mub*sqrt(sum(B.^2,3))*(gjg*mjg - gje*(mjg-1))/h;      
% 
% 
%         Fscatt = zeros(size(X,1),size(X,2),3);
%         for bi=1:size(BEAMS,1)
%             Fscatt = Fscatt + repmat((abs(sum(bsxfun(@times,BeamPol(bi,1,:),conj(T_sm)),3)).^2),[1,1,3]).*bsxfun(@times,BEAMS(bi,1,:),(h/lambda)*(Gamma/2)*so(bi)./(1 + so(bi) + 4*((detuning + Z_sm)/Gamma).^2));
%             Fscatt = Fscatt + repmat((abs(sum(bsxfun(@times,BeamPol(bi,1,:),conj(T_pi)),3)).^2),[1,1,3]).*bsxfun(@times,BEAMS(bi,1,:),(h/lambda)*(Gamma/2)*so(bi)./(1 + so(bi) + 4*((detuning + Z_pi)/Gamma).^2));
%             Fscatt = Fscatt + repmat((abs(sum(bsxfun(@times,BeamPol(bi,1,:),conj(T_sp)),3)).^2),[1,1,3]).*bsxfun(@times,BEAMS(bi,1,:),(h/lambda)*(Gamma/2)*so(bi)./(1 + so(bi) + 4*((detuning + Z_sp)/Gamma).^2));
%         end
% 
%         Fg = m*repmat(g,[size(X,1),size(X,2),1]);
%         F = Fscatt + Fg;
%         Fproj = sum(F.*nR,3);
%         
%         MaxNormalForce = max(Fproj(:));
%     end

    function PlotResults()
        figure(fh);
            spproj = subplot(1,2,1);
                surf(X*1e3,Y*1e3,Z*1e3,Fproj);
                xlabel('x (mm)');
                ylabel('y (mm)');
                zlabel('z (mm)');
                colorbar;
                set(gca,'clim',[min(Fproj(:)),max(Fproj(:))]);
                title('Radial projection of the total force $(\vec{F}\cdot\hat{r})$ (N)','Interpreter','Latex','FontSize',11);
                axis image;
                rotate3d on;
            subplot(2,2,4);
                [fpc,fpb] = hist(Fproj(:),100);
                bar(fpb/(m*sqrt(sum(g.^2,3))),fpc.*double(fpb<0),1,'FaceColor',[0.6 0.6 1.0]);
                hold on;
                bar(fpb/(m*sqrt(sum(g.^2,3))),fpc.*double(fpb>=0),1,'FaceColor',[1.0 0.6 0.6]);
                hold off;
                xlim([-max(abs(fpb)),max(abs(fpb))]/(m*sqrt(sum(g.^2,3))));
                ylim([0 1.05*max(fpc(:))]);
                xlabel('Radial projection of the total force $(\vec{F}\cdot\hat{r})$ (units of $\left| m\vec{g} \right|$)','Interpreter','Latex','FontSize',11);
                ylabel('Frequency');
            sppos = subplot(2,2,2);
                surf(X*1e3,Y*1e3,Z*1e3,2*(double(Fproj>0)-0.5),'linestyle','none');
                xlabel('x (mm)');
                ylabel('y (mm)');
                zlabel('z (mm)');
                set(gca,'clim',[-1,1]);
                axis image;
                rotate3d on;
                
            link = linkprop([spproj,sppos],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
            setappdata(gcf, 'StoreTheLink', link);
                
        set(txtXC,'string',sprintf('xc = %+f mm',xc*1e3));
        set(txtYC,'string',sprintf('yc = %+f mm',yc*1e3));
        set(txtZC,'string',sprintf('zc = %+f mm',zc*1e3));
    end

    function recalculate(~,~)
        xc = get(sldXC,'value');
        yc = get(sldYC,'value');
        zc = get(sldZC,'value');
        
        CalculateForce();
        PlotResults();
    end

    function resetxcyczc(~,~)
        set(sldXC,'value',0);
        set(sldYC,'value',0);
        set(sldZC,'value',0);
        xc = 0;
        yc = 0;
        zc = 0;
        
        CalculateForce();
        PlotResults();
    end

%     function r0 = AutomaticCentering(~,~)
%         options = optimset('TolX',1e-6,'MaxIter',1e3,'MaxFunEvals',1e4,'display','on');
%         [r0,~,~,~] = fminsearch(@(r) CalculateMaxNormalForce(r),[xc yc zc],options);
%         xc = r0(1);
%         yc = r0(2);
%         zc = r0(3);
%         
%         set(sldXC,'value',xc);
%         set(sldYC,'value',yc);
%         set(sldZC,'value',zc);
%         
%         CalculateForce();
%         PlotResults();
%     end
%     function r0 = AutomaticCenteringX(~,~)
%         options = optimset('TolX',1e-6,'MaxIter',1e3,'MaxFunEvals',1e4,'display','on');
%         [r0,~,~,~] = fminsearch(@(r) CalculateMaxNormalForce([r yc zc]),xc,options);
%         xc = r0;
%         set(sldXC,'value',xc);        
%         CalculateForce();
%         PlotResults();
%     end
%     function r0 = AutomaticCenteringY(~,~)
%         options = optimset('TolX',1e-6,'MaxIter',1e3,'MaxFunEvals',1e4,'display','on');
%         [r0,~,~,~] = fminsearch(@(r) CalculateMaxNormalForce([xc r zc]),yc,options);
%         yc = r0;
%         set(sldYC,'value',yc);        
%         CalculateForce();
%         PlotResults();
%     end
%     function r0 = AutomaticCenteringZ(~,~)
%         options = optimset('TolX',1e-6,'MaxIter',1e3,'MaxFunEvals',1e4,'display','on');
%         [r0,~,~,~] = fminsearch(@(r) CalculateMaxNormalForce([xc yc r]),zc,options);
%         zc = r0;
%         set(sldZC,'value',zc);        
%         CalculateForce();
%         PlotResults();
%     end
end