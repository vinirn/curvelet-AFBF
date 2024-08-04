commandwindow; close all; del_temp_images;
imgdir='artigo_fbm_2d_curvelets';

epsfilename=['artigo fbm 2d aniso.eps']

% Hangle=linspace(.1,.9,16);
% Hangle=ones(1,16)*.5;
% Hangle(1)=1;
% y=sin((1:16)*pi/(16/4));
% Hangle=(y+1)/2;

% Hangle=repmat([.2 .2 .2 .2 .8 .8 .8 .8],1,2)
% Hangle=repmat([.5 .5 .5 .5 .8 .8 .8 .8],1,2)

Hangle={};

% Hangle{1}=repmat([.5],1,16)

x=linspace(pi/8,pi-pi/8,8);
Hdist=sin(x).^2;
Hdist=circshift(Hdist',-2)'
Hangle{1}=repmat(Hdist,1,2)

% x=linspace(pi/8,pi,8);
x=linspace(pi/8,pi-pi/8,8);
Hdist=sin(x).^2;
Hangle{2}=repmat(Hdist,1,2)
% Hangle{2}=repmat([.3 .6 .8 .9 .8 .6 .3 .1],1,2)

Hangle{3}=repmat([.2 .2 .2 .2 .9 .9 .9 .9],1,2)


% Hangle{4}=repmat([.2 .2 .2 .2 .7 .7 .7 .7],1,2)

datafactor=[1 1 1 1];

% %TESTES
% Hangle={}
% Hangle{1}=repmat([.1],1,16)
% Hangle{2}=repmat([.5],1,16)
% Hangle{3}=repmat([.9],1,16)
% datafactor=[.1 .5 .9];


select_angles=repmat([1 1 1 1 1 1 1 1],1,2);
% select_angles=repmat([0 1 1 0 0 1 1 0],1,2);



img={};

letters = 'a':'z';

ncols=length(Hangle);

plotfbm=true;

for cont=1:length(Hangle),

    figure(1)
    set(gcf,'color','w')

    RadialGrid=[.5 .9];
    subplot(3,ncols,cont)
%     max(Hangle{cont})
    
    hline=polar_discrete_angles(Hangle{cont}.*select_angles, zeros(1,16),'FontSize',8,'marker','o','MarkerSize',4,...
                      'TextFactor',1.15,'DataFactor',datafactor(cont),'MarkerThreshold',.1,...
                      'AngleStep',1,'GridStep',1,'RadialGrid',RadialGrid,'RadialGridTextAngle',pi/4,...
                      'RadialGridTextFactor',.9,'RadialGridFontSize',7,'RadialGridTextAngle',-(4/5)*(pi/2),...
                      'GridSupportRange',true,'GridLineWidth',.5,...
                      'GridLineStyle',':','GridColor',[.8 .8 .8]);
                      
    delete(findobj(gca,'marker','o'));                  
    delete(hline);
    hold on
    
    %interpolacao
    numrep=16;
    Hangle2=interp_rep(Hangle{cont},numrep);
    select_angles2=interp_rep(select_angles,numrep);
    logic_angles=interp_rep(zeros(1,16),numrep);
   
    
    hline=polar_discrete_angles(Hangle2.*select_angles2, logic_angles,'FontSize',8,'MarkerSize',4,...
                      'TextFactor',1.15,'DataFactor',datafactor(cont),'LinkLastFirst',1);
    
    set(hline,'linestyle',':','marker','.','markersize',3,'MarkerFaceColor','w','Color','k','LineWidth',.5);
   

                  
    title(['(' letters(cont) ')'])
      
%      if cont>1,
%            deslocx=.33;
%            fator_posx=.18;
%            gcapos=get(gca,'Position');
%            gcapos(1)=deslocx+(cont-2)*fator_posx;
%            set(gca,'Position',gcapos)
%       end

       deslocx=.1;
       fator_posx=.25;
       gcapos=get(gca,'Position');
       gcapos(1)=deslocx+(cont-1)*fator_posx;
       set(gca,'Position',gcapos)
              


%     deslocx=.31;
%     fator_posx=.18;
%     if cont>1,
%           set_gca_posx(deslocx+(cont-1)*fator_posx)         
%     end
    

    
    if (plotfbm==true),

        n=1000;
        % img = real(ifdct_wrapping(Ct,1,m,n));
        img{cont}=fbm2d_curvelets(n,n,Hangle{cont},'select_angles',select_angles);

        subplot(3,ncols,cont+ncols)
        imagesc(img{cont});
        colormap gray
        daspect([1 1 1])
    %     set(gca,'visible','off')
        set(gca,'xtick',[]);set(gca,'ytick',[])

        gcapos=get(gca,'Position');
        gcapos(1)=deslocx+(cont-1)*fator_posx;
        set(gca,'Position',gcapos)

        increase_height(.05)

        %Spectral Density
        A=img{cont};
        F = ifftshift(fft2(fftshift(A)));
        S=abs(F).^2;

        M=size(A,1);
        N=size(A,2);
        deltax=1;
        deltay=1;
        kx1=mod(1/2+(0:(M-1))/M,1)-1/2;
        kx=kx1*(2*pi/deltax);
        ky1=mod(1/2+(0:(N-1))/N,1)-1/2;
        ky=ky1*(2*pi/deltay);
        kxshift=[kx((M/2+1):end) kx(1:(M/2))];
        kyshift=[ky((M/2+1):end) ky(1:(M/2))];

        cut=260;
        cutx1=length(kxshift)/2-cut+1
        cutx2=length(kxshift)/2+cut
        cuty1=length(kyshift)/2-cut+1
        cuty2=length(kyshift)/2+cut
        S=S(cutx1:cutx2,cuty1:cuty2);

        subplot(3,ncols,cont+2*ncols)
        imagesc(kxshift,-kyshift,real(log(S).^(1.3)))
        colormap gray
        axis square
        set(gca,'fontsize',7); set(gca,'xtick',-3:3);set(gca,'ytick',-3:3)

        gcapos=get(gca,'Position');
        gcapos(1)=deslocx+(cont-1)*fator_posx;
        set(gca,'Position',gcapos)
        set(gca,'ydir','normal')
        increase_height(.05)  

    %     gcapos=get(gca,'Position');
    %     set_gca_posy(gcapos(2)-.1)

    %     increase_width(.3)

    %     if cont>1,
    %           set_gca_posx(deslocx+(cont-1)*fator_posx)                   
    %     end
    end
    
    savefig(['fbm sintetico 2d'],'dir',imgdir)

end

set(gcf,'color','w')
export_fig(epsfilename)

% %% Surf
% figure('Color','w')
% surf(img{2}, 'LineStyle','none'); colormap gray
% % colorbar
% daspect([1 1 .004])    
% set(gca,'CameraPosition',[5200 -2300 10],'CameraViewAngle',11)
% % set(gca,'visible','off')
% ytick=get(gca,'ytick');
% set(gca,'ytick',ytick(2:end))
% savefig(['fbm sintetico 2d surf'],'dir',imgdir)
% print('-dpng','-f2','-r600','artigo surf');
% maxwin

%% Compute FFT
% A=img{2};
% F = ifftshift(fft2(fftshift(A)));
% S=abs(F).^2;
% figure
% M=size(A,1);
% N=size(A,2);
% deltax=1;
% deltay=1;
% kx1=mod(1/2+(0:(M-1))/M,1)-1/2;
% kx=kx1*(2*pi/deltax);
% ky1=mod(1/2+(0:(N-1))/N,1)-1/2;
% ky=ky1*(2*pi/deltay);
% kxshift=[kx((M/2+1):end) kx(1:(M/2))];
% kyshift=[ky((M/2+1):end) ky(1:(M/2))];
% 
% cut=260;
% cutx1=length(kxshift)/2-cut+1
% cutx2=length(kxshift)/2+cut
% cuty1=length(kyshift)/2-cut+1
% cuty2=length(kyshift)/2+cut
% S=S(cutx1:cutx2,cuty1:cuty2);
% 
% imagesc(kxshift,-kyshift,real(log(S).^(1.3)))
% set(gcf,'color','w')
% colormap gray
% axis square
% savefig('imagesc logS','dir',imgdir)
% export_fig('logS.eps')
% 
% %% Sample Horizontal and vertical lines from FFT
% Shorz=S(round(end/2),(round(end/2)+1):end);
% Svert=S((round(end/2)+1):end,round(end/2))';
% k=1:round(size(S,1)/2);
% %% Fit Horizontal line
% figure
% plot(log10(k),log10(Shorz))
% hold all
% f = fittype('log10((sigma^2)/(abs(x)^(2*H+2)))');
% [cf,gof,output] = fit(k(10:end)',log10(Shorz(10:end))',f,...
%                         'Lower',[0 0],'Upper',[1 10000000],'Startpoint',[.5 .5])
%                     
% est_ps=(cf.sigma)^2./((abs(k)).^(2*cf.H+2));
% plot(log10(k),log10(est_ps),'LineWidth',3);
% title(['Horz H=' num2str(cf.H) ' sigma=' num2str(cf.sigma) ' sse=' num2str(gof.sse) ' rmse=' num2str(gof.rmse)])
% 
% savefig('fft2 horz fit','dir',imgdir)
% %% Fir Vertical line
% figure
% plot(log10(k),log10(Svert))
% hold all
% f = fittype('log10((sigma^2)/(abs(x)^(2*H+2)))');
% [cf,gof,output] = fit(k(10:end)',log10(Svert(10:end))',f,...
%                         'Lower',[0 0],'Upper',[1 10000000],'Startpoint',[.5 .5])
% est_ps=(cf.sigma)^2./((abs(k)).^(2*cf.H+2));
% plot(log10(k),log10(est_ps),'LineWidth',3);
% title(['Vert H=' num2str(cf.H) ' sigma=' num2str(cf.sigma) ' sse=' num2str(gof.sse) ' rmse=' num2str(gof.rmse)])
% 
% savefig('fft2 vert fit','dir',imgdir)

%%
% close all
show_images(imgdir)



