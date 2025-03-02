imgdir='artigo_fbm_2d_curvelets_kurtosis';
del_images(imgdir);

directions={'HORIZONTAL','VERTICAL','DIAGONAL'};
HORZ_IDX=1;
VERT_IDX=2;
DIAG_IDX=3;

N=1000;

kurts_cumul={};
logsigmas2_cumul={};
logvariances_cumul={};
% kurts_cumul_means_cumul={};

multiWaitbar( 'CloseAll' );

iterations=1;

show_plots=false;
    
commandwindow;
close all; 

% del_temp_images;



Hangle={};

%% Minhas distribuições
% Hangle{1}=repmat([.5],1,16);
% Hangle{2}=repmat([.1],1,16);
% Hangle{3}=repmat([.9],1,16);
% Hangle{4}=repmat([.3 .6 .8 .9 .8 .6 .3 .1],1,2);
% Hangle{5}=repmat([.1 .3 .6 .8 .9 .8 .6 .3],1,2);
% Hangle{6}=repmat([.3 .1 .3 .6 .8 .9 .8 .6],1,2);
% Hangle{7}=repmat([.6 .3 .1 .3 .6 .8 .9 .8],1,2);
% Hangle{8}=repmat([.2 .2 .2 .2 .7 .7 .7 .7 ],1,2)
% Hangle{9}=repmat([.7 .7 .7 .7 .2 .2 .2 .2 ],1,2)

% x=linspace(pi/8,pi-pi/8,8);
% Hdist=sin(x).^2;
% Hdist=circshift(Hdist',-2)'
% Hangle{3}=repmat(Hdist,1,2)


%% Distribuições de Deilson
% 
% Hangle={};
% Hangle{1}=repmat(.9,1,16);
% Hangle{2}=repmat([.9 .9 .9 .9 .9 0.9*cos(pi/8) 0.9*cos(pi/8) .9],1,2);
% 
% angles=[11:-2:1 31:-2:13]*pi/16;
% 
% Hmin=0.9*cos(pi/4);
% Hangle{3}=repmat([[.9 .9 .9 .9] Hmin./cos(angles(5:8))],1,2);
% 
% Hmin=0.9*cos(3*pi/8);
% Hangle{4}=repmat([[abs(Hmin./cos(angles(1))) .9 .9] Hmin./cos(angles(4:8))],1,2);

%% Mais Distribuições de Deilson
Hangle={};
angles=[11:-2:1 31:-2:13]*pi/16;

Hmin=0.9*cos(3*pi/8);
Hangle_or=repmat([[abs(Hmin./cos(angles(1))) .9 .9] Hmin./cos(angles(4:8))],1,2);

Hangle{1}=circshift(Hangle_or,[0,4]);   
% Hangle{2}=circshift(Hangle_or,[0,3]);   
% Hangle{3}=circshift(Hangle_or,[0,2]);   
% Hangle{4}=circshift(Hangle_or,[0,1]);   
Hangle{2}=circshift(Hangle_or,[0,0]);   
% Hangle{6}=circshift(Hangle_or,[0,-1]);   
% Hangle{7}=circshift(Hangle_or,[0,-2]);   
% Hangle{8}=circshift(Hangle_or,[0,-3]);   

%%
select_angles=repmat([1 1 1 1 1 1 1 1],1,2);

gen_fbm_2D=true;

nsurfaces=length(Hangle);

for surface=1:nsurfaces,
    for direction=1:3,
        kurts_cumul{surface}{direction}=[];
        logsigmas2_cumul{surface}{direction}=[];
        logvariances_cumul{surface}{direction}=[];
    end
%     kurts_cumul_means_cumul{surface}=[];
end

for iteration=1:iterations,
    multiWaitbar( 'Iterações', iteration/iterations);


    if gen_fbm_2D,
        B={};
        for surface=1:nsurfaces,
            multiWaitbar( 'Superficie', surface/nsurfaces);
            angles=11.25:22.5:360;

            curvelet_std_types={};
            curvelet_std_types{1}='standard';
            % curvelet_std_types{2}='corrected';

            B{surface}=fbm2d_curvelets(N,N,Hangle{surface},'select_angles',select_angles,'std',curvelet_std_types{1});
        end
    end


    for surface=1:nsurfaces,

        X=B{surface};

        nscales=7;
        % Perform decomposition
        [c,s] = wavedec2(X,nscales,'db2');
        sizex = size(X);
        sizec = size(c);
        val_s = s ;
        % Extract details coefficients at levels
        % in each orientation, from wavelet decomposition 
        % structure [c,s]. 
        coefs_hor_ver_diag={}

        for scale=1:nscales,
            [chd,cvd,cdd] = detcoef2('all',c,s,scale); 
            coefs_hor_ver_diag{scale}={chd,cvd,cdd};
        end
        % [chd3,cvd3,cdd3] = detcoef2('all',c,s,3); 
        % [chd2,cvd2,cdd2] = detcoef2('all',c,s,2); 
        % [chd1,cvd1,cdd1] = detcoef2('all',c,s,1); 
        for direction=1:3,
            %histogramas

            if show_plots || (iteration==iterations),
                figure
                subplot(ceil(nscales/2),2,1);
                RadialGrid=[.5 .9];
                hline=polar_discrete_angles(Hangle{surface}.*select_angles, zeros(1,16),'FontSize',8,'MarkerSize',4,...
                      'ShowAnglesNumbering',false,'TextFactor',1.12,'DataFactor',1,'MarkerThreshold',1,...
                      'AngleStep',1,'GridStep',1,'RadialGrid',RadialGrid,'RadialGridTextAngle',pi/4,...
                      'RadialGridTextFactor',1,'RadialGridFontSize',8,'PlotRadius',1,'GridSupportRange',true);
                title('H - Distribuição Angular')
        %         set(hline,'Marker','o','MarkerSize',4,'color',[.5 .5 .5],'LineWidth',2)
            end

            kurts=[];
            logsigmas2=[];
            logvariances=[];
            for scale=1:nscales,  
                if show_plots || iteration==iterations,
                    if (scale/nscales)<0.5,
                        subplot(ceil(nscales/2),2,2*(scale+1)-1);
                    else
                        subplot(ceil(nscales/2),2,2*(scale-floor(nscales/2)));
                    end
                end
                coefs=coefs_hor_ver_diag{nscales-scale+1}{direction}(:);
                kurt=kurtosis(coefs);
                kurts=[kurts kurt];
                variance=var(coefs);
                logvariances=[logvariances log2(variance)];
                
                if show_plots || (iteration==iterations),
                    [n,xout]=hist(coefs,50);

                    plot(xout,n,'o');
                    f=ezfit('gauss');
                %     showfit(f','dispfitlegend','off')
                    legend off;
                    titulo=['Direction ' directions{direction} ' s=' num2str(scale)];
                    title(titulo)

                    a=f.m(1);
                    sigma=f.m(2);
%                     logsigmas2=[logsigmas2 log2(sigma^2)];

                    x0=f.m(3);
                    anot={};
                    anot(1)={['\sigma^{2}=' num2str(sigma^2)]};
                    anot(2)={['log_{2} \sigma^{2}=' num2str(log2(sigma^2))]};
                    anot(3)={['variance=' num2str(variance)]};
                    anot(4)={['kurtosis=' num2str(kurt)]};
                    xlims=get(gca,'xlim');
                    ylims=get(gca,'ylim');
                    posx=xlims(2)-(xlims(2)-xlims(1))/6;
                    posy=ylims(2)-(ylims(2)-ylims(1))/2;
                    text(posx,posy,anot)
                    disp(['direction ' directions(direction) ' s=' num2str(scale)])
                    for cont=1:length(anot)
                        disp(anot{cont}),
                    end
                    disp(' ')

                    x=linspace(xlims(1),xlims(2),50);
                    func=a*exp(-((x-x0).^2)/(2*sigma^2));
                    hold on
                    plot(x,func,'-r');
                end
            end
            kurts_cumul{surface}{direction}=[kurts_cumul{surface}{direction};kurts];
%             logsigmas2_cumul{surface}{direction}=[logsigmas2_cumul{surface}{direction};logsigmas2];
            logvariances_cumul{surface}{direction}=[logvariances_cumul{surface}{direction};logvariances];

            if show_plots || (iteration==iterations),
                maxwin
                imgname=['distribuicao coeficientes db2 surface ' num2str(surface) ' direction ' directions{direction}]
                savefig(imgname,'dir',imgdir);
            end

            if show_plots || (iteration==iterations),
                figure

                subplot(2,2,1)
                RadialGrid=[.5 .9];
                hline=polar_discrete_angles(Hangle{surface}.*select_angles, zeros(1,16),'FontSize',8,'MarkerSize',4,...
                      'ShowAnglesNumbering',false,'TextFactor',1.12,'DataFactor',1,'MarkerThreshold',1,...
                      'AngleStep',1,'GridStep',1,'RadialGrid',RadialGrid,'RadialGridTextAngle',pi/4,...
                      'RadialGridTextFactor',1,'RadialGridFontSize',8,'PlotRadius',1,'GridSupportRange',true);
                % hlines=[hlines hline];
                % legends{1}=['Hurst Sinthesys'];

                set(hline,'Marker','o','MarkerSize',4,'color',[.5 .5 .5],'LineWidth',2)
                title('H - distribuição angular')

                angles_ranges=[3/4 5/8 1/2 3/8 1/4 1/8 0 15/8 7/4 13/8 3/2 11/8 5/4 9/8 1 7/8]*pi;
                angles_ranges_text={'3\pi/4', '5\pi/8', '\pi/2', '3\pi/8', '\pi/4', '\pi/8', '0', '15\pi/8', '7\pi/4', '13\pi/8', '3\pi/2', '11\pi/8', '5\pi/4', '9\pi/8', '\pi', '7\pi/8'};

                for cont_angle=1:length(angles_ranges),
                    text(cos(angles_ranges(cont_angle))*1.15,sin(angles_ranges(cont_angle))*1.15,angles_ranges_text{cont_angle},'HorizontalAlignment','Center','Fontsize',10)
                end

                subplot(2,2,2);
                imagesc(X)
                axis square
                colormap gray
                set(gca,'xtick',[],'ytick',[]);
                title('superfície')

                subplot(2,2,3)
                if iteration==1,
                    plot(1:nscales,kurts,'-o')
                else
                    plot(1:nscales,mean(kurts_cumul{surface}{direction}),'-o')
                end
                ylabel('kurtosis')
                xlabel('scale')
                set(gca,'xtick',1:nscales);
                grid on
                ylim([2.5 6])
                title(['curtose coeficientes db2 ' directions{direction} ' iterações = ' num2str(iteration)])

                subplot(2,2,4)
    %             if iteration==1,
    %                 plot(1:nscales,[logsigmas2;logvariances],'-o')
    %             else
                if iteration>1,
    %                 plot(1:nscales,mean(logsigmas2_cumul{surface}{direction}),'-o')
    %                 hold on
                    plot(1:nscales,mean(logvariances_cumul{surface}{direction}),'-o')
                end
                ylabel('log_{2}\sigma^{2}')
                xlabel('scale')
                set(gca,'xtick',1:nscales);
    %             set(gca,'ytick',sort(logsigmas2));
                grid on
                title(['log_{2}\sigma^{2} coeficientes db2 ' directions{direction} ' iterações = ' num2str(iteration)])

                maxwin
                set(gcf,'color','w')
                imgname=['curtose coeficientes db2 surface ' num2str(surface) ' direction ' directions{direction}]
                savefig(imgname,'dir',imgdir);
            end % if show_plots || (iteration==iterations),
        
        end % end for direction
        
      
        
                                    
                                    
    end %end for surface

    close all
end %end for iteration

%%
% Diferenças das médias das curtoses entre duas superfícies horizontal e
% vertical

figure
iterations=10:10:100;
for iteration=iterations,
    dif_kurt_means=abs(mean(kurts_cumul{1}{VERT_IDX}(1:iteration,:))-mean(kurts_cumul{2}{HORZ_IDX}(1:iteration,:)));

    plot(1:nscales,dif_kurt_means,'-o');
    for scale=1:nscales,
        text(scale,dif_kurt_means(scale),num2str(iteration),'fontsize',6)
    end
    
%     set(gcf,'DefaultAxesColorOrder',[.1:.1:1;.1:.1:1;.1:.1:1]')
%     set(gca,'LineStyleOrder',{'-*',':','o'})
    title('diferença entre as curvas de curtoses medias na vertical (superficie theta=0) da horizontal (superficie theta=pi/2)')
    
    hold all
    
    
end

legend(num2str(iterations'));

maxwin
set(gcf,'color','w');
imgname=['diferenca medias curtoses horz vert '];
savefig(imgname,'dir',imgdir);


%%

close all
%%
show_images(imgdir);