clear

imgdir='artigo_fbm_2d_curvelets_analise_wavelets';
del_images(imgdir);

directions={'HORIZONTAL','VERTICAL','DIAGONAL','3 DIREÇÕES'};
ndirections=3;
HORZ_IDX=1;
VERT_IDX=2;
DIAG_IDX=3;

niterations=100;

N=1000;
multiWaitbar( 'CloseAll' );

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
Hangle{1}=repmat(.5,1,16);
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

angles=[11:-2:1 31:-2:13]*pi/16;

Hmin=0.9*cos(3*pi/8);

Hangle_or=repmat([[abs(Hmin./cos(angles(1))) .9 .9] Hmin./cos(angles(4:8))],1,2);

Hangle{1}=circshift(Hangle_or,[0,0]);   
% Hangle{1}=repmat(.5,1,16);

% Hangle{1}=circshift(Hangle_or,[0,4]);   
% Hangle{2}=circshift(Hangle_or,[0,3]);   
% Hangle{3}=circshift(Hangle_or,[0,2]);   
% Hangle{4}=circshift(Hangle_or,[0,1]);   
% Hangle{2}=circshift(Hangle_or,[0,0]);   
% Hangle{6}=circshift(Hangle_or,[0,-1]);   
% Hangle{7}=circshift(Hangle_or,[0,-2]);   
% Hangle{8}=circshift(Hangle_or,[0,-3]);   

RadialGrid=[.5 .9];

hline=polar_discrete_angles(Hangle{1}, zeros(1,16),'FontSize',8,'ShowAnglesNumbering',false,'MarkerSize',4,...
     'TextFactor',1.12,'DataFactor',1,'MarkerThreshold',1,...
     'AngleStep',1,'GridStep',1,'RadialGrid',RadialGrid,'RadialGridTextAngle',pi/4,...
     'RadialGridTextFactor',1,'RadialGridFontSize',8,'PlotRadius',1,'GridSupportRange',true);
title('H - Distribuição Angular')


%%
select_angles=repmat([1 1 1 1 1 1 1 1],1,2);

gen_fbm_2D=true;

nsurfaces=length(Hangle);
nscales=7;


kurts_scale={};
kurtscale_cumul={};
for scale=1:nscales
    kurts_scale{scale}=[];
    kurtscale_cumul{scale}=[];
end
xouts={};
ns={};
for surface=1:nsurfaces,
    for direction=1:3,
        for scale=1:nscales
            xouts{surface}{direction}{scale}=[];
            ns{surface}{direction}{scale}=[];
        end
    end
%     kurts_cumul_means_cumul{surface}=[];
end


for iteration=1:niterations
                

    multiWaitbar( 'Iteração', iteration/niterations);
    for surface=1:nsurfaces,
        for direction=1:3,
            kurts_cumul{surface}{direction}=[];
            logsigmas2_cumul{surface}{direction}=[];
            logvariances_cumul{surface}{direction}=[];
        end
    %     kurts_cumul_means_cumul{surface}=[];
    end


    if gen_fbm_2D,
        B={};
        for surface=1:nsurfaces,
            multiWaitbar( 'Superficie', surface/nsurfaces);
            angles=11.25:22.5:360;

            curvelet_std_types={};
            curvelet_std_types{1}='standard';
            % curvelet_std_types{2}='corrected';

             B{surface}=fbm2d_curvelets(N,N,Hangle{surface},'select_angles',select_angles,'std',curvelet_std_types{1});
%              B{surface} = fastfBm2D([1024 1024],0.5);
%             B{surface} = synth2(1024,.5);
%             B{surface} = rand(1024,1024);
        end
    end


    for surface=1:nsurfaces,

        X=B{surface};

        % Perform decomposition
        [c,s] = wavedec2(X,nscales,'db2');
        sizex = size(X);
        sizec = size(c);
        val_s = s ;
        % Extract details coefficients at levels
        % in each orientation, from wavelet decomposition 
        % structure [c,s]. 
        coefs_hor_ver_diag={};

        for scale=1:nscales,
            [chd,cvd,cdd] = detcoef2('all',c,s,scale); 
            coefs_hor_ver_diag{scale}={chd,cvd,cdd};
        end
        
        kurtscale={};
        for direction=1:ndirections,
            kurtscale{direction}=[];
        end
            
        for direction=1:ndirections,
            
            logsigmas2=[];
            logvariances=[];

            figure(1)
            for scale=2:nscales,  

                %HORIZONTAL
                coefs=coefs_hor_ver_diag{nscales-scale+1}{direction}(:);

    %                 %VERTICAL
    %                 coefs2=coefs_hor_ver_diag{nscales-scale+1}{2}(:);
    %                 %DIAGONAL
    %                 coefs3=coefs_hor_ver_diag{nscales-scale+1}{3}(:);
    %                 %Efeito cumulativo
    %                 coefs=[coefs;coefs2;coefs3];

                kurt=kurtosis(coefs);
                kurts_scale{scale}=[kurts_scale{scale} kurt];

                variance=var(coefs);
                logvariances=[logvariances log2(variance)];

                [n,xout]=hist(coefs,100);
                xouts{surface}{direction}{scale}=[xouts{surface}{direction}{scale};xout];
                ns{surface}{direction}{scale}=[ns{surface}{direction}{scale};n];

                if iteration>1,
                    xout=mean(xouts{surface}{direction}{scale});
                    n=mean(ns{surface}{direction}{scale});
                end

                logn=log(n);
                logn=logn/max(logn);

                %escalonamento
                xout=xout*(2^(scale-2));

    %                 middle=floor(length(xout)/2);
    %                 sample=15;
    %                 stdx=std(xout(middle-sample:middle+sample));
    %                 xout=xout*(.2/stdx);

                figure(direction+1);
                plot(xout,logn,'o','MarkerSize',5);
    %                 xlim([-2 2])  
                hold all

                %fazer o fit apenas da regiao central
                middle=floor(length(xout)/2);
                sample=15;
                f=ezfit(xout(middle-sample:middle+sample),logn(middle-sample:middle+sample),'f(x)=1-x^2/(2*s^2)');

                xlims=get(gca,'xlim');
                ylims=get(gca,'ylim');

                x=linspace(xlims(1),xlims(2),200);

    %                         a=[1 1 1 1 .99 .98 .97];  
    %                         sigmas=[1 1 .021 .045 .1 .22 .45];
    %                         func=a(scale)-((x).^2)/(2*sigmas(scale)^2);
    %                         hold on
    %                         plot(x,func,'-r');

                sigma=f.m(1);
                xlims=get(gca,'xlim');
                ylims=get(gca,'ylim');

                x=linspace(xlims(1),xlims(2),200);
                func=1-(x.^2)/(2*sigma^2);

                plot(x,func,'-k','color',[.6 .6 .6])

    %                 xlim([-2,2])
    %                 xlim([-10,10])% fastfbm2d
    %                 xlim([-300 300])%synth2
    %                 xlim([-50 50])%rand
                maxx=1.5*max(xout);
                xlim([-maxx maxx])
                ylim([0,1])


            end

            legend off;
            titulo=['Direction ' directions{direction} ' iteration=' num2str(iteration)];
            title(titulo)

            maxwin
            if iteration==niterations
                imgname=['log distribuicao coeficientes db2 surface ' num2str(surface) ' direction ' directions{direction} 'iteration' num2str(iteration)]
                savefig(imgname,'dir',imgdir);
            end
            hold off
            if iteration>1,
%                 figure
%                 kurtscale=[];
                for scale=1:nscales,
                    kurtscale{direction}(scale)=mean(kurts_scale{scale});
                    kurtscale_cumul{scale}=[kurtscale_cumul{scale} kurtscale{direction}(scale)];
%                     plot(2:iteration,kurtscale_cumul{scale},'-');
%                     hold all
                end
%                 hold off
%                 ylim([2 4])
%                 grid on
%                 titulo=['Direction ' directions{direction} ' iteration=' num2str(iteration)];
%                 title(titulo)
%                 hold off
                
               
                figure(5);
                if direction==1,
                    clf
                end
              
                plot(1:nscales,kurtscale{direction},'-o');
                hold all
                if direction==3,
                    legend(directions{1},directions{2},directions{3})
                end
                ylim([0 6])
                grid on
                titulo=[' iteration=' num2str(iteration)];
%                 title(titulo)
                xlabel('Scale');ylabel('Kurtosis')
                imgname=['kurtosis coeficientes db2 surface ' num2str(surface) ' direction ' directions{direction} 'iteration' num2str(iteration)]
                savefig(imgname,'dir',imgdir);
                end


        end % end direction
    end
    

end

% close all
%%
show_images(imgdir);