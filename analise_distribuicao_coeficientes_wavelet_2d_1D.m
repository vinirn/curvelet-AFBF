clear

imgdir='artigo_fbm_2d_curvelets_analise_wavelets';
del_images(imgdir);

niterations=5;

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
% Hangle{1}=repmat(.5,1,16);
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
savefig('H-Distribuicao Angular','dir',imgdir);


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
        for scale=1:nscales
            xouts{surface}{scale}=[];
            ns{surface}{scale}=[];
        end
%     kurts_cumul_means_cumul{surface}=[];
end

for angle=0:10:90,
    
    close all

    for iteration=1:niterations
        angle
        iteration


        multiWaitbar( 'Iteração', iteration/niterations);

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
            rotation_method='bicubic';

            B2=rotate_square_matrix(B{surface},angle,rotation_method);

            X=B2;


            X1D=tr2d1d(X);

            % Perform decomposition
            [c,s] = wavedec(X1D,nscales,'db2');
            sizex = size(X1D);
            sizec = size(c);
            val_s = s ;
            % Extract details coefficients at levels
            % in each orientation, from wavelet decomposition 
            % structure [c,s]. 
            coefs_scale={};

            for scale=1:nscales,
                cd = detcoef(c,s,scale); 
                coefs_scale{scale}=cd;
            end

            kurtscale=[];

            logsigmas2=[];
            logvariances=[];

            figure(1)
            for scale=2:nscales,  

            %HORIZONTAL
            coefs=coefs_scale{nscales-scale+1}(:);
    %         kurt=kurtosis(coefs);
    %         kurts_scale{scale}=[kurts_scale{scale} kurt];

    %         variance=var(coefs);
    %         logvariances=[logvariances log2(variance)];

            [n,xout]=hist(coefs,100);
            xouts{surface}{scale}=[xouts{surface}{scale};xout];
            ns{surface}{scale}=[ns{surface}{scale};n];

            if iteration>1,
                xout=mean(xouts{surface}{scale});
                n=mean(ns{surface}{scale});
            end

            logn=log(n);
            logn=logn/max(logn);

            %escalonamento
            xout=xout*(2^(scale-2));

            figure(2);
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


            sigma=f.m(1);
            xlims=get(gca,'xlim');
            ylims=get(gca,'ylim');

            x=linspace(xlims(1),xlims(2),200);
            func=1-(x.^2)/(2*sigma^2);

            plot(x,func,'-k','color',[.6 .6 .6])

            maxx=1.5*max(xout);
            xlim([-maxx maxx])
            xlim([-.8 .8])
            ylim([0,1])


        end

        legend off;
        titulo=['Iteration=' num2str(iteration) ' Angle=' num2str(angle)];
        title(titulo)

        maxwin
        if iteration==niterations
            imgname=['log distribuicao coeficientes db2 surface ' num2str(surface) 'iteration' num2str(iteration)  ' angle=' num2str(angle)]
            savefig(imgname,'dir',imgdir);
        end
        hold off

        end


    end
end

% close all
%%
show_images(imgdir);