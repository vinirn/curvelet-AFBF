imgdir='artigo_fbm_2d_curvelets_analise_wavelets';
del_images(imgdir);

niterations=1000;

N=2^14;
multiWaitbar( 'CloseAll' );

commandwindow;
close all; 

% del_temp_images;

%%

gen_fbm=true;

nsurfaces=length(Hangle);
nscales=7;


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
                

    multiWaitbar( 'Itera��o', iteration/niterations);
    for surface=1:nsurfaces,
        for direction=1:3,
            kurts_cumul{surface}{direction}=[];
            logsigmas2_cumul{surface}{direction}=[];
            logvariances_cumul{surface}{direction}=[];
        end
    %     kurts_cumul_means_cumul{surface}=[];
    end


    if gen_fbm,
        B={};
        for surface=1:nsurfaces,
            multiWaitbar( 'Superficie', surface/nsurfaces);
            
%             B{surface}= fbmlevinson(N,.5) ;
            B{surface}= fbmwoodchan(N,.5) ;
        end
    end


    for surface=1:nsurfaces,

        X=B{surface};

        % Perform decomposition
        [c,s] = wavedec(X,nscales,'db2');
        sizex = size(X);
        sizec = size(c);
        val_s = s ;
        
        coefs_all={}

        for scale=1:nscales,
            cd = detcoef(c,s,scale); 
            coefs_all{scale}=cd;
        end
        
%         for direction=1:ndirections,
        direction=1;
            

            kurts=[];
            logsigmas2=[];
            logvariances=[];

%                 figure
            for scale=2:nscales,  

                
                coefs=coefs_all{nscales-scale+1};
                kurt=kurtosis(coefs);
                kurts=[kurts kurt];
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
                
                plot(xout,logn,'o','MarkerSize',5);
%                 xlim([-2 2])  
                hold all
                
                %fazer o fit apenas da regiao central
                middle=floor(length(xout)/2);
                sample=25;
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
            
            direction=4;
            legend off;
            titulo=['fbm1D iteration=' num2str(iteration)];
            title(titulo)
        
            maxwin
            if iteration==niterations
                imgname=['log distribuicao coeficientes db2 surface ' num2str(surface) ' direction ' directions{direction} 'iteration' num2str(iteration)]
%                 savefig(imgname,'dir',imgdir);
            end
            hold off
            

%         end
    end
end

close all
%%
show_images(imgdir);