commandwindow;
close all; 
% del_temp_images;
multiWaitbar( 'CloseAll' );

imgdir='artigo_fbm_2d_curvelets2';


N=2000;




angles=11.25:22.5:360;

% rotation_method='nearest';
% rotation_method='bilinear';
rotation_method='bicubic';
curvelet_std_types={};
curvelet_std_types{1}='standard';
% curvelet_std_types{2}='corrected';



Hangle={};
Hangle{1}=repmat([.5],1,16)
% Hangle{2}=repmat([.3 .6 .8 .9 .8 .6 .3 .1 ],1,2)

Hangle{2}=repmat([.2 .2 .2 .2 .7 .7 .7 .7 ],1,2)

x=linspace(pi/8,pi-pi/8,8);
Hdist=.9*sin(x).^2;
Hangle{3}=repmat(Hdist,1,2)




select_angles=repmat([1 1 1 1 1 1 1 1],1,2);
% select_angles=repmat([0 1 1 0 0 1 1 0],1,2);

iterations=5;

compute_Hs=true;

for cont_curvelet_std_types=1:length(curvelet_std_types),
    figure
    if compute_Hs
        save_hursts_angles={};
    end
    nangdists=length(Hangle);
    for cont_angdist=1:nangdists, %para cada distribuicao de Hs (coluna no grafico)

        if compute_Hs,    
            multiWaitbar( 'Distribui��o de �ngulos', cont_angdist/nangdists);
            multiWaitbar( 'Itera��es', 0);
            multiWaitbar( '�ngulos', 0);
            %%
            hursts_angles={};
            hursts_angles2={};

            for cont_iteration=1:iterations, %quantas iteracoes
                multiWaitbar( 'Itera��es', cont_iteration/iterations);
                cont_angdist
                cont_iteration
                B=fbm2d_curvelets(N,N,Hangle{cont_angdist},'select_angles',select_angles,'std',curvelet_std_types{cont_curvelet_std_types});

            %         B = fastfBm2D([1400 1400],0.5);
            %         hurst=load('hurst_N2048_H05.dat');
            %         B = reshape(hurst(:,3),2048,2048)';

            %     B=fbm_2d_spectral_directional(N);
            %     savefig('hurst direcao')

                %Por algum motivo, tenho de flipar a imagem left-right para se adequar a distrib angular de H
                B = fliplr(B);


                for angleidx=1:length(angles); %Para cada angulo
                   multiWaitbar( '�ngulos', angleidx/length(angles));
                   if cont_iteration==1,
                     hursts_angles{angleidx}=[];
                   end
                    angleidx

                    B2=rotate_square_matrix(B,angles(angleidx),rotation_method);

                    % Tirar a media das linhas paralelas e depois calcular H
                    s=mean(B2);
                    [hursts,methods]=estimate_hurst_methods(s);
                    hursts_angles{angleidx}=[hursts_angles{angleidx}; hursts];


            %         %Calcular H para cada linha paralela e depois tirar a media
            %         hursts_lines=[];
            %         for row=1:size(B2,1),
            %             row
            %             [hursts,methods]=estimate_hurst_methods(B2(row,:));
            %             hursts_lines=[hursts_lines; hursts];
            %         end
            %         hursts_angles2{angleidx}=[hursts_angles{angleidx}; mean(hursts_lines)];    
                end
            end
            save_hursts_angles{cont_angdist}=hursts_angles;    

        else
            hursts_angles=save_hursts_angles{cont_angdist};

        end




        %%

        subplot(1,length(Hangle),cont_angdist);

        legends={};
        hlines=[];

        RadialGrid=[.5 .9];
        Hvec=Hangle{cont_angdist}.*select_angles;
        hline=polar_discrete_angles(Hvec, zeros(1,16),'FontSize',8,'MarkerSize',4,...
              'TextFactor',1.12,'DataFactor',1,'MarkerThreshold',1,...
              'AngleStep',1,'GridStep',1,'RadialGrid',RadialGrid,'RadialGridTextAngle',pi/4,...
              'RadialGridTextFactor',1,'RadialGridFontSize',8,'PlotRadius',1,...
                'GridSupportRange',true);
        hlines=[hlines hline];
        legends{1}=['Hurst Synthesis'];

        set(hline,'Marker','o','MarkerSize',4,'color',[.5 .5 .5],'LineWidth',2)
        
        delete(findobj(gca,'marker','o'));                  

        
        hold on
        numrep=16;
        Hvec2=interp_rep(Hvec,numrep);
        logic_angles=interp_rep(zeros(1,16),numrep);


        hline=polar_discrete_angles(Hvec2, logic_angles,'FontSize',8,'MarkerSize',4,...
                              'TextFactor',1.15,'DataFactor',1,'LinkLastFirst',1);

        color_synthesis=[.0 .0 .0];                          
        set(hline,'linestyle',':','marker','.','markersize',5,'MarkerFaceColor','w','Color',color_synthesis,'LineWidth',.5);

        
        legends={};
        hlines=[];

        hline=plot([0 0],[0 1],'Color',color_synthesis,'linewidth',2);
        hlines=[hlines hline];
        legends{1}=['Hurst Synthesis'];

        
        
        used_methods_idx=[7 9 10];
        % used_methods_idx=1:12;
        cont_method=1;
        linestyles={'-.','--','-'};
        markers={'d','s','o'};
        for methodidx= used_methods_idx,
             cont_method=cont_method+1;
             means=[];
             for angleidx=1:length(angles),
                    if size(hursts_angles{angleidx},1)>1,
                        mean_hursts=mean(hursts_angles{angleidx});
                    else
                        mean_hursts=hursts_angles{angleidx};
                    end
                    means=[means mean_hursts(methodidx)];
             end

             erro_medio_quadratico=mean((means-.5-Hangle{cont_angdist}).^2);

             means2=[means(6:-1:1) means(16:-1:7)];
             hold on
             
             RadialGrid=[.5 .9];
             Hvec=means2-.5;
             hline=polar_discrete_angles(Hvec, zeros(1,16),'FontSize',8,'marker','o','MarkerSize',4,...
                      'TextFactor',1.15,'DataFactor',1,'MarkerThreshold',.1,...
                      'AngleStep',1,'GridStep',1,'RadialGrid',RadialGrid,'RadialGridTextAngle',pi/4,...
                      'RadialGridTextFactor',.9,'RadialGridFontSize',7,'RadialGridTextAngle',-(4/5)*(pi/2),...
                      'GridSupportRange',true,'GridLineWidth',.5,...
                      'GridLineStyle',':','GridColor',[.8 .8 .8]);
                      
            delete(findobj(gca,'marker','o'));                  
%             delete(hline);
            %interpolacao
            numrep=16;
            Hvec2=interp_rep(Hvec,numrep);
            logic_angles=interp_rep(zeros(1,16),numrep);


            hline=polar_discrete_angles(Hvec2, logic_angles,'FontSize',8,'MarkerSize',4,...
                              'TextFactor',1.15,'DataFactor',1,'LinkLastFirst',1);

            switch methodidx
                case 7
                    color='r';
                case 9
                    color='g';
                case 10
                    color='b';
                otherwise
                    color='k';
            end
                
            set(hline,'linestyle',':','marker','.','markersize',3,'MarkerFaceColor','w','Color',color,'LineWidth',.5);

            
            
%              hline=polar_discrete_angles(means2-.5, zeros(1,16),'DataFactor',1);        
            hline=plot([0 0],[0 1],'Color',color,'linewidth',1);

            hlines=[hlines hline];

            %Formatar nome do metodo
             methods{methodidx}(1)=upper(methods{methodidx}(1));
             switch methods{methodidx}
                 case 'Boxper'
                    methods{methodidx}='Box Per.'; 
                 case 'Per'
                    methods{methodidx}='Periodogram';
             end

             legends{cont_method}=[methods{methodidx} ', MSE=' num2str(erro_medio_quadratico,'%10.3f')];
%              set(hline,'linestyle',linestyles{cont_method},'marker',markers{cont_method},'markersize',4,'Color','k','LineWidth',.5);
                   
        end

        h=legend(hlines,legends,'FontSize',6,'Location','SouthOutside');
        
        for cont_hline=1:length(hlines)
            delete(hlines(cont_hline))
        end
        
        imgname=['Hursts angles ' curvelet_std_types{cont_curvelet_std_types}];
%         savefig(imgname,'dir',imgdir);
    end

    set(gcf,'Color','w')
    % titlename=['Hurst_N=' num2str(N) '_iterations=' num2str(iterations) ' H=' num2str(Hangle.*select_angles) ' rotationmethod=' rotation_method];
    % titlename = strrep(titlename, '0.', '.');
    % titlename = strrep(titlename, ' ', '');
    % title(titlename);

    imgname=['Hursts angles ' curvelet_std_types{cont_curvelet_std_types}];
    savefig(imgname,'dir',imgdir);
end    
export_fig('hursts angles.eps')

% %% Calcular H para cada linha paralela e depois tirar a media
% 
% figure
% 
% % hline=polar(deg2rad(angles),.1+.2*sin(abs(deg2rad(angles))),'-k');
% % set(hline,'LineWidth',2)
% % hold all
% 
% % polar(pi/2,1);
% % hold all
% 
% Hs = Hangle([6:-1:1 16:-1:7]).*select_angles([6:-1:1 16:-1:7]);
% 
% hline=polar(deg2rad([11.25:22.5:360 11.25]),[Hs Hangle(6)*select_angles(6)],'--ok');
% set(hline,'LineWidth',1)
% hold all
% 
% 
% for methodidx=4:length(methods),
%     means=[];
%      for angleidx=1:length(angles),
%             if size(hursts_angles2{angleidx},1)>1,
%                 mean_hursts=mean(hursts_angles2{angleidx});
%             else
%                 mean_hursts=hursts_angles2{angleidx};
%             end
%             means=[means mean_hursts(methodidx)];
%      end
%      
%      erro=sum((means-.5-Hs).^2);
%      
%      hline=polar(deg2rad([angles angles(1)]),[means means(1)]-.5,'-o');
%      set(hline,'DisplayName',[methods{methodidx} ' err=' num2str(erro)],'MarkerSize',3);
% 
% %      title(['method=' methods{methodidx}]);
% %      savefig(['method=' methods{methodidx}]);
%         
%      hold all
%      maxwin
% %      figure
% end
% h=legend('show');
% set(h,'FontSize',5,'Location','NorthEastOutside')
% titlename=['Hurst2_N=' num2str(N) '_iterations=' num2str(iterations) ' H=' num2str(Hangle.*select_angles) ' rotationmethod=' rotation_method];
% titlename = strrep(titlename, '0.', '.');
% titlename = strrep(titlename, ' ', '');
% title(titlename);
% 
% imgname=['Hurst2_N=' num2str(N) '_iterations=' num2str(iterations) ' Hangle=' num2str(Hangle.*select_angles) ' rotationmethod=' rotation_method];
% imgname = strrep(imgname, '0.', '.');
% imgname = strrep(imgname, ' ', '');
% 
% savefig(imgname,'dir',imgdir);

%%

% close all
show_images(imgdir);
