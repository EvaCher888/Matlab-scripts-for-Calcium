% Stanislav Sherepanov Agnès Martin Pierre Fontanaud
% Calcium Analyse Reseau
%function [Result_Detection]=Analyse_Calcium_Data_Network()


addpath('BCT\2019_03_03_BCT\');

%addpath('F:\Calcium_Matlab_code\Analyse stanilav09102020\emd\');
%addpath('C:\Users\Analyse\Documents\MATLAB\Analyse stanilav09102020\emd\');
close all;
clear all;
clc
sr=2;
dt=1/sr;
Seuil_Corr_Local=0;

[fname,chemin]=uigetfile('*.*','MultiSelect','on');%loads file


if isequal([fname chemin],[0,0])
    return
else



    %chemin='F:\Donnees_Experiences_general\Analyse_Agnes_stanislav14012021\Pierre dat-cre male slice\';
    % coord
    nom_coord=[]; nom_image=[]; nom_Raster=[];nom_coords=[];k_fichiercsv=1;
    for i=1:size(fname,2)
        [rad,suff] = strtok(fname{i},'.');
        f_zip=strcmp(suff,'.zip');
        f_csv=strcmp(suff,'.csv');
        f_tif=strcmp(suff,'.tif');
        f_coord=strcmp(suff,'.coord');
        if f_zip==1
            nom_coord=fname{i};
        end
        if f_tif==1
            nom_image=fname{i};

        end
        if f_csv==1
            nom_Raster=fname{i};
            k_fichiercsv=k_fichiercsv+1;

        end
        if f_coord==1
            nom_coords=fname{i};

        end

    end

    if  ~isempty(nom_coord)
        datname=[chemin nom_coord];
        ROI=ReadImageJROI(datname);
        coord=zeros(size(ROI,2)-1,2);
        for i=1:size(ROI,2)-1
            par=split(ROI{1,i}.strName,'-');
            if ROI{1,i}.nPosition==0
                coord(i,1)=str2num(par{2});
                coord(i,2)=str2num(par{1});

            else
                coord(i,1)=str2num(par{3});
                coord(i,2)=str2num(par{2});
            end
        end
    else
        datname=[chemin char(nom_coords)];
        tmp=importdata(datname,',');

        coord=tmp.data(:,2:3);
    end

    %image
    datname=[chemin nom_image];
    im_tmp= imread( datname);
    %im_tmp=mat2gray(im_tmp);
    im=gray2ind(im_tmp,65536);

    datname=[chemin char(nom_Raster)];
    d=importdata(datname,';',1);
    %      tmp=dt:dt:size(d.data,1)*dt;tmp=tmp';
    %     d.data(:,1)=tmp;
    %d=importdata('20200626 Female PRL cre slice3 PVN ctrl prl 30 min 10 min.csv',';',1);
    %d=importdata('20022020female dat cre slice3 ctrl PRL 30min moco.csv',';',1);
    tmp=dt:dt:size(d.data,1)*dt;tmp=tmp';
    d.data(:,1)=tmp;
    lag=d.data(end,1)/2;
    maxlag=round(lag/dt);
    figure('Name','Signaux')
    plot(d.data(:,1),d.data(:,2:end));
    Fin_period=num2str(size(d.data,1)*dt);

    display('Series Temporelles Calcium: Analyse de Réseau ');


    prompt = {'Start Period 1:','Fin Period 1:',...
        'Start Period 2: ','Fin Period 2: ',...
        'Start Period 3: ','Fin Period 3: ',...
        'Start Period 4: ','Fin Period 4: ',...
        'Start Period 5: ','Fin Period 5: ',...
        'Start Period 6: ','Fin Period 6: '...
        'Start Period 7: ','Fin Period 7: ',...
        'Start Period 8: ','Fin Period 8: ',...
        'Start Period 9: ','Fin Period 9: ',...
        'Start Period 10: ','Fin Period 10: ',...
        'Nom Generique Resultat : ',...
        'Nombre Permutation',...
        'In-vivo pas de col temps 1(O)0(N)',...
        'affichage Correlation 1(O)0(N)',...
        'Seuil Significativité %',...
        'Sampling Rate(Hz)',...
        'percentile Seuil ON-OFF',...
        'Utilise Raster ON-OFF 1(O)0(N)',...
        'Correlation CutOff fixe'};




    dlg_title = 'Series Temporelles Calcium: Detection Transitoires et FFT';
    num_lines = [1  20];


    def = {'0',Fin_period,...%1et2
        '0','0',...%3et4
        '0','0',...%5et6
        '0','0',...%7et8
        '0','0',...%9et10
        '0','0',...%11et12
        '0','0',...%13et14
        '0','0'...%15et16
        '0','0',...%17et18
        '0','0',...%19et20
        'Resultats_Network',...%21
        '100',...%22
        '0',...%23
        '0',...%24
        '99',...%25
        '2',...%26
        '85',...%27
        '0'...%28
        '0'};%29



    Aopts.Resize='on';
    Aopts.WindowStyle='normal';
    Aopts.Interpreter='none';
    %answer = inputdlg(prompt,dlg_title,num_lines,def,opts);
    answer=inputdlgcol(prompt,dlg_title,num_lines,def,Aopts,2);

    PvalueS=str2double(answer{25});
    Value_CutOff=str2double(answer{29});
    sr=str2double(answer{26});
    dt=1/sr;

    dir_result=[chemin answer{21} '\'];
    if (exist(dir_result))==0
        mkdir(dir_result)
    end
    dir_result_cor=[chemin answer{21} '_' char(nom_Raster) '_Corr\'];
    if (exist(dir_result_cor))==0
        mkdir(dir_result_cor)
    end

    for i=1:20;tmp_c(i)=str2double(answer{i});end

    Cursors=[tmp_c(1:2:end-1);tmp_c(2:2:end)]';


    Cursors=Cursors( ~(Cursors(:,1)==0 & Cursors(:,2)==0),:);
    close all

    % Cursor = cell2mat({cursor_info.Position}');
    % Cursor =sort(Cursor);
    % Cursor=reshape(Cursor(:,1),2,[]);

    Time_Interval=[];
    k=1;

    for c_limits=1:size(Cursors,1)
        
      start=round(Cursors(c_limits,1)/dt);
        stop=round(Cursors(c_limits,2)/dt);
        if start==0
            start=1;

        end
        if str2double(answer{23})==0
            data1=d.data(start:stop,2:end);
            [x y1]=size(data1);
            xdata=d.data(start:stop,1); % Time data
            data=data1;
            xt=xdata;
        else
            data1=d.data(start:stop,1:end);
            [x y1]=size(data1);
            xdata=1:x; % Time data
            data=data1;
            xt=xdata';
        end
        Time_Interval{c_limits}=xt;
        for i_emd=1:y1
            bmin= movmin(data1(:,i_emd),9);
            Ligne_Corr =smooth(xt,bmin,0.3,'loess');
            %           pow.trend2 = polyfit(xdata,data(:,i_emd),trendDegree);
            %          data(:,i_emd)   =data(:,i_emd) - polyval(pow.trend2,xdata);
            %
            data(:,i_emd)=data1(:,i_emd)- Ligne_Corr(:,1);

            Seuil(k)=2*std(data(:,i_emd));
            k=k+1;
        end

    end
    G_Seuil=prctile( Seuil,str2double(answer{27}));



    for c_limits=1:size(Cursors,1)
        start=round(Cursors(c_limits,1)/dt);
        stop=round(Cursors(c_limits,2)/dt);
        if start==0
            start=1;

        end
        if str2double(answer{23})==0
            data1=d.data(start:stop,2:end);
            [x y1]=size(data1);
            xdata=d.data(start:stop,1);; % Time data
            data=data1;
            xt=xdata;
        else
            data1=d.data(start:stop,1:end);
            [x y1]=size(data1);
            xdata=1:x; % Time data
            data=data1;
            xt=xdata';
        end
        Time_Interval{c_limits}=xt;
        Raster_ONOFF=zeros(x,y1);

        for i_emd=1:y1
            bmin= movmin(data1(:,i_emd),9);
            Ligne_Corr =smooth(xt,bmin,0.3,'loess');
            %           pow.trend2 = polyfit(xdata,data(:,i_emd),trendDegree);
            %          data(:,i_emd)   =data(:,i_emd) - polyval(pow.trend2,xdata);
            %
            data(:,i_emd)=data1(:,i_emd)- Ligne_Corr(:,1);
            Index_filtre_ON_OFF=find(data(:,i_emd)>G_Seuil);
             
            Raster_ONOFF( Index_filtre_ON_OFF,i_emd)=1;

        end

        %*****************************************************************************
        %
        %*****************************************************************************

        Nperm=str2double(answer{22});

        % Result_Detection{i_emd-1,1,c_limits}=xt;

        Adjacency_Matrix=zeros(y1,y1);
        Carte_Correlation_Org=zeros(y1,y1);
        Carte_Correlation_random=zeros(y1,y1,Nperm);
        Carte_Correlation=zeros(y1,y1);
        Distance_Matrix=zeros(y1,y1);
        Distance_Matrix_Org=zeros(y1,y1);
        Seuil_Corr_moyen=[];Correlation_sauvegarde=[];Correlation_sauvegarde1=[];
        %******************************************************************************
        k_toss_total=1;
        f = waitbar(0,'Please wait...');
        index_correlation_Sauvegarde=1;
       if str2double(answer{28})==0
        Tmp_data=data;
          for i_emd=1:y1
             Tmp_data(:,i_emd)=detrend( Tmp_data(:,i_emd));
          end
       else
          Tmp_data =Raster_ONOFF;
       end
        for k1_cor=1:y1
            for k2_cor=1:y1
                if k1_cor~=k2_cor
                    [correlation_cal,lag_cal]=xcorr(Tmp_data(:,k1_cor),Tmp_data(:,k2_cor),'coeff',maxlag);
                    [AmpCorr,I] = max(correlation_cal);
                    [AmpCorrmin,M] = min(correlation_cal);
                    %                     if abs(AmpCorrmin)>abs(AmpCorr)
                    %                         AmpCorr=AmpCorrmin;
                    %                         lagDiff = lag_cal(M);
                    %                         timeDiff = lagDiff*dt;
                    %
                    %                     else
                    %                         lagDiff =lag_cal(I);
                    %
                    %                         timeDiff = lagDiff*dt;

                    %                     end
                    if isnan(AmpCorr)
                       Carte_Correlation_Org(k1_cor,k2_cor)=0;
                    else
                    Carte_Correlation_Org(k1_cor,k2_cor)=AmpCorr;
                    end
                    Correlation_sauvegarde_label{index_correlation_Sauvegarde} =[num2str(k1_cor),'_',num2str(k2_cor)];
                    if index_correlation_Sauvegarde==1

                        Correlation_sauvegarde{1}=lag_cal;
                    end
                    Correlation_sauvegarde1{index_correlation_Sauvegarde}=correlation_cal;

                    %*********************************************************************************************
                    %
                    %*********************************************************************************************
                    if str2double(answer{24})==1
                        figure('Name','Correlation','position',[100 100 900 900])
                        subplot(2,1,1)
                        plot(xt,Tmp_data(:,k1_cor),'b');
                        hold on
                        plot(xt,Tmp_data(:,k2_cor),'r');
                        hold off
                        subplot(2,1,2)
                        plot(lag_cal.*dt,correlation_cal);


                        saveas(gcf,[dir_result_cor,'Corr','_',num2str(k1_cor),'_', num2str(k2_cor),'.png'],'png');
                        %pause(2)
                        close (gcf)
                    end
                    %*********************************************************************************************
                    %
                    %*********************************************************************************************
                    %--------------------------------------------------------------
                    % tirage aléatoire
                    %----------------------------------------------------------------
                   
                    for N_toss=1:Nperm
                        %perm=randperm(length(xdata));
                        data_random=shuffle(Tmp_data(:,k2_cor)','time_shift');

                        %[correlation_cal_ramdom,lag_cal_random]=xcorr(detrend(data(:,k1_cor)),...
                        %    detrend(data(perm,k2_cor)),'coeff',maxlag);
                        [correlation_cal_ramdom,lag_cal_random]=xcorr(Tmp_data(:,k1_cor),...
                            data_random','coeff',maxlag);
                        %************************************************

                        %************************************************
                        [AmpCorr1,I] = max(correlation_cal_ramdom);
                        [AmpCorrmin1,M] = min(correlation_cal_ramdom);
                        %if abs(AmpCorrmin1)>abs(AmpCorr1)
                        %    AmpCorr1=AmpCorrmin1;
                        %end
                        %Corr_random(N_toss)=AmpCorr1;
                        if Seuil_Corr_Local==0
                             if isnan(AmpCorr1)
                                  Corr_random(N_toss)=0;
                             else
                            Corr_random(N_toss)=AmpCorr1;
                             end
                            k_toss_total=k_toss_total+1;
                        else
                               if isnan(AmpCorr1)
                                    Corr_random( N_toss)=0;
                               else
                                Corr_random( N_toss)=AmpCorr1;

                               end

                        end
                      
                    end
                    
                    if Value_CutOff==0
                        Value_CutOff=prctile (Corr_random,PvalueS);
                       
                    end
                    if Seuil_Corr_Local==1
                        if AmpCorr>Value_CutOff
                            Adjacency_Matrix( k1_cor,k2_cor)=1;
                            Carte_Correlation( k1_cor,k2_cor)=AmpCorr;
                        end
                        Seuil_Corr_moyen(k_toss_total)=Value_CutOff;
                        k_toss_total=k_toss_total+1;
                    else
                        Adjacency_Matrix( k1_cor,k2_cor)=1;
                        Carte_Correlation( k1_cor,k2_cor)=AmpCorr;
                         Carte_Correlation_random(k1_cor,k2_cor,:)=Corr_random;
                    end

                    Distance_Matrix(k1_cor,k2_cor)=sqrt((coord(k1_cor,1)-coord(k2_cor,1))^2+( coord(k1_cor,2)-coord(k2_cor,2))^2);
                    Distance_Matrix_Org(k1_cor,k2_cor)=sqrt((coord(k1_cor,1)-coord(k2_cor,1))^2+( coord(k1_cor,2)-coord(k2_cor,2))^2);


                    index_correlation_Sauvegarde=index_correlation_Sauvegarde+1;

                end


            end

            waitbar(k1_cor/y1,f,sprintf('%3.2f',k1_cor))
        end
        if Seuil_Corr_Local==0
            Seuil_Corr=Value_CutOff;%0.5;
            logicalcorrelation=Carte_Correlation>Seuil_Corr;
            Adjacency_Matrix=Adjacency_Matrix.*logicalcorrelation;
            Carte_Correlation=Carte_Correlation.*logicalcorrelation;
            Distance_Matrix=Distance_Matrix.*logicalcorrelation;
        else
            Distance_Matrix=Distance_Matrix.*Adjacency_Matrix;
            Seuil_Corr=mean(Seuil_Corr_moyen);
        end


        Liste_Adjacency{c_limits}=Adjacency_Matrix;
        Liste_Correlation{c_limits}=Carte_Correlation;
        Liste_Distance{c_limits}=Distance_Matrix;

        close(f)
        %%
        %******************************************************************************************

        %******************************************************************************************
        figure('Name','Ensemble','position',[100 100 900 900])
        ax1=axes;
        imshow(im,'DisplayRange',[2022 11086])
        hold on
        %scatter(ax1,coord(:,1),coord(:,2),'b')


        scatter(ax1,coord(:,1),coord(:,2),'k')
        for j=1:size(coord,1)
            text(coord(j,1),coord(j,2), num2str(j),...
                'Color', 'k', 'FontWeight', 'Bold','FontSize',5,'HorizontalAlignment','center');
        end


        GradientColor=0:0.001:1;
        linecolors = jet(size( GradientColor,2));


        for k1_cor=1:size(coord,1)
            for k2_cor=k1_cor:size(coord,1)
                if k1_cor~=k2_cor
                    x11=coord( k1_cor,1);y11=coord( k1_cor,2);
                    x12=coord(k2_cor,1);y12=coord(k2_cor,2);
                    XX_Tmp=linspace(x11,x12,10);
                    YY_Tmp=linspace(y11,y12,10);

                    if Carte_Correlation( k1_cor,k2_cor)>0
                        c_color=find(GradientColor>= Carte_Correlation( k1_cor,k2_cor),1,'first');
                        if isempty(c_color)
                        %if Carte_Correlation( k1_cor,k2_cor)==1
                            c_color=size(linecolors,1);
                        end

                        plot(ax1,XX_Tmp,YY_Tmp,'Color',linecolors(c_color,:),'LineWidth',1)



                    end
                end
            end

        end



        ax2=axes;
        linkaxes([ax1 ax2],'xy')

        colormap(ax2,linecolors)
        colorbar('position',[0.90 0.1 0.04 0.8]);
        %colorbar
        ax2.Visible='off';
        ax2.XTick=[];
        ax2.YTick=[];
        alpha (0.6)
        %             caxis([0 1])




        saveas(gcf,[dir_result,'Reseau_Cellules','_',char(nom_Raster),'_', num2str(c_limits),'.png'],'png');
        close (gcf)
        %%
        % Yt_Labels={'coeff<0.25','0.25<coeff>0.5','0.5<coeff>0.75','0.75<coeff>1'};
        %****************************************************************************************
        %**********************************************************************
        %************
        %**********************************************************************
        distr_degree=degrees_und(Adjacency_Matrix);
        CoeffClus=clustering_coef_bu(Adjacency_Matrix);
        TransCoeff=transitivity_bu(Adjacency_Matrix);
        r_Assor=assortativity_bin(Adjacency_Matrix,0);
        [dens,Nn,k_link]=density_und(Adjacency_Matrix);
        [OCS MaxM]=modularity_und(Adjacency_Matrix);
        Eglob = efficiency_bin(Adjacency_Matrix);
        Eloc = efficiency_bin(Adjacency_Matrix,1);
        Dist=distance_bin(Adjacency_Matrix);
        [Lambda,efficacy,eccentricity,rayons,diametre] = charpath(Dist,1,0);
        Adjacency_Matrix_Random=makerandCIJ_und(y1,sum(sum(Adjacency_Matrix))/2);
        distr_degree_Random=degrees_und(Adjacency_Matrix_Random);
        CoeffClus_random=clustering_coef_bu(Adjacency_Matrix_Random);
        Dist_Random=distance_bin(Adjacency_Matrix_Random);
        [Lambda_Random,efficacy_Random,eccentricity_Random,rayons_Radom,diametre_Random] = charpath(Dist_Random,1,0);
        Sigma=(mean(CoeffClus)/mean(CoeffClus_random))/(Lambda/Lambda_Random);
        Param_Network=[TransCoeff,Lambda,Lambda_Random,...
            mean(CoeffClus),mean(CoeffClus_random),r_Assor,Eglob,MaxM,Sigma,dens];
        %*************************************************************************
        Distance1=Distance_Matrix(:);
        Correlation1=Carte_Correlation(:);
        Correlation1_Org=Carte_Correlation_Org(:);


        Distance1_Org=Distance_Matrix_Org(:);

        Distance1=Distance1(Correlation1>0);
        Correlation1=Correlation1(Correlation1>0);
        %*****************************************************************************************
        figure('Name','DistancevsCorrelation','position',[100 100 900 900])
        classe_nbrelink=1:1:max(distr_degree);
        [N,edges] = histcounts(distr_degree,classe_nbrelink);
        centerbin=edges(2:end)/2;
        subplot(6,1,1)
        loglog(centerbin,N)
        subplot(6,1,2)

        bar(centerbin,N)
        subplot(6,1,3)
        histogram(Corr_random,100)
        hold on
        line([Seuil_Corr Seuil_Corr],[0 10000],'Color','red');
        hold off
        subplot(6,1,4)
        histogram(distr_degree_Random,classe_nbrelink)
        subplot(6,1,5)
        imagesc(Carte_Correlation)
        colorbar
        subplot(6,1,6)
        scatter(Distance1_Org,Correlation1_Org)
        hold on
        yline(Seuil_Corr,'--',['y = ' num2str(Seuil_Corr)],'LineWidth',3);
        hold off
        ylabel('Cross correlation coefficient');
        xlabel('Distance (pixels)');
        saveas(gcf,[dir_result,'DistancevsCorrelation','_',char(nom_Raster),'_', num2str(c_limits),'.png'],'png');
        close (gcf)
        %************************************************************************************************
        figure('Name','Ensemble_Random','position',[100 100 900 900])
        ax1=axes;
        imshow(im,'DisplayRange',[2022 11086])
        hold on
        %scatter(ax1,coord(:,1),coord(:,2),'b')


        scatter(coord(:,1),coord(:,2),'r','filled')
        for j=1:size(coord,1)
            text(coord(j,1),coord(j,2), num2str(j),...
                'Color', 'y', 'FontWeight', 'Bold','FontSize',5);
        end





        for k1_cor=1:size(coord,1)
            for k2_cor=k1_cor:size(coord,1)
                if k1_cor~=k2_cor
                    if Adjacency_Matrix_Random( k1_cor,k2_cor)>0
                        x11=coord( k1_cor,1);y11=coord( k1_cor,2);
                        x12=coord(k2_cor,1);y12=coord(k2_cor,2);
                        XX_Tmp=linspace(x11,x12,10);
                        YY_Tmp=linspace(y11,y12,10);




                        plot(XX_Tmp,YY_Tmp,'Color','red','LineWidth',1)


                    end
                end
            end

        end
         saveas(gcf,[dir_result,'Reseau_CellulesRandom','_',char(nom_Raster),'_',num2str(c_limits),'.png'],'png');
                        %pause(2)
                        close (gcf)

         %**********************************************************************************************************
         %
         %*********************************************************************************************************
      
       
          Carte_Correlation_Random_moy=mean(Carte_Correlation_random,3);
           



        % *************************************************************************************************
        % Param_Network=[TransCoeff,Lambda,Lambda_Random,mean(CoeffClus),mean(CoeffClus_random),r_Assor,Sigma];
        % ***********************************************************************************************
        filename_Network=[dir_result,strtok(char(nom_Raster),'.'),'_Par', num2str(c_limits)];
        outputFid = fopen([filename_Network, '.csv'],'w');



        fprintf(outputFid,'%s\t','Transitivity');
        fprintf(outputFid,'%s\t','Caracteristic Path length');
        fprintf(outputFid,'%s\t','Caracteristic Path length random');
        fprintf(outputFid,'%s\t','Mean Clustering Coef');
        fprintf(outputFid,'%s\t','Mean Clustering Coef random');
        fprintf(outputFid,'%s\t','Assortativity');
        fprintf(outputFid,'%s\t','Global efficiency');
        fprintf(outputFid,'%s\t','Modularity');

        fprintf(outputFid,'%s\t','sigma');
        fprintf(outputFid,'%s\t','densite');



        fprintf(outputFid,'\n');
        for i=1:size(Param_Network,2)
            fprintf(outputFid,'%s\t',num2str(Param_Network(i)));

        end

        fprintf(outputFid,'\n');


        fclose(outputFid);


        %%
        %*********************** Sauvegarde Correlogramme******************
        filename_CC=[dir_result,strtok(char(nom_Raster),'.'),'_CC', num2str(c_limits)];

        outputFid = fopen([ filename_CC, '.csv'],'w');
        lag_cor=cell2mat(Correlation_sauvegarde);
        Correlation_Result=cell2mat(Correlation_sauvegarde1);
        fprintf(outputFid,'%s\t','Couples_Cellules');
        k=1;
        for k1_cor=1:size(coord,1)
            for k2_cor=1:size(coord,1)
                if k1_cor~=k2_cor
                    if Adjacency_Matrix( k1_cor,k2_cor)>0
                        fprintf(outputFid,'%s\t',Correlation_sauvegarde_label{k});
                    end
                    k=k+1;
                end

            end
        end



        fprintf(outputFid,'\n');

        for i=1:size(Correlation_Result,1)

            fprintf(outputFid,'%s\t',num2str(lag_cor(i)));
            k=1;
            for k1_cor=1:size(coord,1)
                for k2_cor=1:size(coord,1)
                    if k1_cor~=k2_cor
                        if Adjacency_Matrix( k1_cor,k2_cor)>0
                            fprintf(outputFid,'%s\t',num2str(Correlation_Result(i,k)));
                        end
                        k=k+1;

                    end

                end
            end



            fprintf(outputFid,'\n');
        end


        %

        fclose(outputFid);
        % sauvegarde Matrices Correlation Totale,Correlation Significative,AdjacencyMatrix,
        % Distance Totale,Distance Significative
        filename_MAD=[dir_result,strtok(char(nom_Raster),'.'),'_MAD', num2str(c_limits),'.csv'];
        writematrix(Adjacency_Matrix,filename_MAD);

        filename_MCCT=[dir_result,strtok(char(nom_Raster),'.'),'_MCCT', num2str(c_limits),'.csv'];
        writematrix(Carte_Correlation_Org,filename_MCCT);

        filename_MCC=[dir_result,strtok(char(nom_Raster),'.'),'_MCC', num2str(c_limits),'.csv'];
        writematrix(Carte_Correlation,filename_MCC);

        filename_DistT=[dir_result,strtok(char(nom_Raster),'.'),'_DistT', num2str(c_limits),'.csv'];
        writematrix(Distance_Matrix_Org,filename_DistT);

        filename_Dist=[dir_result,strtok(char(nom_Raster),'.'),'_Dist', num2str(c_limits),'.csv'];
        writematrix(Distance_Matrix,filename_Dist);

        filename_Random_Moy=[dir_result,strtok(char(nom_Raster),'.'),'_RandomMoy', num2str(c_limits),'.csv'];
        writematrix( Carte_Correlation_Random_moy,filename_Random_Moy);
        
   

        %%
        %******************************************************************
    end

 filename_wkspace=[dir_result,strtok(char(nom_Raster),'.'),'_wkspace.mat'];
   save( filename_wkspace,'-v7.3')
end

%***********************************************************************



%*****************************************************************




