close all

clear
[fname,chemin]=uigetfile('*.*','MultiSelect','on');%loads file


if isequal([fname chemin],[0,0])
    return
else
    dir_result=[chemin,'Resultats_Ensemble\'];
    
      if (exist(dir_result))==0
        mkdir(dir_result)
    end
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
            nom_Raster{k_fichiercsv}=fname{i};
            k_fichiercsv= k_fichiercsv+1;
            
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
        datname=[chemin nom_coords];
        tmp=importdata(datname,',');
        
        coord=tmp.data(:,2:3);
    end
    
    %image
    datname=[chemin nom_image];
    im_tmp= imread( datname);
    %im_tmp=mat2gray(im_tmp);
    im=gray2ind(im_tmp,65536);
    NbrePerm=1000;
        P_value_Coact=95; 
        datname=[chemin char(nom_Raster{1})];
        EvtRaster=importdata(datname,',');
        Coactive_profile_random=zeros(size(EvtRaster,1),NbrePerm*size(nom_Raster,2));
    for k_fichiercsv=1:size(nom_Raster,2)
    datname=[chemin char(nom_Raster{k_fichiercsv})];
        EvtRaster=importdata(datname,',');
        Profils_Coactif=zeros(size( EvtRaster,1),2);
       
        Coactive_profile=sum(EvtRaster,2);
       
        for N_toss=1:NbrePerm
            
            Raster_Random_CellperCell=zeros(size(EvtRaster,1),size(EvtRaster,2));
            
            for cell_random=1:size(EvtRaster,2)
                
                Raster_Random_CellperCell(:,cell_random)=EvtRaster(randperm(size(EvtRaster,1)),cell_random);
            end
            Coactive_profile_random(:,N_toss*k_fichiercsv)= sum(Raster_Random_CellperCell,2);
        end
    end   
        P_value=prctile(Coactive_profile_random(:),P_value_Coact);
    
    
    for k_fichiercsv=1:size(nom_Raster,2)
        
        datname=[chemin char(nom_Raster{k_fichiercsv})];
        EvtRaster=importdata(datname,',');
        Profils_Coactif=zeros(size( EvtRaster,1),2);
        
        
        Coactive_profile=sum(EvtRaster,2);
       
        
        
        
        figure('Name','Assembly_Raster&Coactive profile','position',[100 100 1500 900])
        subplot(1,2,1)
        imagesc(EvtRaster')
        
        subplot(1,2,2)
        
        plot(Coactive_profile)
        pks=[];
        pks_S=[];
        Frame_Assembly=[];
        
        %[data_h_S,pks_frame_S,pks_S] = findHighactFrames_2(EvtRaster',pks);
        if ~isempty(pks_S)
            hold on
            plot([0 size(Coactive_profile,1)],[pks_S pks_S],'r')
            hold off
        else
            %[pks_S,dep,pks_frame_S,didx] = peakdet_1(Coactive_profile(:,1), P_value,'threshold') ;
            pks_frame_S=find(Coactive_profile(:,1)>P_value);
            for i=1:length(pks_frame_S)
                Frame_Assembly{i}= find(EvtRaster(pks_frame_S(i),:)>0);
            end
            
            hold on
            plot([0 size(Coactive_profile,1)],[P_value P_value],'r')
            hold off
            
        end
        
        saveas(gcf,[dir_result,'Assembly_',strtok(char(nom_Raster{k_fichiercsv}),'.'),'&Coactive profile','.png'],'png');
        Profils_Coactif(:,1)=Coactive_profile;
        Profils_Coactif(:,2)=P_value*ones(size(Coactive_profile,1),1);
        filename_data=[dir_result,'Coactive_',strtok(char(nom_Raster{k_fichiercsv}),'.'),'_','Coactive_profile.csv'];
        writematrix(Profils_Coactif,filename_data);
        
        
        
        for i=1:size(Frame_Assembly,2)
            figure('Name','Assembly_Coactive','position',[100 100 900 900])
            imshow(im,'DisplayRange',[2022 11086])
            hold on
            scatter(coord(:,1),coord(:,2),'b')
            
            num_cell=Frame_Assembly{i};
            scatter(coord(num_cell,1),coord(num_cell,2),'r','filled')
            for j=1:size(coord,1)
                text(coord(j,1),coord(j,2), num2str(j),...
                    'Color', 'y', 'FontWeight', 'Bold','FontSize',10);
            end
            
            hold off
            datname=[chemin nom_coord];
            saveas(gcf,[dir_result,'Coactive_',strtok(char(nom_Raster{k_fichiercsv}),'.'),'_', num2str(pks_frame_S(i)),'.png'],'png');
             close (gcf)
            
        end
        tmp_Assembly=zeros(size(EvtRaster,2),size(Frame_Assembly,2));
        file_name_Assembly=[dir_result,'LCparEpisodeCoAct_',strtok(char(nom_Raster{k_fichiercsv}),'.'),'.csv'];
        
        outputFid = fopen(file_name_Assembly,'w');
        
        
        
        
        for i=1:size(Frame_Assembly,2)
            
            fprintf(outputFid,'%s\t',['Ensemble' num2str(k_fichiercsv) '_' num2str(pks_frame_S(i))]);
            
        end
        fprintf(outputFid,'\n');
        for i=1:size(Frame_Assembly,2)
            a_a=Frame_Assembly{i};
            tmp_Assembly(1:length(a_a),i)=a_a';
            
        end
        for i=1:size(EvtRaster,2)
            for j=1:size(Frame_Assembly,2)
                if tmp_Assembly(i,j)>0
                    fprintf(outputFid,'%s\t',num2str(tmp_Assembly(i,j)));
                else
                    fprintf(outputFid,'%s\t',' ');
                end
                
            end
            fprintf(outputFid,'\n');
        end
        fclose(outputFid);
    end
end
