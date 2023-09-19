close all;clear







%addpath(genpath('..\Entropy_measures\'))


display('Calcul duree activite-silence par cellule');


prompt = {'frequence echantillonage(Herz):'};

dlg_title = 'Calcul duree activite-silence par cellule';
num_lines = [1  30];

def = {'2'};
opts.Resize='on';
opts.WindowStyle='normal';
answer = inputdlg(prompt,dlg_title,num_lines,def,opts);
sr=str2double(answer{1});
dt=1/sr;

[fname,chemin]=uigetfile('*.*','MultiSelect','on');%loads file

if isequal([fname chemin],[0,0])
    return
else
     dir_result=[chemin,'Resultats_duree\'];
     if (exist(dir_result))==0
         mkdir(dir_result)
     end
    if iscell(fname)
        Nbrefichiers=size(fname,2);
    else
       Nbrefichiers=1; 
    end
    for i=1:Nbrefichiers
        if Nbrefichiers==1
            datname=[chemin char(fname)];
            datname_r=[dir_result char(strtrim(strtok(fname,'.'))) '_duree'];
            
        else
            datname=[chemin char(fname{i})];
            datname_r=[dir_result char(strtrim(strtok(fname{i},'.'))) '_duree'];
        end
        Ronoff=importdata(datname,',');
        Resultats_A=zeros(size(Ronoff,2),5000);
        Resultats_S=zeros(size(Ronoff,2),5000);
         dur = waitbar(0,'detection durée activite-silence');
        
        for cell=1:size(Ronoff,2)    % boucle pour la recherche des debuts d'activite ON par cellules
            activite_duree_en_secondes=[];silence_duree_en_secondes=[];
            
            Ronoff(1,cell)=0;
            Ronoff(end,cell)=0;
            loc = (Ronoff(:,cell) >= 0.5);
            d_loc = [false ; diff(loc)];
            rising_edges_index = find(d_loc == 1);
            loc = (Ronoff(:,cell)  <= 0.5);
            d_loc = [false ; diff(loc)];
            falling_edges_index = find(d_loc == 1);
            if ~isempty(rising_edges_index) && ~isempty(falling_edges_index)
                for j = 1:numel(rising_edges_index)
                    if (falling_edges_index(j) - rising_edges_index(j))==0
                        activite_duree_en_secondes(j) =dt;
                    else
                        activite_duree_en_secondes(j)=(falling_edges_index(j) - rising_edges_index(j))*dt;
                    end
                end
                for j = 1:numel(falling_edges_index)-1
                    silence_duree_en_secondes(j)=(rising_edges_index(j+1)-falling_edges_index(j))*dt;
                    
                end
                if    rising_edges_index(1)>1
                    prem_silence_Interval=(rising_edges_index(1)-1)*dt;
                    
                    silence_duree_en_secondes=  [prem_silence_Interval silence_duree_en_secondes];
                end
                if    falling_edges_index(end)<size(Ronoff,1)
                    der_silence_Interval=(size(Ronoff,1)-falling_edges_index(end))*dt;
                    silence_duree_en_secondes=  [silence_duree_en_secondes der_silence_Interval];
                end
                
                %P_entropy=PermEn( Ronoff(:,1)',2);
                
                Resultats_A(cell,1:length(activite_duree_en_secondes))= activite_duree_en_secondes;
                Resultats_S(cell,1:length(silence_duree_en_secondes))= silence_duree_en_secondes;
                %Resultats{cell,3}= P_entropy;
            end
            waitbar(cell/size(Ronoff,2) ,dur,sprintf('%3.2f',cell))
        end
        close(dur)
        
        outputFid = fopen([datname_r, '.csv'],'w');
        
        for cell=1:size(Ronoff,2)
            fprintf(outputFid,'%s\t',['DureeActivite_Cell_' num2str(cell)]);
            fprintf(outputFid,'%s\t',['DureeSilence_Cell_' num2str(cell)]);
        end
        fprintf(outputFid,'\n');
        
        for lignes=1:size(Resultats_A,2)
            for cell=1:size(Resultats_A,1)
                
                    if Resultats_A(cell,lignes)==0
                        fprintf(outputFid,'%s\t','  '); 
                    else
                    fprintf(outputFid,'%s\t',num2str(Resultats_A(cell,lignes)));
                    end
                     if Resultats_S(cell,lignes)==0
                        fprintf(outputFid,'%s\t','  '); 
                    else
                    fprintf(outputFid,'%s\t',num2str(Resultats_S(cell,lignes)));
                    end  
                
            end
            fprintf(outputFid,'\n');
        end
        fprintf(outputFid,'\n');
        
        
        fclose(outputFid);
        
        
        
        
        
        
        
    end
    
    
    
    
end


