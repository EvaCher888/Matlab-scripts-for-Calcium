
function [Result_Detection]=Analyse_Calcium_Data()



%addpath('F:\Calcium_Matlab_code\Analyse stanilav09102020\emd\');
%addpath('C:\Users\Analyse\Documents\MATLAB\Analyse stanilav09102020\emd\');
close all;
clear all;
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - - - - - - Variables that can be changed - - - - - -  %

% Set the accepted time error in seconds for each sample
% This variable checks that the time steps are evenly spaced
timeError = 0.25;

% Set the lowest accepted frequency in periods to be analyzed
% This variable excludes frequencies lower than lowPeriod periods
lowPeriod = 1.0;

% Set the level for auto-storage of the strongest peaks
% This variable stores peaks with powers greater than autoLevel
% of the highest/strongest peak.
autoLevel = 0.75;

% Set increment of the power spectrum
% This variable determines how many values that will be saved in
% the power spectrum file. If set to 1 all values are saved
pow.increment = 1;

sr=2;
dt=1/sr;
%*********************************************************

[fname,pathname]=uigetfile('*.csv');%loads file


if isequal([fname,pathname],[0,0])
    return
else
    datname=[pathname fname];
    d=importdata(datname,';',1);
    %      tmp=dt:dt:size(d.data,1)*dt;tmp=tmp';
    %     d.data(:,1)=tmp;
    %d=importdata('20200626 Female PRL cre slice3 PVN ctrl prl 30 min 10 min.csv',';',1);
    %d=importdata('20022020female dat cre slice3 ctrl PRL 30min moco.csv',';',1);
    tmp=dt:dt:size(d.data,1)*dt;tmp=tmp';
    d.data(:,1)=tmp;
    
    figure('Name','Signaux')
    plot(d.data(:,1),d.data(:,2:end));
    
end


display('Series Temporelles Calcium: Detection Transitoires et FFT ');


prompt = {'Start Period 1:','Fin Period 1:',...
    'Start Period 2: ','Fin Period 2: ',...
    'Start Period 3: ','Fin Period 3: ',...
    'Start Period 4: ','Fin Period 4: ',...
    'Start Period 5: ','Fin Period 5: ',...
    'Start Period 6: ','Fin Period 6: '...
    'Start Period 7: ','Fin Period 7: ',...
    'Start Period 8: ','Fin Period 8: ',...
    'Start Period 9: ','Fin Period 9: ',...
    'Start Period 10: ','Fin Period 10: '},
    

prompt1={'Detection ou FTT ou FFT methode Ulhen(0/1/2)',...
    'Nom Generique Resultat : ',...
    'Nombre Permutation',...
    'Display Oui:1/Non : 0',...
    'percentile Seuil'};



dlg_title = 'Series Temporelles Calcium: Detection Transitoires et FFT';
dlg_title1 = 'Series Temporelles Calcium: Detection Transitoires et FFT';
num_lines = [0.5  20];
num_lines1 = [0.5  20];

def = {'100','300',...
    '400','700',...
    '1200','1500',...
    '0','0',...
    '0','0',...
    '0','0',...
    '0','0',...
    '0','0'...
    '0','0',...
    '0','0',...
    };
def1 = { '0',...
    'Resultat_Date_Slice',...
    '1',...
    '1',...
    '85'};

opts.Resize='on';
opts.WindowStyle='normal';
answer = inputdlg(prompt,dlg_title,num_lines,def,opts);
opts1.Resize='on';
opts1.WindowStyle='normal';
answer1 = inputdlg(prompt1,dlg_title1,num_lines1,def1,opts1);

if sum(size(answer))~=0 
    
    for i=1:size(answer,1);tmp_c(i)=str2double(answer{i});end
    
    Cursors=[tmp_c(1:2:end-1);tmp_c(2:2:end)]';
    
    
    Cursors=Cursors( ~(Cursors(:,1)==0 & Cursors(:,2)==0),:);
    close all
    
    % Cursor = cell2mat({cursor_info.Position}');
    % Cursor =sort(Cursor);
    % Cursor=reshape(Cursor(:,1),2,[]);
    k=1;
    Time_Interval=[];
    for c_limits=1:size(Cursors,1)
        
        data1=d.data(round(Cursors(c_limits,1)/dt):round(Cursors(c_limits,2)/dt),:);
        [x y1]=size(data1);
        xdata=data1(:,1); % Time data
        data=data1;
        xt=xdata;
         Time_Interval{c_limits}=xt;
        for i_emd=2:y1
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
    G_Seuil=prctile( Seuil,str2double(answer1{5}));
    dir_result=[pathname answer1{2} '\'];
    if (exist(dir_result))==0
        mkdir(dir_result)
    end
    
    for c_limits=1:size(Cursors,1)
        
        
        data1=d.data(round(Cursors(c_limits,1)/dt):round(Cursors(c_limits,2)/dt),:);
        
        
        
        
        
        
        [x y1]=size(data1);
        pos_min=1;
        pos_max=x;
        xdata=data1(:,1); % Time data
        data=data1;
        xt=xdata;
        [m,n] = size(data);
        
        
        
        
        if str2double(answer1{1})==0
            Raster=zeros(m,n-1);
            
            Raster_ONOFF=zeros(m,n-1);
            
            RasterAmp=zeros(m,n-1);
             Profils_Coactif=zeros(m,2);
             
            % Result_Detection{i_emd-1,1,c_limits}=xt;
            for i_emd=2:y1
                
                
                
                
                %
                bmin= movmin(data1(:,i_emd),9);
                Ligne_Corr =smooth(xt,bmin,0.3,'loess');
                %           pow.trend2 = polyfit(xdata,data(:,i_emd),trendDegree);
                %          data(:,i_emd)   =data(:,i_emd) - polyval(pow.trend2,xdata);
                %
                data(:,i_emd)=data1(:,i_emd)- Ligne_Corr(:,1);
                % data(data(:,i_emd)<0,i_emd)=0;
                [pidx,criterion] = pickpeaks( data(:,i_emd),0,0);
                
                [didx,criterion] = pickpeaks( data(:,i_emd).*-1,0,0);
                
                
                %  pks - peaks values that are greater than th
                %  dep - depressions values that are less than th
                %  pid - peaks indices
                %  did - depressions indices
                [Amp_Max,Amp_Min,i_pid,i_did] = peakdet(data(:,i_emd), prctile(data(:,i_emd),99), 'threshold', 50);
                
                
                
                
                %
                tmp_didx=[didx(1:end-1) didx(2:end)];
                
                trouve_p=1;
                trouve_c=1;
                for peak=1:length(pidx)
                    % this finds the index of he rows(2) that have x in between
                    x=pidx(peak);
                    idx = find(x < tmp_didx(:,2) & x > tmp_didx(:,1));
                    if ~isempty(idx)
                        tmp_pidx(trouve_p)=x;
                        
                        
                        trouve_p=trouve_p+1;
                        
                        tmpp_didx1(trouve_c)=tmp_didx(idx,1);
                        tmpp_didx2(trouve_c)=tmp_didx(idx,2);
                        
                        trouve_c=trouve_c+1;
                        
                    end
                end
                %
                pidx=[];didx=[];
                pidx=tmp_pidx;
                didx=[tmpp_didx1'  tmpp_didx2'];
                duree_entre_creux=didx(:,2)-didx(:,1);
                
                %AMP=abs(data(didx(:,1),i_emd)-data(pidx(1,:),i_emd)) ;
                
                Index_filtre_AMP_D=find(data(pidx,i_emd)>G_Seuil);
                
                Index_filtre_ON_OFF=find(data(:,i_emd)>G_Seuil);
                
                
                
                pidx=pidx(Index_filtre_AMP_D);
                didx=[didx(Index_filtre_AMP_D,1) didx(Index_filtre_AMP_D,2)];
                
                
                
                
                
                Raster(pidx,i_emd-1)=1;
                
                Raster_ONOFF( Index_filtre_ON_OFF,i_emd-1)=1;
                
                
                
                RasterAmp(pidx,i_emd-1)=  data(pidx,i_emd);
                
                min_level=prctile(data(:,i_emd),2);
                %         echo off
                
                
                
                
                Result_Detection{i_emd-1,1,c_limits}=data1(:,i_emd); % segment brut
                Result_Detection{i_emd-1,2,c_limits}=data(:,i_emd);
                Result_Detection{i_emd-1,3,c_limits}= Ligne_Corr;
                if ~isempty(Index_filtre_AMP_D)
                    Result_Detection{i_emd-1,4,c_limits}=didx(:,1).*dt; %evt
                    Result_Detection{i_emd-1,5,c_limits}=(pidx.*dt)'; % max
                    Result_Detection{i_emd-1,6,c_limits}=abs(xt(pidx(:,1))-xt(didx(:,1)));
                    Result_Detection{i_emd-1,7,c_limits}=abs(xt(pidx(:,1))-xt(didx(:,2)));
                    Result_Detection{i_emd-1,8,c_limits}=diff(didx(:,1)).*dt;
                    Result_Detection{i_emd-1,9,c_limits}= data(pidx,i_emd);
                    
                end
                disp(i_emd);
                
                
                
                %******************************************************************************
                pow.dt      = mean(diff(xdata));
                pow.timeVec=xdata;
                pow.filter='No';
                trendDegree=2;
                
                % Determine sample frequency, frequency resolution and Nyquist frequency
                % mHz is the unit
                pow.mHz            = sr*1e3;
                pow.N           = 2.^(ceil(log(length( xdata))/log(2)));
                pow.sampleFreq =pow.mHz*1/pow.dt;
                pow.nyquist    =pow.mHz*1/(2*pow.dt);
                pow.freqVec    =pow.sampleFreq*(0:pow.N/2)/pow.N;
                pow.B          =pow.freqVec';
                pow.noCells=y1-1;
                % The most dominant frequencies are stored in pow.freq.vec
                %par.freq.vec = [];
                Nperm=str2double(answer1{3});
                %--------------------------------------------------------------
                % tirage aléatoire
                %----------------------------------------------------------------
                f = waitbar(0,'Please wait...');
                
                Power_Spect_Random=zeros(size(pow.freqVec,2),Nperm);
                for N_toss=1:Nperm
                    perm=randperm(length(xdata));
                    
                    tm_spect=Analyse_FFT(xt,detrend(data(perm,i_emd)),pow,0.2);
                    
                    Power_Spect_Random(:,N_toss) =tm_spect.B(:,2);
                    waitbar(N_toss/Nperm,f,sprintf('%3.2f',N_toss))
                end
                
                
                
                close(f)
                
                Seuil_Power=prctile(Power_Spect_Random,99,'all');
                par_Resultat=Analyse_FFT(xdata,detrend(data(:,i_emd)),pow,Seuil_Power);
                
                if str2double(answer1{4})==1
                    figure(i_emd-1)
                    subplot(2,2*size(Cursors,1),c_limits)
                    hold on
                    %
                    plot(xt,data1(:,i_emd),'r')
                    plot(xt,data(:,i_emd),'b')
                    plot(xt,Ligne_Corr,'k')
                    plot(xt(pidx),data(pidx,i_emd),'ro','MarkerSize',8)
                    plot([xt(1) xt(end)],[G_Seuil G_Seuil])
                    plot([xt(1) xt(end)],[ min_level  min_level],'LineWidth',3,'Color','yellow')
                    plot(xt(didx(:,1)),data(didx(:,1),i_emd),'gs','MarkerSize',8)
                    plot(xt(didx(:,2)),data(didx(:,2),i_emd),'cd','MarkerSize',8)
                    
                    hold off
                    
                    subplot(2,2*size(Cursors,1),c_limits+size(Cursors,1))
                    hold on
                    plot(par_Resultat.B(:,1),par_Resultat.B(:,2))
                    plot([par_Resultat.B(1,1) par_Resultat.B(end,1)],[ Seuil_Power  Seuil_Power],'LineWidth',3,'Color','yellow')
                    set(gca,'Title',text('String','Power Spectrum','Color','r'))
                    xlabel('Frequency (mHz)')
                    ylabel('Power Spectral Density')
                    hold off
                    %ylim([0 10^-3])
                    subplot(2,2*size(Cursors,1),c_limits+2*size(Cursors,1))
                    hold on
                    %
                    plot(xt,data1(:,i_emd),'r')
                    plot(xt,data(:,i_emd),'b')
                    plot(xt,Ligne_Corr,'k')
                    plot(xt(i_pid),data(i_pid,i_emd),'ro','MarkerSize',8)
                    % plot([xt(1) xt(end)],[G_Seuil G_Seuil])
                    % plot([xt(1) xt(end)],[ min_level  min_level],'LineWidth',3,'Color','yellow')
                    plot(xt(i_did),data(i_did,i_emd),'gs','MarkerSize',8)
                    hold off
                end
                
                
                Result_Detection{i_emd-1,10,c_limits}=par_Resultat.B(:,1); % spectre data
                Result_Detection{i_emd-1,11,c_limits}=par_Resultat.B(:,2); % spectre data
                Result_Detection{i_emd-1,12,c_limits}=par_Resultat.peaks(:,1); %freq
                Result_Detection{i_emd-1,13,c_limits}=par_Resultat.peaks(:,2); %peak power
                Result_Detection{i_emd-1,14,c_limits}=par_Resultat.peaks(:,3); % percentage Spectrum
                
            end
            figure(3000+c_limits)
            subplot(4,1,1)
            imagesc( RasterAmp')
            xlim([1 m])
            subplot(4,1,2)
            imagesc(Raster');
            xlim([1 m])
            subplot(4,1,3)
            
            plot(sum(Raster,2));
            ylim([0 n-1])
            xlim([1 m])
            
            subplot(4,1,4)
            
            imagesc(Raster_ONOFF')
            xlim([1 m])
            Profils_Coactif(:,1)=sum(Raster,2);
            Profils_Coactif(:,2)=sum(Raster_ONOFF,2);
            
                 filename_data=[pathname answer1{2} '\' 'AmplitudeRaster_' num2str(c_limits) '.csv'];
                 writematrix(RasterAmp,filename_data);
                 filename_data=[pathname answer1{2} '\' 'RasterONFF_' num2str(c_limits) '.csv'];
                 writematrix(Raster_ONOFF,filename_data);
                 filename_data=[pathname answer1{2} '\' 'EvtRaster_' num2str(c_limits) '.csv'];
                 writematrix(Raster,filename_data);
                 filename_data=[pathname answer1{2} '\' 'ProfilRetONOFF_' num2str(c_limits) '.csv'];
                 writematrix(Profils_Coactif,filename_data);
                
                
                %************************************************************************************************
                
                
                
                
                % pause ;
                %  close all
                %
        end
    
    
    
    
    
    
    
    
    %***********************************************************************
    if str2double(answer1{1})==1
        
        for i_emd=2:y1
            
            
            
            
            bmin= movmin(data1(:,i_emd),9);
            Ligne_Corr =smooth(xt,bmin,0.3,'loess');
            
            data(:,i_emd)=data1(:,i_emd)- Ligne_Corr(:,1);
            [res1 res2 res3]=apply(xt,detrend(data(:,i_emd)),10,[1 2 3],i_emd); % This asks forup to 5 candidate periods
            
        end
        
        
        
        
    end
    
    %*****************************************************************
    if str2double(answer1{1})==2
        
        pow.dt      = mean(diff(xdata));
        pow.timeVec=xdata;
        pow.filter='No';
        trendDegree=2;
        
        % Determine sample frequency, frequency resolution and Nyquist frequency
        % mHz is the unit
        pow.mHz            = sr*1e3;
        pow.N           = 2.^(ceil(log(length( xdata))/log(2)));
        pow.sampleFreq =pow.mHz*1/pow.dt;
        pow.nyquist    =pow.mHz*1/(2*pow.dt);
        pow.freqVec    =pow.sampleFreq*(0:pow.N/2)/pow.N;
        pow.B          =pow.freqVec';
        pow.noCells=y1-1;
        % The most dominant frequencies are stored in pow.freq.vec
        %par.freq.vec = [];
        Nperm=str2double(answer1{3});
        %--------------------------------------------------------------
        % tirage aléatoire
        %----------------------------------------------------------------
        f = waitbar(0,'Please wait...');
        for i=2:y1
            bmin= movmin(data1(:,i),9);
            Ligne_Corr =smooth(xt,bmin,0.3,'loess');
            
            data(:,i)=data1(:,i)- Ligne_Corr(:,1);
            
            for N_toss=1:Nperm
                perm=randperm(length(xdata));
                
                tm_spect=Analyse_FFT(xt,detrend(data(perm,i)),pow,0.2);
                
                par_Resultat_Random{i-1,N_toss} =tm_spect.B(:,2);
                
            end
            waitbar(i-1/y1,f,sprintf('%3.2f',i-1))
        end
        
        close(f)
        f = waitbar(0,'Please wait...');
        
        
        
        
        for i=2:y1
            
            Power_Spect_Random=zeros(size(pow.freqVec,2),Nperm);
            for N_toss=1:Nperm
                
                Power_Spect_Random(:,N_toss)=  par_Resultat_Random{ i-1,N_toss};
            end
            Seuil=prctile(Power_Spect_Random,99,'all');
            par_Resultat=Analyse_FFT(xdata,detrend(data(:,i)),pow,Seuil);
            
            Result_Detection{i-1,10,c_limits}=par_Resultat.B(:,2); % spectre data
            Result_Detection{i-1,11,c_limits}=par_Resultat.peaks;
            waitbar(i-1/y1,f,sprintf('%3.2f',i-1))
        end
        close(f)
        
    end
    
  
    end
    
    
    if str2double(answer1{1})==0
    Sauvegarde(dir_result,answer1{2},fname,Result_Detection, Time_Interval);
    end
    
    
    % Statistics
    % pow.freq.meanValue = mean(pow.freq.vec);
    % if pow.noCells > 1
    %   pow.freq.semValue = std(pow.freq.vec)/sqrt(pow.noCells-1);
    % else
    %   pow.freq.semValue = NaN;
    %
    % % Start spectrum analysis cell by cell
    %
    %  % Write info on the screen
    % fprintf('\n');
    % fprintf('Mean frequency: %f mHz (SEM %f) \n',pow.freq.meanValue,pow.freq.semValue);
    %
    %
    %
    % end
    %
    
    
    
    
    
end


