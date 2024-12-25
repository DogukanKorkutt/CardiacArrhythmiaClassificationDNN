clear()
PATH= '/Users/dogukankorkut/Desktop/01_Database_ECGSignal_Unprocessed_DK/';
OutPutFilePath = '/Users/dogukankorkut/Desktop/01_Database_ECGSignal_Rhythm_Denoised_DK/';
%FileTable is the files you want to denoise
FileTable = (readtable('/Users/dogukankorkut/Desktop/ECGSignal_Info_rhythm_DK.csv','ReadVariableNames',true));
[LengthFileList,~] = size(FileTable);
LeadNames =["I";"II";"III";"aVR";"aVL";"aVF";"V1";"V2";"V3";"V4";"V5";"V6"];

DenFileList = dir(strcat(OutPutFilePath,'*.csv'));
DenFileList = struct2table(DenFileList);
DenFileList = table2array(DenFileList(:,1));

parfor i=1:LengthFileList
   DenFileName =strcat(FileTable.ID{i},'.csv');
   DenDataFile = 0;
   if(~ismember(DenFileName,DenFileList))
        DataFileName = strcat(PATH,FileTable.ID{i},'.csv');
        DataFile = table2array(readtable(DataFileName,'ReadVariableNames',true));
        [rows,~] = size(DataFile);
        DenoisingData= zeros(rows,12);
        
        for j=1:12
            OrigECG  = DataFile(:,j);   
            Fs=500;        
            fp=50;fs=60;                    
            rp=1;rs=2.5;                   
            wp=fp/(Fs/2);ws=fs/(Fs/2);     
            [n,wn]=buttord(wp,ws,rp,rs);     
            [bz,az] = butter(n,wn);
            LPassDataFile=filtfilt(bz,az,OrigECG); 
            
            t = 1:length(LPassDataFile);
            yy2 = smooth(t,LPassDataFile,0.1,'rloess');
            BWRemoveDataFile = (LPassDataFile-yy2);
            Dl1=BWRemoveDataFile;
            for k=2:length(Dl1)-1
                Dl1(k)=(2*Dl1(k)-Dl1(k-1)-Dl1(k+1))/sqrt(6);
            end
            NoisSTD = 1.4826*median(abs(Dl1-median(Dl1)));
            DenoisingData(:,j)= NLM_1dDarbon(BWRemoveDataFile,(1.5)*(NoisSTD),5000,10);  
        end

        OutputfileName =strcat(OutPutFilePath, FileTable.ID{i},'.csv');
        csvwrite(OutputfileName,DenoisingData);
        fprintf('Finished File: %s\n',FileTable.ID{i});
  end
end

