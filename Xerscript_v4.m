%Questo programma fa cose di Manu
clear variables
clc
close all

%apri database in lettura
%DBid=fopen('TESTrotateleft.fa');
%DBid=fopen('nucleotideseq.fa');
DBid=fopen('nucleotideseq.fa');
    
    
    
    maxlengthSEQ=-10000;
    
    % Reading Database one time to collect the length of all the sequences
    
    %% Database MUST have this structure:
    %
    % >NX_000000  <- first line: name
    % AAAA .... GGGG      |
    %      ....           } Structure
    %TTTTT .... UUUU      |
    % >MP_11111   <- Nth-Line other name
    % ....
    %
    % It (should) not matter:
    % i. where the EOF is (on the same line or on a  newline
    % ii. If one of the structure is the last one (i.e. it ends with EOF)
    
    
    %% If the structure to be read is the last one of the database:
    %
    % ....
    % >NX_000000  <- Last structure
    % AAAA .... GGGG      |
    %      ....           } Structure
    %TTTTT .... UUUU  EOF |
    %
    % Then the while statement does not allow the last line to be added to SEQ
    % if we use:
    %   while ~feof(DBid).
    %       SEQ=[SEQ linea];
    %       linea=fgetl(DBid);
    %   end
    % Then linea contains the last line of letters AND the EOF. When the
    % condition is checked again the while loop ends and the statement
    % SEQ=[SEQ linea];  is not executed, i.e. the last line is not added to the
    % sequence.
    % If we don't put the condition ~feof(DBid)
    % the reading fails because it tries to read over the EOF.
    %
    % breakval allows to fix this problem:
    %
    %         while ~contains(linea,'>') %while line does not contain >
    %             SEQ=[SEQ linea];
    %             linea=fgetl(DBid);
    % % If linea contains the EOF then, we set breakval =1, when the while starts
    % SEQ is updated with linea we read the previous iteration. Then the loop
    % is escaped with the break
    %             if breakval
    %                 break
    %             end
    %             if feof(DBid)
    %                 breakval=1
    %             end
    %         end
    
    
    
    fprintf('Reading Database!\n')
    linea=fgetl(DBid); %read first line after name
    summaryLengSEQ = [];
    while ~feof(DBid)  %while is not the end of file
        SEQ=[]; %init SEQUENCE
        if contains(linea,'>') %find name of sequence
            linea=fgetl(DBid); %read first line after name
            breakval=0;
            while ~contains(linea,'>') %while line does not contain >
                SEQ=[SEQ linea]; %accocchiare tutte le linee
                if breakval
                    break
                end
                linea=fgetl(DBid); %read next line
                if feof(DBid)
                    breakval=1;
                end
            end
            summaryLengSEQ = [ summaryLengSEQ length(SEQ)];
            if maxlengthSEQ < length(SEQ)
                maxlengthSEQ = length(SEQ);
            end
        end
    end
    
    
    
    fprintf("Max size found: %d",maxlengthSEQ);
    fprintf('\n');
    frewind(DBid);
    
    pp=0;
    noRead=0;
    linea=fgetl(DBid); %take linea from database


filestruc=dir('./inputfiles');
filenames={filestruc(3:end).name}';

for jj=1:length(filenames)
    fprintf(['processing file ',filenames{jj},'\n']);
    inB=importB(['./inputfiles/',filenames{jj}]);
    Name=inB.Name;
    startpos=inB.startpos;
    
    %apri file di output in scrittura (crea se inesistente)
    OUTid=fopen(['./outputfiles/output_',filenames{jj}],'w+');
    
    if ~isempty(strfind(filenames{jj},'FOR'))
        direction=1; %FOR
    else
        direction=0; %REV
    end
    xersite=str2double(string(regexp(filenames{jj},'\d\d','match'))); %28 sono quelle che devi trovare
    cersite=218-xersite; %prima (o dopo) devi recuperarne altre 190
    
    
   % Derivate quantities
   % Max dimension of subSEQ is xersite+cersite
    subSEQ = num2str(zeros(xersite+cersite,1));
 
    
    % indexstruct = if there is more than one struc with the same name this
    % vector contains the index of these structure in the input file (it works
    % with non-consecutive homonim structure as well)
    
    while ~feof(DBid)  %while is not the end of file
        pp = pp+1;
        %fprintf("Analysing %d %s \n",pp,linea);
        ReadSEQ=0;
        indexStruc = [];
        for ii=1:size(Name,1)
            %    SEQ=[]; %init SEQUENCE
            if contains(linea,Name(ii)) %find name of sequence
                fprintf(strcat(num2str(ii),". FOUND! ",Name(ii)," [",num2str(ii/size(Name,1)*100,'%.2f'),"%%] ","\n"));
                indexStruc = [ indexStruc ii ];
                ReadSEQ=1;
            end
        end
        
        linea=fgetl(DBid); %read first line after name
        if ReadSEQ
            ll=0;
            SEQ = num2str(zeros(1,summaryLengSEQ(pp)));
            breakval=0;
            while ~contains(linea,'>')  %while line does not contain >
                for kk = 1:length(linea)
                    ll=ll+1;
                    SEQ(ll)=linea(kk); %accocchiare tutte le linee
                    %                  fprintf("AAA %d \n",ll)
                end
                if breakval
                    break
                end
                linea=fgetl(DBid); %read next line
                if feof(DBid)
                    breakval=1;
                end
                
            end
            
            for ii=1:length(indexStruc)
                indx=indexStruc(ii);
                n=startpos(indx);
                if contains(SEQ,'0') || contains(SEQ,' ')
                    SEQ=replace(SEQ,'0','');
                    SEQ=replace(SEQ,' ','');
                end
                %if contains(SEQ,'0')
                %    SEQ=SEQ( 1 : (min(strfind(SEQ,'0'))-1) ); %truncate sequence
                %end
                if isempty(SEQ)
                    fprintf(strcat(" EMPTY SEQ! ",Name(ii),"\n"));
                    continue
                end
                [ subSEQ  ] = extractSUB( direction, cersite, xersite, n , SEQ);
                writeSEQ(OUTid, subSEQ, Name(indx), n)
            end
        else
            while ~contains(linea,'>') && ~feof(DBid) %while line does not contain >
                linea=fgetl(DBid); %read next line
            end
        end
        
    end
    
    
% Things to do at the end of each reading of the database    
    frewind(DBid);
    
    pp=0;
    noRead=0;
    linea=fgetl(DBid); %take linea from database

    fclose(OUTid);
end
fclose(DBid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                            FUNCTIONS                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ subSEQ  ] = extractSUB( direction, cersite, xersite, startpos, SEQ)


%FIND AND APPEND SUBSEQUENCE
if direction %a seguire FOR
    initialpos=startpos-cersite;
    finalpos=startpos+xersite-1;
    if initialpos<1 %se a cavallo dell'inizio
        fprintf('must rotate right (for)!\n')
        LoShift=cersite-startpos+1; %lo shift da fare a dx
        shiftSEQ=circshift(SEQ,LoShift);
        subSEQ=shiftSEQ(1:cersite+xersite);
    elseif finalpos>length(SEQ) %TODO: DA CHECKARE SE FUNZIA ROTATE LEFT
        fprintf('must rotate left (for)!\n')
        LoShift=initialpos-2; %lo shift da fare a sx %TODO: quel -2 puzzadiculo.
        shiftSEQ=circshift(SEQ,-LoShift);
        subSEQ=shiftSEQ(1:cersite+xersite);
    else
        subSEQ=SEQ(initialpos:finalpos);
    end
else %a precedere REV
    initialpos=startpos;
    finalpos=startpos+xersite+cersite-1;
    if finalpos>length(SEQ) %se a cavallo della fine
        fprintf('must rotate left (rev)!\n')
        LoShift=length(SEQ)-startpos+1; %lo shift da fare a sx
        shiftSEQ=circshift(SEQ,LoShift);
        subSEQ=shiftSEQ(1:cersite+xersite);
        subSEQ=subSEQ(end:-1:1); %invert
        indexA=strfind(subSEQ,'A');
        indexT=strfind(subSEQ,'T');
        indexC=strfind(subSEQ,'C');
        indexG=strfind(subSEQ,'G');
        subSEQ(indexA)='T';
        subSEQ(indexT)='A';
        subSEQ(indexC)='G';
        subSEQ(indexG)='C';
    else
        subSEQ=SEQ(initialpos:finalpos);
        subSEQ=subSEQ(end:-1:1); %invert
        indexA=strfind(subSEQ,'A');
        indexT=strfind(subSEQ,'T');
        indexC=strfind(subSEQ,'C');
        indexG=strfind(subSEQ,'G');
        subSEQ(indexA)='T';
        subSEQ(indexT)='A';
        subSEQ(indexC)='G';
        subSEQ(indexG)='C';
    end
end

return

end


function  writeSEQ(OUTid, subSEQ, nameSEQ, stpos)

fprintf(OUTid,'>');
fprintf(OUTid,nameSEQ);
fprintf(OUTid,'_');
fprintf(OUTid,int2str(stpos));
fprintf(OUTid,'\n');
fprintf(OUTid,subSEQ);
fprintf(OUTid,'\n');

return

end
