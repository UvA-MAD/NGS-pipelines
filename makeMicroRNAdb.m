function [ ] = makeMicroRNAdb( fastafile,gfffile, pathoutput )
%makeMicroRNAdb 
%   
%%

disp('makeMicroRNAdb');

mirna_gff = GFFAnnotation(gfffile);
mirna_seq = fastaread(fastafile);

mirnaAnnotTable = table(mirna_gff.Feature, mirna_gff.Start, mirna_gff.Stop, mirna_gff.Attributes);

seqannot = cell(size(mirnaAnnotTable,1),6);
countmature = 0;

for i = 1:size(mirnaAnnotTable,1)
    
    tmp = strsplit(mirnaAnnotTable.Var4{i},';');
    
    seqannot{i,1} = tmp{1}(4:end);
    seqannot{i,2} = mirnaAnnotTable.Var1{i};
    seqannot{i,3} = mirnaAnnotTable.Var2(i);
    seqannot{i,4} = mirnaAnnotTable.Var3(i);
    
    if length(tmp) == 4
        countmature = countmature +1;
        seqannot{i,6} =  tmp{4}(14:end);
        if tmp{3}(end)=='p';
            seqannot{i,5} = tmp{3}(end-1:end);
        else
            seqannot{i,5} = '?';
        end
       
    elseif length(tmp) == 3
        seqannot{i,5} = '-';
        seqannot{i,6} = '-';      
    end   
end

%%

Hseq = cell(size(mirna_seq,1),1);
for i = 1: size(mirna_seq,1);
    tmp = strsplit(mirna_seq(i).Header, ' ');
    Hseq{i} = tmp{2};
end

%%

newseq = mirna_seq;
for i = 1: size(mirna_seq,1);   
    newseq(i).Header =  Hseq{i};
    newseq(i).Sequence = rna2dna(mirna_seq(i).Sequence);
end

%%

derived = char(seqannot(:,6));
idxpremature = find(derived == '-');
hpnames = seqannot(idxpremature,1);
derivedmature = seqannot(:,6);
derivedmature(idxpremature) = [];

DBpremature = struct();
DBmaturehalhairpin = struct();
DBmatureputative = struct();

seqannotupdated = seqannot;

countputative2build = 0;
countcompletemature = 0;

idxiszero = 0;
idxisone = 0;
idxistwo = 0;

idxDBmp = 1;
idxDBm = 1;
idxDBmp = 1;

for i = 1: size(newseq,1);
    
    headseq = newseq(i).Header;
    seq = newseq(i).Sequence;  
    seqlength = length(seq);
    halfseqlength = round(seqlength/2);
    idx = find(strcmp(headseq,seqannot(:,6)));
    
    % both putative
    if isempty(idx)
        
        idxiszero = idxiszero+1;
        countputative2build = countputative2build + 2; 
            
        DBmatureputative(idxDBmp).Header = strcat('putativemature:',headseq,':5p');
        DBmatureputative(idxDBmp).Sequence = seq(1:halfseqlength);
                
        DBmatureputative(idxDBmp+1).Header = strcat('putativemature:',headseq,':3p');
        DBmatureputative(idxDBmp+1).Sequence = seq(halfseqlength+1:seqlength);
 
        idxDBmp = 2 + idxDBmp;  
        
    %putative and annotated
    elseif length(idx)  == 1  
        
        indexref = find(strcmp(seqannot{idx,6},seqannot(:,1)));
        startref = seqannot{indexref,3};  
        seqannotupdated{idx,3}=(seqannot{idx,3}-startref)+1;
        seqannotupdated{idx,4}=seqannot{idx,4}-startref;  
        seqannotupdated{indexref,5}=seq;  
        
        idxisone = idxisone +1;
        countcompletemature=countcompletemature+1;
        countputative2build = countputative2build + 1; 
        % annotated is a 5' prime
        if seqannot{idx,3}<startref+halfseqlength
            seqannotupdated{idx,5}='5p';
            % annotated  
            DBmaturehalhairpin(idxDBm).Header = strcat('mature:',seqannot{idx,6},':',seqannot{idx,1},':5p');
            DBmaturehalhairpin(idxDBm).Sequence = seq(1:halfseqlength);
            idxDBm = 1 + idxDBm;
            % putative
            DBmatureputative(idxDBmp).Header = strcat('putativemature:',seqannot{idx,6},':3p');
            DBmatureputative(idxDBmp).Sequence = seq(halfseqlength+1:seqlength);
            idxDBmp = 1+idxDBmp;
        % annotated is a 3' prime
        elseif startref+halfseqlength<seqannot{idx,3}
            seqannotupdated{idx,5}='3p';
            % annotated       
            DBmaturehalhairpin(idxDBm).Header = strcat('mature:',seqannot{idx,6},':',seqannot{idx,1},':3p');
            DBmaturehalhairpin(idxDBm).Sequence = seq(halfseqlength+1:seqlength);
            idxDBm = 1 + idxDBm;
            % putative
            DBmatureputative(idxDBmp).Header = strcat('putativemature:',seqannot{idx,6},':5p');
            DBmatureputative(idxDBmp).Sequence = seq(1:halfseqlength);
            idxDBmp = 1+idxDBmp;       
        elseif or(seqannot{idx,3}<startref,seqannot{idx,3}>startref+seqlength)
            error('we have a problem with seq length annotation')            
        end
        
               
    % annotated
    elseif length(idx) == 2
        
        idxistwo = idxistwo + 1;
        for p = 1:2
        indexref = find(strcmp(seqannot{idx(p),6},seqannot(:,1)));
        startref = seqannot{indexref,3};  
        seqannotupdated{idx(p),3}=(seqannot{idx(p),3}-startref)+1;
        seqannotupdated{idx(p),4}=seqannot{idx(p),4}-startref; 
        seqannotupdated{indexref,5}=seq;  
            
        DBmaturehalhairpin(idxDBm).Header = strcat('mature:',seqannot{idx(p),6},':',seqannot{idx(p),1},':',seqannot{idx(p),5});
            if strcmp(seqannot{idx(p),5}, '5p');
                DBmaturehalhairpin(idxDBm).Sequence = seq(1:halfseqlength);
            elseif strcmp(seqannot{idx(p),5}, '3p');
                DBmaturehalhairpin(idxDBm).Sequence = seq(halfseqlength+1:seqlength);
            else
                error('the direction is missing in the annotation!')
            end
        idxDBm = 1 + idxDBm;  
        countcompletemature=countcompletemature+2;
        end
        
    end
    
end

%% make small derived db and mature db
flanking = 3;
DBmature = struct();
DBsmallderived = struct();
idxsd = 1;
idxmat = 1;


for i = 1: size(newseq,1);
 
    headseq = newseq(i).Header;
    seq = newseq(i).Sequence;  
    idx = find(strcmp(headseq,seqannotupdated(:,6)));
    
    if length(idx)==1
        
        matureseq = seq(seqannotupdated{idx,3}:seqannotupdated{idx,4});
        halflenmatureseq = round(length(matureseq)/2);
        
        DBmature(idxmat).Header = strcat('mature:',seqannotupdated{idx,1},':',seqannotupdated{idx,5},':',seqannotupdated{idx,6});
        
        startmature = seqannotupdated{idx,3}-flanking;
        endmature = seqannotupdated{idx,4}+flanking;
        
        if startmature<1
            startmature = 1;
        end
        
        if endmature > length(seq)
            endmature = length(seq);
        end
        
        DBmature(idxmat).Sequence = seq(startmature:endmature);

        DBsmallderived(idxsd).Header = strcat('smallderived:',seqannotupdated{idx,6},':',seqannotupdated{idx,5},'5p');
        DBsmallderived(idxsd).Sequence = matureseq(1:halflenmatureseq);

        DBsmallderived(idxsd+1).Header = strcat('smallderived:',seqannotupdated{idx,6},':',seqannotupdated{idx,5},'3p');
        DBsmallderived(idxsd+1).Sequence = matureseq(halflenmatureseq+1:end);

        idxsd = 2 + idxsd;
        idxmat = idxmat + 1;
    
    elseif length(idx)==2
        for p = 1:2;
            
        matureseq = seq(seqannotupdated{idx(p),3}:seqannotupdated{idx(p),4});
        halflenmatureseq = round(length(matureseq)/2);
        
        DBmature(idxmat).Header = strcat('mature:',seqannotupdated{idx(p),1},':',seqannotupdated{idx(p),5},':',seqannotupdated{idx(p),6});
        
        startmature = seqannotupdated{idx(p),3}-flanking;
        endmature = seqannotupdated{idx(p),4}+flanking;
        
        if startmature<1
            startmature = 1;
        end
        
        if endmature > length(seq)
            endmature = length(seq);
        end
        
        DBmature(idxmat).Sequence = seq(startmature:endmature);

        DBsmallderived(idxsd).Header = strcat('smallderived:',seqannotupdated{idx(p),6},':',seqannotupdated{idx(p),5},'5p');
        DBsmallderived(idxsd).Sequence = matureseq(1:halflenmatureseq);

        DBsmallderived(idxsd+1).Header = strcat('smallderived:',seqannotupdated{idx(p),6},':',seqannotupdated{idx(p),5},'3p');
        DBsmallderived(idxsd+1).Sequence = matureseq(halflenmatureseq+1:end);

        idxsd = 2 + idxsd;
        idxmat = idxmat + 1;
        end
        
    elseif isempty(idx)
        disp (headseq)
        disp('no mature found')
        
    else
        disp (headseq)
        error('idx > 2')
    end 
end

%% premature

DBpremature = newseq;
for i = 1: size(newseq,1);   
    DBpremature(i).Header =  strcat('premature:',newseq(i).Header);
    
end

%% write new version
fastawrite(strcat(pathoutput,'/smallderived.fa'), DBsmallderived);
%%
fastawrite(strcat(pathoutput,'/mature.fa'), DBmature);
%% 
fastawrite(strcat(pathoutput,'/matureputative.fa'), DBmatureputative);
%% 
fastawrite(strcat(pathoutput,'/premature.fa'),DBpremature);

disp(strcat('micrornaDB in', pathoutput));
end

