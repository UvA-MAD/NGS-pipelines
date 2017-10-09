suppressPackageStartupMessages(require(Biostrings))


MirBaseFTP <-function (mir.version) {
   paste("ftp://mirbase.org/pub/mirbase/",mir.version,"/",sep="")
} 

substring.last <- function(x,n) { substring(x,nchar(x)-n+1,nchar(x)) } 
EmptyFile <- function(f) { cat("",file=f); f }

ReadMirBaseGff3<- function(org,mir.version,db.dir="../../Scratch/")
{
     annotation.file  = paste(db.dir,"/Mirbase",mir.version,"_",org,".gff3",sep="")
     if (!file.exists(annotation.file)) {
        cat("Try download ",paste(MirBaseFTP(mir.version),"genomes/",org,".gff3",sep=""),"\n")
         fgff <- try(read.table(paste(MirBaseFTP(mir.version),"genomes/",org,".gff3",sep="")),silent=T)
         cat(class(fgff),"is the class","\n")
         if (class(fgff) == "data.frame") {
              write.table(fgff,annotation.file,quote=F,sep="\t",row.names=F,col.names=F);
	 } 
     }
     cat(annotation.file,"\n");
     if (file.exists(annotation.file))
     {  
        gff <- read.table(annotation.file,comment="#",stringsAsFactors=F)[,c(1,3,4,5,7,9)] #skip non-informative columns
        colnames(gff)=c("reference","type","start","end","direction","comments")
        gff$ID = gsub("ID=","",sapply(as.character(gff$comments),function(x) { strsplit(x,";")[[1]][1] }))
        gff$Alias = gsub("Alias=","",sapply(as.character(gff$comments),function(x) { strsplit(x,";")[[1]][2] }))
        gff$Name = gsub("Name=","",sapply(as.character(gff$comments),function(x) { strsplit(x,";")[[1]][3] }))
        gff$derives.from = gsub("Derives_from=","",sapply(as.character(gff$comments),function(x) { strsplit(x,";")[[1]][4] }))
        rownames(gff) <- gff$ID
        gff
     }
   
} 

WriteMaturesPutativesAndSmallDerived <- function(l,matfile="matures.fa",putfile="putativematures.fa",sdfile="smallderviced.fa")
{
   disp <- c(mature=matfile,putati=putfile,smalld=sdfile)
   for (i in names(l) ) { 
         mf = disp[substr(i,1,6)] 
	 if (mf =="") { cat("no file ",i,":",substr(i,1,6),"\n") } 
         cat(">",i,"\n",l[i],"\n",sep="",file=mf,append=T) 
   }
}

WritePrematures <- function(pmt,fname) 
{ 
   names(pmt) <- paste("premature:",names(pmt),sep="")
   writeXStringSet(BStringSet(pmt),fname, append=TRUE)
}

GetSmallDerived <-function(seq,id,part,dir)
{
    length=nchar(seq)
    smd1=substr(seq,1,length/2)
    names(smd1) = paste("smallderived:",id,":",part,ifelse(dir=="+","5p","3p"),sep="")
    smd2=substr(seq,length/2+1,length)
    names(smd2) = paste("smallderived:",id,":",part,ifelse(dir=="+","3p","5p"),sep="")
    return=c(smd1,smd2)
}


GetMature <- function(id,part,part.start,part.end,precursor.start,sin,mattab,precursor.dir,mat.flank=3,n.flank)
{
         prec.loc.p<-function (x,ps=precursor.start) { x-ps+1 } 
         prec.loc.n<-function (x,ps=precursor.start+nchar(sin)) { ps-x }

         getsubseq <-function(s,st,ed,ps=precursor.start,pd=precursor.dir) {
              if (pd=="+") {
                  substr(s,prec.loc.p(st),prec.loc.p(ed))
              } else {
                  substr(s,prec.loc.n(ed),prec.loc.n(st)) 
              } 
	 
         } 
         m = mattab[mattab$derives.from == id & mattab$start >= part.start & mattab$end <= part.end,]
         if (nrow(m) > 0) {
                 for (cnt in 1:nrow(m) ) { 
                      plmseq=paste(paste(rep("N",n.flank),collapse=""),getsubseq(sin,m[cnt,"start"]-mat.flank,m[cnt,"end"]+mat.flank),paste(rep("N",n.flank),collapse=""),sep="")
                      
                      names(plmseq) = paste("mature:",m[cnt,"ID"],":",part,":",id,sep="");
                      if (grepl("-[35]p$",m[cnt,"Name"])  & (part != substring.last(m[cnt,"Name"],2)))  {
                               cat("WARNING!: part disagreement for ",m$Name,"/",m$ID," in ",id,"\n");
                      } 
                      if (cnt==1) { lmseq=plmseq; } else { lmseq=c(lmseq,plmseq); }
                      lmseq=c(lmseq,GetSmallDerived(getsubseq(sin,m[cnt,"start"],m[cnt,"end"]),id,part,precursor.dir))
                 }
         } else {
                 lmseq=getsubseq(sin,part.start,part.end)
                 names(lmseq) = paste("putativemature:",id,":",part,sep="")
         } 
         lmseq
}


ProcessPrecursor <- function(id,precursor,sequence,matures,flank=3,n.flank=0) { 
   l1=precursor[id,"start"]
   l2=precursor[id,"end"]
   dir=precursor[id,"direction"]
# Split  the precursor in two equal parts.  Check both parts if there is a mature (or even more than one:see GetMature)
# If no mature is found in the annotation then give back the entire half as putative mature. 
   lh <- (l1+l2)/2
   mseq5=GetMature(id,"5p",ifelse(dir=="+",l1,lh+1),ifelse(dir=="+",lh,l2),l1,sequence,matures,dir,flank,n.flank);
   mseq3=GetMature(id,"3p",ifelse(dir=="+",lh+1,l1),ifelse(dir=="+",l2,lh),l1,sequence,matures,dir,flank,n.flank);
   res = c(mseq5,mseq3)
# Check if all matures are used in the approach above. If not try to  repair. 
# The most likely cause is that a mature spans the both halves. First check if we can identify location based on Name.
# Names ending with  -3p or -5p  are assumed to be 3' and 5' matures respectively. No putatives are output even if only 
# one mature is annotated.  
   matures.used = sum(grepl("^mature",names(res)))
   if (matures.used != nrow(matures))  { 
        if (sum(grepl("-[35]p$",matures$Name))== nrow(matures)) 
        { 
            li = 0;
            for (i in 1:nrow(matures)) { 
                if ((matures[i,"start"] >= l1) & (matures[i,"end"] <=l2))
                {
                   part=substring.last(matures[i,"Name"],2)
                   plmseq=paste(paste(rep("N",n.flank),collapse=""),substr(sequence,matures[i,"start"]-l1+1-flank,matures[i,"end"]-l1+1+flank),paste(rep("N",n.flank),collapse=""),sep="") 
                   names(plmseq) = paste("mature:",matures[i,"ID"],":",part,":",id,sep="");
                   if (li==0) { res=plmseq; li=1} else { res=c(res,plmseq); }
                }
            }
        }
   }
# Check again if all matures are used in the approaches  above. If not try to  repair. 
# A possible cause is that a mature without 3p or 5p indication spans the both halves. Output all matures in the premature. Do
# not output putatives. It is not possible to decide which half to regard as putative.
   matures.used = sum(grepl("^mature",names(res)))
   if (matures.used != nrow(matures))  { 
      li = 0;
      for (i in 1:nrow(matures)) { 
         if ((matures[i,"start"] >= l1) & (matures[i,"end"] <=l2))
         {
             part="uu"
             plmseq=paste(paste(rep("N",n.flank),collapse=""),substr(sequence,matures[i,"start"]-l1+1-flank,matures[i,"end"]-l1+1+flank),paste(rep("N",n.flank),collapse=""),sep="")
             names(plmseq) = paste("mature:",matures[i,"ID"],":",part,":",id,sep="");
             if (li==0) { res=plmseq; li=1} else { res=c(res,plmseq); }
         }
      }
   }
   matures.used = sum(grepl("^mature",names(res)))
   if (matures.used != nrow(matures))  { 
     
            cat("Error : Real matures output : ",matures.used,".  Matures found in annotation: ",nrow(matures)," for precursor ",id ,"\n"); 
            cat("premature :  ",id," " ,l1," - ",lh," - ",l2,"\n") 
            for (i in 1:nrow(matures)) { 
                cat("mature ",i,": ",matures[i,"ID"],":",matures[i,"Name"]," " ,matures[i,"start"]," - ",matures[i,"end"],"\n") 
            } 
   } 
   res 
} 

get.mir.sequences <- function(org,mir.version) 
{ 
   hairpin.file=paste("/mad/MAD-RBAB/05_Reference-db/external/Mirbase/Mirbase",mir.version,"_hairpin.fa",sep="")
   if (!file.exists(hairpin.file)) {
       cmd = paste("curl ",MirBaseFTP(mir.version),"hairpin.fa.gz | gunzip > ",hairpin.file,sep="")
       system(cmd)
       if (!file.exists(hairpin.file)) {
             stop("Could not get find or download ",hairpin.file);
       }
   } 
   mir.sequences <- as.character(readAAStringSet(hairpin.file))
   mir.sequences <- mir.sequences[grepl(paste("^",org,sep=""),names(mir.sequences))]
#fasta file comment lines contain name and mirbase ID. Replace fasta id with mirbase ID only
   names(mir.sequences) <- sapply(names(mir.sequences),function(x) { strsplit(x," ")[[1]][2] } )
   mir.sequences
}



ProcessMirBaseData <- function(mir.version,db.dir,orgs)
{
   of.prematures=EmptyFile(paste(db.dir,"/premature.fa",sep=""))
   of.matures=EmptyFile(paste(db.dir,"/mature.fa",sep=""))
   of.putatives=EmptyFile(paste(db.dir,"/matureputative.fa",sep=""))
   of.smallderived=EmptyFile(paste(db.dir,"/smallderived.fa",sep=""))
   for (org in orgs) { 

      mir.ann <- ReadMirBaseGff3(org,mir.version,db.dir)
      if (class(mir.ann) == "data.frame") 
      { 
      mir.precursors <-  mir.ann[is.na(mir.ann$derives.from),]
      sel1 = grepl("_",mir.precursors$ID)
      sel2 = mir.precursors$Alias != mir.precursors$ID
      if ((sum(sel1) != sum(sel2))  |  (sum(sel1) != sum(sel1&sel2))) {  
            cat("Error reading Mirbase annotation : duplicate and alias assumption failed"); 
      }  
      id.not.alias = !mir.precursors$Alias %in% mir.precursors$ID
      if (sum(id.not.alias) != 0)  { 
          cat ("Warning: while reading Mirbase annotation : assumption that alias refers to existing id failed for ",sum(id.not.alias)," cases : ID replaced by Alias\n");
          mir.precursors[!mir.precursors$Alias %in% mir.precursors$ID,"ID"] <- mir.precursors[!mir.precursors$Alias %in% mir.precursors$ID,"Alias"]
      }
      mir.precursors <- mir.precursors[mir.precursors$Alias == mir.precursors$ID,]
      mir.mature <-  mir.ann[!is.na(mir.ann$derives.from),]
   #assumption: all matures derive from a precursor in the precursor list
      if (nrow(mir.mature[!mir.mature$derives.from %in% mir.precursors$ID,]) != 0)  {   cat("Error: mature references unknown precursor");  } 

      mir.sequences <- get.mir.sequences(org,mir.version) 
      if (sum(!rownames(mir.precursors) %in% names(mir.sequences)) != 0) {
              cat("Error: not all precursors in gff3 file (",ann.file,")can be linked to a sequence in fasta file\n") 
      } 
      if (sum(!rownames(mir.sequences) %in% names(mir.precursors)) != 0) {
              cat("Warning: not all secuences in fasta file are included in gff3 file (",ann.file,") no matures wil be generated \n") 
      } 

#Assumptions have been checked. Now start the real work. 
#All sequences in mirbase which refer the organism 
#are written as precursors even if the annotation does not mention them
      WritePrematures(mir.sequences,of.prematures)
#Now loop over all precursors extract matures and others and append to relevant files 
      for ( id in rownames(mir.precursors)) { 
      
         WriteMaturesPutativesAndSmallDerived(ProcessPrecursor(id,mir.precursors[id,],mir.sequences[id],mir.mature[mir.mature$derives.from==id,]),
           of.matures,of.putatives,of.smallderived)
      } 
      } else {
          cat(class(mir.ann),"  skipping ",org," no gff3 file on mirbase \n");
      }
   }
}

#args <- commandArgs(trailingOnly = TRUE)
#if (length(args) < 3) {
#        print(args)
#	stop("make MicroRNADB needs at least 3 args 1.mir.version 2. location to store result 3.... species  ")
#}

#ProcessMirBaseData(as.numeric(args[1]),args[2],args[3:length(args)])
