strReverse <- function(x)
        sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

hash<-c(
'I','T','N','S','L','P','H','R','V','A','D','G','S','F','Y','C','I','T','N','S','L','P','H','R','V','A','D','G','S','F','Y','C',
'I','T','K','R','L','P','Q','R','V','A','E','G','S','L','-Cterm','-Cterm','M','T','K','R','L','P','Q','R','V','A','E','G','S','L','-Cterm','W'	
)
codons=c(
'AUA','ACA','AAC','AGC','CUA','CCA','CAC','CGA','GUA','GCA','GAC','GGA','UCA','UUC','UAC','UGC','AUC','ACC','AAU','AGU','CUC','CCC','CAU','CGC','GUC','GCC','GAU','GGC','UCC','UUU','UAU','UGU',
'AUU','ACG','AAA','AGA','CUG','CCG','CAA','CGG','GUG','GCG','GAA','GGG','UCG','UUA','UAA','UGA','AUG','ACU','AAG','AGG','CUU','CCU','CAG','CGU','GUU','GCU','GAG','GGU','UCU','UUG','UAG','UGG'

)
names(hash)<-c(
'AUA','ACA','AAC','AGC','CUA','CCA','CAC','CGA','GUA','GCA','GAC','GGA','UCA','UUC','UAC','UGC','AUC','ACC','AAU','AGU','CUC','CCC','CAU','CGC','GUC','GCC','GAU','GGC','UCC','UUU','UAU','UGU',
'AUU','ACG','AAA','AGA','CUG','CCG','CAA','CGG','GUG','GCG','GAA','GGG','UCG','UUA','UAA','UGA','AUG','ACU','AAG','AGG','CUU','CCU','CAG','CGU','GUU','GCU','GAG','GGU','UCU','UUG','UAG','UGG'
)

dna<-readline(prompt="enter dna template sequence")

# dna="3` TACGCCATGAATTTTAAAGGGATCGATAG 5`"
# dna="sdfjskdjbfjdbf"


if(substr(dna,1,1)=='5'){
    dna=strReverse(dna);
}
rna=gsub('T','1',dna)
rna=gsub('A','U',rna)
rna=gsub('G','2',rna)
rna=gsub('C','4',rna)

rna=gsub('1','A',rna)
rna=gsub('2','C',rna)
rna=gsub('4','G',rna)

rna=gsub(' ','',rna)
rna=gsub('3','',rna)
rna=gsub('5','',rna)
rna=gsub('`','',rna)

finRNA = paste('5`',rna,'3`');
i=1
stopCodons=c('UAA','UGA','UAG')
codon=substr(rna,i,i+2)
prot<-"Nterm-M"
f=0

while((codon %in% stopCodons)==F && (codon %in% codons)==T){
    if(f==1){
        prot=paste(prot,hash[codon],sep="");
    }
    if(codon =='AUG')
        f=1
    i=i+3
    codon=substr(rna,i,i+2)
}

if(f==1){
cat("\nRNA:\n", finRNA,"\n")
prot=paste(prot,'-Cterm',sep="")
cat("\nPolypeptide:\n",prot)
} else {
    cat("INVALID SEQ")
}
