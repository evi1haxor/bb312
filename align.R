gap_start = 0
gap = -2
match = 2
mismatch = -3

matchchar<-function(a,b){
    if(a==b){
        return (match);
    } else{
        return (mismatch);
    }
}
len<-function(x){
    return(nchar(x))
}
print_mat<-function(x,y,A){
    print(A)
}

printall<-function(){
    sprintf("match = %d; mismatch = %d; gap_start = %d; gap_extend = %d",match, mismatch, gap_start, gap)
}

ok<-function(x=NULL){
    if(is.null(x)){
        print("\nok\n")
    }
    else{
        print(x)
    }
}

global_align<-function(x,y){
    M = matrix(0, nchar(x) + 1, nchar(y) + 1)
    S = matrix(0, nchar(x) + 1, nchar(y) + 1)
    X = matrix(0, nchar(x) + 1, nchar(y) + 1)
    Y = matrix(0, nchar(x) + 1, nchar(y) + 1)
    ax=character()
    ay=character()
    for(i in 2:(nchar(x)+1)){
        M[i,1] = -Inf
        X[i,1] = -Inf
        Y[i,1] = gap_start + (i-1) * gap
        S[i,1]=max(M[i,1],X[i,1],Y[i,1])
    }

    for(i in 2:(nchar(y)+1)){
        M[1,i] = -Inf
        X[1,i] = gap_start + (i-1) * gap
        Y[1,i] = -Inf    
        S[1,i]=max(M[1,i],X[1,i],Y[1,i])
    }
    for(i in 2:(nchar(x)+1)){
        for(j in 2:(nchar(y)+1)){


            M[i,j] = matchchar(substr(x,i-1,i-1), substr(y,j-1,j-1)) + max(
                    M[i-1,j-1],
                    X[i-1,j-1],
                    Y[i-1,j-1]
            )

            X[i,j] = max(
                    gap_start + gap + M[i,j-1],
                    gap + X[i,j-1],
                    gap_start + gap + Y[i,j-1]
            )

            Y[i,j] = max(
                    gap_start + gap + M[i-1,j],
                    gap_start + gap + X[i-1,j],
                    gap + Y[i-1,j]
            )
        S[i,j]=max(M[i,j],X[i,j],Y[i,j])
        }
    }

    seq.x <- unlist(strsplit(x, ''))
    seq.y <- unlist(strsplit(y, ''))

    seq.x <- c(0,seq.x)
    seq.y <- c(0,seq.y)
    i <- length(seq.x)
    j <- length(seq.y)
    while (i > 1 && j >1){
        sc <- S[i-1,j-1]
        
        if (seq.x[i] == seq.y[j]) {
            sc <- sc + match
        } else {
            sc <- sc + mismatch
        }
        
        if (sc == S[i,j]) {
            ax <- c(seq.x[i], ax)
            ay <- c(seq.y[j], ay)
            i <- i -1
            j <- j-1
            next
        }
        
        if ((S[i-1,j] + gap) == S[i,j]) {
            ax <- c(seq.x[i], ax)
            ay <- c("-", ay)
            i <- i-1
            next
        }
        
        if ((S[i,j-1] + gap) == S[i,j]) {
            ax <- c("-", ax)
            ay <- c(seq.y[j], ay)
            j <- j-1
            next
        }
    }
    opt = max(M[nchar(x)+1,nchar(y)+1], X[nchar(x)+1,nchar(y)+1], Y[nchar(x)+1,nchar(y)+1])
    cat ("Optimal alignment:\n",ax,"\n",ay,"\n","Score of Optimal Alignment = ", opt)
}

# global_align("ACTGATTCA","ACGCATCA")

seq1<-readline(prompt="enter sequence 1")
seq2<-readline(prompt="enter sequence 2")
global_align(seq1,seq2)
