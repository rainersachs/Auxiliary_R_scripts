V=1:3
syn_check=barplot(V, col=c('red','green','blue'))
W=rep(.2,3)
barplot(2*V/3, add=TRUE,col='black')
barplot(2*V/3 -.005, col=c('red','green','blue'),add=TRUE)
barplot(W,add=TRUE,col='orange')
