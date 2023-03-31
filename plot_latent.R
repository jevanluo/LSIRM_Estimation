library(ggplot2)
plot_latent <- function(zz,ww,ItemGroup=NULL,ItemName=NULL,title=NULL){
        M.keep <- length(zz)
        N <- nrow(zz[[1]])
        I <- nrow(ww[[1]])
        # calculate the means of latent positions
        arr.z <- array(unlist(zz), c(N, 2, M.keep))
        z = apply(arr.z ,1:2 ,mean)
        arr.w <- array(unlist(ww), c(I, 2, M.keep))
        w = apply(arr.w, 1:2, mean)
        if(!is.null(ItemName)){
                rownames(w) <- ItemName
        }else{
                rownames(w) <- paste0("I", 1:I)
        }
        
        if(is.null(ItemGroup)==TRUE){
                z=data.frame(z)
                w=data.frame(w)
                names(z)=c('coord1','coord2')
                names(w)=c('coord1','coord2')
                p0<-ggplot()
                p0<-p0+geom_point(aes(x=coord1,y=coord2),data=z,cex=1)
                p0<-p0+geom_text(aes(x=coord1,y=coord2,label=rownames(w), colour="w"),data=w ,cex=12, fontface = "bold")+
                        theme(axis.title.x=element_blank(),
                              axis.title.y=element_blank(),
                              axis.text.x=element_blank(),
                              axis.text.y=element_blank(),
                              axis.ticks.x=element_blank(),
                              axis.ticks.y=element_blank(),
                              legend.position = "none")
                print(p0)   
        }else{
                z=data.frame(z)
                w=data.frame(w)
                names(z)=c('coord1','coord2')
                names(w)=c('coord1','coord2')
                p0<-ggplot()
                p0<-p0+geom_point(aes(x=coord1,y=coord2),data=z,cex=1)
                p0<-p0+geom_text(aes(x=coord1,y=coord2,label=rownames(w), colour=ItemGroup),data=w ,cex=12, fontface = "bold")+
                        theme(axis.title.x=element_blank(),
                              axis.title.y=element_blank(),
                              axis.text.x=element_blank(),
                              axis.text.y=element_blank(),
                              axis.ticks.x=element_blank(),
                              axis.ticks.y=element_blank(),
                              legend.position = "none")
                print(p0)  
        }
        
}