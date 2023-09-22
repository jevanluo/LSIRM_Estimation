##################################################
## Project: LSIRM postsamples matching & visualization
## Created: Feb 28th, 2023
## Last Edit: May 8, 2023
## Author: Jinwen Luo, Minjeong Jeon
##################################################
library(coda)
library(ggplot2)

procrustes_trans <-function(z,w,ll,stan=FALSE) {
        # take the first chain
        ll0 <- ll[[1]]
        zz0 <- z[[1]]
        ww0 <- w[[1]]
        
        # define sample sizes
        M = nrow(ll0)
        nz = ncol(zz0)/2
        nw = ncol(ww0)/2
        
        # create a list of M configurations
        if(stan==TRUE){
                zz <- vector("list", M)
                for (i in 1:M) {
                        zz[[i]] <- cbind(zz0[i,2*(1:nz-1)+1], zz0[i,2*(1:nz)])
                }
                
                ww <- vector("list", M)
                for (i in 1:M) {
                        ww[[i]] <- cbind(ww0[i,2*(1:nw-1)+1], ww0[i,2*(1:nw)])
                }
                
        }else{
                zz <- vector("list", M)
                for (i in 1:M) {
                        zz[[i]] <- cbind(zz0[i,1:nz], zz0[i,(nz+1):(2*nz)])
                }
                
                ww <- vector("list", M)
                for (i in 1:M) {
                        ww[[i]] <- cbind(ww0[i,1:nw], ww0[i,(nw+1):(2*nw)])
                }
                
        }
        
        ll <- apply(ll0, 1,  sum)
        ll=ll[1:M]
        # find the iteration with largest ll
        imax=which(ll==max(ll,na.rm=T))
        
        # center the wz matrix that achieve the highest ll 
        wzmax=rbind(zz[[imax]],ww[[imax]]) # add a mu 
        mean_pos=apply(wzmax,2,mean)
        wz0=sweep(wzmax,2,mean_pos)
        # procrustes rotations through MCMCpack
        matched_zz = zz
        matched_ww = ww
        
        # R package 
        require(MCMCpack)
        for (ii in 1: M) {
                
                # target matrix: wz0 (centered)
                #center first
                wz=rbind(zz[[ii]],ww[[ii]])
                mean_pos=apply(wz,2,mean)
                wz_centered=sweep(wz,2,mean_pos)
                
                proc=procrustes(wz_centered, wz0, translation=TRUE, dilation=TRUE)
                
                wz_matched0  <- proc$X.new
                
                # results 
                matched_zz[[ii]]=wz_matched0[1:nz,]
                matched_ww[[ii]]=wz_matched0[(nz+1):(nz+nw),]    
        }
        
        return( list(imax=imax, zz=matched_zz, ww=matched_ww))
}  

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

plot_strength <- function(zz,ww,PersonIndex,ItemGroup=NULL,ItemName=NULL,title=NULL){
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
        ### calculate the distances and inverse distances
        require(pdist)
        strength<-pdist(z[PersonIndex,],w)
        
        if(is.null(ItemGroup)==TRUE){
                p1strength <- data.frame(cbind(item=rownames(w),
                                               strength=1/as.numeric(strength@dist)))
                p1strength$strength <- as.numeric(p1strength$strength)
                z=data.frame(z)
                w=data.frame(w)
                names(z)=c('coord1','coord2')
                names(w)=c('coord1','coord2')
                p0<-ggplot(p1strength, aes(x = reorder(item,-strength), y = strength)) + 
                        geom_col(aes(fill = "grey")) +
                        theme(axis.ticks.y =element_blank(),
                              axis.text.y =element_blank(),
                              axis.title.y=element_blank(),
                              axis.title.x = element_blank(),
                              axis.text=element_text(size=20,angle = 90,face="bold"),
                              legend.position = "none")
                print(p0)   
        }else{
                p1strength <- data.frame(cbind(item=rownames(w),
                                               strength=1/as.numeric(strength@dist),
                                               behavior=ItemGroup))
                p1strength$strength <- as.numeric(p1strength$strength)
                z=data.frame(z)
                w=data.frame(w)
                names(z)=c('coord1','coord2')
                names(w)=c('coord1','coord2')
                p0<-ggplot(p1strength, aes(x = reorder(item,-strength), y = strength)) + 
                        geom_col(aes(fill = behavior)) +
                        theme(axis.ticks.y =element_blank(),
                              axis.text.y =element_blank(),
                              axis.title.y=element_blank(),
                              axis.title.x = element_blank(),
                              axis.text=element_text(size=20,angle = 90,face="bold"),
                              legend.position = "none")
                print(p0)  
        }
        
}
