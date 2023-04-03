procrustes_trans <-function(z,w,ll) {
        # take the first chain
        ll0 <- ll[[1]]
        zz0 <- z[[1]]
        ww0 <- w[[1]]
        
        # define sample sizes
        M = nrow(ll0)
        nz = ncol(zz0)/2
        nw = ncol(ww0)/2
        
        # create a list of M configurations
        zz <- vector("list", M)
        for (i in 1:M) {
                zz[[i]] <- cbind(zz0[i,1:N], zz0[i,(N+1):(2*N)])
        }
        
        ww <- vector("list", M)
        for (i in 1:M) {
                ww[[i]] <- cbind(ww0[i,1:I], ww0[i,(I+1):(2*I)])
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
