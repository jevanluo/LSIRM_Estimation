procrustes_trans <-function(z,w,ll) {
        # find the iteration m with max likelihood 
        nz=dim(z[[1]])[1]
        nw=dim(w[[1]])[1]
        
        M=length(ll)#sum(!is.na(stored_likelihoods))
        ll=ll[1:M]
        imax=which(ll==max(ll,na.rm=T))
        
        # center the wz matrix that achieve the highest ll 
        wzmax=rbind(z[[imax]],w[[imax]]) # add a mu 
        mean_pos=apply(wzmax,2,mean)
        wz0=sweep(wzmax,2,mean_pos)
        # procrustes rotations through MCMCpack
        matched_zz = z
        matched_wz = w 
        # procrustes rotations manually
        matched_zm = z
        matched_wm = w
        
        for (ii in 1: M) {
                
                # target matrix: wz0 (centered)
                #center first
                wz=rbind(z[[ii]],w[[ii]])
                mean_pos=apply(wz,2,mean)
                wz_centered=sweep(wz,2,mean_pos)
                
                # R package 
                require(MCMCpack)
                proc=procrustes(wz_centered, wz0, translation=TRUE, dilation=TRUE)
                
                wz_matched0  <- proc$X.new
                mat_rot  <- proc$R
                mat_t    <- proc$t
                mat_s    <- proc$s
                
                # results 
                matched_zz[[ii]]=wz_matched0[1:nz,]
                matched_ww[[ii]]=wz_matched0[(nz+1):(nz+nw),]    
                
                # manual : translation
                require(expm)
                rot_mat=t(wz0)%*%wz_centered%*%solve(sqrtm(t(wz_centered)%*%wz0%*%t(wz0)%*%wz_centered))
                wz_matched=wz_centered%*%rot_mat
                
                # results 
                matched_zm[[ii]]=wz_matched[1:nz,]
                matched_wm[[ii]]=wz_matched[(nz+1):(nz+nw),]
                
        }
        
        return( list(imax=imax, zz=matched_zz, ww=matched_ww, zm=matched_zm, wm=matched_wm))
}  