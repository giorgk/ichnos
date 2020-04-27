library(rgdal)
library(gwtools)
library(interp)

# Read the polygon domain
c2vsim_outline <- readOGR(dsn = "../../C2VSIM_FG_OR/C2Vsim_FG_v2/gis_data", layer = "C2VsimFG_Outline_26910")

{# Write the domain 
    filename <- "examples/c2vsimDomain.ich"
    con <-  file(filename, open = "w")
    apnd = F
    for (i in 1:length(c2vsim_outline@polygons)) {
        for (j in 1:length(c2vsim_outline@polygons[[i]]@Polygons)) {
            np <- dim(c2vsim_outline@polygons[[i]]@Polygons[[j]]@coords)[1]
            if (c2vsim_outline@polygons[[i]]@Polygons[[j]]@hole == F){
                orien = 1
            }
            else{
                orien = 0
            }
            write(c(np, orien), file = con, append = apnd)
            if (apnd == F){
                apnd = T
            }
            write.table(c2vsim_outline@polygons[[i]]@Polygons[[j]]@coords, file = con, append = apnd, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
        }
    }
    close(con)
}

{# Write the top and bottom elevation
    nd <- gwtools::c2vsim.readNodes(filename = "../../C2VSIM_FG_OR/C2Vsim_FG_v2/c2vsimfg_beta2_publicrelease/C2VSimFG_BETA2_PublicRelease/Preprocessor/C2VSimFG_Nodes.dat",
                                    Nskip = 90,ND = 30179)
    strat <- gwtools::c2vsim.readStrat(filename = "../../C2VSIM_FG_OR/C2Vsim_FG_v2/c2vsimfg_beta2_publicrelease/C2VSimFG_BETA2_PublicRelease/Preprocessor/C2VSimFG_Stratigraphy.dat",
                                       Nnodes = 30179, Nlay = 4, nSkip = 105)
    bottom <- strat$ELV - apply(strat[,3:10],1,sum)
    bottom <- bottom*0.3048
    bottomFile <- "examples/c2vsimBottom.ich"
    r <-  2500
    p <-  4
    t  <- 0.1
    write(c(r,p,t), file = bottomFile, append = FALSE)
    write.table(cbind(nd[,2:3], bottom), file = bottomFile, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
    
    # For the top we will use the ground surface elevation which is easier to setup.
    # However we should use the minimum of the water table and ground surface elevation
    
    top <- strat$ELV*0.3048
    topFile <- "examples/c2vsimTop.ich"
    r <-  2500
    p <-  4
    t  <- 0.1
    write(c(r,p,t), file = topFile, append = FALSE)
    write.table(cbind(nd[,2:3], top), file = topFile, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
}

{ # Read velocity
    velfile <- readLines("../../C2VSIM_FG_OR/C2Vsim_FG_v2/c2vsimfg_beta2_publicrelease/C2VSimFG_BETA2_PublicRelease/Results/Vel_out.out")
    bc_xy <- matrix(data = NA, nrow = 32537, ncol = 2)
    for (i in 1:32537) {
       t <- strsplit(substr(velfile[i+5],1,100), split = " ")[[1]]
       tt <- as.numeric((t[which(t != "")]))
       bc_xy[i,] <- tt[2:3]
    }
    Ntimesteps <- 504
    VX <- array(data = NA, dim = c(32537,504,4))
    VY <- array(data = NA, dim = c(32537,504,4))
    VZ <- array(data = NA, dim = c(32537,504,4))
    TimeMat <- matrix(data = NA, nrow = 504, ncol = 3)
    ix <- seq(2,13,3)
    iy <- seq(3,13,3)
    iz <- seq(4,13,3)
    iline <- 32551
    for (i in 1:504){
        for (j in 1:32537) {
            t <- strsplit(substr(velfile[iline],1,1000),split = " ")[[1]]
            if (j == 1){
                tv <- t[-1]
                tm <- t[1]
                print(tm)
                tm <- as.numeric(strsplit(strsplit(tm,split = "_")[[1]][1], split = "/")[[1]])
                TimeMat[i,] <- tm
            }
            else{
                tv = t
            }
            tv <- as.numeric((tv[which(tv != "")]))
            for (k in 1:4) {
                VX[tv[1],i,k] <- tv[ix[k]]
                VY[tv[1],i,k] <- tv[iy[k]]
                VZ[tv[1],i,k] <- tv[iz[k]]
            }
            iline = iline + 1
        }
    }
    save(bc_xy, VX,VY,VZ,TimeMat, file = "examples/C2Vsim_VelocityField.RData")
}

{# Write x y barycenter locations
    write.table(bc_xy*0.3048, file = "examples/XYnodes.ich", sep = " ", row.names = F, col.names = F,append = F)
}


{# interpolate the stratigraphy values on the element barycenters
    # Top and bottom have been calculated previously
    L1 <- (strat$ELV - apply(strat[,3:4],1,sum))*0.3048
    L2 <- (strat$ELV - apply(strat[,3:6],1,sum))*0.3048
    L3 <- (strat$ELV - apply(strat[,3:8],1,sum))*0.3048
    bc_top <- interp(x = nd$X, y = nd$Y, z = top, xo = bc_xy[,1]*0.3048, yo = bc_xy[,2]*0.3048, input = "points", output = "points" )
    bc_L1 <- interp(x = nd$X, y = nd$Y, z = L1, xo = bc_xy[,1]*0.3048, yo = bc_xy[,2]*0.3048, input = "points", output = "points" )
    bc_L2 <- interp(x = nd$X, y = nd$Y, z = L2, xo = bc_xy[,1]*0.3048, yo = bc_xy[,2]*0.3048, input = "points", output = "points" )
    bc_L3 <- interp(x = nd$X, y = nd$Y, z = L3, xo = bc_xy[,1]*0.3048, yo = bc_xy[,2]*0.3048, input = "points", output = "points" )
    bc_bot <- interp(x = nd$X, y = nd$Y, z = bottom, xo = bc_xy[,1]*0.3048, yo = bc_xy[,2]*0.3048, input = "points", output = "points" )
    
    # write.table(cbind(nd$X,nd$Y,top), file = "examples/fileA.qwer",sep = " ", row.names = F, col.names = F,append = F)
    # write.table(bc_xy*0.3048, file = "examples/fileB.qwer",sep = " ", row.names = F, col.names = F,append = F)
    write.table(bc_xy*0.3048, file = "examples/XY.dat",sep = " ", row.names = F, col.names = F,append = F)
    write.table(cbind(bc_top, bc_L1, bc_L2, bc_L3, bc_bot), file = "examples/LayElev.dat",sep = " ", row.names = F, col.names = F,append = F)
}

{# Write elevation per processor
    prefix <- "examples/iwfm_"
    nproc <- 1
    ElevMat <- cbind(XYid$id,bc_top$z, bc_L1$z, bc_L2$z, bc_L3$z, bc_bot$z)
    for (i in 1:nproc) {
        xyfl <- paste0(prefix,"XY_",i-1,".ich")
        XYid <- read.table(xyfl, header = F, col.names = c("id", "X", "Y", "Proc"), skip = 1)
        elfl <- paste0(prefix,"Elev_",i-1,".ich")
        write.table(ElevMat[XYid$id+1,], file = elfl, sep = " ", row.names = F, col.names = F,append = F)
        
    }
}

{# Write the velocity per processor
    prefix <- "examples/iwfm_"
    nproc <- 1
    istart = which(TimeMat[,3] == 2005 & TimeMat[,1] == 10)
    iend = which(TimeMat[,3] == 2015 & TimeMat[,1] == 9)
    print(paste("NTsteps:", length(istart:iend)))
    for (i in 1:nproc) {
        xyfl <- paste0(prefix,"XY_",i-1,".ich")
        XYid <- read.table(xyfl, header = F, col.names = c("id", "X", "Y", "Proc"), skip = 1)
        for (lay in 1:4) {
            print(paste("Proc:",i,"Lay:", lay))
            fl <- paste0(prefix,"VX_LAY", lay-1, "_", i-1,".ich")
            write.table(cbind(XYid$id, VX[XYid$id+1,istart:iend,lay]), file = fl, sep = " ", row.names = F, col.names = F, append = F)
            fl <- paste0(prefix,"VY_LAY", lay-1, "_", i-1,".ich")
            write.table(cbind(XYid$id, VY[XYid$id+1,istart:iend,lay]), file = fl, sep = " ", row.names = F, col.names = F, append = F)
            fl <- paste0(prefix,"VZ_LAY", lay-1, "_", i-1,".ich")
            write.table(cbind(XYid$id, VZ[XYid$id+1,istart:iend,lay]), file = fl, sep = " ", row.names = F, col.names = F, append = F)
        }
    }
}

{
    for (i in 1:4) {
        print(i)
        flvx <- paste0("examples/VX_Lay",i,".dat")
        write.table(1233.48*VX[,,i], file = flvx, sep = " ", row.names = F, col.names = F, append = F)
        flvy <- paste0("examples/VY_Lay",i,".dat")
        write.table(1233.48*VY[,,i], file = flvy, sep = " ", row.names = F, col.names = F, append = F)
        flvz <- paste0("examples/VZ_Lay",i,".dat")
        write.table(1233.48*VZ[,,i], file = flvz, sep = " ", row.names = F, col.names = F, append = F)
    }
}

