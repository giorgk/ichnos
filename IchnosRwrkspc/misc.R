library(stringr)

convertVelocityFiles <- function(filename, radius = 310, power = 3, threshold = 0.1,interpOut = 1, nmax = 100, nmin = 10, nproc = 6 ){
    for (i in 0:(nproc-1)) {
        print(paste0(filename,str_pad(i, 4, pad = "0"),".vel"))
        tmp <- read.table(file = paste0(filename,str_pad(i, 4, pad = "0"),".vel"), sep = " ")
        newfile <- paste0(filename,str_pad(i, 4, pad = "0"),".ich")
        write(c(radius, power, threshold, interpOut, nmax, nmin), file = newfile, sep = " ", append = F, ncolumns = 6)
        write.table(tmp[, c(7,1:6)],file = newfile,
                    sep = " ", col.names = F, row.names = F, append = T)
    }
}

convertVelocityFiles(filename = "f:/UCDAVIS/ichnos/NPSAT_example/example1__")
