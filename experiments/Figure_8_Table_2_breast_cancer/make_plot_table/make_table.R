genes = NULL
count = rep(0,100)
N = 100



for (arg in 1:N){
    filename = sprintf("output_real/output_%i.Rdata",arg)

    if (file.exists(filename)) {
        load(filename)
        for (gene in set_crt_seqstep){
            if (gene %in% genes){
                count[which(genes == gene)] = count[which(genes == gene)] + 1
            } else{
                genes = c(genes, gene)
                count[length(genes)] = 1
            }
        }
    }
    
}

count = count[1:length(genes)]
order = order(count, decreasing = TRUE)
genes[order]
count[order]




