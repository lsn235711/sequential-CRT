num_disco = NULL
Method = NULL
N = 100

for (arg in 1:N){
    filename = sprintf("output_real/output_%i.Rdata",arg)

    if (file.exists(filename)) {
        load(filename)
        Method = c(Method, "d0CRT", "d1CRT", "HRT", "knockoffs", "CRT_seqstep+")
        dic = c(length(set_d0), length(set_d1), length(set_hrt), length(set_knockoffs), length(set_crt_seqstep))
        num_disco = c(num_disco, dic)
    }
    
}

dat_toplot = data.frame(Method = Method, num_disco = num_disco)
dat_toplot$Method = factor(Method, levels = c("d0CRT", "d1CRT", "HRT", "knockoffs", "CRT_seqstep+"))
p1 = ggplot(dat_toplot, aes(x=Method, y=num_disco, fill= Method, color = Method)) +
    geom_boxplot(alpha = 0.3) + 
    ylab("Number of Discoveries") +
    theme_bw() +
    scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") 
ggsave("real_data.pdf",p1, width=6, height=2.5, units="in")

