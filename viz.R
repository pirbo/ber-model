#library(scales)
library(ggplot2)
#library(plyr)
#library(pracma)
require(splines)
#library(reshape)

setwd("/Users/ramdiv/Documents/work/ptb-vs-rep-paper/")
loadr <- function(dir, grp, fun = identity) {
  setwd(dir)
  temp = list.files(pattern="*.out")
  i = 1
  data <- data.frame()
  while (i <= length(temp)) {
    print(sprintf("reading %d-th of %d\n", i, length(temp)))
    raw <- readChar(temp[i], file.info(temp[i])$size)
    raw = substr(raw, 2, nchar(raw))
    con <- textConnection(raw)
    tmp <- read.table(con, header=T, quote="'")
    n <- length(tmp$time)
    tmp$group <- rep(grp, n)
    tmp$serial <- rep(i, n)
#     tmp$fap_integ <- unlist(lapply(1:length(tmp$time), 
#                                    function(i) trapz(tmp$time[1:i], tmp$free_AP_.Gaped._dgless[1:i])))
    model <- lm(healed~ns(time,3),data=tmp)
    Y <- predict(model,newdata=tmp)
    dY <- diff(Y)/diff(tmp$time)
    tmp$healed_deriv <- unlist(c(0, dY))
    tmp = fun(tmp)
    data = rbind(tmp, data)
    close(con)
    i = i + 1
  }
  
  setwd("/Users/ramdiv/Documents/work/ptb-vs-rep-paper/experiments")
  return(data)
}

setwd("/Users/ramdiv/Documents/work/ptb-vs-rep-paper/experiments")
data <- rbind(loadr("wo_default", "no XRCC1"), loadr("w_default", "with XRCC1"))


p <- ggplot(data, aes(x=time, y=full_cytotoxic, group=group)) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth", aes(color=group)) + 
  xlab("time") + ylab("count") + theme_bw() + ggtitle("Free cytotoxic sites")
ggsave(plot=p, filename="full_cytotoxic.png")

p <- ggplot(data, aes(x=time, y=healed, group=group)) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth", aes(color=group)) + 
  xlab("time") + ylab("count") + theme_bw() + ggtitle("Number of healed bp")
ggsave(plot=p, filename="healed.png")

p <- ggplot(data, aes(x=time, y=healed_deriv, group=group)) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth", aes(color=group)) + 
  xlab("time") + ylab("y") + theme_bw() + ggtitle("Derivative of (smoothed) number of healed bp")
ggsave(plot=p, filename="healed_deriv.png")

p <- ggplot(data, aes(x=time, y=ape_perf, group=group)) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth", aes(color=group)) + 
  xlab("time") + ylab("y") + theme_bw() + ggtitle("APE1 performance")
ggsave(plot=p, filename="ape_perf.png")

p <- ggplot(data, aes(x=time, y=POLb_act_lyase, group=group)) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth", aes(color=group)) + 
  xlab("time") + ylab("y") + theme_bw() + ggtitle("POLb lyase performance")
ggsave(plot=p, filename="polb_lyase_perf.png")

p <- ggplot(data, aes(x=time, y=POLb_act_poly, group=group)) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth", aes(color=group)) + 
  xlab("time") + ylab("y") + theme_bw() + ggtitle("POLb nuclease performance")
ggsave(plot=p, filename="polb_nucl_perf.png") 

p <- ggplot(data, aes(x=time, y=nicked_DNA, group=group)) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth", aes(color=group)) + 
  xlab("time") + ylab("y") + theme_bw() + ggtitle("Nicks in the DNA")
ggsave(plot=p, filename="nicked_dna.png") 

p <- ggplot(data, aes(x=time, y=free_nicked_P, group=group)) +
  stat_summary(fun.data="mean_cl_normal", geom="smooth", aes(color=group)) + 
  xlab("time") + ylab("y") + theme_bw() + ggtitle("Free (P-)nicks in the DNA")
ggsave(plot=p, filename="free_pnicked_dna.png") 

# data2 <- loadr("wo_polb_2", "no XRCC1", fun=function(tmp) {tmp$polb_count <- rep("/2", length(tmp$time)); tmp})
# data2 = rbind(loadr("wo_polb_neut", "no XRCC1", fun=function(tmp) {tmp$polb_count <- rep("o", length(tmp$time)); tmp}), data2)
# data2 = rbind(loadr("wo_polb_4", "no XRCC1", fun=function(tmp) {tmp$polb_count <- rep("/4", length(tmp$time)); tmp}), data2)
# data2 = rbind(loadr("wo_polb_8", "no XRCC1", fun=function(tmp) {tmp$polb_count <- rep("/8", length(tmp$time)); tmp}), data2)
# f <- c();
# data2 = rbind(loadr("w_polb_2", "with XRCC1", fun=function(tmp) {tmp$polb_count <- rep("/2", length(tmp$time)); tmp}), data2)
# data2 = rbind(loadr("w_polb_neut", "with XRCC1", fun=function(tmp) {tmp$polb_count <- rep("o", length(tmp$time)); 
#                                                                     tmp[setdiff(colnames(tmp), colnames(data2))] <- list(NULL);
#                                                                     tmp}), data2)
# data2 = rbind(loadr("w_polb_4", "with XRCC1", fun=function(tmp) {tmp$polb_count <- rep("/4", length(tmp$time)); tmp}), data2)
# data2 = rbind(loadr("w_polb_8", "with XRCC1", fun=function(tmp) {tmp$polb_count <- rep("/8", length(tmp$time)); tmp}), data2)
# 
# data2$polb_count <- factor(data2$polb_count, levels=c("o", "/2", "/4", "/8"))
# 
# print("var_healed")
# p <- ggplot(data2, aes(x=time, y=healed, group=group)) +
#   facet_grid(polb_count ~ ., scales = "free_y") + 
#   stat_summary(fun.data="mean_cl_normal", geom="smooth", aes(color=group)) + 
#   xlab("time") + ylab("count") + theme_bw() + ggtitle("Number of healed bp")
# png(filename="var_healed.png", width= 6, height = 5, units="in", res=300)
# show(p)
# dev.off()
# 
# print("var_healed_deriv")
# p <- ggplot(data2, aes(x=time, y=healed_deriv, group=group)) +
#   facet_grid(polb_count ~ ., scales = "free_y") + 
#   stat_summary(fun.data="mean_cl_normal", geom="smooth", aes(color=group)) + 
#   xlab("time") + ylab("count") + theme_bw() + ggtitle("Derivative of (smoothed) number of healed bp")
# png(filename="var_healed_deriv.png", width= 6, height = 5, units="in", res=300)
# show(p)
# dev.off()
# 
# print("var_nicked")
# p <- ggplot(data2, aes(x=time, y=nicked_DNA, group=group)) +
#   facet_grid(polb_count ~ ., scales = "free_y") + 
#   stat_summary(fun.data="mean_cl_normal", geom="smooth", aes(color=group)) + 
#   xlab("time") + ylab("count") + theme_bw() + ggtitle("Number of nicks in DNA")
# png(filename="var_nicked.png", width= 6, height = 5, units="in", res=300)
# show(p)
# dev.off()
# 
# print("var_free_nicked")
# data2$free_cyto <- data2$nicked_DNA - data2$LIG3_on_ligatable
# p <- ggplot(data2, aes(x=time, y=free_cyto, group=group)) +
#   facet_grid(polb_count ~ ., scales = "free_y") + 
#   stat_summary(fun.data="mean_cl_normal", geom="smooth", aes(color=group)) + 
#   xlab("time") + ylab("count") + theme_bw() + ggtitle("Approx. number of cytotoxic nicks in DNA")
# png(filename="var_free_nicked.png", width= 6, height = 5, units="in", res=300)
# show(p)
# dev.off()