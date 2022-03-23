library(dplyr)
library(stringr)
library(ggplot2)

std <- read.csv("STD_SIG.csv")

sample <- read.csv("ALL_SIG_METABS.csv") %>% mutate(code=str_sub(FileName,-7,-5))


std2 <- std %>% filter(Name=="D-(+)-Cellobiose, (isomer 2), 8TMS derivative" | 
                         Name=="Silanol, trimethyl-, phosphate (3:1)")

methanol <- std2[c(1,2,42,43),]
methanol2 <- methanol %>% group_by(FileName,Name) %>% summarise(median=median(Area))

water <- std2[c(3,4,44,45),]
water2 <- water %>% group_by(FileName,Name) %>% summarise(median=median(Area))


blank2 <- sample %>% filter(Name=="D-(+)-Cellobiose, (isomer 2), 8TMS derivative" | 
                               Name=="Silanol, trimethyl-, phosphate (3:1)") %>%
  filter(code=="003") %>%
  group_by(Name) %>%
  summarise(median=median(Area))


sample2 <- sample %>% filter(Name=="D-(+)-Cellobiose, (isomer 2), 8TMS derivative" | 
                               Name=="Silanol, trimethyl-, phosphate (3:1)") %>%
  filter(!(code %in% c("001","002","003","004","005","012","019","026","027",
                       "028","029","030","031","032"))) %>%
  group_by(Name) %>%
  summarise(median=median(Area))

check <- sample %>% filter(!(code %in% c("001","002","028","029","030","031","032")))

length(unique(check$FileName))

check2 <- sample %>% filter(code=="003") 
length(unique(check2$FileName))
########################################

blank3 <- sample %>% filter(Name=="D-(+)-Cellobiose, (isomer 2), 8TMS derivative" | 
                              Name=="Silanol, trimethyl-, phosphate (3:1)") %>%
  filter(code=="003")

sample3 <- sample %>% filter(Name=="D-(+)-Cellobiose, (isomer 2), 8TMS derivative" | 
                               Name=="Silanol, trimethyl-, phosphate (3:1)") %>%
  filter(!(code %in% c("003","004","005","012","019","026","027")))


ggplot(blank3[blank3$Name=="D-(+)-Cellobiose, (isomer 2), 8TMS derivative",], aes(x=1:55,y=Area)) + 
  geom_point() + geom_hline(yintercept = 6753184) + geom_hline(yintercept = 7639531) + 
  ggtitle ("D-Cellobiose") + theme_minimal() + 
  scale_y_continuous(name="Area", limits=c(1e+05, 6e+07)) 


ggplot(blank3[blank3$Name=="Silanol, trimethyl-, phosphate (3:1)",], aes(x=1:55,y=Area)) +
  geom_point() + geom_hline(yintercept = 208778927) + geom_hline(yintercept = 301517027) +
  ggtitle ("Silanol") + theme_minimal()


ggplot(sample3[sample3$Name=="D-(+)-Cellobiose, (isomer 2), 8TMS derivative",], aes(x=1:966,y=Area)) + 
  geom_point(colour=I("Red")) + geom_hline(yintercept = 6753184) + geom_hline(yintercept = 7639531) + 
  ggtitle ("D-Cellobiose") + theme_minimal()


ggplot(sample3[sample3$Name=="Silanol, trimethyl-, phosphate (3:1)",], aes(x=1:1029,y=Area)) + 
  geom_point() + geom_hline(yintercept = 6753184) + geom_hline(yintercept = 7639531) + 
  ggtitle ("D-Cellobiose") + theme_minimal()
