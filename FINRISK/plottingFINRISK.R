
# Importing Packages ------------------------------------------------------

library(ggplot2)
library(plyr)
library(reshape2)
library(dplyr)
library(tidyr)
library(data.table)
library(purrr)


# Reading in Data ---------------------------------------------------------

setwd("~/Dropbox (Partners HealthCare)/Susan-Mir Shared/Eicosanoids-FS-CODE/eicdata/FINRISK")

# SLIGHTLY DIFFERENT FOR FINRISK
clean_df <- function(df){
  df <- t(df)
  df <- df[,-1] %>% as.data.frame()
  rownames(df) <- NULL
  return(df)
}

df <- read.csv('./NormalizationAK/FR02_eic_data.csv')[,-c(1:2)]
#mzids <- colnames(df)

df_a <- read.csv('./NormalizationAK/FR02_centr_a.csv', row.names = 1) %>% clean_df()
df_b <- read.csv('./NormalizationAK/FR02_centr_b.csv', row.names = 1) %>% clean_df()
df_cb <- read.csv('./NormalizationAK/FR02_centr_cb.csv', row.names = 1) %>% clean_df()
df_EigenMS <- read.csv('./NormalizationAK/FR02_centr_EigenMS.csv', row.names = 1) %>% clean_df()
df_q <- read.csv('./NormalizationAK/FR02_centr_q.csv', row.names = 1) %>% clean_df()
df_qcb <- read.csv('./NormalizationAK/FR02_centr_qcb.csv', row.names = 1) %>% clean_df()
df_ruv2 <- read.csv('./NormalizationAK/FR02_centr_ruv2.csv', row.names = 1) %>% clean_df()
df_ruv4 <- read.csv('./NormalizationAK/FR02_centr_ruv4.csv', row.names = 1) %>% clean_df()

setwd("~/Dropbox (Partners HealthCare)/2018 Applied Bioinformatics Work/Personal folders/Andy Kim/NormalizationYY/FINRISK")

testdf <- df %>%
  keep(is.numeric) %>%                     # Keep only numeric columns
  gather() %>%                             # Convert to key-value pairs
  mutate(index = rep_len(1:8293, length.out=nrow(df)*ncol(df)))

# Plotting ----------------------------------------------------------------


df %>%
  keep(is.numeric) %>%                     # Keep only numeric columns
  gather() %>%                             # Convert to key-value pairs
  mutate(index = rep_len(1:nrow(df), length.out=nrow(df)*ncol(df))) %>%               # Create index column
  ggplot(aes(index, value)) +                     # Plot the values
  facet_wrap(~ key, scales = "free") +   # In separate panels
  geom_point(size = 0.1) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"cm"))

ggsave("FR02_eicdata.png",height=160,width=160,dpi=100,limitsize = F)

df_a %>%
  keep(is.numeric) %>%                     # Keep only numeric columns
  gather() %>%                             # Convert to key-value pairs
  mutate(index = rep_len(1:nrow(df), length.out=nrow(df)*ncol(df))) %>%               # Create index column
  ggplot(aes(index, value)) +                     # Plot the values
  facet_wrap(~ key, scales = "free") +   # In separate panels
  geom_point(size = 0.1) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"cm"))

ggsave("FR02_norm_mean.png",height=160,width=160,dpi=72,limitsize = F)

df_b %>%
  keep(is.numeric) %>%                     # Keep only numeric columns
  gather() %>%                             # Convert to key-value pairs
  mutate(index = rep_len(1:nrow(df), length.out=nrow(df)*ncol(df))) %>%               # Create index column
  ggplot(aes(index, value)) +                     # Plot the values
  facet_wrap(~ key, scales = "free") +   # In separate panels
  geom_point(size = 0.1) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"cm"))

ggsave("FR02_norm_median.png",height=160,width=160,dpi=72,limitsize = F)

df_q %>%
  keep(is.numeric) %>%                     # Keep only numeric columns
  gather() %>%                             # Convert to key-value pairs
  mutate(index = rep_len(1:nrow(df), length.out=nrow(df)*ncol(df))) %>%               # Create index column
  ggplot(aes(index, value)) +                     # Plot the values
  facet_wrap(~ key, scales = "free") +   # In separate panels
  geom_point(size = 0.1) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"cm"))

ggsave("FR02_norm_quantile.png",height=160,width=160,dpi=72,limitsize = F)

df_cb %>%
  keep(is.numeric) %>%                     # Keep only numeric columns
  gather() %>%                             # Convert to key-value pairs
  mutate(index = rep_len(1:nrow(df), length.out=nrow(df)*ncol(df))) %>%               # Create index column
  ggplot(aes(index, value)) +                     # Plot the values
  facet_wrap(~ key, scales = "free") +   # In separate panels
  geom_point(size = 0.1) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"cm"))

ggsave("FR02_norm_combat.png",height=160,width=160,dpi=72,limitsize = F)

df_qcb %>%
  keep(is.numeric) %>%                     # Keep only numeric columns
  gather() %>%                             # Convert to key-value pairs
  mutate(index = rep_len(1:nrow(df), length.out=nrow(df)*ncol(df))) %>%               # Create index column
  ggplot(aes(index, value)) +                     # Plot the values
  facet_wrap(~ key, scales = "free") +   # In separate panels
  geom_point(size = 0.1) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"cm"))

ggsave("FR02_norm_quantilecombat.png",height=160,width=160,dpi=72,limitsize = F)

df_ruv2 %>%
  keep(is.numeric) %>%                     # Keep only numeric columns
  gather() %>%                             # Convert to key-value pairs
  mutate(index = rep_len(1:nrow(df), length.out=nrow(df)*ncol(df))) %>%               # Create index column
  ggplot(aes(index, value)) +                     # Plot the values
  facet_wrap(~ key, scales = "free") +   # In separate panels
  geom_point(size = 0.1) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"cm"))

ggsave("FR02_norm_ruv2.png",height=160,width=160,dpi=72,limitsize = F)

df_ruv4 %>%
  keep(is.numeric) %>%                     # Keep only numeric columns
  gather() %>%                             # Convert to key-value pairs
  mutate(index = rep_len(1:nrow(df), length.out=nrow(df)*ncol(df))) %>%               # Create index column
  ggplot(aes(index, value)) +                     # Plot the values
  facet_wrap(~ key, scales = "free") +   # In separate panels
  geom_point(size = 0.1) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"cm"))

ggsave("FR02_norm_ruv4.png",height=160,width=160,dpi=72,limitsize = F)

df_EigenMS %>%
  keep(is.numeric) %>%                     # Keep only numeric columns
  gather() %>%                             # Convert to key-value pairs
  mutate(index = rep_len(1:nrow(df), length.out=nrow(df)*ncol(df))) %>%               # Create index column
  ggplot(aes(index, value)) +                     # Plot the values
  facet_wrap(~ key, scales = "free") +   # In separate panels
  geom_point(size = 0.1) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"cm"))

ggsave("FR02_norm_EigenMS.png",height=160,width=160,dpi=72,limitsize = F)
