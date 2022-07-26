library(ggplot2)
library(dplyr)

df = read.csv("C:\\Users\\klnag\\Downloads\\submission_data_by_primary_lin_7-24.csv")

numbers_only <- df[df$Primary.Lineage.Family %in% c(1, 2, 3, 4),]
numbers_only$Primary.Lineage.Family = as.numeric(numbers_only$Primary.Lineage.Family)

## PROPORTION RESISTANT BY PRIMARY LINEAGE PLOT
ggplot(data = numbers_only, aes(x = Drug, y = Proportion.Resistant, group = Primary.Lineage.Family, fill = Primary.Lineage.Family)) +
  geom_col(position = 'dodge') +
  scale_fill_manual(values = c('red','orange','green','blue','purple','black')) +
  theme_bw()


#
## MDR or not?
#
df2 = read.csv("C:\\Users\\klnag\\Downloads\\full_df_7-18 (1).csv")[, c(2,3,4,5)]

resistances = data.frame(matrix(nrow = length(unique(df2$ID)), ncol = 4))
colnames(resistances) = c('ID', 'Lineage', 'MDR', 'Fluoroquinolones')
resistances$ID = unique(df2$ID)
resistances$Lineage = sapply(resistances$ID, function(ID) unique(df2$Lineage[df2$ID==ID]))

is_mdr <- function(ID) {
  temp = df2[df2$ID == ID & df2$Drug %in% c('rif', 'inh'),]
  return(temp$Resistant[temp$Drug=='rif'] == 1 & temp$Resistant[temp$Drug=='inh'] == 1)
}

resistances$MDR = sapply(resistances$ID, is_mdr)
resistances$Primary_Lineage = sapply(resistances$Lineage, function(x) substr(as.character(x), 1, 1))

labels <- resistances %>%
  group_by(Primary_Lineage, MDR) %>%
  summarise(n=n()) %>%
  mutate(Label = n / sum(n))

only_1234 <- labels[labels$Primary_Lineage %in% c(1, 2, 3, 4),]
only_1234$Label <- round(only_1234$Label, 3)

## MULTI-DRUG RESISTANCE BARPLOT
ggplot(data = only_1234, aes(x = Primary_Lineage, y = Label, fill = MDR)) +
  geom_col(position = 'Stack') +
  scale_fill_manual(values = c('lightgrey', ' chartreuse3'), labels = c('Not MDR', 'MDR')) +
  theme_bw() +
  theme(axis.text.x = element_text(size=14)) +
  xlab('Primary Lineage') +
  ylab('Proportion') +
  geom_text(aes(label = Label), position = position_stack(vjust = .5), size = 15)


# Fluoroquinolones
is_fluoro <- function(ID) {
  temp = df2[df2$ID == ID & df2$Drug %in% c('oflx', 'levo'),]
  if (temp$Resistant[temp$Drug == 'oflx'] == 1 & temp$Resistant[temp$Drug == 'levo'] == 1) {
    return('Both')
  }
  else if (temp$Resistant[temp$Drug == 'oflx'] == 1) {
    return('Oflx Only')
  }
  else if (temp$Resistant[temp$Drug == 'levo'] == 1) {
    return('Levo Only')
  }
  return('Neither')
}

resistances$Fluoroquinolones <- sapply(resistances$ID, is_fluoro)
resistances$Fluoroquinolones <- factor(resistances$Fluoroquinolones, levels = c('Oflx Only', 'Levo Only', 'Both', 'Neither'))

labels1 <- resistances[resistances$Primary_Lineage %in% c(1, 2, 3, 4),] %>%
  group_by(Primary_Lineage, Fluoroquinolones) %>%
  summarise(n=n()) %>%
  mutate(Label = n / sum(n))

labels1$Label = round(labels1$Label, 3)

ggplot(data = labels1, aes(x = Primary_Lineage, y = Label, fill = Fluoroquinolones)) +
  geom_col(position = 'Stack') +
  # scale_fill_manual(values = c('lightgrey', ' chartreuse3'), labels = c('Not MDR', 'MDR')) +
  theme_bw() +
  xlab('Primary Lineage') +
  ylab('Proportion') +
  geom_text(aes(label = Label), position = position_stack(vjust = .5))

