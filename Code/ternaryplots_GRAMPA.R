# Illustrate GRAMPA dataset in ternary plots

## Working directory
getwd()
setwd('./Desktop/MODELS/')

## packages
install.packages("ggplot2")
install.packages('ggtern')
library(ggplot2)
library(ggtern)

## Load datasets
grampa_db <- read.csv('./Data/GRAMPA/grampa_pep2d.csv', header=TRUE)
grampa_db

ecoli_db <- read.csv('./Data/GRAMPA/Strain_datasets/ecoli_structure.csv', header=TRUE)
ecoli_db <- ecoli_db[!duplicated(ecoli_db[c('id')]), ]
ecoli_db

saureus_db <- read.csv('./Data/GRAMPA/Strain_datasets/saureus_structure.csv', header=TRUE)
saureus_db <- saureus_db[!duplicated(saureus_db[c('id')]), ]
saureus_db

paeru_db <- read.csv('./Data/GRAMPA/Strain_datasets/paeruginosa_structure.csv', header=TRUE)
paeru_db <- paeru_db[!duplicated(paeru_db[c('id')]), ]
paeru_db

pep2D_db <- read.csv("./Data/GRAMPA/pep2d_one.csv", header=TRUE)
pep2D_db

pdbe_db <- read.csv("./Data/GRAMPA/oneid_pdbe.csv", header=TRUE)
pdbe_db

#Subsets
grampa_sub <- grampa_db[,4:6]
colnames(grampa_sub) <- c('helix_H', 'sheet_E', 'coil_C')

ecoli_sub <- ecoli_db[,18:20]
colnames(ecoli_sub) <- c('helix_H', 'sheet_E', 'coil_C')
ecoli_sub

saureus_sub <- saureus_db[,17:19]
colnames(saureus_sub) <- c('helix_H', 'sheet_E', 'coil_C')
saureus_sub 

paeru_sub <- paeru_db[,17:19]
colnames(paeru_sub) <- c('helix_H', 'sheet_E', 'coil_C')
paeru_sub 

pep2D_sub <- pep2D_db[,4:6]
colnames(pep2D_sub) <- c('helix_H', 'sheet_E', 'coil_C')
pep2D_sub

pdbe_sub <- pdbe_db[,7:9]
colnames(pdbe_sub) <- c('helix_H', 'sheet_E', 'coil_C')
pdbe_sub

# Change null values from 0.0 to 0.5
grampa_sub[grampa_sub == 0.0] <- 0.5
grampa_sub

ecoli_sub[ecoli_sub == 0.0] <- 0.5
ecoli_sub

saureus_sub[saureus_sub == 0.0] <- 0.5
saureus_sub

paeru_sub[paeru_sub == 0.0] <- 0.5
paeru_sub

pep2D_sub[pep2D_sub == 0.0] <- 0.5
pep2D_sub

pdbe_sub[pdbe_sub == 0.0] <- 0.5
pdbe_sub

TernBasic <- ggtern::ggtern(data = pep2D_sub,
                            aes(x = sheet_E,
                                y = helix_H, 
                                z = coil_C)) +  
  geom_point(size=0.6, colour='dodgerblue3', alpha=0.8) + 
  #scale_color_viridis_c() +
  #scale_fill_gradient(low='green', high='red') +
  theme_bw() +
  theme_showarrows() +
  theme_anticlockwise() +
  ggtitle('Distribution PEP2D N=261')

TernBasic
ggsave('./Figures//TernaryPlots/Pep2D_points.png')

#Overlap geom points
Overlay <- ggtern(NULL,aes(x = sheet_E,y = helix_H,z = coil_C)) +  
  geom_point(data = pdbe_sub, size=1.4, colour='dodgerblue3', alpha=0.7) +
  geom_point(data = pep2D_sub, size=1.4, colour='darkslategray3', alpha=0.7) +
  theme_bw() +
  theme_showarrows() +
  theme_anticlockwise() +
  ggtitle('Distributions PDBe and PEP2D (N=261)')

Overlay
ggsave('./Figures//TernaryPlots/PDBe&Pep2D_points.png')

#Create a common dataframe to PDBe and Pep2D dataframes
##Make sure both dataframes share a column (pdb_id)
pep2D_sub$pdb_id <- pep2D_db$peptide_ID
pdbe_sub$pdb_id <- pdbe_db$pdb_id

## Merge their subsets into a new dataframe, rename columns
struc_db <- merge(x=pdbe_sub,y=pep2D_sub, by = 'pdb_id')
colnames(struc_db) <- c('pdb_id', 'pdbe_helix_H', 'pdbe_sheet_E', 'pdbe_coil_C', 'pep2d_helix_H', 'pep2d_sheet_E', 'pep2d_coil_C' )
write.csv(struc_db, file='./Data/GRAMPA/structures_pdbe_pep2d_HEC_values.csv')


TernDens <- ggtern::ggtern(data = paeru_sub,
                           aes(x = sheet_E,
                               y = helix_H, 
                               z = coil_C),
                           aes(x,y,z))  + 
  #geom_point(size=0.4, colour='darkgrey', alpha=0.8) +
  stat_density_tern(geom='polygon', 
                    #color='grey',
                    #n=300,
                    bins=50,
                    expand = 1, h=0.1,
                    base='identity',
                    aes(fill   = ..level.., alpha = 0.5),
                    na.rm = TRUE) +
  #geom_point(size=0.4, colour='darkgrey', alpha=0.8) +
  #scale_fill_distiller(palette = 'RdYlBu') +
  scale_fill_continuous(low = 'darkolivegreen1', high='red', name = "density") +
  #scale_fill_viridis(option= "B", direction = -1) +
  theme_bw() +
  theme_showarrows() +
  theme_anticlockwise() +
  ggtitle('Density Map P.aeruginosa N=2499')

TernDens
ggsave('./Figures/TernaryPlots/P_aeruginosa_density.png')

TernDens2 <- ggtern::ggtern(data = ecoli_sub,
                           aes(x = sheet_E,
                               y = helix_H, 
                               z = coil_C),
                           aes(x,y,z))  + 
  #geom_point(size=0.4, colour='darkgrey', alpha=0.8) +
  stat_density_tern(geom='polygon', 
                    #color='black',
                    #n=300, 
                    bins=50,
                    expand = 1, h=0.15,
                    base='identity',
                    aes(fill   = ..level.., alpha = 0.5),
                    na.rm = TRUE) +
  #geom_point(size=0.8, colour='darkgrey', alpha=0.8) +
  #scale_fill_distiller(palette = 'RdYlBu') +
  scale_fill_continuous(low = 'darkolivegreen1', high='red', name = "density", limits=c(0,25)) +
  #scale_fill_viridis(option= "B", direction = -1) +
  theme_bw() +
  theme_showarrows() +
  theme_anticlockwise() +
  ggtitle('Density Map E.coli N=4540')

TernDens2
ggsave('./Figures/TernaryPlots/E_coli_density.png')



#Display Mahalanobis distance
pep2D_sub <- pep2D_db[,4:6]
colnames(pep2D_sub) <- c('helix_H', 'sheet_E', 'coil_C')

#add reference points 
A <- c(90.0, 5.0, 5.0)
B <- c(5.0, 90.0, 5.0)
C <- c(5.0, 5.0, 90.0)
D <- c(30.0, 30.0, 40.0)

pep2D_sub <- rbind(pep2D_sub, A)
pep2D_sub <- rbind(pep2D_sub, B)
pep2D_sub <- rbind(pep2D_sub, C)
pep2D_sub <- rbind(pep2D_sub, D)
pep2D_sub

mahalanobis(pep2D_sub, colMeans(pep2D_sub), cov(pep2D_sub))
pep2D_sub$MD <- mahalanobis(pep2D_sub, colMeans(pep2D_sub), cov(pep2D_sub))

MD_density <- ggplot(pep2D_sub, aes(x=MD)) +
             geom_freqpoly() +
             geom_area(stat = "bin", bins = 30, color = "black", fill = 'azure4')
MD_density

#MD is not separating coiled and helical structures between 0 and 5, beta-sheets MD~40.


