seq_counts <- read.table("//wsl$/Ubuntu/home/rdecourc/bioinformatics/group-projects-marine_microbiomes/seqcounts_mod.txt", header=TRUE)

head(seq_counts[,1])

sea_vec <- c(1, 2, 3, 4, 10, 11, 12, 13, 14)
pyr_vec <- c(1, 5, 6, 7, 8, 9, 15, 16)

sea_counts <- (seq_counts[,sea_vec])
pyr_counts <- (seq_counts[,pyr_vec])

# I have these separated at this point
# I want to loop through and find median/mean (summary()) of each column, to pick cutoff

# Should also look into producing an abundance table of the counts


## Gathering means and medians of counts for individualized environments
## After removing counts of 0

dim(pyr_counts)

pyr_medians <- c()
pyr_means <- c()

for (x in c(2:8)){
  work_vec <- pyr_counts[,x]
  work_vec <- work_vec[work_vec>=1]
  pyr_medians <- c(pyr_medians, summary(work_vec)[3])
  pyr_means <- c(pyr_means, summary(work_vec)[4])
}

pyr_medians
pyr_means


# Now for the seawater

dim(sea_counts)

sea_medians <- c()
sea_means <- c()

for (x in c(2:9)){
  work_vec <- sea_counts[,x]
  work_vec <- work_vec[work_vec>=1]
  sea_medians <- c(sea_medians, summary(work_vec)[3])
  sea_means <- c(sea_means, summary(work_vec)[4])
}

sea_medians
sea_means

# pyrosome median 8-38
# seawater median 48 - 75

dim(pyr_counts)

atest <- (pyr_counts[2,2:8])
new <- as.numeric(as.vector(atest))
new

# Holding on to only cutoffs above 99, when taking the mean across all samples
pyr_asvlist <- c()
cutoff <- 99

for (row in c(1:5799)){
  current <- as.numeric(as.vector(pyr_counts[row,2:8]))
  if (mean(current) > cutoff){        # <<--- This line is where the avg. cutoff is selected
    pyr_asvlist <- c(pyr_asvlist, (as.vector(pyr_counts[row,1])))
  }
}

length(pyr_asvlist)

# Doing same for seawater

sea_asvlist <- c()

for (row in c(1:5799)){
  current <- as.numeric(as.vector(sea_counts[row,2:9]))
  if (mean(current) > cutoff){
    sea_asvlist <- c(sea_asvlist, (as.vector(sea_counts[row,1])))
  }
}

sea_asvlist
length(sea_asvlist)


# Hmmm - maybe make a histogram to see where a good spot to make the cutoff would be?

avec <- c()

dim(pyr_counts)

for (x in c(2:8)){
  avec <- c(avec, pyr_counts[,x])
}

## ^ Here (above) grabbed all the counts in the pyrosome only ASV table

# Removed all the 0 counts from the pyrosome ASV list
pyr_nozero <- (avec[avec>=1])


hist(pyr_nozero, breaks=5000, ylim=c(0,200), xlim=c(0,2000))
print(pyr_nozero)
# ^ Plotted as a histogram...

# Decided to do bins, to make it more easy to visualize
# Binned by several large chunks, 10, 20, 50, 100
bin_pnz <- (pyr_nozero%/%10)*10
hist(bin_pnz, breaks=10000, xlim=c(0,1000), ylim=c(0,200),
main='Occurrences - Pyrosome Sample ASV counts
Binned by blocks of 10
all 0 counts removed')

# Let's do the same for the seawater, just for consistency's sake

dim(sea_counts)
bvec <- c()
for (x in c(2:9)){
  bvec <- c(bvec, sea_counts[,x])
}
sea_nozero <- bvec[bvec>=1]
bin_snz <- (sea_nozero%/%10)*10
hist(bin_snz, breaks=10000, xlim=c(0,500), main="Occurrences - Sea Sample ASV counts
Binned by blocks of 10
all 0 counts removed")


### Testing - at least one of the samples must have it's ASV count above 100 to be
# included

