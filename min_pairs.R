# Filename: min_pairs.R
# Author: Michael Dow, 2023
# Email: michael.dow@umontreal.ca
# Description:
#   Extracts all minimal pairs in a corpus, as defined as a distance of one operation (substitution, deletion or insertion).
# Inputs:
#   x (data frame): corpus with columns `phon` (phonetic transcriptions) and `ortho` (orthographic forms, optional).
# Notes:
#   - I suggest removing any unnecessary characters from transcriptions (including //, [], etc.), especially if they're Unicode instead of ASCII.
#   - Distinctions represented by combining characters or additional modifiers (e.g., aspiration) can be handled by this script but are treated as insertions. It may be useful to remove the modifying or combining character and transform the main segment in question.
#   - In order to reduce runtime, I suggest removing duplicate transcriptions and/or filtering by relevant context(s) at the beginning.
#   - For most corpora it'll be prohibitive to get all pairs at once. This script works in chunks, whose size is defined below as `n`. Between 1000 and 10000 should work, depending on your machine. Keep in mind that before filtering, the number of rows in `y` will be (n+1)*nrow(x).

# Load required libraries
library(tidyverse)
library(stringdist)
library(gtools)

# Define chunk size and generate steps
n = 1000

vec = c(seq(1, nrow(x), by=n), nrow(x))

all = data.frame()

# Loop over chunks, keeping only pairs with a Levenshtein distance of 1
for (i in 1:(length(vec)-1)) {
  temp = x[vec[i]:vec[i+1],]$phon
  
  y = expand.grid(temp, x$phon, stringsAsFactors = FALSE) %>% 
    filter(
      (Var1 != Var2) & 
        stringdist(Var1, Var2, method = "lv") == 1)
  
  all = bind_rows(all, y)
  
  rm(y, temp)
  gc()
  
  message(i, " of ", length(vec)-1, " done")
}

# Get differences in terms of operations
all = all %>% 
  rowwise() %>% 
  mutate(diffc = diag(
    attr(
      adist(Var1, Var2, counts = TRUE),
      "trafos")
  ),
  longer = ifelse(
    (nchar(Var1) >= nchar(Var2)), 
    Var1, 
    Var2))

diffc = all$diffc

# Get differences between pairs. `diff` and `diff2` are the differing segments between the two, while `diff.pre` and `diff.post` are the preceding and following contexts (blank if word boundary)
all = transform(
  all,
  diff = ifelse(!grepl("[ID]", diffc),
                regmatches(all$Var1,
                           regexpr("[^M]", diffc)),
                ""),
  diff2 = ifelse(!grepl("[ID]", diffc),
                 regmatches(all$Var2,
                            regexpr("[^M]", diffc)),
                 regmatches(all$longer,
                            regexpr("[^M]", diffc))
  ),
  diff.pre = regmatches(all$longer,
                        regexpr("[^M]", diffc) - 1),
  diff.post = regmatches(all$longer,
                         regexpr("[^M]", diffc) + 1)
)

# Run below if you want to remove pairs which differ in zero vs. something. You can also filter at this point by pre or post if desired.
all = all %>% filter(nchar(diff) != 0 & 
                       nchar(diff2) != 0)

# If interested in only pairs for a BINARY distinction, define in the following fashion, as illustrated below for oral vs. nasal vowel pairs. Each vector is imagined as a natural class, and the difference between the two vectors at each index should differ only in the desired feature. For instance, "a" is the oral equivalent of "@", and so on. 
oral = c("a", "A", "e", "E", "o", "O", "9", "2")
nasal = c("@", "@", "5", "5", "§", "§", "1", "1")
pairs = c(paste0(oral, nasal), paste0(nasal, oral))

# If interested in pairs with MORE THAN TWO values, use the following code below to define your pairs instead. Vowel height used as an example below. If any elements belong to more than one group, you may need to modify `permutations` such that `n = length(unique(c(high, mid, low)))`.
high = c("i", "u")
mid = c("e", "o")
low = c("a")

pairs = apply(
  permutations(
    n = length(c(high, mid, low)), 
    r = 2, v = c(high, mid, low)), 
  1, paste0, collapse = "")

# Get only those pairs which differ along the above lines
minimals = all %>% 
  unite(vv, diff, diff2, sep = "") %>% 
  filter(vv %in% pairs) %>% 
  distinct(Var1, Var2, .keep_all = T)

# Depending on your parameters, you may notice that many duplicate pairs exist, differing only in the order of the words between the two lines. Run the following block to reduce these pairs to a single instance.
m = vector()
for (i in 1:nrow(minimals)){
  if (length(unique(unlist(minimals[i:(i+1),c("Var1","Var2")]))) == 2) {
    m[i] = 0
  } else {
    m[i] = 1
  }
}

minimals$match = m

minimals = minimals %>% filter(match == 1)

# Assuming your original database in `x` has a column of orthographic forms called `ortho`, you can get these forms for each word.
minimals$ortho1 = x[match(minimals$Var1, x$phon),]$ortho
minimals$ortho2 = x[match(minimals$Var2, x$phon),]$ortho

# That's all!
