---
  title: "HW Class 6"
author: "Victor Yu"
format: pdf
editor: 
  markdown: 
  wrap: 72
---
  
  ```{r}
df <- data.frame(a=1:10, b=seq(200,400,length=10),c=11:20,d=NA)
df$a <- (df$a - min(df$a)) / (max(df$a) - min(df$a))
df$b <- (df$b - min(df$a)) / (max(df$b) - min(df$b))
df$c <- (df$c - min(df$c)) / (max(df$c) - min(df$c))
df$d <- (df$d - min(df$d)) / (max(df$a) - min(df$d))
```

```{r}
install.packages("bio3d")
```

```{r}
library(bio3d)
s1 <- read.pdb("4AKE") # kinase with drug
s2 <- read.pdb("1AKE") # kinase no drug
s3 <- read.pdb("1E4Y") # kinase with drug
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s1, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b
plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")
```

Q1. What type of object is returned from the read.pdb() function?
  
  ```{r}
library("bio3d")
s1 <- read.pdb("4AKE") # kinase with drug
s2 <- read.pdb("1AKE") # kinase no drug
s3 <- read.pdb("1E4Y") # kinase with drug
```

Answer: The read.pdb() function returns the pdb file of interest as a
class.

Q2. What does the trim.pdb() function do?
  
  ```{r}
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s1, chain="A", elety="CA")
```

```{r}
head(s1)
```

Answer: The trim.pdb() pretty much returns a smaller object of the pdb
file we are looking at.

Q3. What input parameter would turn off the marginal black and grey
rectangles in the plots and what do they represent in this case?
  
  ```{r}

s1 <- read.pdb("4AKE") # kinase with drug
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
plotb3(s1.b, typ="l", ylab="Bfactor")
?plotb3()

```

Answer: The input parameter that turns off the black & grey rectangles
are sse=s#.chainA. The black and grey rectangles in the plots represent
secondary structure objects.

Q4. What would be a better plot to compare across the different
proteins?
  
  Answer: maybe a scatter plot? See the spread and varaition of each protein?
  
  Q5. Which proteins are more similar to each other in their B-factor
trends. How could you quantify this? HINT: try the rbind(), dist() and
hclust() functions together with a resulting dendrogram plot. Look up
the documentation to see what each of these functions does. \## Quarto

```{r}
hc <-  hclust( dist(rbind(s1.b, s2.b, s3.b) ))
plot(hc)
cor(s1.b, s3.b)
cor(s1.b, s2.b)
cor(s2.b, s3.b)
```

Answer: According to the dendrogram plot, it seems like s1.b and s3.b
are the most similar to each oter in their B-factor trends. We could
quantify there correlation / similarity with each using the cor ()
function.


Q6  How would you generalize the original code above to work with any set of input
protein structures?
  
  ```{r}
s1 <- read.pdb("4AKE") # kinase with drug

s1.chainA <- trim.pdb(s1, chain="A", elety="CA")

s1.b <- s1.chainA$atom$b

plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
```

```{r}
s2 <- read.pdb("1AKE") # kinase no drug

s2.chainA <- trim.pdb(s2, chain="A", elety="CA")

s2.b <- s2.chainA$atom$b

plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
```

```{r}
s3 <- read.pdb("1E4Y") # kinase with drug

s3.chainA <- trim.pdb(s1, chain="A", elety="CA")

s3.b <- s3.chainA$atom$b

plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")
```

Notes: The only things that changes in all 3 executions are
x <- read.pdb(______)
x_chainA <- trim.pdb (x, chain = "A", elety = "CA")
x.b <- x.chainA$atom$b
plotb3 (x.b, sse= x.chainA, typ = "l", ylab = "BFactor")

Pretty much variables names & file name. Variable names don't affect output (can literally name it anything I want. So, probably file name is important then.)

```{r}
#My function! To plot the pdb files wiithout having to retype everything. Only thing you need to input into the function is the file name since that was the only differentiating aspect of all 3 chunks of example code.
filereader <- function (file) {

x <- read.pdb(file)

x.chainA <- trim.pdb (x, chain = "A", elety = "CA")

x.b <- x.chainA$atom$b

plotb3(x.b, sse = x.chainA, typ="l", ylab = "Bfactor")
}
#Output of function is a Bfactor vs Residue graph of the proteins?

#Test Example
filereader("1E4Y") #s3

#WHAAAA IT WORKED ... ain't no shot
```

```{r}
#All 3 files / s# running through my function. Results produced the same graph are the previous above.
#Just input file name into function
filereader("4AKE") #s1
filereader("1AKE") #s2
filereader("1E4Y") #s3

#3 outputs of the code given to us.
```
Just using the original code as a visual comparison for myself here

```{r}
s1 <- read.pdb("4AKE") # kinase with drug
s2 <- read.pdb("1AKE") # kinase no drug
s3 <- read.pdb("1E4Y") # kinase with drug
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s1, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b
plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")
```
