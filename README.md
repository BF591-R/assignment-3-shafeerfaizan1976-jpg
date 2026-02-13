# Assignment Bioinformatics Basics

## Problem Statement
Arranging the structure of our data inputs is vital to using R, as is creating
the correct environment for analysis by installing packages. It is also
important to share data and results easily using R markdown notebooks. We will
talk in more detail about microarrays in-class. 

## Required Readings
Please read sections 9.2 through 9.7.4 in the textbook. Section 9.2 starts [here](https://bu-bioinfo.github.io/r-for-biological-sciences/biology-bioinformatics.html#biological-data-overview)

## Learning Objectives
- Install various packages needed for analysis
- Load data, filter that data, and retrieve HGNC ids for the data
- Create a small plot to display results
- Utilize R Markdown to create an attractive format for sharing data.

## Skill List
- How to utilize R markdown to create a report of data analysis
- Installing and loading packages in R
- Utilizing Bioconductor to equate affy ids to HGNC ids (gene names)

## Instructions
Our main focus for this assignment are installing packages, manipulating data,
and plotting our manipulated data.

Please accept the GitHub classroom link provided on blackboard for this assignment.

The project is laid out as such:  
```
main.R
test_main.R
report.Rmd
```
A skeleton of the functions you need to complete is in `main.R`. Tests have been
pre-written to test your code and help you ensure it is running correctly, these
are in `test_main.R`. Finally, we are also introducing the concept of R
Markdown, which for this assignment is `report.Rmd`. The document itself goes
into greater detail, but you will:  

1. Complete the functions in `main.R` and use
`testthat::test_file('test_main.R')` to ensure they work correctly.  
2. Read the R Markdown file and complete the section called **"Assignment"**. To
do this, you can `source('main.R')` to bring over the functions you wrote in
step one.  
3. Finally, annotate the functions you wrote and **Knit** the R Markdown report,
complete with your additional comments and code execution.  

This page will go into detail on how the functions and their associated tests
should work.

## Function Details

### 1. Bioconductor
While many useful R packages can be loaded through
[CRAN](https://cran.r-project.org/) using the `install.packages()` syntax, a lot
of specifically bioinformatics packages are exclusively released on
[Bioconductor](https://bioconductor.org/install/). For this assignment we only
need the package called
[bioma**R**t](https://bioconductor.org/packages/release/bioc/html/biomaRt.html).
R programmers fancy themselves very clever, so **R**s show up a lot.

Naturally, we want to load our packages at the beginning of our script so all of
the code we write beneath can access it as it runs. However, if a user or
ourselves _already_ has this package installed we don't want to waste their time
installing it again. The function `require()` can help us avoid unnecessary
installation time and will help us develop faster. The Bioconductor link above
has an example of this method.

This section is untested.

### 2. `load_expression()`

Perhaps the most integral part of using the many data wrangling abilities of R
is actually entering your data into the R environment. While there are many ways
to do this in R, we ultimately want this data to be in a **tibble**, which means
the current form of the CSV will make R very angry if you attempt to load it in.
This is because tibbles don't support row names very well, and the first column
of our data doesn't have a name. Try to use the `load_expression()` data to load
data from a `filename` parameter and **return** a tibble of that information. I
called my first column "probeids".

#### Tests

The tests for this function are :  

```
describe("load_expression()", {
  result_tib <- load_expression('data/example_intensity_data_subset.csv')
  
  it("returns a tibble, not a dataframe", {
    expect_true(is_tibble(result_tib))
  })
  it("returns a tibble with the right dimensions", {
    expect_equal(dim(result_tib), c(1000, 36))
  })
})
```

This test uses the `load_expression()` you write to store the returned tibble in `result_tib`.  


While this test is using the same file for data input as you are, this may not always be the case.


The test then compares the dimensions of that result, and expects 1000 rows
and 36 columns. These are the dimensions of the input CSV. It also checks to
confirm it is a tibble object (because tibbles are better than dataframes).

### 3. `filter_15()`

In order to filter the numerous rows we have for this data, we introduce a
function that filters the probe IDs in the tibble our data is stored in. We want
to capture probes that have a suitably high level of expression, so we are
setting `log2(15)` as the cutoff for an expression level. We will keep a row if
15% of the values in that row exceed `log2(15)` (about 3.9). Since we may want
to examine the probe IDs we find, the function simple returns the values of the
probe IDs (column 1) instead of returning the entire tibble.

This function presents an important concept in R: using built-ins to speed up
our code. Built-ins are functions and packages that are optimized to process
data in a certain way. Since we're looking at each row of a table, we could
simply use a for loop to iterate one row at a time. This is slow, though, and
for this function might take 5-10 seconds to run (a long time for a program like
this!). Instead, you could use a function like `apply()` or `lapply()` to filter
every row at once. This solution takes mere moments.

#### Tests
```
describe("filter_15()", {
  test_tib <- tibble(probe=c('1_s_at', '2_s_at', '3_s_at', '4_s_at'),
                     GSM1=c(1.0, 3.95, 4.05, 0.5),
                     GSM2=rep(1.6, 4),
                     GSM3=rep(2.5, 4),
                     GSM4=rep(3.99, 4),
                     GSM5=rep(3.0, 4),
                     GSM6=rep(1.0, 4),
                     GSM7=rep(0.5, 4))
  
  it("correctly filters the rows meeting the specified criteria", {
    expect_equal(c(filter_15(test_tib)$probe), c("2_s_at", "3_s_at"))
  })
})
```

In order to test this function, we create a small sample tibble of expression
data containing only seven samples and four IDs. Two of the rows _do_ have more
than 15% of their values exceeding `log2(15)`, the other two do not. This test
ensures that `filter_15()` selects the correct rows. Creating a small sample
table like this can be very useful when testing your own code since you don't
need to look at a large amount of data to see if it's working correctly or not.
It also avoids dependency issues by avoiding requiring the prior function to 
work correctly. We have structured as many tests as possible in this fashion, 
which will allow you to work out of order if you'd like. Not _all_ tests were
able to be structured this way, so please read the instructions for each 
assignment for further guidance. 

### 4. `affy_to_hgnc()`
This is an important, but sometimes painful, part of using R. There is a great
built-in package for connecting to Ensembl (a database of genomic information
for many species) called `biomaRt`. We will use `biomaRt` to connect the
affymetrix probe IDs to more recognizable HGNC gene IDs. The problem is that
`biomaRt` depends on an external API (application program interface) to retrieve
data, and this connection sometimes (oftentimes) doesn't work. While their may
be more nuanced approaches to an unstable resource like this like automatically
retrying failed connections, the best advice for the time being is to try
running this function a few times if it doesn't work at first. The errors are
clear when it comes to a failed connection, so know that when you get to this
stage it likely isn't your code's fault.

To build a `biomaRt` query, read the documentation in section 3
[here](https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html#how-to-build-a-biomart-query).
The `biomart` you should use is `ENSEMBL_MART_ENSEMBL`, the data set
`hsapiens_gene_ensembl`, and you want to find the attributes
`c("affy_hg_u133_plus_2", "hgnc_symbol")`. The data you filter using
`filter_15()` returns a list of `affy_hg_u133_plus_2` probe IDs, and the gene
names we're interested in are stored in `hgnc_symbol`.

This function should return a tibble, but `biomaRt`'s `getBM()` will only accept
and return a `data.frame`. You can use `dplyr::pull()` to turn a tibble into a
simple character vector, and `dplyr::as_tibble()` to go from a data frame to a
tibble.

#### Tests
```
describe("affy_to_hgnc()", {
  success <- FALSE
  attempt <- 0
  while (!success ) {
    try({
      attempt <- attempt + 1
      tryCatch({
        response <- affy_to_hgnc(tibble('1553551_s_at'))},
        warning = function(w) stop("Could not connect to ENSEMBL. Reattempting connection")
      )
      Sys.sleep(1)
      success <- TRUE
    })
    if (attempt == 3) stop("Process halted. Function faced non-connection errors. Try again later.")
  }
  
  it("returns the right identifiers for the chosen probe, 1553551_s_at", {
    expect_equal(response$hgnc_symbol, c("MT-ND1", "MT-TI", "MT-TM", "MT-ND2"))
  })
})
```
As fun as it is to try to get `biomaRt` function to connect correctly, it is
even more fun to test them. Since a failure to connect doesn't indicate an
actual failure in _our_ code we must use a `try()` block in order to capture if
there is a connection error.


**Note** that the `try()` function is a part of programming called
***error-handling** which extends to many other languages. We often expect
*errors when running our programs (such as right now) but don't want to shut
*down our entire operation if it's an error we can expect. Using `try`,
*`except`, and `finally` (the latter two not appearing here) we can account for
*issues outside of our control and adapt our code to change the outcome.  

Using `try except` is **_not_** a replacement for writing code that doesn't
generate errors. If you can avoid an error in the first place, that is far
better than using error-handling.

In this case we `try()` to use `affy_to_hgnc()` to connect to Ensembl and store
the resulting error in `response`. We then check: if there is an "Error" then we
throw a `warning()` to our testing output. This doesn't stop further testing
from happening, but it ensures we know that something isn't quite right. If the
`response` does not contain an error, we simply test that it returned the
correct gene symbols for our random affy probe ID of choice.

### 5. `reduce_data()`

We have one final step in manipulating our data before we plot it. We have our
original data, the probe IDs and their associated HGNC symbols, and a list of
good gene names and bad gene names. `reduce_data` takes these four inputs and
returns a tibble that has reduced our expression data to only the genes of
interest and has a column describing which set of genes it belongs to ("good" or
"bad"). Changing the shape of the data is incredibly useful for `ggplot`, the
tidyverse package we will use for plotting. While there is flexibility when
using `ggplot` to plot data, having data in
[long](https://www.datacamp.com/community/tutorials/long-wide-data-R) format is
typically ideal.  

Once again, there are multiple ways to reorganize data in this way. We used the
base function `match()` to connect our probe IDs to with our HGNC IDs in
`name_ids`. We then used `tibble::add_column()` to insert the new data in the
correct location. Finally, we created two tibbles of good and bad genes using
`which()` and the `%in%` modifier. Which evaluates true conditions across a
range of data, so we can pass the list of genes we want and select the correct
ones. For instance:  

```
library(tibble)
tib <- tibble(gene = c("gene1", "gene2", "gene3", "gene4"),
              affy = c("a_s_1", "a_s_2", "a_s_3", "a_s_4"))
which(tib$gene %in% c("gene2", "gene3")) # returns the index of the TRUE rows
tib$affy[which(tib$gene %in% c("gene2", "gene3"))] # use [] to get data back
```

Once again, there are many ways to reshape this data (some maybe more elegant
than this!) and all we need is the data to be correctly shaped when it is
returned.

#### Tests
```
describe("reduce_data()", {
  t_tibble <- tibble(probe = c("1_s_at", "2_s_at", "3_s_at"),
                     GSM1 = c(9.5, 7.6, 5.5),
                     GSM2 = c(9.7, 7.2, 2.9),
                     GSM3 = c(6.9, 4.3, 6.8))
  names <- tibble(affy_hg_u133_plus_2 = c("1_s_at", "2_s_at", "3_s_at"),
                  hgnc_symbol = c("A-REAL-GENE", "SONIC", "UTWPU"))
  good <- c("A-REAL-GENE")
  bad <- c("SONIC")
  
  reduce_test <- reduce_data(t_tibble, names, good, bad)
  expected_result <- tibble(probe = c("1_s_at", "2_s_at"),
                   hgnc_symbol = c("A-REAL-GENE", "SONIC"),
                   gene_set = c("good", "bad"),
                   GSM1 = c(9.5, 7.6),
                   GSM2 = c(9.7, 7.2),
                   GSM3 = c(6.9, 4.3))
  
  it("returns a reduced tibble with the specified columns", {
    expect_equal(reduce_test, expected_result)
  })
})
```

In order to test a function that changes a tibble, we need to use a tibble. We
only create a test table that is three rows by four columns, but that is enough
to get the gist of the function. We simply pass the four parameters to
`reduce_data()` and we expect it to create a tibble like `result`. **Ensure your
output column names are the same here or your tests my fail.** While not crucial
to the success of your assignment, maintaining correct column names across
multiple data transformations is an important skill.

### 6. `convert_to_long`

In order to make plotting this data in the format desired, it is easier to manipulate
it to a "long" format as you learned in class. Tidyverse has a paired set of functions
`pivot_longer` and `pivot_wider` that alter the shape, but not the underlying content
of the tibble. Use the reduced tibble you generated in the last function, and create
a function that successfully converts it to a "long format" in order to facilitate 
creating a boxplot using ggplot2. 

#### Tests
```
describe("convert_to_long()", {
  # wide format tibble
  wide_tib <- tibble(probe = c("202274_at", "202541_at", "202542_s_at", "203919_at"),
                     hgnc = c("ACTG2", "AIMP1", "AIMP1", "TCEA2"),
                     gene_set = rep("good", 4),
                     GSM1 = c(8.05, 8.40, 9.55, 4.44),
                     GSM2 = c(7.74, 7.11, 8.48, 5.39))
  
  # expected long format tibble
  expected_tib <- tibble(probe = rep(c("202274_at", "202541_at", "202542_s_at", "203919_at"), each = 2),
                         hgnc = rep(c("ACTG2", "AIMP1", "AIMP1", "TCEA2"), each = 2),
                         gene_set = rep(rep("good", 4), each = 2),
                         sample = rep(c("GSM1", "GSM2"), 4),
                         value = c(8.05, 7.74, 8.40, 7.11, 9.55, 8.48, 4.44, 5.39))
  
  long_format <- convert_to_long(wide_tib)
  
  #this tests for shape
  it("correctly converts a wide format tibble to long format", {
    expect_equal(long_format, expected_tib)
  })
  #this tests for column named properly
  it("contains a column named sample", {
    expect_true("sample" %in% colnames(long_format))
  })
})
```

### 7. Plotting a boxplot

This is the first time we have asked you to fill in the actual Rmarkdown.
Please attempt to recreate the boxplot of samples using the data / variables
you have generated in the course of this assignment. Do your best to capture
the same general formatting, but don't worry too much about matching it exactly.
