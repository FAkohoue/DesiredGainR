# Summarise a fitted selection index

Returns one tidy table per question a reader asks of an index: what the
coefficients are, what response is expected, and what the index is
worth. The tables are data frames rather than printed text so that they
can be written into a report without reformatting.

## Usage

``` r
# S3 method for class 'desiredgainr_index'
summary(object, ...)
```

## Arguments

- object:

  A fitted index from
  [`selection_index()`](https://FAkohoue.github.io/DesiredGainR/reference/selection_index.md)
  or
  [`restricted_index()`](https://FAkohoue.github.io/DesiredGainR/reference/restricted_index.md).

- ...:

  Unused.

## Value

An object of class `desiredgainr_index_summary`.

## Examples

``` r
set.seed(1)
traits <- c("yield", "protein")
G <- matrix(c(1.0, 0.2, 0.2, 0.5), 2, dimnames = list(traits, traits))
P <- matrix(c(2.5, 0.4, 0.4, 1.2), 2, dimnames = list(traits, traits))
values <- as.data.frame(matrix(
  stats::rnorm(40),
  ncol = 2, dimnames = list(paste0("g", 1:20), traits)
))
fit <- selection_index(
  values, traits,
  method = "smith_hazel", G = G, P = P,
  economic_weights = c(yield = 2, protein = 1), n_select = 5
)
summary(fit)
#> <desiredgainr_index_summary>
#>   Method: smith_hazel 
#>   Candidates: 20   Selected: 5   Intensity: 1.271
#>   Traits centred: TRUE  scaled: TRUE (sample) 
#> 
#>   Coefficients:
#>      Trait Coefficient Effective_weight Aggregate_weight
#>     <char>       <num>            <num>            <num>
#> 1:   yield   0.8034940               NA                2
#> 2: protein   0.4791643               NA                1
#> 
#>   Response:
#>      Trait Expected_response Observed_differential
#>     <char>             <num>                 <num>
#> 1:   yield         0.8407274             0.7112215
#> 2: protein         0.4014156             0.8755802
#> 
#>   Criteria:
#>         Criterion     Value
#>            <char>     <num>
#> 1:           R_HI 0.6447193
#> 2:        Delta_H 2.0828704
#> 3:             RE 0.9550691
#> 4:       h2_index 0.4166670
#> 5: accuracy_index 0.6454975
#> 6:           CV_I        NA
#>                                                     Meaning
#>                                                      <char>
#> 1:              Correlation between the index and net merit
#> 2:                         Expected gain in aggregate merit
#> 3:         Efficiency relative to direct selection on yield
#> 4:             Heritability of the index treated as a trait
#> 5: Correlation between the index and the value it estimates
#> 6:          Coefficient of variation of the index (percent)
```
