---
title: "R for Data Science"
output:
  html_notebook: default
  pdf_document: default
---


If we need to be explicit about where a function (or dataset) comes from, we’ll use the special form 
``` package::function()```. For example, ```ggplot2::ggplot()``` tells you explicitly that we’re using the ```ggplot()``` function from the ggplot2 package.

# EXPLORE

## INSTALLATION

```{r}
library(tidyverse)
library(nycflights13)
library(stringr)
library(forcats)
library(lubridate)
library(magrittr)
library(modelr)
```


## DATA VISUALIZATION

A data frame is a rectangular collection of variables (in the columns) and observations (in the rows).

```{r}
ggplot2::mpg
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy))
```

With ggplot2, you begin a plot with the function ```ggplot()```. ```ggplot()``` creates a coordinate system that you can add layers to. The first argument of ```ggplot()``` is the dataset to use in the graph.

You complete your graph by adding one or more layers to ```ggplot()```. 

The function ```geom_point()``` adds a layer of points to your plot, which creates a scatterplot. ggplot2 comes with many geom functions that each add a different type of layer to a plot.


Each geom function in ggplot2 takes a ```mapping``` argument. This defines how variables in your dataset are mapped to visual properties. 
The ```mapping``` argument is always paired with ```aes()```, and the ```x``` and ```y``` arguments of ```aes()``` specify which variables to map to the x and y axes. 
ggplot2 looks for the mapped variable in the ```data``` argument, in this case, ```mpg```.

### A Graphing Template

```
ggplot(data = <DATA>) + 
  <GEOM_FUNCTION>(mapping = aes(<MAPPINGS>))
```

### Aesthetic Mappings

An aesthetic is a visual property of the objects in your plot. Aesthetics include things like the size, the shape, or the color of your points. You can display a point in different ways by changing the values of its aesthetic properties.

The alpha aesthetic controls the transparency of the points, or the shape of the points.

ggplot2 will only use six shapes at a time. By default, additional groups will go unplotted when you use the shape aesthetic

```{r}
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy, color = class))

ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy, alpha = class))

ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy, shape = class))
```

R has 25 built in shapes that are identified by numbers.

One common problem when creating ggplot2 graphics is to put the ```+``` in the wrong place: it has to come at the end of the line, not the start

### Facets

One way to add additional variables is with aesthetics. Another way, particularly useful for categorical variables, is to split your plot into **facets**, subplots that each display one subset of the data.

To facet your plot by a ***single variable***, use ```facet_wrap()```. 
The first argument of ```facet_wrap()``` should be a formula, which you create with ```~``` followed by a variable name (here “formula” is the name of a data structure in R, not a synonym for “equation”). 
The variable that you pass to ```facet_wrap()``` should be discrete.

```{r}
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy)) + 
  facet_wrap(~ class, nrow = 2)
```

To facet your plot on the combination of two variables, add ```facet_grid()``` to your plot call. 
The first argument of ```facet_grid()``` is also a formula. This time the formula should contain two variable names separated by a ```~```.

```{r}
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy)) + 
  facet_grid(drv ~ cyl)
```

If you prefer to not facet in the rows or columns dimension, use a ```.``` instead of a variable name, e.g. ```+ facet_grid(. ~ cyl)```.


### Geometric Objects

A ***geom*** is the geometrical object that a plot uses to represent data.

```{r}
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy))
ggplot(data = mpg) + 
  geom_smooth(mapping = aes(x = displ, y = hwy))
```


Every geom function in ggplot2 takes a mapping argument. However, not every aesthetic works with every geom. 

ggplot2 provides over 30 geoms, and extension packages provide even more

To display multiple geoms in the same plot, add multiple geom functions to ggplot():

```{r}
ggplot(data = mpg, mapping = aes(x = displ, y = hwy)) + 
  geom_point() + 
  geom_smooth()
```


### Statistical Transformations

Consider a basic bar chart, as drawn with ```geom_bar()```:

```{r}
ggplot(data = diamonds) + 
  geom_bar(mapping = aes(x = cut))
```

On the x-axis, the chart displays ```cut```, a variable from ```diamonds```. On the y-axis, it displays ```count```, but count is not a variable in diamonds! Where does count come from? 

Many graphs, like scatterplots, plot the raw values of your dataset. Other graphs, like bar charts, calculate new values to plot:

* bar charts, histograms, and frequency polygons bin your data and then plot bin counts, the number of points that fall in each bin.

* smoothers fit a model to your data and then plot predictions from the model.

* boxplots compute a robust summary of the distribution and then display a specially formatted box.

The algorithm used to calculate new values for a graph is called a **stat**, short for statistical transformation. 

You can learn which stat a geom uses by inspecting the default value for the ```stat``` argument. For example, ```?geom_bar``` shows that the default value for ```stat``` is “count”, which means that ```geom_bar()``` uses ```stat_count()```. 

### Position Adjustments

You can colour a bar chart using either the ```colour``` aesthetic, or, more usefully,```fill```:

```{r}
ggplot(data = diamonds) + 
  geom_bar(mapping = aes(x = cut, colour = cut))
ggplot(data = diamonds) + 
  geom_bar(mapping = aes(x = cut, fill = cut))
```

If you map the fill aesthetic to another variable, like ```clarity``` the bars are automatically stacked. 
Each colored rectangle represents a combination of ```cut``` and ```clarity```.

The stacking is performed automatically by the ***position adjustment*** specified by the position argument. 
If you don’t want a stacked bar chart, you can use one of three other options: ```"identity"```, ```"dodge"``` or ```"fill"```.

* ```position = "identity"``` will place each object exactly where it falls in the context of the graph. This is not very useful for bars, because it overlaps them. To see that overlapping we either need to make the bars slightly transparent by setting ```alpha``` to a small value, or completely transparent by setting ```fill = NA```.

```{r}
ggplot(data = diamonds, mapping = aes(x = cut, fill = clarity)) + 
  geom_bar(alpha = 1/5, position = "identity")
ggplot(data = diamonds, mapping = aes(x = cut, colour = clarity)) + 
  geom_bar(fill = NA, position = "identity")
```

* ```position = "fill"``` works like stacking, but makes each set of stacked bars the same height. This makes it easier to compare proportions across groups.

```{r}
ggplot(data = diamonds) + 
  geom_bar(mapping = aes(x = cut, fill = clarity), position = "fill")
```


* ```position = "dodge"``` places overlapping objects directly beside one another. This makes it easier to compare individual values.

```{r}
ggplot(data = diamonds) + 
  geom_bar(mapping = aes(x = cut, fill = clarity), position = "dodge")
```

There’s one other type of adjustment that’s not useful for bar charts, but it can be very useful for scatterplots. The values of ```hwy``` and ```displ``` are rounded so the points appear on a grid and many points overlap each other. This problem is known as **overplotting**.

You can avoid this gridding by setting the position adjustment to **“jitter”**. ```position = "jitter"``` adds a small amount of random noise to each point. This spreads the points out because no two points are likely to receive the same amount of random noise:

```{r}
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy), position = "jitter")
```


Because this is such a useful operation, ggplot2 comes with a shorthand for ```geom_point(position = "jitter")```: ```geom_jitter()```.


### Coordinate Systems

The default coordinate system is the Cartesian coordinate system where the x and y positions act independently to determine the location of each point. There are a number of other coordinate systems that are occasionally helpful:

* ```coord_flip()``` switches the x and y axes. 
This is useful (for example), if you want horizontal boxplots. It’s also useful for long labels: it’s hard to get them to fit without overlapping on the x-axis.

```{r}
ggplot(data = mpg, mapping = aes(x = class, y = hwy)) + 
  geom_boxplot()
ggplot(data = mpg, mapping = aes(x = class, y = hwy)) + 
  geom_boxplot() +
  coord_flip()
```

* ```coord_quickmap()``` sets the aspect ratio correctly for maps. This is very important if you’re plotting spatial data with ggplot2

* ```coord_polar()``` uses polar coordinates. Polar coordinates reveal an interesting connection between a bar chart and a Coxcomb chart.

### The Layered Grammar of Graphics

```
ggplot(data = <DATA>) + 
  <GEOM_FUNCTION>(
     mapping = aes(<MAPPINGS>),
     stat = <STAT>, 
     position = <POSITION>
  ) +
  <COORDINATE_FUNCTION> +
  <FACET_FUNCTION>
```

The seven parameters in the template compose the grammar of graphics, a formal system for building plots. The grammar of graphics is based on the insight that you can uniquely describe any plot as a combination of a dataset, a geom, a set of mappings, a stat, a position adjustment, a coordinate system, and a faceting scheme.


## DATA TRANSFORMATION

```{r}
library(nycflights13)
```

```{r}
flights
```


 The dataset ```flights``` prints differently because it’s a **tibble**. Tibbles are data frames, but slightly tweaked to work better in the tidyverse.
 
These describe the type of each variable:

* *int* stands for integers.

* *dbl* stands for doubles, or real numbers.

* *chr* stands for character vectors, or strings.

* *dttm* stands for date-times (a date + a time).

* *lgl* stands for logical, vectors that contain only TRUE or FALSE.

* *fctr* stands for factors, which R uses to represent categorical variables with fixed possible values.

* *date* stands for dates.

There are five key **dply**r functions that allow you to solve the vast majority of your data manipulation challenges:

* Pick observations by their values: ```filter()```
* Reorder the rows ```arrange()```
* Pick variables by their names ```select()```
* Create new variables with functions of existing variables ```mutate()```
* Collapse many values down to a single summary ```summarise()```

These can all be used in conjunction with ```group_by()``` which changes the scope of each function from operating on the entire dataset to operating on it group-by-group. 

These six functions provide the verbs for a language of data manipulation.

### Filter rows with ```filter```

```filter()``` allows you to subset observations based on their values. The first argument is the name of the data frame. The second and subsequent arguments are the expressions that filter the data frame.

For example, we can select all flights on January 1st with:

```{r}
filter(flights, month == 1, day ==1)
```

R either prints out the results, or saves them to a variable. **If you want to do both, you can wrap the assignment in parentheses**:

```{r}
(dec25 <- filter(flights, month == 12, day == 25))
```

#### Comparisons

R provides the standard suite of comparison operators: ```>```, ```>=```, ```<```, ```<=```, ```!=``` (not equal), and ```==``` (equal).

#### Logical operators

Multiple arguments to ```filter()``` are combined with “and”: every expression must be true in order for a row to be included in the output. 

For other types of combinations, you’ll need to use Boolean operators yourself: ```&``` is “and”, ```|``` is “or”, and ```!``` is “not”.

A useful short-hand for this problem is ```x %in% y```. **This will select every row where x is one of the values in y**.

#### Missing values

One important feature of R that can make comparison tricky are missing values, or ```NA```s (“not availables”).

If you want to determine if a value is missing, use ```is.na()```

```filter()``` only includes rows where the condition is ```TRUE```; it excludes both ```FALSE``` and ```NA``` values. If you want to preserve missing values, ask for them explicitly:

```
filter(df, is.na(x) | x > 1)
```

### Arrange rows with ```arrange()```

```arrange()``` works similarly to ```filter()``` except that instead of selecting rows, it changes their order. It takes a data frame and a set of column names (or more complicated expressions) to order by.

```{r}
arrange(flights, year, month, day)
```

Use ```desc()``` to re-order by a column in descending order:

```{r}
arrange(flights, desc(arr_delay))
```

Missing values are always sorted at the end


### Select columns with ```select()```

```select()``` allows you to rapidly zoom in on a useful subset using operations based on the names of the variables.

```{r}
# Select columns by name
select(flights, year, month, day)
# Select all columns between year and day (inclusive)
select(flights, year:day)
# Select all columns except those from year to day (inclusive)
select(flights, -(year:day))
```

There are a number of helper functions you can use within ```select()```:

* ```starts_with("abc")```: matches names that begin with “abc”.

* ```ends_with("xyz")```: matches names that end with “xyz”.

* ```contains("ijk")```: matches names that contain “ijk”.

* ```matches("(.)\\1")```: selects variables that match a regular expression. This one matches any variables that contain repeated characters. You’ll learn more about regular expressions in strings.

* ```num_range("x", 1:3)```: matches x1, x2 and x3.

See ```?select``` for more details.

To rename variables, use ```rename()```:

```{r}
rename(flights, tail_num = tailnum)
```


### Add new variables with ```mutate()```

```mutate()``` always adds new columns at the end of your dataset.

In RStudio, the easiest way to see all the columns is ```View()```.

```{r}
flights_sml <- select(flights, 
  year:day, 
  ends_with("delay"), 
  distance, 
  air_time
)
mutate(flights_sml,
  gain = arr_delay - dep_delay,
  speed = distance / air_time * 60
)
```

```{r}
mutate(flights_sml,
  gain = arr_delay - dep_delay,
  hours = air_time / 60,
  gain_per_hour = gain / hours
)
```

If you only want to keep the new variables, use ```transmute()```:

```{r}
transmute(flights,
  gain = arr_delay - dep_delay,
  hours = air_time / 60,
  gain_per_hour = gain / hours
)
```

There are many functions for creating new variables that you can use with ```mutate()```. The key property is that the function must be vectorised: it must take a vector of values as input, return a vector with the same number of values as output.

* Arithmetic operators: ```+```,```-```, ```*```, ```/```, ```^```.

* Modular arithmetic: ```%/%``` (integer division) and ```%%``` (remainder)

* Logs: ```log()```, ```log2()```, ```log10()```.

* Offsets: ```lead()``` and ```lag()```

* Cumulative and rolling aggregates: R provides functions for running sums, products, mins and maxes: ```cumsum()```, ```cumprod()```, ```cummin()```, ```cummax()```; and dplyr provides ```cummean()``` for cumulative means.

* Logical comparisons, ```<```, ```<=```,```>```, ```>=```, ```!=```

* Ranking: there are a number of ranking functions, but you should start with ```min_rank()```.

* If ```min_rank()``` doesn’t do what you need, look at the variants ```row_number()```, ```dense_rank()```, ```percent_rank()```, ```cume_dist()```, ```ntile()```.

### Grouped summaries with ```summarise()```

It collapses a data frame to a single row:

```summarise()``` is not terribly useful unless we pair it with ```group_by()```. This changes the unit of analysis from the complete dataset to individual groups. Then, when you use the dplyr verbs on a grouped data frame they’ll be automatically applied "by group".

```{r}
by_day <- group_by(flights, year, month, day)
summarise(by_day, delay = mean(dep_delay, na.rm = TRUE))
```

**Together ```group_by()``` and ```summarise()``` provide one of the tools that you’ll use most commonly when working with dplyr: grouped summaries**.

#### Combining operations with pipe ```%>%```

```{r}
delays <- flights %>% 
  group_by(dest) %>% 
  summarise(
    count = n(),
    dist = mean(distance, na.rm = TRUE),
    delay = mean(arr_delay, na.rm = TRUE)
  ) %>% 
  filter(count > 20, dest != "HNL")
```

You can read it as a series of imperative statements: group, then summarise, then filter. As suggested by this reading, a good way to pronounce ```%>%``` when reading code is **“then”**.

Behind the scenes, ```x %>% f(y)``` turns into ```f(x, y)```, and ```x %>% f(y) %>% g(z)``` turns into ```g(f(x, y), z)``` and so on. 

The only exception is ggplot2: it was written before the pipe was discovered.

#### Missing values

Aggregation functions obey the usual rule of missing values: if there’s any missing value in the input, the output will be a missing value. Fortunately, all aggregation functions have an ```na.rm``` argument which removes the missing values prior to computation:

```{r}
flights %>% 
  group_by(year, month, day) %>% 
  summarise(mean = mean(dep_delay, na.rm = TRUE))
```

Where missing values represent cancelled flights, we could also tackle the problem by first removing the cancelled flights.

```{r}
not_cancelled <- flights %>% 
  filter(!is.na(dep_delay), !is.na(arr_delay))

not_cancelled %>% 
  group_by(year, month, day) %>% 
  summarise(mean = mean(dep_delay))
```

#### Counts

Whenever you do any aggregation, it’s always a good idea to include either a count (```n()```), or a count of non-missing values (```sum(!is.na(x))```). That way you can check that you’re not drawing conclusions based on very small amounts of data.


```{r}
delays <- not_cancelled %>% 
  group_by(tailnum) %>% 
  summarise(
    delay = mean(arr_delay, na.rm = TRUE),
    n = n()
  )

ggplot(data = delays, mapping = aes(x = n, y = delay)) + 
  geom_point(alpha = 1/10)
```


The following code shows you a handy pattern for integrating ggplot2 into dplyr flows:

```{r}
delays %>% 
  filter(n > 25) %>% 
  ggplot(mapping = aes(x = n, y = delay)) + 
    geom_point(alpha = 1/10)
```

#### Summary functions

* Measures of location: we’ve used ```mean(x)```, but ```median(x)``` is also useful. The mean is the sum divided by the length; the median is a value where 50% of x is above it, and 50% is below it.


* Measures of spread: ```sd(x)```, ```IQR(x)```, ```mad(x)```. The mean squared deviation, or standard deviation or sd for short, is the standard measure of spread. The interquartile range ```IQR()``` and median absolute deviation ```mad(x)``` are robust equivalents that may be more useful if you have outliers.

* Measures of rank: ```min(x)```, ```quantile(x, 0.25)```, ```max(x)```. Quantiles are a generalisation of the median. For example, ```quantile(x, 0.25)``` will find a value of ```x``` that is greater than 25% of the values, and less than the remaining 75%.

* Measures of position:```first(x)```, ```nth(x, 2)```, ```last(x)```. These work similarly to ```x[1]```, ```x[2]```, and ```x[length(x)]``` but let you set a default value if that position does not exist (i.e. you’re trying to get the 3rd element from a group that only has two elements).

* Counts: You’ve seen ```n()```, which takes no arguments, and returns the size of the current group. To count the number of non-missing values, use ```sum(!is.na(x))```. To count the number of distinct (unique) values, use ```n_distinct(x)```.

* Counts and proportions of logical values: ```sum(x > 10)```, ```mean(y == 0)```. When used with numeric functions, ```TRUE``` is converted to 1 and ```FALSE``` to 0. This makes ```sum()``` and ```mean()``` very useful: ```sum(x)``` gives the number of ```TRUE```s in x, and ```mean(x)``` gives the proportion.


#### Ungrouping

If you need to remove grouping, and return to operations on ungrouped data, use ```ungroup()```.


### Grouped Mutates

Grouping is most useful in conjunction with ```summarise()```, but you can also do convenient operations with ```mutate()``` and ```filter()```:

* Find the worst members of each group:

```{r}
flights_sml %>% 
  group_by(year, month, day) %>%
  filter(rank(desc(arr_delay)) < 10)
```

* Find all groups bigger than a threshold:

```{r}
popular_dests <- flights %>% 
  group_by(dest) %>% 
  filter(n() > 365)
popular_dests
```

* Standardise to compute per group metrics:

```{r}
popular_dests %>% 
  filter(arr_delay > 0) %>% 
  mutate(prop_delay = arr_delay / sum(arr_delay)) %>% 
  select(year:day, dest, arr_delay, prop_delay)
```

A grouped filter is a grouped mutate followed by an ungrouped filter. I generally avoid them except for quick and dirty manipulations: otherwise it’s hard to check that you’ve done the manipulation correctly.

You should always start your script with the packages that you need. That way, if you share your code with others, they can easily see what packages they need to install. 

**You should never include ```install.packages()``` or ```setwd()``` in a script that you share. It’s very antisocial to change settings on someone else’s computer!**

You can also execute the complete script in one step: ```Cmd/Ctrl + Shift + S```. Doing this regularly is a great way to check that you’ve captured all the important parts of your code in the script.


## EXPLORATORY DATA ANALYSIS

To do data cleaning, you’ll need to deploy all the tools of EDA: visualisation, transformation, and modelling.

There is no rule about which questions you should ask to guide your research. However, two types of questions will always be useful for making discoveries within your data. You can loosely word these questions as:

1. What type of variation occurs within my variables?

2. What type of covariation occurs between my variables?

To make the discussion easier, let’s define some terms:

* A variable is a quantity, quality, or property that you can measure.

* A value is the state of a variable when you measure it. The value of a variable may change from measurement to measurement.

* An observation is a set of measurements made under similar conditions (you usually make all of the measurements in an observation at the same time and on the same object). An observation will contain several values, each associated with a different variable. I’ll sometimes refer to an observation as a data point.

* Tabular data is a set of values, each associated with a variable and an observation. Tabular data is tidy if each value is placed in its own “cell”, each variable in its own column, and each observation in its own row.


### Variation

Variation is the tendency of the values of a variable to change from measurement to measurement.

#### Visualising distributions

How you visualise the distribution of a variable will depend on whether the variable is **categorical** or **continuous**. 

A variable is **categorical** if it can only take one of a small set of values. In R, categorical variables are usually saved as factors or character vectors. 
To examine the distribution of a categorical variable, use a *bar chart*:

```{r}
ggplot(data = diamonds) +
  geom_bar(mapping = aes(x = cut))
```

The height of the bars displays how many observations occurred with each x value. You can compute these values manually with ```dplyr::count()```:

```{r}
diamonds %>% 
  count(cut)
```

A variable is **continuous** if it can take any of an infinite set of ordered values. Numbers and date-times are two examples of continuous variables. To examine the distribution of a continuous variable, use a **histogram**:

```{r}
ggplot(data = diamonds) +
  geom_histogram(mapping = aes(x = carat), binwidth = 0.5)
```

You can compute this by hand by combining ```dplyr::count()``` and ```ggplot2::cut_width()```:

```{r}
diamonds %>% 
  count(cut_width(carat, 0.5))
```

A histogram divides the x-axis into equally spaced bins and then uses the height of a bar to display the number of observations that fall in each bin. 

You can set the width of the intervals in a histogram with the ```binwidth``` argument, which is measured in the units of the ```x``` variable. You should always explore a variety of binwidths when working with histograms, as **different binwidths can reveal different patterns**.

If you wish to overlay multiple histograms in the same plot, use ```geom_freqpoly()```.
```geom_freqpoly()``` performs the same calculation as ```geom_histogram()```, but instead of displaying the counts with bars, uses lines instead. It’s much easier to understand overlapping lines than bars.

```{r}
smaller <- diamonds %>% 
  filter(carat < 3)
ggplot(data = smaller, mapping = aes(x = carat, colour = cut)) +
  geom_freqpoly(binwidth = 0.1)
```


#### Typical values

In both bar charts and histograms, tall bars show the common values of a variable, and shorter bars show less-common values. Places that do not have bars reveal values that were not seen in your data. 

To turn this information into useful questions, look for anything unexpected


#### Unusual values

Outliers are observations that are unusual; data points that don’t seem to fit the pattern. Sometimes outliers are data entry errors; other times outliers suggest important new science. 

To make it easy to see the unusual values, we need to zoom to small values of the y-axis with ```coord_cartesian()```:

```{r}
ggplot(diamonds) + 
  geom_histogram(mapping = aes(x = y), binwidth = 0.5) +
  coord_cartesian(ylim = c(0, 50))
```


```coord_cartesian()``` also has an ```xlim()``` argument for when you need to zoom into the x-axis. ggplot2 also has ```xlim()``` and ```ylim()``` functions that work slightly differently: they throw away the data outside the limits.

```{r}
unusual <- diamonds %>% 
  filter(y < 3 | y > 20) %>% 
  select(price, x, y, z) %>%
  arrange(y)
unusual
```


### Missing Values

If you’ve encountered unusual values in your dataset, and simply want to move on to the rest of your analysis, you have two options.

1. Drop the entire row with the strange values:

```{r}
diamonds2 <- diamonds %>% 
  filter(between(y, 3, 20))
```

I don’t recommend this option because just because one measurement is invalid, doesn’t mean all the measurements are.


2. I recommend replacing the unusual values with missing values. The easiest way to do this is to use ```mutate()``` to replace the variable with a modified copy. You can use the ```ifelse()``` function to replace unusual values with ```NA```:


```{r}
diamonds2 <- diamonds %>% 
  mutate(y = ifelse(y < 3 | y > 20, NA, y))
```


```ifelse()``` has **three arguments**. 
The first argument ```test``` should be a logical vector. The result will contain the value of the second argument, yes, when test is TRUE, and the value of the third argument, no, when it is false.

Like R, ggplot2 subscribes to the philosophy that missing values should never silently go missing. It’s not obvious where you should plot missing values, so ggplot2 doesn’t include them in the plot, but it does warn that they’ve been removed:

```{r}
ggplot(data = diamonds2, mapping = aes(x = x, y = y)) + 
  geom_point()
```

To suppress that warning, set ```na.rm = TRUE```:

```{r}
ggplot(data = diamonds2, mapping = aes(x = x, y = y)) + 
  geom_point(na.rm = TRUE)
```

Other times you want to understand what makes observations with missing values different to observations with recorded values. For example, in ```nycflights13::flights```, missing values in the dep_time variable indicate that the flight was cancelled. So you might want to compare the scheduled departure times for cancelled and non-cancelled times. You can do this by making a new variable with ```is.na()```.

```{r}
nycflights13::flights %>% 
  mutate(
    cancelled = is.na(dep_time),
    sched_hour = sched_dep_time %/% 100,
    sched_min = sched_dep_time %% 100,
    sched_dep_time = sched_hour + sched_min / 60
  ) %>% 
  ggplot(mapping = aes(sched_dep_time)) + 
    geom_freqpoly(mapping = aes(colour = cancelled), binwidth = 1/4)
```



### Covariation

If variation describes the behavior within a variable, **covariation** describes the behavior between variables. 

**Covariation** is the tendency for the values of two or more variables to vary together in a related way. 

The best way to spot covariation is to visualise the relationship between two or more variables.


#### A categorical and continuous variable

To make the comparison easier we need to swap what is displayed on the y-axis. Instead of displaying count, we’ll display ```density```, which is the count standardised so that the area under each frequency polygon is one.

```{r}
ggplot(data = diamonds, mapping = aes(x = price, y = ..density..)) + 
  geom_freqpoly(mapping = aes(colour = cut), binwidth = 500)
```


Another alternative to display the distribution of a continuous variable broken down by a categorical variable is the **boxplot** ; each boxplot consists of:

* A box that stretches from the 25th percentile of the distribution to the 75th percentile, a distance known as the interquartile range (IQR). In the middle of the box is a line that displays the median, i.e. 50th percentile, of the distribution. These three lines give you a sense of the spread of the distribution and whether or not the distribution is symmetric about the median or skewed to one side.

* Visual points that display observations that fall more than 1.5 times the IQR from either edge of the box. These outlying points are unusual so are plotted individually.

* A line (or whisker) that extends from each end of the box and goes to the
farthest non-outlier point in the distribution.

Many categorical variables don’t have such an intrinsic order, so you might want to reorder them to make a more informative display. One way to do that is with the ```reorder()``` function.

```{r}
ggplot(data = mpg) +
  geom_boxplot(mapping = aes(x = reorder(class, hwy, FUN = median), y = hwy))
```

We just reordered ```class``` based on the median value of ```hwy```

If you have long variable names, ```geom_boxplot()``` will work better if you flip it 90°. You can do that with ```coord_flip()```.

```{r}
ggplot(data = mpg) +
  geom_boxplot(mapping = aes(x = reorder(class, hwy, FUN = median), y = hwy)) +
  coord_flip()
```


#### Two categorical variables

To visualise the covariation between categorical variables, you’ll need to count the number of observations for each combination. One way to do that is to rely on the built-in ```geom_count()```:

```{r}
ggplot(data = diamonds) +
  geom_count(mapping = aes(x = cut, y = color))
diamonds %>% 
  count(color, cut)
```

Then visualise with geom_tile() and the fill aesthetic:

```{r}
diamonds %>% 
  count(color, cut) %>%  
  ggplot(mapping = aes(x = color, y = cut)) +
    geom_tile(mapping = aes(fill = n))
```


#### Two categorical variables

Scatterplots become less useful as the size of your dataset grows, because points begin to overplot, and pile up into areas of uniform black (as above). We can use the ```alpha``` aesthetic to add transparency.

```{r}
ggplot(data = diamonds) + 
  geom_point(mapping = aes(x = carat, y = price), alpha = 1 / 100)
```

```geom_bin2d()``` and ```geom_hex()``` divide the coordinate plane into 2d bins and then use a fill color to display how many points fall into each bin. 

```geom_bin2d()``` creates rectangular bins. 

```geom_hex()``` creates hexagonal bins. You will need to install the hexbin package to use ```geom_hex()```.

```{r}
ggplot(data = smaller) +
  geom_bin2d(mapping = aes(x = carat, y = price))

# install.packages("hexbin")
ggplot(data = smaller) +
  geom_hex(mapping = aes(x = carat, y = price))
```


### Patterns and Models

Patterns provide one of the most useful tools for data scientists because they reveal covariation. 

If you think of variation as a phenomenon that creates uncertainty, covariation is a phenomenon that reduces it. If two variables covary, you can use the values of one variable to make better predictions about the values of the second. 

If the covariation is due to a causal relationship (a special case), then you can use the value of one variable to control the value of the second.

Models are a tool for extracting patterns out of data.

The following code fits a model that predicts ```price``` from ```carat``` and then computes the **residuals** (the difference between the predicted value and the actual value). The residuals give us a view of the price of the diamond, once the effect of carat has been removed.

```{r}
library(modelr)

mod <- lm(log(price) ~ log(carat), data = diamonds)

diamonds2 <- diamonds %>% 
  add_residuals(mod) %>% 
  mutate(resid = exp(resid))

ggplot(data = diamonds2) + 
  geom_point(mapping = aes(x = carat, y = resid))
```


### ggplot2 calls

Typically, the first one or two arguments to a function are so important that you should know them by heart. The first two arguments to ```ggplot()``` are ```data``` and ```mapping```, and the first two arguments to ```aes()``` are ```x``` and ```y```.

Sometimes we’ll turn the end of a pipeline of data transformation into a plot. Watch for the transition from ```%>%``` to ```+```.

```{r}
diamonds %>% 
  count(cut, clarity) %>% 
  ggplot(aes(clarity, cut, fill = n)) + 
    geom_tile()
```


### RStudio projects

R experts keep all the files associated with a project together — input data, R scripts, analytical results, figures.


# WRANGLE

Data wrangling is the art of getting your data into R in a useful form for visualisation and modelling. There are three main parts to data wrangling:

* Import

* Tidy

* Transform

## TIBBLES

Tibbles *are* data frames, but they tweak some older behaviours to make life a little easier. 


### Creating Tibbles

**Tibbles** are one of the unifying features of the tidyverse. Most other R packages use regular data frames, so you might want to coerce a data frame to a tibble. You can do that with ```as_tibble()```:

You can create a new tibble from individual vectors with ```tibble()```.

```{r}
tibble(
  x = 1:5, 
  y = 1, 
  z = x ^ 2 + y
)
```


**`tibble()`never changes the type of the inputs (e.g. it never converts strings to factors!), it never changes the names of variables, and it never creates row names.**

It’s possible for a tibble to have column names that are not valid R variable names, aka **non-syntactic names**. For example, they might not start with a letter, or they might contain unusual characters like a space. To refer to these variables, you need to surround them with backticks, **`**:

Another way to create a tibble is with `tribble()`, short for **transposed tibble**. `tribble()` is customised for data entry in code: column headings are defined by formulas (i.e. they start with `~`), and entries are separated by commas. This makes it possible to lay out small amounts of data in easy to read form.

```{r}
tribble(
  ~x, ~y, ~z,
  #--|--|----
  "a", 2, 3.6,
  "b", 1, 8.5
)
```


### Tibbles vs. data.frame

**1. Printing**

Tibbles have a refined print method that shows only the first 10 rows, and all the columns that fit on screen. This makes it much easier to work with large data.


You can explicitly `print()` the data frame and control the number of rows (`n`) and the `width` of the display. `width = Inf` will display all columns:

```{r}
nycflights13::flights %>% 
  print(n = 10, width = Inf)
```


**2. Subsetting**

`[[` can extract by name or position; `$` only extracts by name but is a little less typing.

To use these in a pipe, you’ll need to use the special placeholder `.`:

```{r}
df <- tibble(
  x = runif(5),
  y = rnorm(5)
)
df %>% .$x
df %>% .[["x"]]
```

## Data Import

* `read_csv()` reads comma delimited files, `read_csv2()` reads semicolon separated files (common in countries where , is used as the decimal place), `read_tsv()`reads tab delimited files, and `read_delim()` reads in files with any delimiter.


* `read_fwf()` reads fixed width files. You can specify fields either by their widths with `fwf_widths()` or their position with `fwf_positions()`. `read_table()` reads a common variation of fixed width files where columns are separated by white space.


* `read_log()` reads Apache style log files


The first argument to `read_csv()` is the most important: **it’s the path to the file to read**.

`heights <- read_csv("<path_to_fileName>")`


`read_csv()` uses the first line of the data for the column names, which is a very common convention. 

**There are two cases where you might want to tweak this behaviour:**

**1. Sometimes there are a few lines of metadata at the top of the file.** You can use `skip = n` to skip the first `n` lines; or use `comment = "#"` to drop all lines that start with (e.g.) `#`.

**2. The data might not have column names** You can use `col_names = FALSE` to tell `read_csv()` not to treat the first row as headings, and instead label them sequentially from `X1` to `Xn`

Alternatively you can pass `col_names` a character vector which will be used as the column names:


```{r}
read_csv("1,2,3\n4,5,6", col_names = c("x", "y", "z"))
```

Another option that commonly needs tweaking is `na`: this specifies the value (or values) that are used to represent missing values in your file:

```{r}
read_csv("a,b,c\n1,2,.", na = ".")
```

You might wonder **why we’re not using `read.csv()`**. There are a few good reasons to favour readr functions over the base equivalents:

* They are typically much faster (~10x) than their base equivalents.

* They produce tibbles, they don’t convert character vectors to factors, use row names, or munge the column names.

* They are more reproducible. 


## Parsing a vector

The `parse_*()` functions take a character vector and return a more specialised vector like a logical, integer, or date.


The first argument is a character vector to parse, and the na argument specifies which strings should be treated as missing:

```{r}
parse_integer(c("1", "231", ".", "456"), na = ".")
```

1. `parse_logical()` and `parse_integer()` parse logicals and integers respectively. 

2. `parse_double()` is a strict numeric parser, and `parse_number()` is a flexible numeric parser. 

When parsing numbers, the most important option is the character you use for the decimal mark. You can override the default value of `.` by creating a new locale and setting the decimal_mark argument:

```{r}
parse_double("1,23", locale = locale(decimal_mark = ","))
```

`parse_number()` ignores non-numeric characters before and after the number. This is particularly useful for currencies and percentages, but also works to extract numbers embedded in text.

```{r}
parse_number("It cost $123.45")
```


3. `parse_character()`

UTF-8 can encode just about every character used by humans today, as well as many extra symbols.

```{r}
x1 <- "El Ni\xf1o was particularly bad this year"
parse_character(x1, locale = locale(encoding = "Latin1"))
```

readr provides `guess_encoding()` to help you figure out the right encoding. It’s not foolproof, and it works better when you have lots of text 

```{r}
guess_encoding(charToRaw(x1))
```

4. `parse_factor()` create factors, the data structure that R uses to represent categorical variables with fixed and known values.

Give `parse_factor()` a vector of known levels to generate a warning whenever an unexpected value is present:

```{r}
fruit <- c("apple", "banana")
parse_factor(c("apple", "banana"), levels = fruit)
```

5. `parse_datetime()`, `parse_date()`, and `parse_time()` allow you to parse various date & time specifications.

* `parse_datetime()` expects an ISO8601 date-time.

```{r}
parse_datetime("2010-10-01T2010")
```


* `parse_date()` expects a four digit year, a - or /, the month, a - or /, then the day:

```{r}
parse_date("2010-10-01")
```


* `parse_time()` expects the hour, :, minutes, optionally : and seconds, and an optional am/pm specifier



## Parsing a File

readr uses a heuristic to figure out the type of each column: it reads the first 1000 rows and uses some (moderately conservative) heuristics to figure out the type of each column. 

You can emulate this process with a character vector using `guess_parser()`, which returns readr’s best guess, and `parse_guess()` which uses that guess to parse the column.

* If you’re reading a very large file, you might want to set `n_max` to a smallish number like 10,000 or 100,000. That will accelerate your iterations while you eliminate common problems.

* If you’re having major parsing problems, sometimes it’s easier to just read into a character vector of lines with `read_lines()`, or even a character vector of length 1 with `read_file()`. 


## Writing to a File

readr also comes with two useful functions for writing data back to disk: `write_csv()` and `write_tsv()`. Both functions increase the chances of the output file being read back in correctly by:

* Always encoding strings in UTF-8.

* Saving dates and date-times in ISO8601 format so they are easily parsed elsewhere.

If you want to export a csv file to Excel, use `write_excel_csv()`

`write_csv(fileName, "filePath")`

Note that the type information is lost when you save to csv.

We can use `write_rds()` and `read_rds()` are uniform wrappers around the base functions `readRDS()` and `saveRDS()`. These store data in R’s custom binary format called RDS


## Other types of Data

* **haven** reads SPSS, Stata, and SAS files.

* **readxl** reads excel files (both .xls and .xlsx).

* **DBI**, along with a database specific backend (e.g. RMySQL, RSQLite, RPostgreSQL etc) allows you to run SQL queries against a database and return a data frame.

* **jsonlite** for json, and xml2 for XML. 


## Tidy Data

There are three interrelated rules which make a dataset tidy:

1. Each variable must have its own column.

2. Each observation must have its own row.

3. Each value must have its own cell.


## Spreading and Gathering

The first step is always to figure out what the variables and observations are.

The second step is to resolve one of two common problems:

* One variable might be spread across multiple columns.

* One observation might be scattered across multiple rows.

### Gathering

A common problem is a dataset where some of the column names are not names of variables, but values of a variable.

```{r}
table4a
table4a %>% 
  gather(`1999`, `2000`, key = "year", value = "cases")
table4b
table4b %>% 
  gather(`1999`, `2000`, key = "year", value = "population")
```

* Note that “1999” and “2000” are non-syntactic names (because they don’t start with a letter) so we have to surround them in backticks. 

### Spreading

Spreading is the opposite of gathering. You use it when an observation is scattered across multiple rows. 

```{r}
table2
```


We need two parameters:

1. The column that contains variable names, the `key` column. Here, it’s `type`.

2. The column that contains values forms multiple variables, the `value` column. Here it’s `count`.

```{r}
spread(table2, key = type, value = count)
```


**`gather()` makes wide tables narrower and longer; `spread()` makes long tables shorter and wider.**




## Separating and uniting

### Separate

`separate()` pulls apart one column into multiple columns, by splitting wherever a separator character appears.

```{r}
table3
```

`separate()` takes the name of the column to separate, and the names of the columns to separate into

```{r}
table3 %>%
  separate(rate, into = c("cases", "population"))
```

If you wish to use a specific character to separate a column, you can pass the character to the `sep` argument of `separate()`. 

We can ask `separate()` to try and convert to better types using `convert = TRUE`

You can also pass a vector of integers to `sep`. `separate()` will interpret the integers as positions to split at. Positive values start at 1 on the far-left of the strings; negative value start at -1 on the far-right of the strings. When using integers to separate strings, the length of sep should be one less than the number of names in into.

```{r}
table3 %>%
  separate(year, into = c("century", "year"), sep = 2)
```


### Unite

`unite()` is the inverse of `separate()`: it combines multiple columns into a single column.

`unite()` takes a data frame, the name of the new variable to create, and a set of columns to combine.

The default will place an underscore `(_)` between the values from different columns. Here we don’t want any separator so we use `""`:

```{r}
table5 %>%
  unite(new, century, year, sep = "")
```


## Missing Values

An explicit missing value is the presence of an absence; an implicit missing value is the absence of a presence.

You can set `na.rm = TRUE` in `gather()` to turn explicit missing values implicit.

`complete()` takes a set of columns, and finds all unique combinations. It then ensures the original dataset contains all those values, filling in explicit NAs where necessary.

You can also fill in these missing values with `fill()`. It takes a set of columns where you want missing values to be replaced by the most recent non-missing value (sometimes called last observation carried forward).


## Case Study

The best place to start is almost always to gather together the columns that are not variables. 

```{r}
who1 <- who %>% 
  gather(new_sp_m014:newrel_f65, key = "key", value = "cases", na.rm = TRUE)
who1 %>% 
  count(key)
```


# RELATIONAL DATA

Collectively, multiple tables of data are called **relational data** because it is the relations, not just the individual datasets, that are important.

Each relation always concerns a pair of tables. You don’t need to understand the whole thing; you just need to understand the chain of relations between the tables that you are interested in.


## Keys

The variables used to connect each pair of tables are called **keys**. A key is a variable (or set of variables) that uniquely identifies an observation.


* A **primary key** uniquely identifies an observation in its own table. For example, `planes$tailnum` is a primary key because it uniquely identifies each plane in the `planes` table.

* A **foreign key** uniquely identifies an observation in another table. For example, the `flights$tailnum` is a foreign key because it appears in the flights `table` where it matches each flight to a unique plane.

A variable can be both a primary key **and** a foreign key.


Once you’ve identified the primary keys in your tables, it’s good practice to verify that they do indeed uniquely identify each observation. One way to do that is to `count()` the primary keys and look for entries where `n` is greater than one:


```{r}
planes %>% 
  count(tailnum) %>% 
  filter(n > 1)
weather %>% 
  count(year, month, day, hour, origin) %>% 
  filter(n > 1)
```


Sometimes a table doesn’t have an explicit primary key: each row is an observation, but no combination of variables reliably identifies it. For example, what’s the primary key in the `flights` table? You might think it would be the date plus the flight or tail number, but neither of those are unique:

```{r}
flights %>% 
  count(year, month, day, flight) %>% 
  filter(n > 1)
```


If a table lacks a primary key, it’s sometimes useful to add one with `mutate()` and `row_number()`. That makes it easier to match observations if you’ve done some filtering and want to check back in with the original data. This is called a **surrogate key**.


A primary key and the corresponding foreign key in another table form a **relation**. 


## Mutating joints

A **mutating join** allows you to combine variables from two tables. It first matches observations by their keys, then copies across variables from one table to the other.

The join functions add variables to the right:

Imagine you want to add the full airline name to the `flights2` data. You can combine the `airlines` and `flights2` data frames with `left_join()`:

```{r}
flights %>%
  select(-origin, -dest) %>% 
  left_join(airlines, by = "carrier")
```

The result of joining airlines to `flights2` is an additional variable: `name`.


### Inner join

The simplest type of join is the **inner join**. An inner join matches pairs of observations whenever their keys are equal:

The output of an inner join is a new data frame that contains the key, the x values, and the y values. We use by to tell `dplyr` which variable is the key:

```{r}
x <- tribble(
  ~key, ~val_x,
     1, "x1",
     2, "x2",
     3, "x3"
)
y <- tribble(
  ~key, ~val_y,
     1, "y1",
     2, "y2",
     4, "y3"
)
x %>% 
  inner_join(y, by = "key")
```


### Outer joins

An outer join keeps observations that appear in at least one of the tables. There are three types of outer joins:

* A **left join** keeps all observations in x.

* A **right join** keeps all observations in y.

* A **full join** keeps all observations in x and y.

These joins work by adding an additional “virtual” observation to each table. This observation has a key that always matches (if no other key matches), and a value filled with `NA`.

The most commonly used join is the left join: you use this whenever you look up additional data from another table, because it preserves the original observations even when there isn’t a match. 


### Duplicate keys

**1.One table has duplicate keys.** This is useful when you want to add in additional information as there is typically a one-to-many relationship.

**2.Both tables have duplicate keys.** This is usually an error because in neither table do the keys uniquely identify an observation. 


### Defining the key columns

* The default, `by = NULL`, uses all variables that appear in both tables, the so called natural join. 


* A character vector, `by = "x"`. This is like a natural join, but uses only some of the common variables.


* A named character vector: `by = c("a" = "b")`. This will match variable `a` in table `x` to variable `b` in table `y`. The variables from `x` will be used in the output.


## Filtering Joins

Filtering joins match observations in the same way as mutating joins, but **affect the observations, not the variables**. There are two types:

* `semi_join(x, y)` keeps all observations in x that have a match in y.

* `anti_join(x, y)` drops all observations in x that have a match in y.

```{r}
top_dest <- flights %>%
  count(dest, sort = TRUE) %>%
  head(10)
flights %>% 
  semi_join(top_dest)
```

Anti-joins are useful for diagnosing join mismatches:


```{r}
flights %>%
  anti_join(planes, by = "tailnum") %>%
  count(tailnum, sort = TRUE)
```



## Join Problems

1. Start by identifying the variables that form the primary key in each table.

2. Check that none of the variables in the primary key are missing.

3. Check that your foreign keys match primary keys in another table. The best way to do this is with an `anti_join()`. 


## Set Operations

Useful when you want to break a single complex filter into simpler pieces. 

* `intersect(x, y)`: return only observations in both `x` and `y`.

*`union(x, y)`: return unique observations in `x` and `y`.

* `setdiff(x, y)`: return observations in `x`, but not in `y`.



# STRINGS


## String Basics

You can create strings with either **single quotes or double quotes**. Unlike other languages, there is no difference in behaviour. I recommend always using ", unless you want to create a string that contains multiple ".

`string2 <- 'If I want to include a "quote" inside a string, I use single quotes'`

To include a literal single or double quote in a string you can use `\` to “escape” it:

```
double_quote <- "\"" # or '"'
single_quote <- '\'' # or "'"
```

There are a handful of other special characters. The most common are `"\n"`, newline, and `"\t"`, tab.

### String length

We’ll use functions from **stringr**. These have more intuitive names, and all start with `str_`.

```{r}
str_length(c("a", "R for data science", NA))
```


### Combining strings

To combine two or more strings, use str_c() and use `sep` TO CONTROL THE SEPARATION:

```{r}
str_c("x", "y", sep = ", ")
```


### Subsetting strings

As well as the string, `str_sub()` takes `start` and end `arguments` which give the (inclusive) position of the substring:

```{r}
x <- c("Apple", "Banana", "Pear")
str_sub(x, 1, 3)
# negative numbers count backwards from end
str_sub(x, -3, -1)
```


### Locales

You can use `str_to_lower()` to change the text to lower case. You can also use `str_to_upper()` or `str_to_title()`.

The locale is specified as a ISO 639 language code, which is a two or three letter abbreviation. If you don’t already know the code for your language, Wikipedia has a good list. If you leave the locale blank, it will use the current locale, as provided by your operating system.

Another important operation that’s affected by the locale is sorting. The base R `order()` and `sort()` functions sort strings using the current locale. 



## Matching Patterns with Regular Expressions

To learn regular expressions, we’ll use `str_view()` and `str_view_all()`. **These functions take a character vector and a regular expression, and show you how they match**.


### Basic matches

The simplest patterns match exact strings:

```{r}
x <- c("apple", "banana", "pear")
str_view(x, "an")
```


* `.` matches any character (except a newline):

```{r}
str_view(x, ".a.")
```


To match an `.`, you need the regexp `\.`. Unfortunately this creates a problem. **We use strings to represent regular expressions, and `\` is also used as an escape symbol in strings. So to create the regular expression `\.` we need the string `"\\."`.**

```{r}
# To create the regular expression, we need \\
dot <- "\\."

# But the expression itself only contains one:
writeLines(dot)
#> \.

# And this tells R to look for an explicit .
str_view(c("abc", "a.c", "bef"), "a\\.c")
```


To match a literal `\` you need to escape it, creating the regular expression `\\`. To create that regular expression, you need to use a string, which also needs to escape `\`. That means to match a literal `\` you need to write `"\\\\"` — **you need four backslashes to match one**


### Anchors

By default, regular expressions will match any part of a string. It’s often useful to **anchor** the regular expression so that it matches from the start or end of the string. You can use:

* `^` to match the start of the string.

* `$` to match the end of the string.


```{r}
x <- c("apple", "banana", "pear")
str_view(x, "^a")
str_view(x, "a$")
str_view(x, "^ana$")
```


**If you begin with power (`^`), you end up with money (`$`).**


### Character classes and alternatives

There are four other useful tools:

* `\d`: matches any digit.

* `\s`: matches any whitespace (e.g. space, tab, newline).

* `[abc]`: matches a, b, or c.

* `[^abc]`: matches anything except a, b, or c.


### Repetition

The next step up in power involves controlling **how many times a pattern matches**:


* `?`: 0 or 1

* `+`: 1 or more

* `*`: 0 or more


```{r}
x <- "1888 is the longest year in Roman numerals: MDCCCLXXXVIII"
str_view(x, "CC?")
str_view(x, "CC+")
str_view(x, 'C[LX]+')
```


You can also specify the number of matches precisely:

* `{n}`: exactly n

* `{n,}`: n or more

* `{,m}`: at most m

* `{n,m}`: between n and m

```{r}
str_view(x, "C{2}")
str_view(x, "C{2,}")
str_view(x, "C{2,3}")
```

By default these matches are “greedy”: they will match the longest string possible. You can make them “lazy”, matching the shortest string possible by putting a ? after them. 


```{r}
str_view(x, 'C{2,3}?')
str_view(x, 'C[LX]+?')
```


### Grouping and backreferences

Parentheses are a way to disambiguate complex expressions. 

They also define **groups** that you can refer to with **backreferences**, like `\1`, `\2` etc. 

For example, the following regular expression finds all fruits that have a repeated pair of letters:

```{r}
str_view(fruit, "(..)\\1", match = TRUE)
```


## Tools

### Detect matches

To determine if a character vector matches a pattern, use `str_detect()`. It returns a logical vector the same length as the input:

```{r}
x <- c("apple", "banana", "pear")
str_detect(x, "e")
```


When you use a logical vector in a numeric context, `FALSE` becomes `0` and `TRUE` becomes `1`. That makes `sum()` and `mean()` useful if you want to answer questions about matches across a larger vector:

```{r}
# How many common words start with t?
sum(str_detect(words, "^t"))
# What proportion of common words end with a vowel?
mean(str_detect(words, "[aeiou]$"))
```

It’s often easier to combine multiple `str_detect()` calls with logical operators, rather than trying to create a single regular expression.

```{r}
# Find all words containing at least one vowel, and negate
no_vowels_1 <- !str_detect(words, "[aeiou]")
# Find all words consisting only of consonants (non-vowels)
no_vowels_2 <- str_detect(words, "^[^aeiou]+$")
identical(no_vowels_1, no_vowels_2)
```


A common use of `str_detect()` is to select the elements that match a pattern. You can do this with logical subsetting, or the convenient `str_subset()` wrapper:

```{r}
str_subset(words, "x$")
```

Typically, however, your strings will be one column of a data frame, and you’ll want to use filter instead:

```{r}
df <- tibble(
  word = words, 
  i = seq_along(word)
)
df %>% 
  filter(str_detect(words, "x$"))
```


A variation on `str_detect()` is `str_count()`: rather than a simple yes or no, it tells you how many matches there are in a string:

```{r}
df %>% 
  mutate(
    vowels = str_count(word, "[aeiou]"),
    consonants = str_count(word, "[^aeiou]")
  )
```


Note that **matches never overlap**. For example, in `"abababa"`, how many times will the pattern `"aba"` match? Regular expressions say two, not three:

```{r}
str_count("abababa", "aba")
str_view_all("abababa", "aba")
```


Many stringr functions come in **pairs**: one function works with a single match, and the other works with all matches. The second function will have the suffix `_all`.


### Extract matches

To extract the actual text of a match, use `str_extract()`.

Imagine we want to find all sentences that contain a colour. We first create a vector of colour names, and then turn it into a single regular expression:

```{r}
colours <- c("red", "orange", "yellow", "green", "blue", "purple")
colour_match <- str_c(colours, collapse = "|")
colour_match
```


Now we can select the sentences that contain a colour, and then extract the colour to figure out which one it is:


```{r}
has_colour <- str_subset(sentences, colour_match)
head(has_colour)
matches <- str_extract(has_colour, colour_match)
head(matches)
```


Note that `str_extract()` only extracts the first match. We can see that most easily by first selecting all the sentences that have more than 1 match:


```{r}
more <- sentences[str_count(sentences, colour_match) > 1]
str_view_all(more, colour_match)
```


```{r}
str_extract(more, colour_match)
```

To get all matches, use `str_extract_all()`. It returns a list:

```{r}
str_extract_all(more, colour_match)
```

If you use `simplify = TRUE`, `str_extract_all()` will return a matrix with short matches expanded to the same length as the longest:

```{r}
str_extract_all(more, colour_match, simplify = TRUE)
```



### Grouped matches

You can also use parentheses to extract parts of a complex match. 


```{r}
noun <- "(a|the) ([^ ]+)"

has_noun <- sentences %>%
  str_subset(noun) %>%
  head(10)
has_noun %>% 
  str_extract(noun)
```

`str_extract()` gives us the complete match; `str_match() gi`ves each individual component. Instead of a character vector, it returns a matrix, with one column for the complete match followed by one column for each group:

```{r}
has_noun %>% 
  str_match(noun)
```

If your data is in a tibble, it’s often easier to use `tidyr::extract()`. It works like `str_match()` but requires you to name the matches, which are then placed in new columns:

```{r}
tibble(sentence = sentences) %>% 
  tidyr::extract(
    sentence, c("article", "noun"), "(a|the) ([^ ]+)", 
    remove = FALSE
  )
```

Like `str_extract()`, if you want all matches for each string, you’ll need `str_match_all()`.


### Replacing matches

`str_replace()` and `str_replace_all()` allow you to replace matches with new strings. The simplest use is to replace a pattern with a fixed string:

```{r}
x <- c("apple", "pear", "banana")
str_replace(x, "[aeiou]", "-")
str_replace_all(x, "[aeiou]", "-")
```


With `str_replace_all()` you can perform multiple replacements by supplying a named vector:


```{r}
x <- c("1 house", "2 cars", "3 people")
str_replace_all(x, c("1" = "one", "2" = "two", "3" = "three"))
```


### Splitting

Use `str_split()` to split a string up into pieces. For example, we could split sentences into words:

```{r}
sentences %>%
  head(5) %>% 
  str_split(" ")
```


Like the other stringr functions that return a list, you can use `simplify = TRUE` to return a matrix:

```{r}
sentences %>%
  head(5) %>% 
  str_split(" ", simplify = TRUE)
```


Instead of splitting up strings by patterns, you can also split up by character, line, sentence and word `boundary()`s:

```{r}
x <- "This is a sentence.  This is another sentence."
str_view_all(x, boundary("word"))
str_split(x, boundary("word"))[[1]]
```


### Find matches

`str_locate()` and `str_locate_all()` give you the starting and ending positions of each match. These are particularly useful when none of the other functions does exactly what you want. 

You can use `str_locate()` to find the matching pattern, `str_sub()` to extract and/or modify them.


## Other Types of Pattern

When you use a pattern that’s a string, it’s automatically wrapped into a call to `regex()`:

```{r}
# The regular call:
str_view(fruit, "nana")
# Is shorthand for
str_view(fruit, regex("nana"))
```


You can use the other arguments of `regex()` to control details of the match:


* `ignore_case = TRUE` allows characters to match either their uppercase or lowercase forms. This always uses the current locale.

```{r}
bananas <- c("banana", "Banana", "BANANA")
str_view(bananas, "banana")
str_view(bananas, regex("banana", ignore_case = TRUE))
```


* `multiline = TRUE` allows `^` and `$` to match the start and end of each line rather than the start and end of the complete string.

```{r}
x <- "Line 1\nLine 2\nLine 3"
x
str_extract_all(x, "^Line")[[1]]
str_extract_all(x, regex("^Line", multiline = TRUE))[[1]]
```


* `comments = TRUE` allows you to use comments and white space to make complex regular expressions more understandable. Spaces are ignored, as is everything after `#`. To match a literal space, you’ll need to escape it: `"\\ "`.

```{r}
phone <- regex("
  \\(?     # optional opening parens
  (\\d{3}) # area code
  [)- ]?   # optional closing parens, dash, or space
  (\\d{3}) # another three numbers
  [ -]?    # optional space or dash
  (\\d{3}) # three more numbers
  ", comments = TRUE)
str_match("514-791-8141", phone)
```


* `dotall = TRUE` allows `.` to match everything, including `\n`.


## Other Uses of Regular Expressions

* `apropos()` searches all objects available from the global environment. This is useful if you can’t quite remember the name of the function.

```{r}
apropos("replace")
```


* `dir()` lists all the files in a directory. The pattern argument takes a regular expression and only returns file names that match the pattern. 

```{r}
head(dir(pattern = "\\.Rmd$"))
```


# FACTORS

In R, factors are used to work with categorical variables, variables that have a fixed and known set of possible values. They are also useful when you want to display character vectors in a non-alphabetical order.

Historically, factors were much easier to work with than characters. As a result, many of the functions in base R automatically convert characters to factors. 


## Creating Factors

To create a factor you must start by creating a list of the valid **levels**:

```{r}
month_levels <- c(
  "Jan", "Feb", "Mar", "Apr", "May", "Jun", 
  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
)
x1 <- c("Dec", "Apr", "Jan", "Mar")
y1 <- factor(x1, levels = month_levels)
sort(y1)
```

If you omit the levels, they’ll be taken from the data in alphabetical order.

Sometimes you’d prefer that the order of the levels match the order of the first appearance in the data. You can do that when creating the factor by setting levels to `unique(x)`


```{r}
f1 <- factor(x1, levels = unique(x1))
f1
levels(f1)
```

When factors are stored in a tibble, you can’t see their levels so easily. One way to see them is with **`count()**`:

```{r}
gss_cat %>%
  count(race)
```


Or with a **bar chart**:

```{r}
ggplot(gss_cat, aes(race)) +
  geom_bar()
```


By default, ggplot2 will drop levels that don’t have any values. You can force them to display with:


```{r}
ggplot(gss_cat, aes(race)) +
  geom_bar() +
  scale_x_discrete(drop = FALSE)
```



## Modifying Factor Order

`fct_reorder()` takes three arguments:

* `f`, the factor whose levels you want to modify.

* `x`, a numeric vector that you want to use to reorder the levels.

* Optionally, `fun`, a function that’s used if there are multiple values of x for each value of f. The default value is median.

```{r}
relig_summary <- gss_cat %>%
  group_by(relig) %>%
  summarise(
    age = mean(age, na.rm = TRUE),
    tvhours = mean(tvhours, na.rm = TRUE),
    n = n()
  )
ggplot(relig_summary, aes(tvhours, fct_reorder(relig, tvhours))) +
  geom_point()
```


`fct_relevel()` takes a factor, `f`, and then any number of levels that you want to move to the front of the line.


When you are colouring the lines on a plot, `fct_reorder2()` reorders the factor by the `y` values associated with the largest `x` values. This makes the plot easier to read because the line colours line up with the legend.

```{r}
by_age <- gss_cat %>%
  filter(!is.na(age)) %>%
  group_by(age, marital) %>%
  count() %>%
  mutate(prop = n / sum(n))

ggplot(by_age, aes(age, prop, colour = marital)) +
  geom_line(na.rm = TRUE)

ggplot(by_age, aes(age, prop, colour = fct_reorder2(marital, age, prop))) +
  geom_line() +
  labs(colour = "marital")
```


For bar plots, you can use `fct_infreq()` to **order levels in increasing frequency**: this is the simplest type of reordering because it doesn’t need any extra variables. You may want to combine with `fct_rev()`

```{r}
gss_cat %>%
  mutate(marital = marital %>% fct_infreq() %>% fct_rev()) %>%
  ggplot(aes(marital)) +
    geom_bar()
```


## Modifying Factor Levels

The most general and powerful tool is `fct_recode()`. It allows you to recode, or change, the value of each level. 

```{r}
gss_cat %>%
  mutate(partyid = fct_recode(partyid,
    "Republican, strong"    = "Strong republican",
    "Republican, weak"      = "Not str republican",
    "Independent, near rep" = "Ind,near rep",
    "Independent, near dem" = "Ind,near dem",
    "Democrat, weak"        = "Not str democrat",
    "Democrat, strong"      = "Strong democrat"
  )) %>%
  count(partyid)
```

`fct_recode()` will leave levels that aren’t explicitly mentioned as is, and will warn you if you accidentally refer to a level that doesn’t exist.

To combine groups, you can assign multiple old levels to the same new level:


```{r}
gss_cat %>%
  mutate(partyid = fct_recode(partyid,
    "Republican, strong"    = "Strong republican",
    "Republican, weak"      = "Not str republican",
    "Independent, near rep" = "Ind,near rep",
    "Independent, near dem" = "Ind,near dem",
    "Democrat, weak"        = "Not str democrat",
    "Democrat, strong"      = "Strong democrat",
    "Other"                 = "No answer",
    "Other"                 = "Don't know",
    "Other"                 = "Other party"
  )) %>%
  count(partyid)
```


If you want to collapse a lot of levels, `fct_collapse()` is a useful variant of `fct_recode()`. For each new variable, you can provide a vector of old levels:

```{r}
gss_cat %>%
  mutate(partyid = fct_collapse(partyid,
    other = c("No answer", "Don't know", "Other party"),
    rep = c("Strong republican", "Not str republican"),
    ind = c("Ind,near rep", "Independent", "Ind,near dem"),
    dem = c("Not str democrat", "Strong democrat")
  )) %>%
  count(partyid)
```


Sometimes you just want to **lump** together all the small groups to make a plot or table simpler. That’s the job of `fct_lump()`:

```{r}
gss_cat %>%
  mutate(relig = fct_lump(relig)) %>%
  count(relig)
```

We can use the `n` parameter to specify how many groups (excluding other) we want to keep:

```{r}
gss_cat %>%
  mutate(relig = fct_lump(relig, n = 10)) %>%
  count(relig, sort = TRUE) %>%
  print(n = Inf)
```



# DATES AND TIMES


## Creating date/times

You should always use the simplest possible data type that works for your needs. That means if you can use a date instead of a date-time, you should.

To get the current date or date-time you can use `today()` or `now()`.


### 1. From strings

`ymd()` and friends create dates.

```{r}
ymd(20170131)
ymd_hms("2017-01-31 20:11:59")
mdy_hm("01/31/2017 08:01")
ymd(20170131, tz = "UTC")
```


### 2. From individual components

Use `make_date()` for dates, or `make_datetime()` for date-times:

```{r}
flights %>% 
  select(year, month, day, hour, minute) %>% 
  mutate(departure = make_datetime(year, month, day, hour, minute))
```


### 3. From other types

* `as_datetime()` and `as_date()`:

```{r}
as_datetime(today())
as_date(now())
as_datetime(60 * 60 * 10)
as_date(365 * 10 + 2)
```


## Date-time Components


### Getting components

You can pull out individual parts of the date with the accessor functions `year()`, `month()`, `mday()` (day of the month), `yday()` (day of the year), `wday()` (day of the week), `hour()`, `minute()`, and `second()`.

```{r}
datetime <- ymd_hms("2016-07-08 12:34:56")
yday(datetime)
```


For `month()` and `wday()` you can set `label = TRUE` to return the abbreviated name of the month or day of the week. Set `abbr = FALSE` to return the full name.

```{r}
wday(datetime, label = TRUE, abbr = FALSE)
```


### Rounding

An alternative approach to plotting individual components is to round the date to a nearby unit of time, with `floor_date()`, `round_date()`, and `ceiling_date()`. 

Each function takes a vector of dates to adjust and then the name of the unit round down (floor), round up (ceiling), or round to.



## Time Spans


### Durations

In R, when you subtract two dates, you get a difftime object:


### Periods

Periods are time spans but don’t have a fixed length in seconds, instead they work with “human” times, like days and months.



# PIPES

The pipe, `%>%`, comes from the `magrittr` package. Packages in the tidyverse load %>% for you automatically, so you don’t usually load magrittr explicitly.


This means that the pipe won’t work for two classes of functions:

* 1. Functions that use the current environment. For example, `assign()` will create a new variable with the given name in the current environment.

```{r}
assign("x", 10)
x
"x" %>% assign(100)
x
```

The use of `assign` with the pipe does not work because it assigns it to a temporary environment used by `%>%`.


* 2. Functions that use lazy evaluation.


## When not to Use the Pipe

* Your pipes are longer than (say) ten steps. In that case, create intermediate objects with meaningful names.

* You have multiple inputs or outputs. If there isn’t one primary object being transformed, but two or more objects being combined together, don’t use the pipe.

* You are starting to think about a directed graph with a complex dependency structure. Pipes are fundamentally linear and expressing complex relationships with them will typically yield confusing code.


## Other Tools from `magrittr`

* `%T>%` works like `%>%` except that it returns the left-hand side instead of the right-hand side. 


* `%$%` "Explodes" out the variables in a data frame so that you can refer to them explicitly.


* For **assignment** `magrittr` provides the `%<>%` operator which allows you to replace code like:

```{r}
mtcars <- mtcars %>% 
  transform(cyl = cyl * 2)
```

with:

```{r}
mtcars %<>% transform(cyl = cyl * 2)
```


# FUNCTIONS

You should consider writing a function whenever you’ve copied and pasted a block of code more than twice.

There are three key steps to creating a new function:

1. You need to pick a **name** for the function.

2. You list the inputs, or **arguments**, to the function inside `function`.

3. You place the code you have developed in **body** of the function, a `{` block that immediately follows `function(...)`.

Generally, function names should be verbs, and arguments should be nouns.

If you have a family of functions that do similar things, make sure they have consistent names and arguments.

Where possible, avoid overriding existing functions and variables. 

Use **comments**, lines starting with `#`, to explain the "why"" of your code. 

Another important use of comments is to break up your file into easily readable chunks. Use long lines of `-` and `= `to make it easy to spot the breaks.

## Conditional Execution

An `if` statement allows you to conditionally execute code. It looks like this:

```
if (condition) {
  # code executed when condition is TRUE
} else {
  # code executed when condition is FALSE
}
```

**The standard return rule**: a function returns the last value that it computed. 

### Conditions

The `condition` must evaluate to either `TRUE` or `FALSE`. 

If it’s a vector, you’ll get a warning message; if it’s an `NA`, you’ll get an error. 

You can use `||` (or) and `&&` (and) to combine multiple logical expressions. 

These operators are "short-circuiting": as soon as `||` sees the first `TRUE` it returns `TRUE` without computing anything else. As soon as `&&` sees the first `FALSE` it returns `FALSE`.

### Multiple conditions

```
if (this) {
  # do that
} else if (that) {
  # do something else
} else {
  # 
}
```

If you end up with a very long series of chained if statements, you should consider rewriting.


*  The `switch()` function allows you to evaluate selected code based on position or name.


### Code style

Both `if` and `function` should (almost) always be followed by squiggly brackets (`{}`), and the contents should be indented by two spaces. 

An opening curly brace should never go on its own line and should always be followed by a new line. A closing curly brace should always go on its own line, unless it’s followed by `else`. Always indent the code inside curly braces.


```
# Good
if (y < 0 && debug) {
  message("Y is negative")
}

if (y == 0) {
  log(x)
} else {
  y ^ x
}
```

It’s ok to drop the curly braces if you have a very short if statement that can fit on one line:

```{r}
y <- 10
x <- if (y < 20) "Too low" else "Too high"
```


## Function Arguments

The arguments to a function typically fall into two broad sets: one set supplies the **data** to compute on, and the other supplies arguments that control the **details** of the computation. 

Generally, data arguments should come first. Detail arguments should go on the end, and usually should have default values.


### Choosing names

Generally you should prefer longer, more descriptive names, but there are a handful of very common, very short names:

* `x`, `y`, `z`: vectors.

* `w`: a vector of weights.

* `df`: a data frame.

* `i`, `j`: numeric indices (typically rows and columns).

* `n`: length, or number of rows.

* `p`: number of columns.


### Checking values

It’s good practice to check important preconditions, and throw an error (with `stop()`), if they are not true:

```{r}
wt_mean <- function(x, w) {
  if (length(x) != length(w)) {
    stop("`x` and `w` must be the same length", call. = FALSE)
  }
  sum(w * x) / sum(w)
}
wt_mean(c(1,1),3)
```

Be careful not to take this too far. There’s a tradeoff between how much time you spend making your function robust, versus how long you spend writing it. 


A useful compromise is the built-in `stopifnot()`: **it checks that each argument is TRUE, and produces a generic error message if not**.

```{r}
wt_mean <- function(x, w, na.rm = FALSE) {
  stopifnot(is.logical(na.rm), length(na.rm) == 1)
  stopifnot(length(x) == length(w))
  
  if (na.rm) {
    miss <- is.na(x) | is.na(w)
    x <- x[!miss]
    w <- w[!miss]
  }
  sum(w * x) / sum(w)
}
wt_mean(1:6, 6:1, na.rm = "foo")
```


Note that when using `stopifnot()` you assert what should be true rather than checking for what might be wrong.


###  Dot-dot-dot `(...)`

Many functions in R take an arbitrary number of inputs, they rely on a special argument: `...` (pronounced dot-dot-dot). This special argument captures any number of arguments that aren’t otherwise matched.

It’s useful because you can then send those `...` on to another function. This is a useful catch-all if your function primarily wraps another function. 

If you just want to capture the values of the `...`, use `list(...)`.


### Lazy evaluation

Arguments in R are lazily evaluated: they’re not computed until they’re needed. That means if they’re never used, they’re never called. 


## Return Values

There are two things you should consider when returning a value:

* Does returning early make your function easier to read?

* Can you make your function pipeable?

### Explicit return statements

The value returned by the function is usually the last statement it evaluates, but you can choose to return early by using `return()`.

It’s best to save the use of `return()` to signal that you can return early with a simpler solution.

```{r}
f <- function() {
  if (!x) {
    return(something_short)
  }

  # Do 
  # something
  # that
  # takes
  # many
  # lines
  # to
  # express
}
```


### Writing pipeable functions

If you want to write your own pipeable functions, it’s important to think about the return value. Knowing the return value’s object type will mean that your pipeline will "just work".

There are two basic types of pipeable functions:

* With **transformations**, an object is passed to the function’s first argument and a modified object is returned. 

* With **side-effects**, the passed object is not transformed. Instead, the function performs an action on the object.


## Environment

The **environment** of a function controls how R finds the value associated with a name.

```
f <- function(x) {
  x + y
} 
```

In many programming languages, this would be an error, because `y` is not defined inside the function. 

In R, this is valid code because R uses rules called **lexical scoping** to find the value associated with a name. 

Since `y` is not defined inside the function, R will look in the **environment** where the function was defined:

This behaviour seems like a recipe for bugs, and indeed you should avoid creating functions like this deliberately, but by and large it doesn’t cause too many problems.


# VECTORS

## Vector Basics

There are two types of vectors:

**1.** **Atomic vectors**, of which there are six types: **logical**, **integer**, **double**, **character**, **complex**, and **raw**. Integer and double vectors are collectively known as **numeric** vectors.

**2.** **Lists**, which are sometimes called recursive vectors because lists can contain other lists.


The chief difference between atomic vectors and lists is that atomic vectors are **homogeneous**, while lists can be **heterogeneous**.

Every vector has two key properties:

1. Its **type**, which you can determine with `typeof()`.

2. Its **length**, which you can determine with `length()`.

Vectors can also contain arbitrary additional metadata in the form of attributes. These attributes are used to create **augmented vectors** which build on additional behaviour. 

## Atomic Vectors

### Logical

Logical vectors are the simplest type of atomic vector because they can take only three possible values: `FALSE`, `TRUE`, and `NA`.

```{r}
1:10 %% 3 ==0
```


### Numeric

Integer and double vectors are known collectively as numeric vectors. In R, numbers are doubles by default. To make an integer, place an `L` after the number:

```{r}
typeof(1)
typeof(1L)
```

The distinction between integers and doubles is not usually important, but there are two important differences that you should be aware of:

*1* Doubles are approximations. 

*2* Integers have one special value: `NA`, while doubles have four: `NA`, `NaN`, `Inf` and `-Inf`.


### Character

Character vectors are the most complex type of atomic vector, because **each element of a character vector is a string, and a string can contain an arbitrary amount of data.**


## Using Atomic Vectors

### Coercion

*1* **Explicit coercion** happens when you call a function like `as.logical()`, `as.integer()`, `as.double()`, or `as.character()`.

*2* **Implicit coercion** happens when you use a vector in a specific context that expects a certain type of vector.


### Test functions

Sometimes you want to do different things based on the type of vector. One option is to use `typeof()`.

You can also use the `is_*` functions provided by purrr(e.g `is_numeric`)


### Scalars and recycling rules

R will also implicitly coerce the length of vectors. This is called vector **recycling**, because the shorter vector is repeated, or recycled, to the same length as the longer vector.


### Naming vectors

All types of vectors can be named. You can name them during creation with `c()`:

```{r}
c(x = 1, y = 2, z = 4)
```

The subsetting function is called like `x[a]`. You can subset a vector with:

*1.A numeric vector containing only integers.*

```{r}
x <- c("one", "two", "three", "four", "five")
x[c(3, 2, 5)]
```

*2.Subsetting with a logical vector keeps all values corresponding to a `TRUE` value.* 

```{r}
x <- c(10, 3, NA, 5, 8, 1, NA)
# All non-missing values of x
x[!is.na(x)]
# All even (or missing!) values of x
x[x %% 2 == 0]
```


*3.If you have a named vector, you can subset it with a character vector*


```{r}
x <- c(abc = 1, def = 2, xyz = 5)
x[c("xyz", "def")]
```


*4. The simplest type of subsetting is nothing, `x[]`, which returns the complete x.*


## Lists

You create a list with `list()`:

```{r}
x <- list(1, 2, 3)
x
```

A very useful tool for working with lists is `str()` because it focusses on the **structure**, not the contents.

```{r}
x_named <- list(a = 1, b = 2, c = 3)
str(x_named)
```

**A list can contain a mix of objects**

### Subsetting lists

```{r}
a <- list(a = 1:3, b = "a string", c = pi, d = list(-1, -5))
```

* `[` Extracts a sublist. The result will always be a list.

```{r}
str(a[4])
```

* `[[` Extracts a single component from a list. It removes a level of hierarchy from the list.

```{r}
str(a[[1]])
```


* `$` Shorthand for extracting named elements of a list. It works similarly to `[[` except that you don’t need to use quotes.

```{r}
a$a
a[["a"]]
```

**The distinction between `[` and `[[` is really important for lists, because `[[` drills down into the list while `[` returns a new, smaller list. **



## Attributes

Any vector can contain arbitrary additional metadata through its **attributes**. 

There are three very important attributes that are used to implement fundamental parts of R:

*1. Names* are used to name the elements of a vector.

*2. Dimensions* (dims, for short) make a vector behave like a matrix or array.

*3. Class*


# ITERATION


## For Loops

```{r}
df <- tibble(
  a = rnorm(10),
  b = rnorm(10),
  c = rnorm(10),
  d = rnorm(10)
)
output <- vector("double", ncol(df))  # 1. output
for (i in seq_along(df)) {            # 2. sequence
  output[[i]] <- median(df[[i]])      # 3. body
}
output
```

**Every for loop has three components:**

**1** The **output**: `output <- vector("double", length(x))`. Before you start the loop, you must always allocate sufficient space for the output. This is very important for efficiency: if you grow the for loop at each iteration using `c()` (for example), your for loop will be very slow.

A general way of creating an empty vector of given length is the `vector()` function. It has two arguments: the type of the vector (“logical”, “integer”, “double”, “character”, etc) and the length of the vector.

**2** The **sequence**: `i in seq_along(df)`. This determines what to loop over: each run of the for loop will assign `i` to a different value from `seq_along(df)`. It’s useful to think of `i` as a pronoun, like "it".

`seq_along()` it's a safe version of the familiar `1:length(l)`

**3** The **body**: `output[[i]] <- median(df[[i]])`. This is the code that does the work. It’s run repeatedly, each time with a different value for `i`. The first iteration will run `output[[1]] <- median(df[[1]])`, the second will run `output[[2]] <- median(df[[2]])`, and so on.


## For Loop Variations

### Modifying an existing object (instead of creating a new one)

**Output**: we already have the output — it’s the same as the input.

**Sequence**: we can think about a data frame as a list of columns, so **we can iterate over each column with `seq_along(df)`**.

**Body**: apply `rescale01()`.

```{r}
df <- tibble(
  a = rnorm(10),
  b = rnorm(10),
  c = rnorm(10),
  d = rnorm(10)
)
rescale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

for (i in seq_along(df)) {
  df[[i]] <- rescale01(df[[i]])
}
```

Typically you’ll be modifying a list or data frame with this sort of loop, so remember to use `[[`, not `[`. 

**I used `[[` in all my for loops: I think it’s better to use `[[` even for atomic vectors because it makes it clear that I want to work with a single element.**


### Looping over names or values, instead of indices

* Loop over the elements: `for (x in xs)`. This is most useful if you only care about side-effects, like plotting or saving a file, because it’s difficult to save the output efficiently.

* Loop over the names: `for (nm in names(xs))`. This gives you name, which you can use to access the value with `x[[nm]]`.

Iteration over the numeric indices is the most general form, because given the position you can extract both the name and the value:

```
for (i in seq_along(x)) {
  name <- names(x)[[i]]
  value <- x[[i]]
}
```


### Handling outputs of unknown length

A better solution to save the results in a list, and then combine into a single vector after the loop is done:

```{r}
means <- c(0, 1, 2)
out <- vector("list", length(means))
for (i in seq_along(means)) {
  n <- sample(100, 1)
  out[[i]] <- rnorm(n, means[[i]])
}
str(out)
str(unlist(out))
```

`unlist()` flatten a list of vectors into a single vector.


### Handling sequences of unknown length

Sometimes you don’t even know how long the input sequence should run for. This is common when doing simulations. For example, you might want to loop until you get three heads in a row. 

A **while** loop is simpler than for loop because it only has two components, a condition and a body:

```
while (condition) {
  # body
}
```

```{r}
flip <- function() sample(c("T", "H"), 1)

flips <- 0
nheads <- 0

while (nheads < 3) {
  if (flip() == "H") {
    nheads <- nheads + 1
  } else {
    nheads <- 0
  }
  flips <- flips + 1
}
flips
```


## For Loops vs. Functionals

**For loops are not as important in R as they are in other languages because R is a functional programming language**. 

This means that it’s possible to wrap up for loops in a function, and call that function instead of using the for loop directly.

```{r}
df <- tibble(
  a = rnorm(10),
  b = rnorm(10),
  c = rnorm(10),
  d = rnorm(10)
)

col_summary <- function(df, fun) {
  out <- vector("double", length(df))
  for (i in seq_along(df)) {
    out[i] <- fun(df[[i]])
  }
  out
}
col_summary(df, median)
col_summary(df, mean)
```

**The idea of passing a function to another function is extremely powerful idea, and it’s one of the behaviours that makes R a functional programming language.** 

The goal of using **purrr** functions instead of for loops is to allow you break common list manipulation challenges into independent pieces:

* How can you solve the problem for a single element of the list? Once you’ve solved that problem, purrr takes care of generalising your solution to every element in the list.

* If you’re solving a complex problem, how can you break it down into bite-sized pieces that allow you to advance one small step towards a solution? With purrr, you get lots of small pieces that you can compose together with the pipe.


## The map Functions

The pattern of looping over a vector, doing something to each element and saving the results is so common that the purrr package provides a family of functions to do it for you. 

* `map()` makes a list.

* `map_lgl()` makes a logical vector.

* `map_int()` makes an integer vector.

* `map_dbl()` makes a double vector.

* `map_chr()` makes a character vector.

Each function takes a vector as input, applies a function to each piece, and then returns a new vector that’s the same length (and has the same names) as the input. The type of the vector is determined by the suffix to the map function.

```{r}
df %>% map_dbl(median)
```


##Dealing with Failure

`safely()` always returns a list with two elements:

*1* `result` is the original result. If there was an error, this will be `NULL`.

*2* `error` is an error object. If the operation was successful, this will be `NULL`.

```{r}
safe_log <- safely(log)
str(safe_log(10))
str(safe_log("a"))
```


## Mapping Over Multiple Elements

Often you have multiple related inputs that you need iterate along in parallel. That’s the job of the `map2()` and `pmap()` functions. 

`map2()` iterates over two vectors in parallel:

```{r}
mu <- list(5, 10, -3)
sigma <- list(1, 5, 10)
map2(mu, sigma, rnorm, n = 5) %>% str()
```


`pmap()` takes a list of arguments. It's better to name the arguments:

```{r}
n <- list(1, 3, 5)
args2 <- list(mean = mu, sd = sigma, n = n)
args2 %>% 
  pmap(rnorm) %>% 
  str()
```


## Walk

You typically do this because you want to render output to the screen or save files to disk - the important thing is the action, not the return value. 

```{r}
x <- list(1, "a", 3)

x %>% 
  walk(print)
```

# MODEL BASICS

The goal of a model is to provide a simple low-dimensional summary of a dataset. 

There are two parts to a model:

*1* First, you define a **family of models** that express a precise, but generic, pattern that you want to capture.

*2* Next, you generate a **fitted model** by finding the model from the family that is the closest to your data. 


** A fitted model is just the closest model from a family of models**. That implies that you have the “best” model (according to some criteria); it doesn’t imply that you have a good model and it certainly doesn’t imply that the model is “true”.

**The `modelr` package which wraps around base R’s modelling functions to make them work naturally in a pipe.**






















