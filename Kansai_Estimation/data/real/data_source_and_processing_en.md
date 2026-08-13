# Data Sources and Processing for `transformed_data198_for_R.csv`

## 1. Data Overview

`transformed_data198_for_R.csv` is a two-period panel dataset created from municipal data in the Statistics Bureau of Japan's System of Social and Demographic Statistics. It covers 198 municipalities in the Kansai region. Two variables representing population aging and industrial structure are included as dependent variables, and four variables related to income, foreign residents, population migration, and inbound commuting are included as explanatory variables.

This CSV is not the original data released as official statistics. It is a processed dataset in which variables were selected, log-transformed, standardized, and renamed for analysis.

| Item | Description |
|---|---|
| Output file | `transformed_data198_for_R.csv` |
| Spatial units | All 198 municipalities in Shiga, Kyoto, Osaka, Hyogo, Nara, and Wakayama prefectures |
| Time periods | FY2015 and FY2020 |
| Data format | Long-form balanced panel |
| Number of observations | 396 rows (198 municipalities × 2 time periods) |
| Number of variables | 8 columns (2 region/time columns and 6 analysis-variable columns) |
| Missing values | None |

The numbers of municipalities by prefecture are 19 in Shiga, 26 in Kyoto, 43 in Osaka, 41 in Hyogo, 39 in Nara, and 30 in Wakayama.

## 2. Data Sources

### 2.1 System of Social and Demographic Statistics

The source data are the "Municipal Data: Social Indicators (Adjusted for Municipal Mergers and Boundary Changes)" from the Statistics Bureau of Japan's System of Social and Demographic Statistics.

The System of Social and Demographic Statistics is a systematically organized collection of regional statistical data that have been compiled and processed to describe a wide range of aspects of people's lives, including population and households, the natural environment, the economic base, the administrative base, education, labor, housing, health and medical care, and welfare and social security. For details, see the Statistics Bureau of Japan's [Statistical Observations of Prefectures and Municipalities](https://www.stat.go.jp/data/ssds/index.html) (Japanese).

The following three e-Stat tables were used for this dataset.

| Table ID | Title | Indicators used in this dataset |
|---|---|---|
| [0000020301](https://www.e-stat.go.jp/dbview?sid=0000020301) | A: Population and Households | Percentage of population aged 65 and over; foreign population (per 100,000 population); net in-migration rate (Japanese migrants) |
| [0000020304](https://www.e-stat.go.jp/dbview?sid=0000020304) | D: Administrative Base | Taxable income (per taxpayer) |
| [0000020306](https://www.e-stat.go.jp/dbview?sid=0000020306) | F: Labor | Percentage of workers employed in secondary industry; percentage of commuters from other municipalities |

The source data were retrieved on **July 3, 2026**. Because statistical tables on e-Stat may be updated, values retrieved at a later date may differ. When using e-Stat data, also consult the [Terms of Use for the Portal Site of Official Statistics of Japan (e-Stat)](https://www.e-stat.go.jp/terms-of-use).

### 2.2 Study Area and Region IDs

The study area comprises 198 municipalities in Shiga, Kyoto, Osaka, Hyogo, Nara, and Wakayama prefectures.

The output column `region` is not a municipal code. It is an internal ID sequentially assigned from 1 to 198 for analysis. The correspondence between these IDs, municipal codes, and municipality names is recorded in a separately created region-ID crosswalk. Therefore, when merging this dataset with other data, do not treat `region` as a municipal code; always use the crosswalk.

### 2.3 Time Periods

The e-Stat survey years "FY2015" and "FY2020" are used. In the output CSV, they are replaced with integers for analysis as follows.

| `time` | Survey year |
|---:|---|
| 1 | FY2015 |
| 2 | FY2020 |

## 3. Indicators and Transformations

### 3.1 Variable Mapping

| Output column | Role in the model | e-Stat indicator code | e-Stat indicator name | Original unit | Construction of output value |
|---|---|---|---|---|---|
| `y1` | Dependent variable 1 | `#A03506` | Percentage of population aged 65 and over | % | Standardized using pooled FY2015 and FY2020 observations |
| `y2` | Dependent variable 2 | `#F01202` | Percentage of workers employed in secondary industry | % | Standardized using pooled FY2015 and FY2020 observations |
| `x_common1` | Common explanatory variable 1 included in both equations | `#D02206` | Taxable income (per taxpayer) | JPY 1,000 | Natural logarithm `ln(x)` |
| `x_common2` | Common explanatory variable 2 included in both equations | `#A01601` | Foreign population (per 100,000 population) | Persons | Natural logarithm `ln(x + 1)` |
| `x_specific1_1` | Explanatory variable specific to the equation for dependent variable 1 | `#A05301` | Net in-migration rate (Japanese migrants) | % | No transformation (e-Stat value used as is) |
| `x_specific2_1` | Explanatory variable specific to the equation for dependent variable 2 | `#F02702` | Percentage of commuters from other municipalities | % | No transformation (e-Stat value used as is) |

`#A01601` is not the actual count of foreign residents; it is the foreign population per 100,000 population. `x_common2` is calculated by adding 1 to this value and taking its natural logarithm.

Variables reported as percentages were not converted to proportions by dividing by 100. For example, `x_specific2_1 = 23.5` represents 23.5%, not 0.235.

### 3.2 Standardization

For each of `y1` and `y2`, all 396 observations (198 municipalities × 2 time periods) are pooled and standardized using the following formula.

```text
z(i,t) = {x(i,t) - mean(x)} / sd_pop(x)
```

Here, `mean(x)` and `sd_pop(x)` are the mean and population standard deviation, respectively, calculated over all 396 observations. The denominator used for the standard deviation is `N`, not `N - 1`. Standardization is not performed separately by fiscal year or municipality.

### 3.3 Log Transformations

`x_common1` and `x_common2` are transformed as follows. `ln` denotes the natural logarithm.

```text
x_common1 = ln(taxable income)
x_common2 = ln(foreign population (per 100,000 population) + 1)
```

The value 1 is added to `x_common2` so that the log transformation can also be applied to observations whose original value is 0.

## 4. Data Processing Procedure

The data were processed as follows.

1. Obtain municipal data for FY2015 and FY2020 from e-Stat's "Social Indicators (Adjusted for Municipal Mergers and Boundary Changes)."
2. Restrict the study area to 198 municipalities in the six Kansai prefectures and organize the indicator codes, indicator names, FY2015 values, and FY2020 values.
3. Extract the six indicators and standardize `y1` and `y2` using all 396 observations.
4. Transform `x_common1` using `ln(x)` and `x_common2` using `ln(x + 1)`. Use the e-Stat values without transformation for the remaining two explanatory variables.
5. Replace the survey years with `time = 1` (FY2015) and `time = 2` (FY2020).
6. Order the observations within each municipality as `time = 1`, followed by `time = 2`, and output the long-form file `transformed_data198_for_R.csv`.

## 5. Output File Column Definitions

The columns in the output CSV are ordered as follows.

```text
region,time,y1,y2,x_common1,x_common2,x_specific1_1,x_specific2_1
```

| Column | Data type | Description |
|---|---|---|
| `region` | integer | Internal municipality ID (1–198) |
| `time` | integer | Time-period ID (1 = FY2015; 2 = FY2020) |
| `y1` | numeric | Standardized percentage of population aged 65 and over |
| `y2` | numeric | Standardized percentage of workers employed in secondary industry |
| `x_common1` | numeric | Natural logarithm of taxable income (per taxpayer) |
| `x_common2` | numeric | Natural logarithm of the foreign population (per 100,000 population) plus 1 |
| `x_specific1_1` | numeric | Original value of the net in-migration rate (Japanese migrants) (%) |
| `x_specific2_1` | numeric | Original value of the percentage of commuters from other municipalities (%) |

Each `region` has one row for `time = 1` and one row for `time = 2`. Each combination of `region` and `time` is unique.

## 6. Descriptive Statistics

The descriptive statistics were calculated from the 396 observations in the output CSV. SD is the population standard deviation.

| Column | Minimum | Maximum | Mean | SD |
|---|---:|---:|---:|---:|
| `y1` | -1.8887 | 3.5265 | 0.0000 | 1.0000 |
| `y2` | -2.1456 | 2.9951 | 0.0000 | 1.0000 |
| `x_common1` | 7.7212 | 8.7828 | 8.0081 | 0.1342 |
| `x_common2` | 0.0000 | 8.6415 | 6.6332 | 0.8319 |
| `x_specific1_1` | -5.7900 | 1.4100 | -0.4419 | 0.6902 |
| `x_specific2_1` | 7.8000 | 231.1000 | 37.7076 | 21.9559 |

## 7. Notes on Use

- This dataset is secondary data produced by processing official statistics for research and analysis. To confirm the definitions of the values, consult the e-Stat indicator names, explanatory notes, and original source materials.
- `region` is an ID for analysis and is not a municipal code.
- Because `y1` and `y2` have been standardized, they are not expressed in their original percentage units.
- Because `x_common1` and `x_common2` have been natural-log transformed, they are not expressed in their original units of thousands of yen or persons.
- e-Stat data described as "adjusted for municipal mergers and boundary changes" align municipalities observed at different times with administrative boundaries at a common reference date. Check the documentation for the version published on e-Stat that was used here to identify the reference date for those administrative boundaries.
- This dataset contains no missing values. However, this is the result of checks on the selected study area, time periods, and indicators; it does not imply that the entire e-Stat statistical table contains no missing values.

## 8. Example Source Citation

> Created from the Statistics Bureau of Japan, "System of Social and Demographic Statistics: Municipal Data, Social Indicators (Adjusted for Municipal Mergers and Boundary Changes)," A: Population and Households (Table ID 0000020301), D: Administrative Base (Table ID 0000020304), and F: Labor (Table ID 0000020306), Portal Site of Official Statistics of Japan (e-Stat), retrieved July 3, 2026.

When using this processed dataset in an academic paper, report, or redistributed material, cite the original data above as well as the location of this document, its author, its version or commit ID, and the date accessed.
