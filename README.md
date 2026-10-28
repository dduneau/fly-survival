# Fly Survival

**Version 1.0.** A web app that compares the survival of groups of infected flies with a Cox model, as described in the protocol (section "Survival"). It draws the survival curves and the hazard ratios, runs the likelihood ratio test and the check of proportional hazards, and gives a methods paragraph and the R code of the analysis.

The analysis is R itself, with the packages `survival` and `coxme`. The app is built with shinylive, so R runs in the browser: there is no server, and the data are never uploaded. The first visit downloads R and the packages, which takes a few seconds.

## The data

A table with **one row per fly**, separated by commas, semicolons or tabs. Column names are free: you tell the app which column plays which role. The protocol uses these:

| Column | Content |
|---|---|
| `Day_of_injection` | The injection session. Used as a stratum. |
| `Box_ID` | The box the fly was kept in. Unique across the experiment, with letters in the name (`Box_03`, not `3`). |
| `Fly_ID` | One identifier per fly. |
| `Genotype`, `Sex`, `Treatment`, `Inoculum` | What the fly is and what it received. |
| `Time_post_injection` | Time of death, or of the last observation for flies that did not die. |
| `Censor` | 1 = died at that time, 0 = alive at that time (end of the observation, escape or loss). |

### Example data sets

- **`survival_CS_Prett.csv`**: real data. 779 Canton-S males and females injected with *P. rettgeri*, scored daily for 10 days, over 8 injection days with one box each.
- **`simulated_survival_3lines.csv`**: **simulated data**, to show designs with several genotypes and several boxes per day. 2,400 flies: 3 lines × 2 sexes, 4 injection days, 2 boxes per line per day, 50 flies of each sex per box. The daily mortality of each sex is that of the Canton-S data above; Line_B and Line_C multiply the hazard by 0.55 and 1.8 relative to Line_A, and each box adds a random effect with an SD of 0.3 on the log hazard. These lines do not exist.

## Using the app

1. **Data.** Load a file, or try one of the examples. **Keep only rows where** restricts the analysis to some values of a column (for example some genotypes), without editing the file.
2. **Columns.** Choose the time, the event, and the factor to **compare** (for example `Sex`), with its **reference group**: hazard ratios are given relative to it. **Analyse separately for each** runs the same analysis for each level of another column (for example each `Genotype`), with one panel per level. Choose the **day of injection** (stratum) and the **box**.
3. **Model.**
   - *Cox model with stratum (coxph)*: each day of injection has its own baseline hazard.
   - *Mixed Cox model (coxme)*: the box is also a random effect, as in the protocol. When each day holds a single box, box and day cannot be separated, and the app says so: use the first model.
   - **Analyse deaths up to**: restricts the hazard ratio to the acute phase. Flies still alive then are censored at that time. The window is shaded on the plot and reported everywhere.

The app checks the data first. It stops when a time is not a number, an event is not 0 or 1, or the comparison has fewer than two groups. It warns when a group has no death, when box labels are numbers or repeat across days, when each day holds a single box, and when no box holds more than one group (the comparison is then between boxes, so several boxes per group are needed).

## What it gives

- **Plots**: the proportion of survivors at each scoring, joined by straight lines, and the log hazard ratio of each group with its 95% CI and the numbers of flies and boxes. PDF or PNG.
- **Results**: hazard ratios with 95% CI (Wald, on the log scale), the likelihood ratio test against the model without the comparison, the box SD of the mixed model, flies and deaths per group, and `cox.zph` on the fixed-effects fit, with time scaled by survival (the default) and as recorded. A non-significant `cox.zph` is weak evidence when deaths fall on a few scorings; what matters is whether the curves cross.
- **Methods text**: a paragraph with the model, the window, the filter and the results.
- **R code**: the same analysis as a script, with the settings chosen, to rerun or extend it in R.

## Repository layout and updating

```
app_survival/   the app: app.R and the two example data sets (edit here)
docs/           the published site, made by shinylive from app_survival/ (do not edit)
```

To update the app, edit `app_survival/app.R`, then in R, from the repository folder:

```r
shinylive::export("app_survival", "docs")
file.create("docs/.nojekyll")
httpuv::runStaticServer("docs")   # check it locally
```

The empty file `docs/.nojekyll` lets GitHub Pages serve all of shinylive's files. Then commit and push. GitHub Pages serves `docs/` from the `main` branch.

To run the app in R instead: `shiny::runApp("app_survival")`.

## Citation

If you use this app, please cite the protocol: [authors, title, journal, year, DOI].
