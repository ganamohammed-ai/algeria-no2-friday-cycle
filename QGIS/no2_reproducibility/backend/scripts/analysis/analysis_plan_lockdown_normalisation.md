# Analysis plan: weather normalisation of TROPOMI NO2 and the 2020 lockdown

Written on 2026-11-19, before any lockdown, placebo, gradient-boosting or version-step
estimate was computed. The elasticity results (Tables 3 to 7 of the draft) were already known
when this plan was written; nothing below has been run.

## Question
How much of the apparent change in NO2 during the spring 2020 containment in Algeria is
weather, and does ventilation-based weather adjustment give a result that passes checks at
sites with no local emissions?

## Data
Daily tables already in the project: `no2_ventilation_daily_v2.csv` (six study sites) and
`control/control_no2_ventilation_daily.csv` (two background controls). Outcome: ln of the
10 km disc NO2 column, usable days only.

## Primary model (sample A: 2019-08-07 to 2021-06-30)
Sample A ends before processor version v2.2 became operational (2021-07-01, van Geffen et al.
2022; TEMIS version table). OLS with Newey-West errors, lag 7 (as in the main model).

    ln NO2 = a + [b1 ln W10 + b2 ln BLH + g T2m] + month + weekday + trend + Ramadan + L * lockdown

* Raw model: the bracket is left out.
* Adjusted model: the bracket is included.
* Lockdown change in percent = 100 * (exp(L) - 1).
* Primary lockdown window: 2020-03-23 to 2020-06-13 (partial containment from 23 March in
  Algiers; main measures extended to 13 June in 44 of 48 wilayas).
* Alternative windows: 2020-03-29 to 2020-06-13 (partial containment in Oran, Constantine,
  Sétif and others from 29 March); 2020-03-23 to 2020-07-28.
* Ramadan dummy: 2020-04-24 to 2020-05-23 and 2021-04-13 to 2021-05-12 (approximate, +/- 1 day).

## Outputs per site
1. L raw and L adjusted with 95% CI.
2. Weather contribution = L raw minus L adjusted.
3. Placebo windows: windows of the same length starting every 7 days from 2019-08-07, inside
   sample A and not overlapping the true window. The true lockdown dummy is kept in the model.
   Report the SD of placebo L (raw and adjusted) and the rank of the actual |L| among placebo |L|.

## Negative-control rule (fixed now)
The adjustment is judged valid if (i) at both control sites the adjusted L has |L| below 5%
or a 95% CI that includes 0, and (ii) the SD of the adjusted placebo L is smaller than that of the
raw placebo L at most of the six study sites. If this fails, it is reported as a failure.

## Secondary analyses
* Gradient boosting counterfactual (scikit-learn HistGradientBoostingRegressor): features ln W10,
  ln BLH, T2m, u10, v10, day of year, weekday, trend; trained on sample A days outside the
  primary window; L_ML = mean of ln(observed) - ln(predicted) over lockdown days. Time-blocked
  cross-validated R2 is reported. Interval: block bootstrap (7-day blocks) of the residuals.
* Full-period trend (2019-08-07 to 2025-12-31), exploratory: same model with version steps at
  2021-07-01, 2021-11-14, 2022-07-17, 2023-03-12, 2023-11-26, 2024-09-08. Steps on 2024-11-16 and
  2025-11-22 are not used because the segments would be shorter than 90 and 60 days. Trend in
  percent per year, raw and adjusted, with and without the steps. Version changes before
  2021-07-01 are not modelled.

## Reporting
All sites are reported, in whatever direction. Any departure from this plan is listed in the paper
as a deviation.
