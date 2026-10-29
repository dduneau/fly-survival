# Survival analysis after infection: Cox models for co-housed groups
# One row per fly; the user chooses which columns hold time, event, the comparison,
# the day of injection (stratum) and the box (random effect).
# Runs as a normal Shiny app, or in the browser with shinylive (no server).

library(shiny)
library(bslib)
library(survival)
library(ggplot2)
library(patchwork)

if (FALSE) { library(coxme); library(svglite); library(jsonlite) }   # declares them to shinylive, which bundles the packages it finds in library() calls
has_coxme <- requireNamespace("coxme", quietly = TRUE)
has_svglite <- requireNamespace("svglite", quietly = TRUE)
APP_VERSION <- "1.0"
`%||%` <- function(a, b) if (is.null(a)) b else a
EXAMPLES <- c(real = "survival_CS_Prett.csv",            # Canton-S males and females, real data
              sim  = "simulated_survival_3lines.csv")    # 3 lines, 6 boxes per day, simulated
PALETTE <- c("#1f5fbf", "#e8850c", "#2a9d6f", "#c0398a", "#7a5cc4", "#8c6d1f", "#3a3a3a", "#d1453b")
# Female orange and male blue in every app of the protocol, whatever the order; other groups take the palette in order
group_colours <- function(lv) {
  cols <- setNames(rep(NA_character_, length(lv)), lv)
  low <- tolower(lv)
  cols[low %in% c("female", "females", "f")] <- "#e8850c"
  cols[low %in% c("male", "males", "m")] <- "#1f5fbf"
  free <- setdiff(PALETTE, cols)
  cols[is.na(cols)] <- rep(free, length.out = sum(is.na(cols)))
  cols
}

# ------------------------------------------------------------------ reading
read_any <- function(path) {
  first <- readLines(path, n = 1, warn = FALSE)
  seps <- c(";" = ";", "," = ",", "\t" = "\t")
  counts <- sapply(seps, function(s) lengths(regmatches(first, gregexpr(s, first, fixed = TRUE))))
  sep <- seps[[which.max(counts)]]
  d <- read.table(path, header = TRUE, sep = sep, quote = "\"'", dec = ".",
                  stringsAsFactors = FALSE, check.names = FALSE, comment.char = "",
                  strip.white = TRUE, na.strings = c("NA", ""))
  attr(d, "sep") <- sep
  d
}

guess <- function(cols, candidates) {
  hit <- candidates[candidates %in% cols]
  if (length(hit)) hit[1] else ""
}

# ------------------------------------------------------------------ checks
check_data <- function(d, time, event, group, ref, stratum, box) {
  problems <- character(0); warnings <- character(0)
  if (!nzchar(time) || !nzchar(event) || !nzchar(group)) {
    return(list(problems = "Choose the columns for time, event and comparison.", warnings = warnings))
  }
  tt <- suppressWarnings(as.numeric(d[[time]]))
  ev <- suppressWarnings(as.numeric(d[[event]]))
  if (any(is.na(tt))) problems <- c(problems, sprintf("%d row(s) have a time that is not a number.", sum(is.na(tt))))
  if (any(!is.na(tt) & tt < 0)) problems <- c(problems, "Some times are negative.")
  if (any(is.na(ev)) || any(!ev %in% c(0, 1)))
    problems <- c(problems, sprintf("The event column must hold 1 (died) or 0 (censored); %d row(s) do not.",
                                    sum(is.na(ev) | !ev %in% c(0, 1))))
  g <- d[[group]]
  if (any(is.na(g) | g == "")) problems <- c(problems, sprintf("%d row(s) have no value for %s.", sum(is.na(g) | g == ""), group))
  if (length(unique(na.omit(g))) < 2) problems <- c(problems, sprintf("%s needs at least two groups to compare.", group))
  for (col in c(stratum, box)) if (nzchar(col)) {
    v <- d[[col]]
    if (any(is.na(v) | v == "")) problems <- c(problems, sprintf("%d row(s) have no value for %s.", sum(is.na(v) | v == ""), col))
    if (is.numeric(v)) warnings <- c(warnings, sprintf("%s holds numbers. It is used as a label, but names such as \"Box_1\" avoid confusion with a numeric variable.", col))
  }
  if (!length(problems)) {
    deaths <- tapply(ev, g, sum)
    if (any(deaths == 0)) warnings <- c(warnings, sprintf("No death in group(s): %s. Its hazard ratio cannot be estimated.",
                                                          paste(names(deaths)[deaths == 0], collapse = ", ")))
    if (nzchar(stratum) && nzchar(box)) {
      per_box <- tapply(d[[stratum]], d[[box]], function(x) length(unique(x)))
      if (any(per_box > 1)) warnings <- c(warnings, sprintf("Some boxes appear under more than one %s: box labels should be unique across the experiment.", stratum))
      boxes_per_day <- tapply(d[[box]], d[[stratum]], function(x) length(unique(x)))
      if (all(boxes_per_day == 1))
        warnings <- c(warnings, sprintf("Each %s holds a single box, so box and %s cannot be separated: the stratum already accounts for both. Use the Cox model with stratum only (coxph).", stratum, stratum))
    }
    within <- if (nzchar(box)) tapply(d[[group]], d[[box]], function(x) length(unique(x))) else NULL
    if (!is.null(within) && all(within == 1))
      warnings <- c(warnings, sprintf("No box holds more than one level of %s: the comparison is between boxes, so several boxes per group are needed.", group))
  }
  list(problems = problems, warnings = warnings)
}

# ------------------------------------------------------------------ analysis
prep <- function(d, time, event, group, ref, window) {
  d$.time <- as.numeric(d[[time]])
  d$.event <- as.numeric(d[[event]])
  lv <- unique(as.character(d[[group]]))
  d$.group <- factor(as.character(d[[group]]), levels = c(ref, setdiff(sort(lv), ref)))
  if (is.finite(window)) {
    d$.time_w <- pmin(d$.time, window)
    d$.event_w <- ifelse(d$.time <= window, d$.event, 0)
  } else {
    d$.time_w <- d$.time; d$.event_w <- d$.event
  }
  d
}

fit_models <- function(d, stratum, box, model) {
  st <- if (nzchar(stratum)) sprintf(" + strata(`%s`)", stratum) else ""
  if (model == "coxme") {
    re <- sprintf(" + (1 | `%s`)", box)
    f1 <- as.formula(paste0("Surv(.time_w, .event_w) ~ .group", st, re))
    f0 <- as.formula(paste0("Surv(.time_w, .event_w) ~ 1", st, re))
    m1 <- coxme::coxme(f1, data = d); m0 <- coxme::coxme(f0, data = d)
    beta <- coxme::fixef(m1)
    nfr <- nrow(m1$var) - length(beta)
    se <- sqrt(diag(as.matrix(m1$var))[nfr + seq_along(beta)])
    chisq <- 2 * unname(m1$loglik["Integrated"] - m0$loglik["Integrated"])
    re_sd <- sqrt(unlist(coxme::VarCorr(m1)))
  } else {
    f1 <- as.formula(paste0("Surv(.time_w, .event_w) ~ .group", st))
    f0 <- as.formula(paste0("Surv(.time_w, .event_w) ~ 1", st))
    m1 <- coxph(f1, data = d); m0 <- coxph(f0, data = d)
    beta <- coef(m1); se <- sqrt(diag(m1$var))
    chisq <- 2 * (tail(m1$loglik, 1) - tail(m0$loglik, 1))   # null model has a single log-likelihood
    re_sd <- NA
  }
  df <- length(beta)
  # proportional hazards, on the fixed-effects fit as cox.zph does not take random effects
  fz <- as.formula(paste0("Surv(.time_w, .event_w) ~ .group", st))
  zph <- tryCatch({
    z <- cox.zph(coxph(fz, data = d))
    data.frame(Test = c("Survival-scaled time (default)", "Time as recorded"),
               Chisq = c(z$table["GLOBAL", "chisq"], cox.zph(coxph(fz, data = d), transform = "identity")$table["GLOBAL", "chisq"]),
               df = as.integer(z$table["GLOBAL", "df"]),
               p = c(z$table["GLOBAL", "p"], cox.zph(coxph(fz, data = d), transform = "identity")$table["GLOBAL", "p"]))
  }, error = function(e) NULL)
  lvls <- levels(d$.group)
  hr <- data.frame(Group = lvls[-1], beta = unname(beta), se = unname(se))
  hr$lower <- hr$beta - qnorm(0.975) * hr$se
  hr$upper <- hr$beta + qnorm(0.975) * hr$se
  hr$HR <- exp(hr$beta); hr$HR_lower <- exp(hr$lower); hr$HR_upper <- exp(hr$upper)
  hr$p_wald <- 2 * pnorm(-abs(hr$beta / hr$se))
  hr <- rbind(data.frame(Group = lvls[1], beta = 0, se = 0, lower = 0, upper = 0, HR = 1, HR_lower = 1, HR_upper = 1, p_wald = NA), hr)
  list(hr = hr, lrt = data.frame(Chisq = chisq, df = df, p = pchisq(chisq, df, lower.tail = FALSE)),
       zph = zph, re_sd = re_sd)
}

counts <- function(d, stratum, box) {
  out <- data.frame(Group = levels(d$.group))
  out$Flies <- as.vector(table(d$.group))
  out$Deaths_in_window <- as.integer(tapply(d$.event_w, d$.group, sum))
  out$Deaths_total <- as.integer(tapply(d$.event, d$.group, sum))
  if (nzchar(box)) out$Boxes <- as.vector(tapply(d[[box]], d$.group, function(x) length(unique(x))))
  if (nzchar(stratum)) out$Days <- as.vector(tapply(d[[stratum]], d$.group, function(x) length(unique(x))))
  out
}

analyse <- function(d, s) {
  d <- prep(d, s$time, s$event, s$group, s$ref, s$window)
  parts <- if (nzchar(s$split)) split(d, d[[s$split]], drop = TRUE) else list(All = d)
  res <- lapply(names(parts), function(nm) {
    x <- parts[[nm]]
    x$.group <- droplevels(x$.group)
    fit <- tryCatch(fit_models(x, s$stratum, s$box, s$model), error = function(e) e)
    list(name = nm, data = x, counts = counts(x, s$stratum, s$box), fit = fit)
  })
  names(res) <- names(parts)
  res
}

# ------------------------------------------------------------------ plots
surv_points <- function(x) {
  times <- sort(unique(c(0, x$.time)))
  out <- lapply(levels(x$.group), function(g) {
    xx <- x[x$.group == g, ]
    if (!nrow(xx)) return(NULL)
    f <- survfit(Surv(.time, .event) ~ 1, data = xx)
    s <- summary(f, times = times, extend = TRUE)
    data.frame(Group = g, Time = s$time, Survival = s$surv)
  })
  do.call(rbind, out)
}

theme_surv <- theme_classic(base_size = 12) +
  theme(axis.line = element_line(linewidth = 0.6), axis.ticks = element_line(linewidth = 0.6),
        axis.text = element_text(colour = "black"), legend.position = "top", legend.title = element_blank(),
        strip.background = element_rect(fill = NA, colour = "black"), strip.text = element_text(face = "italic"))

plot_survival <- function(res, s) {
  pts <- do.call(rbind, lapply(res, function(r) { p <- surv_points(r$data); p$Panel <- r$name; p }))
  lv <- unique(unlist(lapply(res, function(r) levels(r$data$.group))))
  pts$Group <- factor(pts$Group, levels = lv)
  cols <- group_colours(lv)
  p <- ggplot(pts, aes(Time, Survival, colour = Group))
  if (is.finite(s$window) && s$window < max(pts$Time))
    p <- p + annotate("rect", xmin = 0, xmax = s$window, ymin = 0, ymax = 1, fill = grey(0.92))
  p <- p + geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
    scale_colour_manual(values = cols) +
    scale_y_continuous("Proportion of survivors", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = function(l) unique(round(pretty(l)))) +
    labs(x = s$xlab) + theme_surv + theme(aspect.ratio = 1)
  if (nzchar(s$split)) p <- p + facet_wrap(~Panel)
  p
}

plot_hr <- function(res, s) {
  rows <- do.call(rbind, lapply(res, function(r) {
    if (inherits(r$fit, "error")) return(NULL)
    h <- r$fit$hr; h$Panel <- r$name
    cn <- r$counts
    h$label <- paste0("(", cn$Flies[match(h$Group, cn$Group)],
                      if (!is.null(cn$Boxes)) paste0(" flies,\n", cn$Boxes[match(h$Group, cn$Group)], " boxes") else " flies", ")")
    h
  }))
  if (is.null(rows) || !nrow(rows)) return(NULL)
  lv <- unique(unlist(lapply(res, function(r) levels(r$data$.group))))
  rows$Group <- factor(rows$Group, levels = lv)
  cols <- group_colours(lv)
  lo <- min(rows$lower, 0, na.rm = TRUE); hi <- max(rows$upper, 0, na.rm = TRUE)
  ylab_y <- lo - 0.15 * (hi - lo + 0.1)
  win <- if (is.finite(s$window)) paste0(", up to ", s$window, " ", s$unit) else ""
  p <- ggplot(rows, aes(Group, beta)) +
    geom_hline(yintercept = 0, colour = "red", linetype = 4) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, linewidth = 0.5) +
    geom_point(aes(colour = Group), size = 2.5, show.legend = FALSE) +
    geom_text(aes(y = ylab_y, label = label), size = 3, vjust = 1, lineheight = 0.9) +
    scale_colour_manual(values = cols) +
    scale_y_continuous(paste0("log(HR) (± 95% CI)\nrelative to ", s$ref, win)) +
    coord_cartesian(ylim = c(ylab_y - 0.2 * (hi - lo + 0.1), hi + 0.05 * (hi - lo + 0.1))) +
    labs(x = NULL) + theme_surv
  if (nzchar(s$split)) p <- p + facet_wrap(~Panel)
  p
}

# ------------------------------------------------------------------ reports
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3)))

methods_text <- function(res, s, n_rows) {
  model_txt <- if (s$model == "coxme")
    sprintf("a mixed-effects Cox model (coxme) with %s as fixed effect%s and %s as random effect",
            s$group, if (nzchar(s$stratum)) sprintf(", %s as stratum", s$stratum) else "", s$box)
  else sprintf("a Cox proportional hazards model (coxph) with %s as fixed effect%s",
               s$group, if (nzchar(s$stratum)) sprintf(" and %s as stratum", s$stratum) else "")
  win <- if (is.finite(s$window)) sprintf("Deaths up to %s %s after injection were analysed; flies still alive then were censored at %s %s.", s$window, s$unit, s$window, s$unit)
         else "The whole observation period was analysed."
  res_txt <- sapply(res, function(r) {
    if (inherits(r$fit, "error")) return(sprintf("%s: the model could not be fitted (%s).", r$name, conditionMessage(r$fit)))
    h <- r$fit$hr[-1, ]; l <- r$fit$lrt
    hr_txt <- paste(sprintf("%s: HR = %.2f (95%% CI %.2f to %.2f)", h$Group, h$HR, h$HR_lower, h$HR_upper), collapse = "; ")
    paste0(if (nzchar(s$split)) paste0(r$name, ". ") else "", hr_txt,
           sprintf(", relative to %s; likelihood ratio test, chi2 = %.1f, df = %d, p = %s.", s$ref, l$Chisq, l$df, fmt_p(l$p)))
  })
  filt <- if (nzchar(s$filter_col) && length(s$filter_vals))
    sprintf(" The analysis was restricted to flies with %s = %s.", s$filter_col, paste(s$filter_vals, collapse = ", ")) else ""
  paste0("Survival was compared with ", model_txt, ", fitted in R (package survival",
         if (s$model == "coxme") " and coxme" else "", "). ", win,
         filt, " The effect was tested by a likelihood ratio test against the model without it, and hazard ratios are given with Wald 95% confidence intervals computed on the log scale. ",
         "Analysis made with the Survival app version ", APP_VERSION, " on ", n_rows, " flies.\n\n",
         paste(res_txt, collapse = "\n"))
}

r_code <- function(s, file_name, sep = ";") {
  q <- function(x) sprintf("`%s`", x)
  st <- if (nzchar(s$stratum)) sprintf(" + strata(%s)", q(s$stratum)) else ""
  re <- if (s$model == "coxme") sprintf(" + (1 | %s)", q(s$box)) else ""
  fun <- if (s$model == "coxme") "coxme" else "coxph"
  win <- if (is.finite(s$window)) s$window else "Inf"
  sub <- if (nzchar(s$split)) sprintf("\n# Analysis made separately for each level of %s; repeat for each:\n# d <- subset(d, %s == \"...\")\n", s$split, q(s$split)) else ""
  paste0(
'library(survival)', if (s$model == "coxme") '\nlibrary(coxme)' else '', '

d <- read.table("', file_name, '", header = TRUE, sep = "', if (sep == "\t") "\\t" else sep, '")
', if (nzchar(s$filter_col) && length(s$filter_vals))
  paste0('d <- subset(d, ', q(s$filter_col), ' %in% c(', paste0('"', s$filter_vals, '"', collapse = ", "), '))\n') else '', 'd$group <- relevel(factor(d$', q(s$group), '), ref = "', s$ref, '")', sub, '

# Window: deaths after ', win, ' are censored at ', win, '
window <- ', win, '
d$time  <- pmin(d$', q(s$time), ', window)
d$event <- ifelse(d$', q(s$time), ' <= window, d$', q(s$event), ', 0)

m1 <- ', fun, '(Surv(time, event) ~ group', st, re, ', data = d)
m0 <- ', fun, '(Surv(time, event) ~ 1', st, re, ', data = d)
',
if (s$model == "coxme") '
# Likelihood ratio test on the integrated log-likelihood
chisq <- 2 * unname(m1$loglik["Integrated"] - m0$loglik["Integrated"])
df <- length(fixef(m1))
c(chisq = chisq, df = df, p = pchisq(chisq, df, lower.tail = FALSE))

# Hazard ratios with 95% CI (log scale, then exponentiated)
beta <- fixef(m1)
se <- sqrt(diag(as.matrix(m1$var))[nrow(m1$var) - length(beta) + seq_along(beta)])
' else '
# Likelihood ratio test
anova(m0, m1)

# Hazard ratios with 95% CI (log scale, then exponentiated)
beta <- coef(m1)
se <- sqrt(diag(m1$var))
', '
data.frame(HR = exp(beta), lower = exp(beta - 1.96 * se), upper = exp(beta + 1.96 * se))

# Proportional hazards, on the fixed-effects fit
cox.zph(coxph(Surv(time, event) ~ group', st, ', data = d))
')
}

# ------------------------------------------------------------------ downloads
# Files are made in R and handed to the browser as a data URL. Shiny's downloadHandler is not used:
# under shinylive it depends on the service worker, and fails (a 404 page saved as .html) whenever the
# worker is not active, for example after a hard reload.
dl_button <- function(id, label) actionButton(id, label, icon = icon("download"), class = "btn-outline-secondary btn-sm me-1")
DL_JS <- "
Shiny.addCustomMessageHandler('send-file', function(m) {
  var a = document.createElement('a');
  a.href = 'data:' + m.mime + ';base64,' + m.data;
  a.download = m.name;
  document.body.appendChild(a); a.click(); a.remove();
});"
send_file <- function(session, name, mime, write_fun) {
  f <- tempfile(fileext = paste0(".", tools::file_ext(name)))
  write_fun(f)
  data <- jsonlite::base64_enc(readBin(f, "raw", file.info(f)$size))
  session$sendCustomMessage("send-file", list(name = name, mime = mime, data = gsub("\n", "", data)))
}

# ------------------------------------------------------------------ UI
ui <- page_sidebar(
  title = "Survival after infection",
  tags$head(tags$script(HTML(DL_JS)), tags$style(HTML("select.shiny-input-select:not(.selectized) { appearance: auto; -webkit-appearance: menulist; }"))),
  sidebar = sidebar(width = 330,
    h6("1. Data"),
    fileInput("file", NULL, accept = c(".csv", ".txt", ".tsv"), placeholder = "One row per fly"),
    div("Or try an example:"),
    actionLink("use_real", "Canton-S males and females (real data)"), br(),
    actionLink("use_sim", "Three lines, several boxes per day (simulated)"),
    div(class = "small text-muted mt-2", textOutput("data_info")),
    selectInput("filter_col", "Keep only rows where (optional)", NULL, selectize = FALSE),
    selectizeInput("filter_vals", NULL, choices = NULL, multiple = TRUE,
                   options = list(placeholder = "values to keep")),
    hr(),
    h6("2. Columns"),
    selectInput("time", "Time of death or censoring", NULL),
    selectInput("event", "Event (1 = died, 0 = censored)", NULL),
    selectInput("group", "Compare (fixed effect)", NULL),
    selectInput("ref", "Reference group: the least susceptible", NULL),
    div(class = "small text-muted", style = "margin-top:-0.6rem; margin-bottom:0.8rem",
        "Take the group that survives best, so that hazard ratios are above one and read in one direction (see the protocol).",
        textOutput("mortality", inline = TRUE)),
    selectInput("split", "Analyse separately for each (optional)", NULL, selectize = FALSE),
    selectInput("stratum", "Day of injection (stratum, optional)", NULL, selectize = FALSE),
    selectInput("box", "Box (random effect)", NULL, selectize = FALSE),
    hr(),
    h6("3. Model"),
    radioButtons("model", NULL,
                 c("Cox model with stratum (coxph)" = "coxph",
                   "Mixed Cox model, box as random effect (coxme)" = "coxme")),
    uiOutput("coxme_note"),
    numericInput("window", "Analyse deaths up to (empty = whole period)", value = NA, min = 0),
    textInput("unit", "Time unit", "days"),
    textInput("xlab", "x-axis label", "Days post-injection")
  ),
  uiOutput("messages"),
  navset_card_underline(
    nav_panel("Plots",
      plotOutput("plots", height = "430px"),
      div(dl_button("dl_pdf", "PDF"), dl_button("dl_png", "PNG"),
          if (has_svglite) dl_button("dl_svg", "SVG"))),
    nav_panel("Results",
      h6("Hazard ratios"), tableOutput("hr_table"),
      h6("Likelihood ratio test"), tableOutput("lrt_table"),
      h6("Flies and deaths"), tableOutput("count_table"),
      h6("Proportional hazards (cox.zph)"), tableOutput("zph_table"),
      p(class = "text-muted small",
        "A non-significant test is weak evidence when deaths fall on a few scorings, as with daily scoring.",
        "What matters in practice is whether the survival curves cross."),
      dl_button("dl_results", "Results (CSV)")),
    nav_panel("Methods text", verbatimTextOutput("methods"), dl_button("dl_methods", "Text")),
    nav_panel("R code", verbatimTextOutput("rcode"), dl_button("dl_code", "R script")),
    nav_panel("Data", tableOutput("preview"))
  )
)

# ------------------------------------------------------------------ server
server <- function(input, output, session) {
  rv <- reactiveValues(d = NULL, name = NULL)

  observeEvent(input$file, {
    d <- tryCatch(read_any(input$file$datapath), error = function(e) NULL)
    if (is.null(d)) { showNotification("This file could not be read as a table.", type = "error"); return() }
    rv$d <- d; rv$name <- input$file$name
  })
  observeEvent(input$use_real, { rv$d <- read_any(EXAMPLES[["real"]]); rv$name <- EXAMPLES[["real"]] })
  observeEvent(input$use_sim,  { rv$d <- read_any(EXAMPLES[["sim"]]);  rv$name <- EXAMPLES[["sim"]] })

  # rows kept by the optional filter
  dat <- reactive({
    req(rv$d); d <- rv$d
    if (nzchar(input$filter_col %||% "") && length(input$filter_vals))
      d <- d[as.character(d[[input$filter_col]]) %in% input$filter_vals, , drop = FALSE]
    d
  })
  observeEvent(input$filter_col, {
    vals <- if (nzchar(input$filter_col)) sort(unique(as.character(rv$d[[input$filter_col]]))) else character(0)
    updateSelectizeInput(session, "filter_vals", choices = vals, selected = vals)
  }, ignoreInit = TRUE)

  observeEvent(rv$d, {
    cols <- names(rv$d); none <- c("(none)" = "")
    updateSelectInput(session, "time", choices = cols, selected = guess(cols, c("Time_to_day", "Time_to_death", "Time", "time", "Time_post_injection")))
    updateSelectInput(session, "event", choices = cols, selected = guess(cols, c("Censor", "Status", "status", "Event", "event")))
    updateSelectInput(session, "group", choices = cols, selected = guess(cols, c("Sex", "Genotype", "Treatment")))
    updateSelectInput(session, "split", choices = c(none, cols), selected = "")
    updateSelectInput(session, "filter_col", choices = c(none, cols), selected = "")
    updateSelectInput(session, "stratum", choices = c(none, cols), selected = guess(cols, c("Day_of_injection", "Day", "Experiment")))
    updateSelectInput(session, "box", choices = c(none, cols), selected = guess(cols, c("Box_ID", "Box", "box")))
    t <- suppressWarnings(as.numeric(rv$d[[guess(cols, c("Time_to_day", "Time_to_death", "Time", "time", "Time_post_injection"))]]))
    if (length(t) && any(is.finite(t))) updateNumericInput(session, "window", value = NA, max = max(t, na.rm = TRUE))
  })

  observeEvent(list(dat(), input$group), {
    req(nzchar(input$group))
    lv <- sort(unique(as.character(dat()[[input$group]])))
    cur <- isolate(input$ref)
    updateSelectInput(session, "ref", choices = c("Choose the reference group" = "", lv),
                      selected = if (!is.null(cur) && cur %in% lv) cur else "")
  })
  # mortality of each group over the analysed window, to help choose the reference
  output$mortality <- renderText({
    req(rv$d, nzchar(input$group), nzchar(input$time), nzchar(input$event))
    d <- dat(); t <- suppressWarnings(as.numeric(d[[input$time]])); e <- suppressWarnings(as.numeric(d[[input$event]]))
    w <- if (is.null(input$window) || is.na(input$window)) Inf else input$window
    died <- e == 1 & t <= w
    m <- tapply(died, as.character(d[[input$group]]), mean, na.rm = TRUE)
    m <- sort(m)
    paste0(" Mortality", if (is.finite(w)) paste0(" up to ", w, " ", input$unit) else "", ": ",
           paste(sprintf("%s %.0f%%", names(m), 100 * m), collapse = ", "), ".")
  })

  output$data_info <- renderText({
    if (is.null(rv$d)) return("No data loaded.")
    n <- nrow(dat())
    sprintf("%s: %d flies%s, %d columns.", rv$name, nrow(rv$d), if (n < nrow(rv$d)) sprintf(" (%d kept by the filter)", n) else "", ncol(rv$d))
  })
  output$coxme_note <- renderUI({
    if (!has_coxme) p(class = "text-danger small", "The coxme package is not available here: only the Cox model with stratum can be fitted.")
  })

  settings <- reactive({
    w <- if (is.null(input$window) || is.na(input$window)) Inf else input$window
    list(time = input$time, event = input$event, group = input$group, ref = input$ref,
         split = input$split, stratum = input$stratum, box = input$box,
         model = if (input$model == "coxme" && has_coxme) "coxme" else "coxph",
         window = w, unit = input$unit, xlab = input$xlab,
         filter_col = input$filter_col %||% "", filter_vals = input$filter_vals %||% character(0))
  })

  checks <- reactive({
    req(rv$d); s <- settings()
    if (!nrow(dat())) return(list(problems = "The filter keeps no row: choose the values to keep.", warnings = character(0)))
    ch <- check_data(dat(), s$time, s$event, s$group, s$ref, s$stratum, s$box)
    if (s$model == "coxme" && !nzchar(s$box)) ch$problems <- c(ch$problems, "The mixed model needs a box column.")
    if (nzchar(s$split) && s$split == s$group)
      ch$problems <- c(ch$problems, sprintf("%s is chosen both to compare and to analyse separately: each analysis would then hold a single group. Set \"Analyse separately for each\" to (none), or to another column.", s$group))
    if (nzchar(s$ref) && nzchar(s$group) && !s$ref %in% as.character(dat()[[s$group]])) ch$problems <- c(ch$problems, "Choose a reference group.")
    ch$no_ref <- !nzchar(s$ref %||% "")
    ch
  })

  results <- reactive({
    req(rv$d); ch <- checks(); req(!length(ch$problems), !ch$no_ref)
    analyse(dat(), settings())
  })

  output$messages <- renderUI({
    if (is.null(rv$d)) return(div(class = "alert alert-info", "Load a file with one row per fly, or try one of the examples."))
    ch <- checks()
    if (!length(ch$problems) && ch$no_ref)
      return(div(class = "alert alert-info", strong("Choose the reference group."),
                 " Take the least susceptible group, the one that survives best: hazard ratios are then above one, and mortality and load read in one direction. The mortality of each group is shown under the menu."))
    fit_err <- if (!length(ch$problems) && !ch$no_ref) {
      errs <- Filter(function(r) inherits(r$fit, "error"), results())
      vapply(errs, function(r) sprintf("%s: the model could not be fitted (%s).", r$name, conditionMessage(r$fit)), "")
    } else character(0)
    tagList(
      if (length(ch$problems)) div(class = "alert alert-danger", tags$ul(lapply(ch$problems, tags$li))),
      if (length(c(ch$warnings, fit_err))) div(class = "alert alert-warning", tags$ul(lapply(c(ch$warnings, fit_err), tags$li))))
  })

  both_plots <- reactive({
    res <- results(); s <- settings()
    p1 <- plot_survival(res, s); p2 <- plot_hr(res, s)
    if (is.null(p2)) p1 else p1 + p2 + plot_layout(widths = c(2, 1))
  })
  output$plots <- renderPlot(both_plots(), res = 96)

  hr_tab <- reactive({
    res <- results()
    do.call(rbind, lapply(res, function(r) {
      if (inherits(r$fit, "error")) return(NULL)
      h <- r$fit$hr
      data.frame(Panel = r$name, Group = h$Group, `log HR` = round(h$beta, 3), SE = round(h$se, 3),
                 HR = round(h$HR, 2), `95% CI` = ifelse(h$se == 0, "reference", sprintf("%.2f to %.2f", h$HR_lower, h$HR_upper)),
                 `p (Wald)` = fmt_p(h$p_wald), check.names = FALSE)
    }))
  })
  output$hr_table <- renderTable({ h <- hr_tab(); if (!nzchar(settings()$split)) h$Panel <- NULL; h })
  output$lrt_table <- renderTable({
    do.call(rbind, lapply(results(), function(r) if (!inherits(r$fit, "error"))
      data.frame(Panel = r$name, Chisq = round(r$fit$lrt$Chisq, 2), df = r$fit$lrt$df, p = fmt_p(r$fit$lrt$p),
                 `Box SD (log hazard)` = if (is.na(r$fit$re_sd[1])) "" else round(r$fit$re_sd[1], 3), check.names = FALSE)))
  })
  output$count_table <- renderTable({
    do.call(rbind, lapply(results(), function(r) cbind(Panel = r$name, r$counts)))
  })
  output$zph_table <- renderTable({
    do.call(rbind, lapply(results(), function(r) if (!inherits(r$fit, "error") && !is.null(r$fit$zph))
      cbind(Panel = r$name, transform(r$fit$zph, Chisq = round(Chisq, 2), p = fmt_p(p)))))
  })

  methods <- reactive(methods_text(results(), settings(), nrow(dat())))
  output$methods <- renderText(methods())
  code <- reactive({ req(rv$d); r_code(settings(), rv$name, attr(rv$d, "sep") %||% ";") })
  output$rcode <- renderText(code())
  output$preview <- renderTable(head(dat(), 20))

  observeEvent(input$dl_pdf, send_file(session, "survival_plots.pdf", "application/pdf",
    function(f) ggsave(f, both_plots(), width = 9, height = 4.5)))
  observeEvent(input$dl_png, send_file(session, "survival_plots.png", "image/png",
    function(f) ggsave(f, both_plots(), width = 9, height = 4.5, dpi = 300)))
  observeEvent(input$dl_svg, send_file(session, "survival_plots.svg", "image/svg+xml",
    function(f) ggsave(f, both_plots(), width = 9, height = 4.5, device = svglite::svglite)))
  observeEvent(input$dl_results, send_file(session, "survival_results.csv", "text/csv",
    function(f) write.table(hr_tab(), f, sep = ";", row.names = FALSE)))
  observeEvent(input$dl_methods, send_file(session, "survival_methods.txt", "text/plain",
    function(f) writeLines(methods(), f)))
  observeEvent(input$dl_code, send_file(session, "survival_analysis.R", "text/plain",
    function(f) writeLines(code(), f)))
}

shinyApp(ui, server)
