#' Interpret Odds Ratio
#'
#' @param OR Value or vector of (log) odds ratio values.
#' @param rules Can be "`chen2010"` (default), `"cohen1988"` (through
#'   transformation to standardized difference, see [oddsratio_to_d()]) or custom set
#'   of [rules()].
#' @param p0 Baseline risk
#' @param log Are the provided values log odds ratio.
#' @inheritParams interpret
#'
#' @section Rules:
#'
#' Rules apply to OR as ratios, so OR of 10 is as extreme as a OR of 0.1 (1/10).
#'
#' - Chen et al. (2010) (`"chen2010"`; default)
#'   - **OR < 1.68** - Very small
#'   - **1.68 <= OR < 3.47** - Small
#'   - **3.47 <= OR < 6.71** - Medium
#'   - **OR >= 6.71 ** - Large
#' - Cohen (1988) (`"cohen1988"`, based on the [oddsratio_to_d()] conversion, see [interpret_cohens_d()])
#'   - **OR < 1.44** - Very small
#'   - **1.44 <= OR < 2.48** - Small
#'   - **2.48 <= OR < 4.27** - Medium
#'   - **OR >= 4.27 ** - Large
#'
#' @examples
#' interpret_oddsratio(1)
#' interpret_oddsratio(c(5, 2))
#'
#' @references
#' - Cohen, J. (1988). Statistical power analysis for the behavioral sciences
#' (2nd Ed.). New York: Routledge.
#'
#' - Chen, H., Cohen, P., & Chen, S. (2010). How big is a big odds ratio?
#' Interpreting the magnitudes of odds ratios in epidemiological studies.
#' Communications in Statistics-Simulation and Computation, 39(4), 860-864.
#'
#' - Sánchez-Meca, J., Marín-Martínez, F., & Chacón-Moscoso, S. (2003).
#' Effect-size indices for dichotomized outcomes in meta-analysis. Psychological
#' methods, 8(4), 448.
#'
#' @keywords interpreters
#' @export
interpret_oddsratio <- function(OR, p0 = 0.01, rules = "chen2010", log = FALSE, ...) {
  if (log) {
    f_transform <- function(.x) exp(abs(.x))
  } else {
    f_transform <- function(.x) exp(abs(log(.x)))
  }


  if (is.character(rules) && rules == "cohen1988") {
    d <- oddsratio_to_d(OR, log = log)
    return(interpret_cohens_d(d, rules = rules))
  }

  if (is.numeric(p0) && p0 >= 0 && p0 <= 1) {
    calculate_thresholds <- function(p0, d = c(0.2,0.5,0.8)) {
      z0 <- qnorm(p0)
      z1 <- z0 + d
      p1 <- pnorm(z1)
      or <- (p1*(1-p0))/(p0*(1-p1))
      return(or)
    }

    thresholds <- calculate_thresholds(p0)
  }

  rules <- .match.rules(
    rules,
    list(
      chen2010 = rules(thresholds, c("very small", "small", "medium", "large"),
        name = paste("chen2010; baseline risk:", p0), right = FALSE
      ),
      cohen1988 = NA # for correct error msg
    )
  )

  interpret(OR, rules, transform = f_transform)
}
