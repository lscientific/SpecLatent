#' American National Election Studies 2022 pilot study data
#'
#' A subset of data from the American National Election Studies 2022 Pilot Study
#'
#' @format ## `anes`
#' A data frame with 1,307 rows and 146 columns:
#' \describe{
#'   \item{party}{party reported as \code{pid3} in the raw data}
#'   \item{follow}{response to "Would you say you follow what’s going on in government and public affairs most of the time, some of the time, only now and then, or hardly at all?"}
#'   \item{reg}{response to "Are you registered to vote, or not?"}
#'   ...
#' }
#' @source <https://electionstudies.org/data-center/>
"anes"


#' 112th Senate Roll Call Data
#'
#' A subset of binary data from the 112th Senate roll call data
#'
#' @format ## `senate`
#' A data frame with 94 rows and 487 columns. Missing responses are imputed with average response rate.
#' \describe{
#'   \item{party}{party information: Dem. or Repub.}
#'   \item{roll1-486}{binary response to 486 roll calls}
#'   ...
#' }
#' @source <https://legacy.voteview.com/senate112.htm>
"senate"