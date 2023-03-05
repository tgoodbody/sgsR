#' Sample existing
#'
#' @description Sub-sample an existing sample. Four sampling methods are available:
#' \code{clhs}, \code{balanced}, \code{srs} and \code{strat}.
#'
#' @family sample functions
#'
#' @inheritParams sample_systematic
#' @inheritParams extract_strata
#' @inheritParams sample_clhs
#'
#' @param raster SpatRaster. Raster to guide the location of the samples. If \code{type = "clhs"} this raster can also
#' be used to define the population distributions to be used for sampling.
#' @param type Character. A string indicating the type of sampling method to use.
#' Possible values are \code{"clhs"}, \code{"balanced"}, \code{"srs"} and \code{"strat"}.
#' @param ... Additional arguments for the sampling method selected.
#'
#' @return An sf object of samples or a list object if \code{details = TRUE}
#'
#' @note When \code{type = "clhs"} or \code{type = "balanced"} all attributes in \code{existing} will be used for sampling.
#' Remove attributes not indented for sampling' prior to using this algorithm.
#'
#' @author Tristan R.H. Goodbody
#'
#' @export

sample_existing <- function(existing,
                            nSamp,
                            raster = NULL,
                            type = "clhs",
                            access = NULL,
                            buff_inner = NULL,
                            buff_outer = NULL,
                            plot = FALSE,
                            details = FALSE,
                            filename = NULL,
                            overwrite = FALSE,
                            ...) {
  #--- error handling ---#

  check_existing(
    existing = existing,
    raster = raster,
    nSamp = nSamp,
    plot = plot,
    details = details
  )

  existing <- prepare_existing(
    existing = existing,
    raster = raster,
    access = access,
    buff_inner = buff_inner,
    buff_outer = buff_outer
  )

  #--- sampling ---#

  if (type == "clhs") {
    samples <- sample_existing_clhs(
      existing = existing,
      nSamp = nSamp,
      filename = filename,
      details = details,
      overwrite = overwrite,
      raster = raster,
      ...
    )
  }

  if (type == "balanced") {
    samples <- sample_existing_balanced(
      existing = existing,
      nSamp = nSamp,
      filename = filename,
      overwrite = overwrite,
      ...
    )
  }

  if (type == "srs") {
    samples <- sample_existing_srs(
      existing = existing,
      nSamp = nSamp,
      filename = filename,
      overwrite = overwrite
    )
  }

  if (type == "strat") {
    toSample <- calculate_allocation_existing(
      existing = existing,
      nSamp = nSamp,
      ...
    )

    samples <- sample_existing_strat(
      existing = existing,
      toSample = toSample,
      filename = filename,
      overwrite = overwrite
    )

    if (isTRUE(details)) {
      samples <- list(samples = samples, details = toSample)
    }
  }

  #   #--- plotting ---#
  #   if (is.null(raster)) {
  #     if (any(anyCat != "numeric")) {
  #       nonNumeric <- oclass[oclass$Class != "numeric", ]
  #
  #       pecdfcat <- all %>%
  #         dplyr::select(dplyr::any_of(nonNumeric$Name)) %>%
  #         tidyr::pivot_longer(c(!type), names_to = "metric") %>%
  #         dplyr::group_by(type, metric, value) %>%
  #         dplyr::summarize(counts = dplyr::n()) %>%
  #         dplyr::mutate(counts = dplyr::case_when(
  #           type == "population" ~ counts / nrow(all %>% dplyr::filter(type == "population")),
  #           type == "sample" ~ counts / nrow(all %>% dplyr::filter(type == "sample")),
  #         )) %>%
  #         ggplot2::ggplot(ggplot2::aes(x = value, y = counts, class = type, fill = type)) +
  #         ggplot2::geom_col(position = "fill") +
  #         ggplot2::facet_grid(. ~ metric, scales = "free") +
  #         ggplot2::ggtitle(label = "Percent stacked barchart for existing samples") +
  #         ggplot2::xlab("Sampling metrics") +
  #         ggplot2::ylab("Representation percent")
  #     }
  #
  #     Numeric <- oclass[oclass$Class == "numeric", ] %>% dplyr::filter(!Name %in% c("X", "Y"))
  #
  #     pecdf <- all %>%
  #       dplyr::select(Numeric$Name, type) %>%
  #       dplyr::select(-dplyr::all_of(names(costLoc)[cost]), type) %>%
  #       tidyr::pivot_longer(c(!type), names_to = "metric") %>%
  #       ggplot2::ggplot(ggplot2::aes(value, class = type, colour = type)) +
  #       ggplot2::stat_ecdf(geom = "step") +
  #       ggplot2::facet_grid(. ~ metric, scales = "free") +
  #       ggplot2::ggtitle(label = "Empirical cumulative distributions for existing samples") +
  #       ggplot2::xlab("Sampling metrics") +
  #       ggplot2::ylab("Cumulative distribution")
  #   } else {
  #     if (!is.null(access)) {
  #       #--- plot sample locations and raster output ---#
  #       terra::plot(raster[[1]])
  #       suppressWarnings(plot(masked$buffer, add = TRUE, alpha = 0.1))
  #       suppressWarnings(plot(samples, add = TRUE, col = "black"))
  #     } else {
  #       #--- plot sample locations and raster output ---#
  #       terra::plot(raster[[1]])
  #       suppressWarnings(plot(samples, add = TRUE, col = "black"))
  #     }
  #
  #     if (any(anyCat != "numeric")) {
  #       nonNumeric <- oclass[oclass$Class != "numeric", ]
  #
  #       pecdfcat <- all %>%
  #         dplyr::select(nonNumeric$Name) %>%
  #         tidyr::pivot_longer(c(!type), names_to = "metric") %>%
  #         dplyr::group_by(type, metric, value) %>%
  #         dplyr::summarize(counts = dplyr::n()) %>%
  #         dplyr::mutate(counts = dplyr::case_when(
  #           type == "population" ~ counts / nrow(all %>% dplyr::filter(type == "population")),
  #           type == "sample" ~ counts / nrow(all %>% dplyr::filter(type == "sample")),
  #         )) %>%
  #         ggplot2::ggplot(ggplot2::aes(x = value, y = counts, class = type, fill = type)) +
  #         ggplot2::geom_col(position = "fill") +
  #         ggplot2::facet_grid(. ~ metric, scales = "free") +
  #         ggplot2::ggtitle(label = "Percent stacked barchart for existing samples") +
  #         ggplot2::xlab("Sampling metrics") +
  #         ggplot2::ylab("Representation percent")
  #     }
  #
  #     Numeric <- oclass[oclass$Class == "numeric", ] %>% dplyr::filter(Name %in% names(raster) & !Name %in% c("X", "Y"))
  #
  #     #--- generate ecdf curves comparing population and sample ---#
  #     pecdf <- all[, c(names(raster), "type")] %>%
  #       dplyr::select(Numeric$Name, type, -dplyr::all_of(names(costLoc)[cost])) %>%
  #       tidyr::pivot_longer(c(!type), names_to = "metric") %>%
  #       ggplot2::ggplot(ggplot2::aes(value, class = type, colour = type)) +
  #       ggplot2::stat_ecdf(geom = "step") +
  #       ggplot2::facet_grid(. ~ metric, scales = "free") +
  #       ggplot2::ggtitle(label = "Empirical cumulative distributions for existing samples") +
  #       ggplot2::xlab("Sampling metrics") +
  #       ggplot2::ylab("Cumulative distribution")
  #   }
  #
  #   if (exists("pecdfcat")) print(pecdfcat)
  #
  #   print(pecdf)
  # }

  #--- details ---#

  return(samples)
}
