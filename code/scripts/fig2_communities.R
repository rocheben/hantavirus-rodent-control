## ---------------------------------------------------------------------------
## Figure 2 -- The two host communities: who is in them, which guild each
## species belongs to, which species are competent for the virus of that
## community, and how strongly the species meet each other.
##
## Everything is read from data/, so the figure cannot drift from the model:
## the species and their competence from data/species_traits.csv, the contacts
## from data/contact_matrix_<community>.csv.
## ---------------------------------------------------------------------------

source("code/R/paths.R"); source("code/R/parameters.R"); source("code/R/theme.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork) })

COMPETENCE_COLOURS <- c("primary"   = "#A50F15",
                        "secondary" = "#FB8072",
                        "none"      = "grey88")
CONTACT_WIDTH <- c("0.001" = 0.25, "0.01" = 0.9, "0.1" = 2.4)

## Species are laid out on a circle, commensal species on the upper half and
## wild species on the lower half, in the order of the data file.
community_layout <- function(p) {
  sp <- p$species
  up   <- which(sp$guild == "commensal")
  down <- which(sp$guild == "wild")
  ang <- numeric(nrow(sp))
  ang[up]   <- seq(150, 30, length.out = length(up))
  ang[down] <- seq(210, 330, length.out = length(down))
  data.frame(code = sp$code, latin = sp$latin, guild = sp$guild,
             competence = sp$competence,
             x = cos(ang * pi / 180), y = sin(ang * pi / 180),
             angle = ang, stringsAsFactors = FALSE)
}

community_panel <- function(community, title) {
  p   <- seov_parameters(community)
  pos <- community_layout(p)
  C   <- p$C

  edges <- expand.grid(i = pos$code, j = pos$code, stringsAsFactors = FALSE) %>%
    filter(i < j) %>%
    mutate(contact = mapply(function(a, b) C[a, b], i, j)) %>%
    filter(contact > 0) %>%
    left_join(pos %>% select(i = code, x, y), by = "i") %>%
    left_join(pos %>% select(j = code, xend = x, yend = y), by = "j") %>%
    mutate(level = factor(c("0.001", "0.01", "0.1")[match(round(contact, 5),
                                                          c(0.001, 0.01, 0.1))],
                          levels = names(CONTACT_WIDTH)))
  stopifnot(!any(is.na(edges$level)))

  pos$hjust <- ifelse(pos$x > 0.1, -0.12, ifelse(pos$x < -0.1, 1.12, 0.5))
  pos$vjust <- ifelse(abs(pos$x) > 0.1, 0.5, ifelse(pos$y > 0, -0.9, 1.9))
  pos$competence <- factor(pos$competence, levels = names(COMPETENCE_COLOURS))

  ggplot() +
    geom_hline(yintercept = 0, colour = "grey85", linetype = "dotted") +
    geom_segment(data = edges,
                 aes(x = x, y = y, xend = xend, yend = yend, linewidth = level),
                 colour = "grey45", alpha = 0.7, lineend = "round") +
    geom_point(data = pos, aes(x = x, y = y, fill = competence),
               shape = 21, size = 7, stroke = 0.6, colour = "grey25") +
    geom_text(data = pos, aes(x = x, y = y, label = latin,
                              hjust = hjust, vjust = vjust),
              fontface = "italic", size = 3.1) +
    annotate("text", x = -1.45, y = 0.92, label = "Commensal",
             size = 3.2, colour = "grey35", fontface = "bold", hjust = 0) +
    annotate("text", x = -1.45, y = -0.92, label = "Wild",
             size = 3.2, colour = "grey35", fontface = "bold", hjust = 0) +
    scale_linewidth_manual(values = CONTACT_WIDTH, drop = FALSE,
                           name = "Contact rate",
                           labels = c("0.001 (low)", "0.01 (moderate)", "0.1 (high)")) +
    scale_fill_manual(values = COMPETENCE_COLOURS, drop = FALSE,
                      name = "Competence for the virus",
                      labels = c("primary reservoir", "secondary reservoir",
                                 "not competent")) +
    coord_equal(xlim = c(-1.5, 1.5), ylim = c(-1.35, 1.35), clip = "off") +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(plot.title   = element_text(size = 10, face = "bold",
                                      margin = margin(b = 6)),
          legend.title = element_text(size = 9),
          legend.text  = element_text(size = 8),
          plot.margin  = margin(6, 10, 6, 10))
}

fig2 <- community_panel("France", "A. Lyon, urban park (SEOV)") +
        community_panel("Mexico", "B. Northeastern Mexico, peri-urban (SNV)") +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom", legend.box = "horizontal")

ggsave("figures/Fig2_communities.png", fig2, width = 11, height = 5.6, dpi = 300)
