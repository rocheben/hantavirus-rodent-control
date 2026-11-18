## ---------------------------------------------------------------------------
## Figure 1 -- Structure of the multi-host SI model: the two compartments of a
## focal species, the demographic and epidemiological flows between them, and
## the coupling to the other species of the community through the contact
## matrix.
##
## The figure is drawn from the equations, not imported, so it stays in step
## with the model.
## ---------------------------------------------------------------------------

source("code/R/paths.R"); source("code/R/theme.R")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

BOX_W <- 1.5; BOX_H <- 0.9
FOCAL <- "#FFF4D8"; OTHER <- "grey96"; INF <- "#A50F15"

box <- function(x, y, label, fill, w = BOX_W, h = BOX_H)
  data.frame(x = x, y = y, label = label, fill = fill, w = w, h = h)

boxes <- bind_rows(
  box(0,   0, "S[i]", FOCAL),
  box(3.4, 0, "I[i]", FOCAL)
)

## the focal species, its flows -----------------------------------------------
arrow_style <- arrow(length = unit(0.16, "cm"), type = "closed")

flows <- data.frame(
  x    = c(1.5 + 0.25,  0 - 0.75,  0,          3.4,        0 - 0.75,   3.4 + 0.75),
  xend = c(3.4 - 0.85, -1.5,       0,          3.4,       -1.5,        4.9),
  y    = c(0,           0,        -0.45,      -0.45,       0,           0),
  yend = c(0,           0,        -1.15,      -1.15,       0,           0)
)

## infection flow, births, deaths and control, drawn one by one so that each
## label sits on its own arrow
seg <- function(x, y, xend, yend) data.frame(x = x, y = y, xend = xend, yend = yend)
arrows_df <- bind_rows(
  seg(0.80,  0.00,  2.60,  0.00),   # S -> I, force of infection
  seg(-1.60, 0.00, -0.80,  0.00),   # births -> S
  seg(0.00, -0.48,  0.00, -1.20),   # S -> death
  seg(3.40, -0.48,  3.40, -1.20),   # I -> death
  seg(0.00,  0.48,  0.00,  1.20),   # S -> control
  seg(3.40,  0.48,  3.40,  1.20)    # I -> control
)

labels_df <- data.frame(
  x     = c(1.70,  -1.75,  0.00,  3.40,  0.00,  3.40),
  y     = c(0.22,   0.00, -1.42, -1.42,  1.42,  1.42),
  label = c("lambda[i](t)", "b[i](t)", "d[i]~N[i]/K[i]", "d[i]~N[i]/K[i]",
            "rho[i](t)", "rho[i](t)"),
  hjust = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
)

## the rest of the community --------------------------------------------------
others <- bind_rows(
  box(-0.9, 3.1, "S[j]", OTHER, w = 1.1, h = 0.7),
  box( 1.2, 3.1, "I[j]", OTHER, w = 1.1, h = 0.7),
  box( 3.5, 3.1, "...",  OTHER, w = 1.1, h = 0.7)
)

coupling <- data.frame(x = 1.2, y = 2.70, xend = 1.68, yend = 0.34)

fig1 <- ggplot() +
  ## community box
  annotate("rect", xmin = -1.9, xmax = 4.5, ymin = 2.45, ymax = 3.75,
           fill = NA, colour = "grey70", linetype = "dashed") +
  annotate("text", x = -2.0, y = 3.98, hjust = 0, size = 3.1, colour = "grey35",
           parse = TRUE,
           label = "'The other species of the community,'~italic(j)==1~'...'~italic(m)") +
  geom_tile(data = others, aes(x = x, y = y, width = w, height = h),
            fill = OTHER, colour = "grey60") +
  geom_text(data = others, aes(x = x, y = y, label = label), parse = TRUE,
            size = 3.6) +
  ## coupling arrow
  geom_segment(data = coupling, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow_style, colour = INF, linewidth = 0.5,
               linetype = "22") +
  annotate("text", x = 4.95, y = 1.55, hjust = 0, size = 3.1, colour = INF,
           label = "Frequency-dependent transmission") +
  annotate("text", x = 4.95, y = 0.75, hjust = 0, size = 3.6, colour = INF,
           parse = TRUE,
           label = "lambda[i](t)==tau[i]~sum(C[i*','*j]~I[j](t)/N[tot](t), j==1, m)") +
  ## focal species
  geom_segment(data = arrows_df, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow_style, linewidth = 0.5, colour = "grey25") +
  geom_tile(data = boxes, aes(x = x, y = y, width = w, height = h),
            fill = FOCAL, colour = "grey25", linewidth = 0.7) +
  geom_text(data = boxes, aes(x = x, y = y, label = label), parse = TRUE,
            size = 5.5) +
  geom_text(data = labels_df, aes(x = x, y = y, label = label, hjust = hjust),
            parse = TRUE, size = 3.6) +
  annotate("text", x = -1.75, y = -0.32, size = 2.9, colour = "grey35",
           label = "births") +
  annotate("text", x = 0.00, y = -1.78, size = 2.9, colour = "grey35",
           label = "deaths") +
  annotate("text", x = 3.40, y = -1.78, size = 2.9, colour = "grey35",
           label = "deaths") +
  annotate("text", x = 0.00, y = 1.78, size = 2.9, colour = "grey35",
           label = "rodent control") +
  annotate("text", x = 3.40, y = 1.78, size = 2.9, colour = "grey35",
           label = "rodent control") +
  annotate("text", x = -2.2, y = -2.25, hjust = 0, size = 3.1, colour = "grey35",
           parse = TRUE,
           label = "'Focal species'~italic(i)*','~'with abundance'~N[i]==S[i]+I[i]") +
  coord_equal(xlim = c(-2.4, 8.9), ylim = c(-2.5, 4.2), clip = "off") +
  theme_void(base_size = 11) +
  theme(plot.margin = margin(6, 6, 6, 6))

ggsave("figures/Fig1_model_diagram.png", fig1, width = 9, height = 5.1, dpi = 300)
