### themes for ggplot2

library(ggplot2)

theme.pres   <- theme(axis.line    = element_blank(),
                      panel.border = element_blank(),
                      axis.title.y = element_text(face = "bold.italic",
                                                 color = "darkred",
                                                 size = 30),
                      axis.title.x = element_text(face = "bold.italic",
                                                 color = "darkred",
                                                 size = 30),
                      axis.text    = element_text(size = 30,
                                                 face = "bold",
                                                 color = "black"),
                      
                      plot.margin  = unit(c(20, 20, 20, 20),
                                          "points"),
                      
                      legend.key.size = unit(22, "mm"),
                      legend.text     = element_text(size = 22),
                      legend.title    = element_text(size = 24,
                                                    face = "bold"),
                      strip.text.x = element_text(size = 18,
                                                  face = "bold"),
                      strip.text.y = element_text(size = 18,
                                                  face = "bold"))

theme.paper <- theme(axis.title.y = element_text(face = "bold.italic",
                                                 color = "darkred",
                                                 size = 32),
                     axis.title.x = element_text(face = "bold.italic",
                                                 color = "darkred",
                                                 size = 32),
                     axis.text    = element_text(size = 32,
                                                 face = "bold",
                                                 color = "black"),
                     
                     plot.margin  = unit(c(20, 20, 20, 20),
                                         "points"),
                     
                     legend.key.size = unit(22, "mm"),
                     legend.text     = element_text(size = 32),
                     legend.title    = element_text(size = 32,
                                                    face = "bold"),
                     strip.text.x = element_text(size = 18,
                                                 face = "bold"),
                     strip.text.y = element_text(size = 18,
                                                 face = "bold"))