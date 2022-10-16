# first run nyc_routine.R and nyc_urgent.R
# then run this file

library(cowplot)
library(ggplot2)

p_urgent <-  ggplot(data = res_urgent,
                    aes(x = method,y = mean, ymin = lower, ymax = upper ))+
  geom_pointrange(aes(col = method, shape = method)) +
  scale_y_continuous(limits = c(0, 0.65)) + 
  xlab('Collection period: Urgent care') + ylab("Seroprevalence")+
  geom_errorbar(aes(ymin=lower, ymax=upper, col=method), width=0.2, cex=1)+ 
  scale_color_manual(name = "Method", values = c("gray25", "gray50", "black", "gray75"),
                     labels = c(expression(hat(rho)),
                                expression(hat(pi)[RG]), 
                                expression(hat(pi)[SRG]),
                                expression(hat(pi)[SRGM]))) + 
  scale_shape_manual(name = "Method", values = 15:18,
                     labels = c(expression(hat(rho)),
                                expression(hat(pi)[RG]), 
                                expression(hat(pi)[SRG]),
                                expression(hat(pi)[SRGM]))) + 
  facet_wrap(~month, strip.position="top", nrow=1, scales = "free_x") +
  theme_bw() + 
  theme(plot.title = element_text(size=16,face="bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.ticks.x = element_blank(),
        legend.position = c(.08,0.76),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.title = element_text(size=20),
        strip.text = element_text(size = 13)) 

p_routine <-  ggplot(data = res_routine,
                     aes(x = method,y = mean, ymin = lower, ymax = upper ))+
  geom_pointrange(aes(col = method, shape = method)) +
  scale_y_continuous(limits = c(0, 0.3)) + 
  xlab('Collection period: Routine care') + ylab("Seroprevalence")+
  geom_errorbar(aes(ymin=lower, ymax=upper, col=method), width=0.2, cex=1)+ 
  scale_color_manual(name = "Method", values = c("gray25", "gray50", "black", "gray75"),
                     labels = c(expression(hat(rho)),
                                expression(hat(pi)[RG]), 
                                expression(hat(pi)[SRG]),
                                expression(hat(pi)[SRGM]))) + 
  scale_shape_manual(name = "Method", values = 15:18,
                     labels = c(expression(hat(rho)),
                                expression(hat(pi)[RG]), 
                                expression(hat(pi)[SRG]),
                                expression(hat(pi)[SRGM]))) + 
  facet_wrap(~month, strip.position="top", nrow=1, scales = "free_x") +
  theme_bw() + 
  theme(plot.title = element_text(size=16,face="bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.ticks.x = element_blank(),
        legend.position = c(.08,0.76),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.title = element_text(size=20),
        strip.text = element_text(size = 13)) 

pdf(here("routine_urgent.pdf"),
    paper = "USr", width = 8.5, height = 11)
plot_grid(p_routine, 
          p_urgent,
          ncol = 1,
          labels = c("A.", "B."),
          label_size = 18,
          hjust = 0
)
dev.off()

