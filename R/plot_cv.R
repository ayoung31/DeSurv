
plot_cv_over_alpha = function(res,metric,errorbar=TRUE){
  ggplot(res,aes(x=alpha,y=!!sym(metric)))+
    geom_point()+
    geom_line()+
    theme_bw()+
    labs
}

plot_cv_over_k = function(metrics){
  
}