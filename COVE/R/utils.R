print.ve_nai <- function(object, what = "ve", ...){
  
  if(what == "ve"){
    for(i in 1:length(object$ve_fit$t0)){
      boot_t0 <- object$boot_est[[paste0("t0_", object$ve_fit$t0[i])]]
      cat(paste0("------- Results for t0 = ", object$ve_fit$t0[i], " -------\n"))
      for(ve in c("ve_i", "ve_s", "ve_ai", "ve_nai")){
        cat(paste0(ve, ": ", round(object$ve_fit[[ve]][i], 4), 
              " (", round(boot_t0[[ve]]$lower,4), ", ", 
              round(boot_t0[[ve]]$upper,4), ") \n"))
      }
    }
  } else if(what == "ps"){
    for(i in 1:length(object$ve_fit$t0)){
      boot_t0 <- object$boot_est[[paste0("t0_", object$ve_fit$t0[i])]]
      cat(paste0("------- Results for t0 = ", object$ve_fit$t0[i], " -------\n"))
      for(p in c("p_immune", "p_doomed", "p_alwaysinf", "p_converted", "p_helped", "p_helpedplus")){
        cat(paste0(p, ": ", round(object$ve_fit[[p]][i], 4), 
                   " (", round(boot_t0[[p]]$lower,4), ", ", 
                   round(boot_t0[[p]]$upper,4), ") \n"))
      }
    }
  }
  
}
