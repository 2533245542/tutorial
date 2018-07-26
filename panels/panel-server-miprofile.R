#output UI for covariates and confoundings

output$miprofile_out_cova <- renderUI({
  selectInput("miprofile_cova", "SELECT COVA", COVA(), selected = "NULL", multiple = TRUE)
})
output$miprofile_out_conf <- renderUI({
  selectInput("miprofile_conf", "SELECT CONF", COVA(), selected = "NULL",multiple = TRUE)
})

#render function for generating .txt
output$test <- renderText({
  miprofile_txt()
  c("")
})

#actual function for generating .txt
miprofile_txt = reactive({
  if(input$action_miprofile < 1){
    return(NULL)
  }
  isolate({
    if(input$action_miprofile != 0){
      fileinfo = file('information.txt')
      cat(file = fileinfo)
      write("## === specify the variable names for the covariate in meta file === ##", fileinfo, append = TRUE)
      write("COVA = Normal", 'information.txt', append = TRUE)
      write("STRATA = triplet", 'information.txt', append = TRUE)

      
      write("\n## === request distances === ##", 'information.txt', append = TRUE)
      out = input$miprofile_outfiles
      vec = c()
      if(input$miprofile_bcd){
        write("DIST = Bray-Curtis", 'information.txt', append = TRUE)
        disnorm = paste("DIST_NORM",input$miprofile_bcd_norm, sep = " = " )
        write(disnorm, 'information.txt', append = TRUE)
        #output dist if true
        if(out){
          write("DIST_OFILE = Dbc.dist", 'information.txt', append = TRUE)
        }
        write("", 'information.txt', append = TRUE)
        
      }
      if(input$miprofile_jcd){
        write("DIST = Jaccard", 'information.txt', append = TRUE)
        disnorm = paste("DIST_NORM",input$miprofile_jcd_norm, sep = " = " )
        write(disnorm, 'information.txt', append = TRUE)
        
        # vec = c(vec, "DIST", "Jaccard Distance", "DIST_NORM", input$miprofile_jcd_norm)
        if(out){
          write("DIST_OFILE = Dj.dist", 'information.txt', append = TRUE)
        }
        write("", 'information.txt', append = TRUE)
        
      }	
      if(input$miprofile_uwu){
        write("DIST = uwUniFrac", 'information.txt', append = TRUE)
        disnorm = paste("DIST_NORM",input$miprofile_uwu_norm, sep = " = " )
        write(disnorm, 'information.txt', append = TRUE)
        
        # vec = c(vec, "DIST", "UnWeighted UniFrac", "DIST_NORM", input$miprofile_uwu_norm)
        if(out){
          write("DIST_OFILE = Duw.dist", 'information.txt', append = TRUE)
        }
        write("", 'information.txt', append = TRUE)
        
      }
      if(input$miprofile_gu){
        write("DIST = wUniFrac", 'information.txt', append = TRUE)
        disnorm = paste("DIST_NORM",input$miprofile_gu_norm, sep = " = " )
        write(disnorm, 'information.txt', append = TRUE)
        
        # vec = c(vec, "DIST", paste("Generalized UniFrac", output$miprofile_gu_alpha), "DIST_NORM", input$miprofile_gu_norm)
        if(out){
          write("DIST_OFILE = Dw.dist", 'information.txt', append = TRUE)
        }
        write("", 'information.txt', append = TRUE)
        
      }	
      if(input$miprofile_pwu){
        write("DIST = pwUniFrac", 'information.txt', append = TRUE)
        disnorm = paste("DIST_NORM",input$miprofile_pwu_alpha, sep = " = " )
        write(disnorm, 'information.txt', append = TRUE)
        
        # vec = c(vec, "DIST", paste("Presence Weighted UniFrac", output$miprofile_pwu_alpha), "DIST_NORM", input$miprofile_pwu_norm)
        if(out){
          write("DIST_OFILE = Dpw.dist", 'information.txt', append = TRUE)
          # vec = c(vec, "DIST_OFILE", "Duw.dist")
        }
        write("", 'information.txt', append = TRUE)
        
      }
      # print(vec)
      # x = matrix( vec, ncol = 2, nrow = length(vec)/2, byrow = TRUE)
      
      # write.table(x,file = "information.txt", sep = " = ", row.names = FALSE, col.names = FALSE, quote = FALSE)
      close(fileinfo)
      print("end generating files")
      }
    
  })
})
# 
# test = reactive({
#   if(input$action_miprofile != 0){
# 	print("start generating files")
# 	out = input$miprofile_outfiles
# 	vec = c()
# 	if(input$miprofile_bcd){
# 		vec = c(vec, "DIST", "Bray-Curtis Distance", "DIST_NORM", input$miprofile_bcd_norm)
# 		if(out){
# 			vec = c(vec, "DIST_OFILE", "Dbc.dist")
# 		}
# 	}
# 	if(input$miprofile_jcd){
# 		vec = c(vec, "DIST", "Jaccard Distance", "DIST_NORM", input$miprofile_jcd_norm)
# 		if(out){
# 			vec = c(vec, "DIST_OFILE", "Dj.dist")
# 		}
# 	}	
# 	if(input$miprofile_uwu){
# 		vec = c(vec, "DIST", "UnWeighted UniFrac", "DIST_NORM", input$miprofile_uwu_norm)
# 		if(out){
# 			vec = c(vec, "DIST_OFILE", "Duw.dist")
# 		}
# 	}
# 	if(input$miprofile_gu){
# 		vec = c(vec, "DIST", paste("Generalized UniFrac", output$miprofile_gu_alpha), "DIST_NORM", input$miprofile_gu_norm)
# 		if(out){
# 			vec = c(vec, "DIST_OFILE", "Dgu.dist")
# 		}
# 	}	
# 	if(input$miprofile_pwu){
# 		vec = c(vec, "DIST", paste("Presence Weighted UniFrac", output$miprofile_pwu_alpha), "DIST_NORM", input$miprofile_pwu_norm)
# 		if(out){
# 			vec = c(vec, "DIST_OFILE", "Duw.dist")
# 		}
# 	}
# 	print(vec)
# 	x = matrix( vec, ncol = 2, nrow = length(vec)/2, byrow = TRUE)
# 
# 	write.table(x,file = "information.txt", sep = " = ", row.names = FALSE, col.names = FALSE, quote = FALSE)
# 	print("end generating files")
# 
# }})
