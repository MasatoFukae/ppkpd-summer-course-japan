# server.R
shinyServer(function(input, output, session) {
  # reactive values
  df = reactiveValues()
  df$df_sim = data.frame(x = 1)
  df$plt = ggplot()
  output$plt = renderPlot(df$plt)
  
  observeEvent(input$file, {
    df$df_sim = read_csv(input$file$datapath)
    print(df$df_sim)
    output$r = renderUI({
      pickerInput("regimen", "Regimen", choices = unique(df$df_sim$Regimen))
    })
    output$up = renderUI({
      actionBttn(
        "update",
        label = "Update plot",
        style = "simple",
        color = "danger",
        size = "xs"
      )
    })
    output$out = renderUI({
      actionBttn(
        "save",
        label = "Save plot as .pdf",
        style = "simple",
        color = "primary",
        size = "xs"
      )
    })
    
  })
  
  observeEvent(input$regimen, {
    if(input$x == "day") {
      plt = df$df_sim %>%
        filter(Regimen == input$regimen) %>%
        mutate(name = ifelse(name == "DV_CP", "CP (mg/L)", "IL6 (pg/mL)")) %>%
        ggplot(aes(day, Median))
    } else {
      plt = df$df_sim %>%
        filter(Regimen == input$regimen) %>%
        mutate(name = ifelse(name == "DV_CP", "CP (mg/L)", "IL6 (pg/mL)")) %>%
        ggplot(aes(time, Median))
    }
    
    if(input$ggtheme == "default") {
      plt = plt +
        theme_grey()
    }
    if(input$ggtheme == "white panel") {
      plt = plt +
        theme_bw()
    }
    if(input$ggtheme == "simple") {
      plt = plt +
        theme_classic()
    }
    
    plt = plt +
      theme(
        text = element_text(size = input$fontsize),
        axis.title = element_text(size = input$fontsize + 2),
        axis.text = element_text(size = input$fontsize),
        aspect.ratio = 1/1
      ) +
      scale_x_continuous(limits = c(0, input$xmax), breaks = seq(0, input$xmax, input$xtick)) +
      labs(x = input$xlabel, y = input$ylabel, title = input$regimen) +
      geom_line() +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.5, fill="darkblue") +
      facet_wrap(~name, scales = "free_y")
    if(input$y == "ylog") {
      plt = plt +
        scale_y_log10()
    }
    
    output$plt = renderPlot(plt)
  })
  
  observeEvent(input$update, {
    if(input$x == "day") {
      plt = df$df_sim %>%
        filter(Regimen == input$regimen) %>%
        mutate(name = ifelse(name == "DV_CP", "CP (mg/L)", "IL6 (pg/mL)")) %>%
        ggplot(aes(day, Median))
    } else {
      plt = df$df_sim %>%
        filter(Regimen == input$regimen) %>%
        mutate(name = ifelse(name == "DV_CP", "CP (mg/L)", "IL6 (pg/mL)")) %>%
        ggplot(aes(time, Median))
    }
    
    if(input$ggtheme == "default") {
      plt = plt +
        theme_grey()
    }
    if(input$ggtheme == "white panel") {
      plt = plt +
        theme_bw()
    }
    if(input$ggtheme == "simple") {
      plt = plt +
        theme_classic()
    }
    
    plt = plt +
      theme(
        text = element_text(size = input$fontsize),
        axis.title = element_text(size = input$fontsize + 2),
        axis.text = element_text(size = input$fontsize),
        aspect.ratio = 1/1
      ) +
      scale_x_continuous(limits = c(0, input$xmax), breaks = seq(0, input$xmax, input$xtick)) +
      labs(x = input$xlabel, y = input$ylabel, title = input$regimen) +
      geom_line() +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.5, fill="darkblue") +
      facet_wrap(~name, scales = "free_y")
    if(input$y == "ylog") {
      plt = plt +
        scale_y_log10()
    }
    
    output$plt = renderPlot(plt)
  })
  
  observeEvent(input$save, {
    filename = str_replace(input$file$name, ".csv", "")
    pdf(paste0(filename, "_shiny.pdf"), paper = "a4r", width = 20)
    regimens = unique(df$df_sim$Regimen)
    for(i in 1:length(regimens)) {
      regimen = regimens[i]
      if(input$x == "day") {
        plt = df$df_sim %>%
          filter(Regimen == regimen) %>%
          mutate(name = ifelse(name == "DV_CP", "CP (mg/L)", "IL6 (pg/mL)")) %>%
          ggplot(aes(day, Median))
      } else {
        plt = df$df_sim %>%
          filter(Regimen == regimen) %>%
          mutate(name = ifelse(name == "DV_CP", "CP (mg/L)", "IL6 (pg/mL)")) %>%
          ggplot(aes(time, Median))
      }
      
      if(input$ggtheme == "default") {
        plt = plt +
          theme_grey()
      }
      if(input$ggtheme == "white panel") {
        plt = plt +
          theme_bw()
      }
      if(input$ggtheme == "simple") {
        plt = plt +
          theme_classic()
      }
      
      plt = plt +
        theme(
          text = element_text(size = input$fontsize),
          axis.title = element_text(size = input$fontsize + 2),
          axis.text = element_text(size = input$fontsize),
          aspect.ratio = 1/1
        ) +
        scale_x_continuous(limits = c(0, input$xmax), breaks = seq(0, input$xmax, input$xtick)) +
        labs(x = input$xlabel, y = input$ylabel, title = input$regimen) +
        geom_line() +
        geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.5, fill="darkblue") +
        facet_wrap(~name, scales = "free_y")
      if(input$y == "ylog") {
        plt = plt +
          scale_y_log10()
      }
      print(plt)
    }
    dev.off()
    
    show_alert(
      title = "message",
      text = "Plot output completed"
    )
  })
})

