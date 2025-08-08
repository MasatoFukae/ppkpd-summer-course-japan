# ui.R
shinyUI(
  fluidPage(
    theme = shinytheme(theme = "cerulean"),
    
    sidebarLayout(
      sidebarPanel(
        fileInput(
          "file",
          "Choose simulation result (.csv)",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv"
          )
        ),
        pickerInput("ggtheme", "Plot theme", choices = c("default", "white panel",  "simple")),
        htmlOutput("r"),
        prettyRadioButtons("x", "day/hour", choices = c("day", "hour")),
        prettyRadioButtons("y", "linear/ylog", choices = c("linear", "ylog")),
        numericInput("xmax", "xmax", value = 70),
        numericInput("xtick", "xtick", value = 21),
        textInput("xlabel", "xlabel", value = "Time (day)"),
        textInput("ylabel", "ylabel", value = "Median (90% PI)"),
        numericInput("fontsize", "font size", value = 10)
      ),
      mainPanel(
        br(),
        htmlOutput("up"),
        br(),
        plotOutput("plt"),
        br(),
        htmlOutput("out")
      )
    )
  )
)