library(shinydashboard)
library(shiny)
library(DT)


ui <- dashboardPage(
  dashboardHeader(title = "ORBIT UI"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("MODEL", tabName = "MODEL", icon = icon("dashboard")),
      menuItem("BURN ITME CALCULATOR", tabName = "BT", icon = icon("th")),
      menuItem("DISCUSSION", tabName = "discussion", icon = icon("th"))        
      
    )#ENDOF sidebarMenu
  ),#ENDOF dashboardSidebar()
  
  
  dashboardBody(
    
    tabItems(
      
      # TAB1
      tabItem( tabName = "MODEL",
               
               fluidRow(
                 box(plotOutput("plotMO", width = "100%",height = "400px")),#ENDOF box plot
                 box(plotOutput("plotMOf", width = "100%",height = "400px")),#ENDOF box plot
                
               ),#ENDF fluid ROW 1 DONT FORGET the COMMA!
               
               
               fluidRow(
                 
                 box(plotOutput("plotRdot0", width = "100%",height = "400px")),#ENDOF box plot
                 box(plotOutput("plotRdotF", width = "100%",height = "400px"))#ENDOF box plot
               ),#ENDF fluid ROW 2 DONT FORGET the COMMA!
               
               
               #fluidPage(box(DT::dataTableOutput("TELECOMM"),width = "200%",height = "400px"))
               
               
               
      ),#tabItem1
      
      
      
      #TAB2
      tabItem( tabName = "BT",
               fluidRow(
                 
                # box(plotOutput("plot2", width = "100%",height = "400px")),
                # box(plotOutput("plot3", width = "100%",height = "400px"))
                 
               ),#ENDOF row1 
               
               fluidRow(
                 
                box(plotOutput("PLOT_DVslider", width = "100%",height = "400px")),#ENDOF box plot1
                box(title = "Velocity",sliderInput("DVslider", "Number of observations:", V_Po, 1500, 50))#ENDOF box plot2
               )#ENDF fluid ROW 2
               
      ),#ENDOF tabItem2
      
      
      
      
      
      
      ##tab3
      tabItem( tabName = "discussion", 
               fluidRow(
                # box(textOutput("discussion")),
                # box(plotOutput("plot4")),box(plotOutput("plot5")),box(plotOutput("plot6")),
                # box(plotOutput("plot7")),box(plotOutput("plot8")),box(plotOutput("plot9"))
               )#ENDOF row
               
               
               
      )###ENDOF tab3
      
      ########## 
    )#ENDOF tabItems
  )#ENDOF body
)#ENDOF ui






server <- function(input, output) {
  #TAB1 PLOTS
  output$plotMO<- renderPlot({PLOT_MODEL_ORBIT()})#ENDOF render plotMO
  output$plotMOf <- renderPlot({PLOT_FINAL_ORBIT_MODEL()})#ENDOF render plotOP
  output$plotRdot0<- renderPlot({PLOT_RDOT_0()})#ENDOF render plotrd0
  output$plotRdotF<- renderPlot({PLOT_RDOT_F()})#ENDOF render plotrdF
  
  
  output$PLOT_DVslider <- renderPlot({PLOT_ELLIPSE_DVslider(input$DVslider)})#ENDOF render plot1
  #output$TELECOMM <-DT::renderDataTable(DT::datatable({data <- TelCo_Customer_Data}))#ENDOF table output
  #output$discussion <- renderText({DISC})
}#ENDOF server()

shinyApp(ui, server)

