inline_ui <- function(tag) {
  div(style = "display: inline-block", tag)
}

fluidPage(
  useShinyjs(),
  titlePanel(paste0("coveffectsplot: ",utils::packageVersion("coveffectsplot"))),
  fluidRow(
    column(
      2,
      tabsetPanel(
        tabPanel(
          "Inputs",
          br(),
          tags$div(
            tags$strong("Choose csv file to upload"),
            "or", actionLink("sample_data_btn", "use sample data")
          ),
          fileInput("datafile", NULL,
                    multiple = FALSE, accept = c("csv")),
          shinyjs::hidden(
            selectizeInput(
              'exposurevariables',
              label = "Parameter(s)",
              choices = c(),
              multiple = TRUE,
              options = list(plugins = list('remove_button', 'drag_drop')),
              width = '800px'
            ),
            checkboxInput('shapebyparamname', 'Change Symbol by Parameter(s) ?', value = TRUE),
            checkboxInput('colourbyparamname', 'Change Color by Parameter(s) ?', value = FALSE),
            sliderInput("vdodgeheight", "Vertical Space Between Parameters(s)",
                        min=0.5, max=2, value=0.8,width = '800px'),
            selectizeInput(
              "covariates",
              "Covariates Top to Bottom (Remove/Drag and Drop to Desired Order):",
              choices = c(),
              multiple = TRUE,
              options = list(
                placeholder = 'Please select one or more variables',
                plugins = list('remove_button', 'drag_drop')
              ),
              width = '800px'
            ),
            selectizeInput(
              'covvalueorder',
              label = paste("Drag and Drop to Desired Order within facets", "values"),
              choices = c(),
              multiple = TRUE,
              options = list(plugins = list('remove_button', 'drag_drop')),
              width = '800px'
            )
          )
        ), # tabPanel
        tabPanel("Facets",
                 selectInput(  "facetformula", "Facet Formula:",
                               choices = c("covname ~ .",
                                           "covname ~ paramname"),
                               selected = c("covname ~ ."),
                               multiple = FALSE),
                 
                 tabsetPanel(
                   tabPanel("Size/Angle/Face",
                 sliderInput("facettexty", "Facet Text Size Y",
                             min = 0, max = 32, step = 1, value = 22),
                 sliderInput("facettextx", "Facet Text Size X",
                             min = 0, max = 32, step = 1, value = 22),
                 sliderInput("facettextyangle", "Facet Text Angle Y",
                             min = 0, max = 180+90, step = 90, value = 0),
                 sliderInput("facettextxangle", "Facet Text Angle X",
                             min = 0, max = 90, step = 90, value = 0),
                 checkboxInput('boldfacettext', "Bold Facet Text", value = TRUE)
                   ),
                 tabPanel("Justification",
                 sliderInput("x_facet_text_vjust", "Facet Text Vertical Justification X",
                                      min = 0, max = 1, step = 0.1, value = 0.5),
                 sliderInput("y_facet_text_vjust", "Facet Text Vertical Justification Y",
                                      min = 0, max = 1, step = 0.1, value = 0.5),
                 sliderInput("x_facet_text_hjust", "Facet Text Horizontal Justification X",
                                      min = 0, max = 1, step = 0.1, value = 0.5),
                 sliderInput("y_facet_text_hjust", "Facet Text Horizontal Justification y",
                                      min = 0, max = 1, step = 0.1, value = 0.5)
                 ),
                 tabPanel("Facet Labels",
                          selectInput('facetlabeller' ,'Facet Label:',c(
                            "Variable(s) Name(s) and Value(s)" ="label_both",
                            "Value(s)"="label_value",
                            "Parsed Expression" ="label_parsed",
                            "Depends on Context" ="label_context",
                            "Wrap lines" ="label_wrap_gen"),
                            selected="label_value"),
                          conditionalPanel(
                            condition = "input.facetlabeller== 'label_wrap_gen'  " ,
                            sliderInput("labelwrapwidth", "N Characters to Wrap Labels:",
                                        min=1, max=100, value=c(25),step=1)
                          ),
                          checkboxInput('facetwrapmultiline',
                                        'Strip labels on multiple lines?',
                                        value = TRUE)
                 )
                 ),
                 selectizeInput(  "stripplacement", "Strip Placement:",
                                  choices = c("inside","outside"),
                                  selected = c("outside"),
                                  options = list(  maxItems = 1 )  ),
                 selectInput(  "facetswitch", "Facet Switch to Near Axis:",
                               choices = c("both","y","x","none"),
                               selected = c("both"),
                               multiple = FALSE),
                 selectInput(  "facetscales", "Facet Scales:",
                               choices = c("free_y","fixed","free_x","free"),
                               selected = c("free_y"),
                               multiple = FALSE),
                 selectInput('facetspace' ,'Facet Spaces:',
                             c("fixed","free_x","free_y","free") ,
                             selected = c("free_y"))
        ),
        tabPanel(
          "X/Y Axes",
          sliderInput("ylablesize", "Y axis labels size", min=1, max=32, value=24,step=0.5),
          checkboxInput('breakylabel', "Break long Y axis labels?", value = FALSE),
          conditionalPanel(
            condition = "input.breakylabel" ,
            sliderInput("ylabeltextwidth", "Y axis Width to break at:", min=1, max=40, value=25, step=1)),
          sliderInput("xlablesize", "X axis labels size", min=1, max=32, value=24,step=0.5),
          checkboxInput('showyaxisgridlines', "Keep Y axis Gridlines", value = TRUE),
          checkboxInput('showxaxisgridlines', "Keep X axis Gridlines", value = TRUE),
          checkboxInput('customxticks', 'Custom X axis Ticks ?', value = FALSE),
          conditionalPanel(
            condition = "input.customxticks" ,
            textInput("xaxisbreaks",label ="X axis major Breaks",
                      value = as.character(paste(
                        0,0.25,0.5,0.8,1,1.25,1.5,2
                        ,sep=",") )
            ),
            checkboxInput('customxlabels', 'Custom X axis Labels ?', value = FALSE),
            conditionalPanel(
              condition = "input.customxticks & input.customxlabels" ,
              textInput("xaxislabels",label ="X axis major Labels",
                        value = as.character(paste(
                          "0","1/4","1/2","0.8","1","1.25","1.5","2"
                          ,sep=",") )
              )
            ),
            textInput("xaxisminorbreaks",label ="X axis minor Breaks",
                      value = as.character(paste(
                        0.75,1.333
                        ,sep=",") )
            ),
            hr()
          ),
          checkboxInput('userxzoom', 'Custom X axis Range ?', value = FALSE),
          conditionalPanel(
            condition = "input.userxzoom" ,
            numericInput("lowerxin",label = "Lower X Limit",value = 0.01,min=NA,max=NA,width='100%'),
            numericInput("upperxin",label = "Upper X Limit",value = 2,min=NA,max=NA,width='100%')
          ),
          checkboxInput('logxscale', 'Log-scale X axis ?', value = FALSE),
          textInput("yaxistitle", label = "Y axis Title", value = ""),
          conditionalPanel(
            condition = "input.yaxistitle!=''" ,
          checkboxInput('parseyaxistitle', 'Parse Y axis Title?', value = FALSE)),
          textInput("xaxistitle", label = "X axis Title", value = ""),
          conditionalPanel(
            condition = "input.xaxistitle!=''" ,
          checkboxInput('parsexaxistitle', 'Parse X axis Title?', value = FALSE))
        ),
        tabPanel(
          "How To",
          hr(),
          includeMarkdown(file.path("text", "howto.md"))
        ) # tabpanel
      ) # tabsetPanel
    ), # column3

    column(
      8,
      plotOutput('plot', height = "auto", width = "100%"),
      shinyjs::hidden(
        actionButton("get_code", "Show Code", icon = icon("code")), br(), br()
      )
    ), # column6

    column(
      2,
      tabsetPanel(
        tabPanel(
          "Table Options",
              numericInput("sigdigits",label = "Significant Digits",value = 2,min=0,max=NA),
              sliderInput("tabletextsize", "Table Text Size", min=1, max=12,step=0.5, value=7),
          checkboxInput('tabletextcoloroverwrite',
                        'Overwrite Table Text Color?', value = FALSE),    
          conditionalPanel(condition = "input.tabletextcoloroverwrite" ,
                           colourpicker::colourInput("tabletextcolour",
                                    "Table Text Color:",
                                    value= "black",
                                    showColour = "both",allowTransparent=TRUE,
                                    returnName = TRUE)),
              sliderInput("plottotableratio", "Plot to Table Ratio", min=1, max=5,
                          value=4,step=0.25,
                          animate = FALSE),
              selectInput('tableposition','Table Position:',
                          c("on the right" = "right", "below" = "below", "none" = "none") ),
              selectInput(  "showtablefacetstrips", "Show Table Facet Strip on:",
                            choices = c("none","both","y","x"),
                            selected = c("none"),
                            multiple = FALSE),
              selectInput(  "tablefacetswitch", "Table Facet Switch to Near Axis:",
                            choices = c("both","y","x","none"),
                            selected = c("both"),
                            multiple = FALSE),
              checkboxInput('showtableyaxisticklabel',
                            'Show Table y axis ticks/labels ?', value = FALSE),
              checkboxInput('reservetablexaxislabelspace', 'Reserve Table x axis space ?',
                            value = FALSE),
              checkboxInput('tablepanelborder',
                            'Draw Table Panel Borders ?',
                            value = TRUE)
        ),#tabpanel
        tabPanel(
          "Reference Options",
          checkboxInput('showrefvalue', 'Show Reference Line?', value = TRUE),
          conditionalPanel(condition = "input.showrefvalue" ,
          numericInput("refvalue","Reference Line",value = 1,step = 0.1)),
          checkboxInput('showrefarea', 'Show Reference Area?', value = TRUE),
          conditionalPanel(condition = "input.showrefarea" ,
                           uiOutput("refarea")),
          
          colourpicker::colourInput("fillrefarea",
                                    "Reference Area Fill:",
                                    value= "#BEBEBE50",
                                    showColour = "both",allowTransparent=TRUE,
                                    returnName = TRUE),
          div( actionButton("fillrefareareset", "Reset Reference Area Fill"),
               style="text-align: right"),
          
          colourpicker::colourInput("colorrefvalue",
                                    "Reference Line Color:",
                                    value= "black",
                                    showColour = "both",allowTransparent=TRUE,
                                    returnName = TRUE),
          div( actionButton("colorrefvaluereset", "Reset Reference Line Color"),
               style="text-align: right"),
          sliderInput("sizerefvalue", "Reference Line Size:",
                      min = 0, max = 10, step = 0.1, value = 1),
          selectInput('linetyperefvalue',
                      'Reference Line Linetype:',c("solid","dashed",
                                     "dotted", "dotdash",
                                     "longdash", "twodash","blank"),selected = "dashed")
          ),#tabpanel
        
        tabPanel(
          "Colour/Legend Options/Theme",
          colourpicker::colourInput("stripbackgroundfill",
                                    "Strip Background Fill:",
                                    value="#E5E5E5",
                                    showColour = "both",allowTransparent=TRUE, returnName = TRUE),
          div( actionButton("stripbackfillreset", "Reset Strip Background Fill"),
               style="text-align: right"),
          colourpicker::colourInput("y_facet_text_col",
                                    "Facet Text Color Y:",
                                    value="black",
                                    showColour = "both",allowTransparent=TRUE, returnName = TRUE),
          colourpicker::colourInput("x_facet_text_col",
                                    "Facet Text Color X:",
                                    value="black",
                                    showColour = "both",allowTransparent=TRUE, returnName = TRUE),
          checkboxInput('removestrip', "Show Strip Background",value = TRUE),
          conditionalPanel(condition = "!input.colourbyparamname",
          colourpicker::colourInput("colourpointrange",
                                    "Point Range Colour:",
                                    value="blue",
                                    showColour = "both",allowTransparent=TRUE, returnName = TRUE),
          div( actionButton("colourpointrangereset", "Reset Point Range Colour"),
               style="text-align: right")
          ),
          conditionalPanel(condition = "!input.colourbyparamname",
                           colourpicker::colourInput("colourbsvrange",
                                                     "BSV Range Colour:",
                                                     value="red",
                                                     showColour = "both",allowTransparent=TRUE, returnName = TRUE),
                           div( actionButton("colourbsvrangereset", "Reset BSV Range Colour"),
                                style="text-align: right")
          ),
          conditionalPanel(condition = "!input.shapebyparamname", 
                           selectInput(inputId = "shapepointrange",
                                       label = "Point Range Shape:",
                                       choices = c("circle small", "triangle" ,"square","plus","square cross","asterisk",
                                                   "square open","cross","diamond open","triangle down open","square cross",
                                                   "diamond plus","circle plus","star","star","square plus","circle cross",
                                                   "square triangle","diamond","circle","bullet",
                                                   "circle filled","square filled","diamond filled","triangle filled",
                                                   "triangle down filled"),
                                       selected = "circle small"
                           ) 
          ),
          conditionalPanel(condition = "!input.shapebyparamname", 
                           selectInput(inputId = "shapebsvrange",
                                       label = "BSV Range Shape:",
                                       choices = c("circle small", "triangle" ,"square","plus","square cross","asterisk",
                                                   "square open","cross","diamond open","triangle down open","square cross",
                                                   "diamond plus","circle plus","star","star","square plus","circle cross",
                                                   "square triangle","diamond","circle","bullet",
                                                   "circle filled","square filled","diamond filled","triangle filled",
                                                   "triangle down filled"),
                                       selected = "circle small"
                           ) 
          ),
          conditionalPanel(condition = "input.colourbyparamname",uiOutput('userdefinedcolorui')
          ),
          conditionalPanel(condition = "input.shapebyparamname",uiOutput('userdefinedshapeui')
          ),
          sliderInput("sizepointrange", "Point range size",
                      min = 0, max = 10, step = 0.1, value = 1),
          sliderInput("fattenpointrange", "Point range fatten",
                      min = 0, max = 10, step = 0.1, value = 4),
          sliderInput("linewidthpointrange", "Point range linewidth",
                      min = 0, max = 10, step = 0.1, value = 1),
          sliderInput("base_size", "Base size for the theme",
                      min = 1, max = 30, step = 0.1, value = 22),
          sliderInput("height", "Plot Height", min=1080/4, max=1080,
                      value=900, animate = FALSE),
          checkboxInput('theme_benrich', "Apply Ben's Theme",value = FALSE),
          conditionalPanel(
            condition = "input.theme_benrich",
          textInput("custom_table_title", label ="Table Title",
                    value=""),
          sliderInput("table_title_size", "Size for Table Title",
                      min = 1, max = 30, step = 0.1, value = 15)
          ) ,
          selectizeInput(
            'legendposition',
            label = "Legend Position",
            choices = c("top","bottom","right","none"),
            selected = c("top"),
            multiple=FALSE)
        ),#tabpanel
        tabPanel(
          "Custom Legend Ordering/Spacing",
          inline_ui(
              numericInput("ncolinterval",label = "Number of columns for the Interval legend",
                           value = 1,min=NA,max=NA,width='120px')),
          inline_ui(
              numericInput("ncolshape",label = "Number of columns for the shape legend",
                           value = 1,min=NA,max=NA,width='120px')),
          
          selectizeInput(
            'legendordering',
            label = paste("Drag/Drop to reorder","Colour, Ref, Area Legends"),
            choices = c("pointinterval","ref","area","shape"),
            selected = c("pointinterval","ref","area","shape"),
            multiple=TRUE,  options = list(
              plugins = list('drag_drop')
            )),
          checkboxInput('legendshapereverse',
                        'Reverse the order of shape legend items ?',value = TRUE),
          checkboxInput('legendcoloreverse',
                        'Reverse the order of colour legend items ?',value = FALSE),
          sliderInput("legendspacex", "Multiplier for Space between Legends",
                      min = 0, max = 1.5, step = 0.1, value = 1),
          numericInput("panelspacing",label = "Strip Panel Spacing",
                       value = 5.5,min=0,step=0.1,
                       max=20,width='100%'),
          inline_ui(
              numericInput("margintop",label = "Plot Top Margin",
                           value = 0,min=0,max=NA,width='80px')),
          inline_ui(
              numericInput("tabletop",label = "Table Top Margin",
                           value = 0,min=0,max=NA,width='80px')),
          inline_ui(
              numericInput("legendtop",label = "Legend Top Margin",
                           value = 0,min=0,max=NA,width='80px')),
          inline_ui(
              numericInput("marginleft",label = "Plot Left Margin",
                           value = 5.5,min=0,max=NA,width='80px')),
          inline_ui(
              numericInput("tableleft",label = "Table Left Margin",
                           value = 5.5,min=0,max=NA,width='80px')),
          inline_ui(
              numericInput("legendleft",label = "Legend Left Margin",
                           value = 5.5,min=0,max=NA,width='80px')),
          inline_ui(
              numericInput("marginright",label = "Plot Right Margin",
                           value = 5.5,min=0,max=NA,width='80px')),
          inline_ui(
              numericInput("tableright",label = "Table Right Margin",
                           value = 5.5,min=0,max=NA,width='80px')),
          inline_ui(
             numericInput("legendright",label = "Legend Right Margin",
                          value = 5.5,min=0,max=NA,width='80px')),
          inline_ui(
            numericInput("marginbottom",label = "Plot Bottom Margin",
                         value = 0,min=0,max=NA,width='80px')),
          inline_ui(
            numericInput("tablebottom",label = "Table Bottom Margin",
                         value = 0,min=0,max=NA,width='80px')),
          inline_ui(
            numericInput("legendbottom",label = "Legend Bottom Margin",
                         value = 0,min=0,max=NA,width='80px'))
          
        ),#tabpanel
        tabPanel(
          "Custom Legend Text",
          textInput("customplottitle", label ="Plot Title",value = ""),
          textInput("customcolourtitle", label ="Pointinterval Legend text",
                    value="Median (points)\\n95% CI (horizontal lines)"),
          textInput("customcolourtitletext", label ="Pointinterval Legend title",
                    value=""),
          textInput("customshapetitletext", label ="Shape Legend title",
                    value=""),
          textInput("custombsvtitle", label ="BSV Legend text",
                    value="BSV (points)\\nPrediction Intervals (horizontal lines)"),
          textInput("customlinetypetitle", label ="Ref Legend text",
                    value="Reference (vertical line)\\nClinically relevant limits (colored area)"),
          textInput("customfilltitle", label ="Area Legend text",
                    value="Reference (vertical line)\\nClinically relevant limits (colored area)"),
          checkboxInput('combineareareflegend',
                        'Combine Ref and Area Legends if they share the same text ?',value = TRUE),
          checkboxInput('combineintervalshapelegend',
                        'Combine Interval and Shape Legends?',value = TRUE),
          sliderInput("legendtitlesize", "Text size of legend(s) title(s) ",
                      min=1, max=32, value= 16, step=0.5)
          
        )#tabpanel
      )  # tabsetpanel
    ) # closes the column 3
  )# fluidrow
)#fluidpage
