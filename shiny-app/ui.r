
shinyUI(pageWithSidebar(
	
	headerPanel("GLEE: differential protein expression test"),
	
	sidebarPanel(
		checkboxInput("show_format_help", label="Show data file formatting instructions", value=FALSE),
		checkboxInput("show_settings_help", label="Show explanation of input/output settings", value=TRUE),
		checkboxInput("show_faq_answers", label="Show FAQ and answers", value=FALSE),
		br(),
		strong("Input settings"),
		inputPanel(
			numericInputMini("num_A", label="# reps(A)", value=3, min=1, max=10, step=1), br(),
			numericInputMini("num_B", label="# reps(B)", value=3, min=1, max=10, step=1), br(),
			numericInputMini("num_iter", label="# iterations", value=1000, min=1000, max=50000, step=1000), br(),
			radioButtons("fit_type", label="fit type", choices=c("linear","cubic"), selected="cubic", inline=TRUE)
		),
		br(),
		strong("Output settings"),
		inputPanel(
			numericInputMini("num_show", label="# display", value=10, min=1, max=100, step=1), br(),
			numericInputMini("num_digits", label="# digits", value=4, min=1, max=10, step=1), br(),
			uiOutput("output_settings")
		)
	),
	
	mainPanel(
		HTML('<strong>Note</strong>'), br(),
		HTML('<ul>
			<li>The analysis begins once the data file is uploaded, so ensure that the input settings are correct before uploading</li>
			<li>The output settings may be modified after the analysis has finished running</li>
			<li>To reset everything and start over, hit the refresh button on your web browser</li>
			</ul>'), br(),
		conditionalPanel(condition="input.show_format_help",
			uiOutput("help_format"), br()
		),
		conditionalPanel(condition="input.show_settings_help",
			uiOutput("help_settings"), br()
		),
		conditionalPanel(condition="input.show_faq_answers",
			uiOutput("faq_answers"), br()
		),
		fileInput("data_file", label="Upload data (text file)", multiple=FALSE), br(),
		conditionalPanel(condition="input.show_model_fit_plots",
			# plotOutput("model_fit_plots", width=plotDispW, height=plotDispHL), br()
			plotOutput("model_fit_plots"), br()
		),
		conditionalPanel(condition="input.show_stn_pval_plots",
			# plotOutput("stn_pval_plots", width=plotDispW, height=plotDispH), br()
			plotOutput("stn_pval_plots"), br()
		),
		conditionalPanel(condition="input.get_diff_exp_table",
			uiOutput("download_button"), br()
		),
		br(),
		uiOutput("diff_exp_title"), br(),
		tableOutput("diff_exp_table"), br()
	)
	
))
