
shinyServer(function(input,output) {

	read_numA = reactive({
		y = tryCatch({ as.numeric(input$num_A) }, error = function(e) { NA }, warning = function(w) { NA })
		if (!is.finite(y)) { return(NULL) }
		return(floor(abs(y)))
	})
	
	read_numB = reactive({
		y = tryCatch({ as.numeric(input$num_B) }, error = function(e) { NA }, warning = function(w) { NA })
		if (!is.finite(y)) { return(NULL) }
		return(floor(abs(y)))
	})
	
	read_numIter = reactive({
		y = tryCatch({ as.numeric(input$num_iter) }, error = function(e) { NA }, warning = function(w) { NA })
		if (!is.finite(y)) { return(NULL) }
		return(floor(abs(y)))
	})
	
	read_dataFile = reactive({
		if (is.null(input$data_file)) { return(NULL) }
		dataFile = input$data_file
		
		D = tryCatch({ read.table(dataFile$datapath, header=TRUE, sep="\t", quote="") }, error = function(e) { NA }, warning = function(w) { NA })
		if (any(is.na(D))) { return(NULL) }
		
		# ensure that the input file has as many columns as claimed, and at least one row
		if (is.null(read_numA()) || is.null(read_numB())) { return(NULL) }
		if (dim(D)[2] != (1 + read_numA() + read_numB())) { return(NULL) }
		if (dim(D)[1]<2) { return(NULL) }
		
		# ensure that the spectral-count data contains only finite positive values
		dat = tryCatch({ matrix(as.numeric(unlist(D[,-1])), nrow=dim(D)[1]) }, error = function(e) { NA }, warning = function(w) { NA })
		if (!all(is.finite(dat)) || any(dat<0)) { return(NULL) }
		
		# return a list containing the replicates for each condition separately (as matrices)
		colnames(dat) = names(D)[-1]
		Prot = as.character(D[,1])
		A = dat[,1:read_numA()]
		B = dat[,read_numA()+(1:read_numB())]
		return(list(Prot=Prot,A=A,B=B))
	})
	
	fit_model = reactive({
		withProgress(message="reading data from file ...", {
			if (is.null(read_dataFile())) { return(NULL) }
			if (class(read_dataFile()) != "list") { return(HTML("error reading the input data file: ensure that the spreadsheet has the specified number of colums, and that the data contains only positive finite values for at least 2 named proteins")) }
			
			A = read_dataFile()$A
			B = read_dataFile()$B
			
			conditions = c("A","B")
			
			xbar = cbind(rowMeans(A), rowMeans(B))
			stdev = cbind(apply(A, FUN=sd, MARGIN=1), apply(B, FUN=sd, MARGIN=1))
			
			setProgress(message="fitting model ...")
			m = vector("list", length(conditions))
			names(m) = conditions
			for (j in 1:length(conditions)) {
				ind_valid = is.finite(xbar[,j]) & xbar[,j]>0 & is.finite(stdev[,j]) & stdev[,j]>0
				xbar_log = log(xbar[ind_valid,j])
				stdev_log = log(stdev[ind_valid,j])
				
				if (input$fit_type=="linear") {
					model = tryCatch({ lm(stdev_log ~ xbar_log) }, error = function(e) { NA }, warning = function(w) { NA })
				} else if (input$fit_type=="cubic") {
					model = tryCatch({ lm(stdev_log ~ poly(xbar_log,3)) }, error = function(e) { NA }, warning = function(w) { NA })
				}
				if (class(model) != "lm") { stop("error fitting model: condition ", conditions[j]) }
				
				adjRsq = summary(model)$adj.r.squared
				coeffs = model$coefficients
				
				# pick 20 values to show model fitted values
				xbar_showfit = seq(min(xbar_log), max(xbar_log), length.out=20)
				stdev_showfit = predict(model, data.frame(xbar_log = xbar_showfit))
				
				# remove outliers, i.e. keep points that are within 5 "sigma" of stdev
				ind_show = (stdev_log>(mean(stdev_log)-5*sd(stdev_log))) & (stdev_log<(mean(stdev_log)+5*sd(stdev_log)))
				xbar_showpts = xbar_log[ind_show]
				stdev_showpts = stdev_log[ind_show]
				
				# collect into a list to be returned
				m[[j]] = list(model = model, adjRsq = adjRsq, coeffs = coeffs, xbar_showfit = xbar_showfit, stdev_showfit = stdev_showfit, xbar_showpts = xbar_showpts, stdev_showpts = stdev_showpts, xbar = xbar[,j], stdev = stdev[,j])
			}
			return(m)
		})
	})
	
	output$output_settings = renderUI({
		if (is.null(calc_stn_pval())) { return(NULL) }
		y = paste( as.character(checkboxInput("show_model_fit_plots", label="Show model-fit plots", value=FALSE)),
					as.character(checkboxInput("show_stn_pval_plots", label="Show STN-pValue plots", value=FALSE)),
					as.character(checkboxInput("get_diff_exp_table", label="Download results", value=FALSE)) )
		HTML(y)
	})
	
	output$model_fit_plots = renderPlot({
		if (is.null(fit_model())) { return(NULL) }
		m = fit_model()
		par(mfrow=c(2,2), mar=c(5,4,4,2)+0.1)
		for (j in 1:length(m)) {
			plot(m[[j]]$xbar_showpts, m[[j]]$stdev_showpts, type="p", pch=20, col="gray", xlab="log(mean)", ylab="log(stdev)", main=paste("Condition ",names(m)[j],sep=""))
			lines(m[[j]]$xbar_showfit, m[[j]]$stdev_showfit, lty=1, col="blue", lwd=2)
			legend("topleft", paste("adj.Rsq = ",round(m[[j]]$adjRsq,3),sep=""), bty="n")
			plot(sort(m[[j]]$xbar_showpts), type="p", pch=20, col="gray", xlab="protein", ylab="signal level (mean)", main=paste("Condition ",names(m)[j],sep=""))
		}
	})
	# , width=plotRendW, height=plotRendHL, res=plotRendRes
	
	calc_stn_pval = reactive({
		if (is.null(fit_model())) { return(NULL) }
		if (is.null(read_numIter())) { return(NULL) }
		
		m = fit_model()
		best_idx = as.numeric(which.max(sapply(m, FUN = function(x) { x$adjRsq })))
		best_model = m[[best_idx]]$model
		
		xbar = c()
		for (j in 1:length(m)) {
			xbar_j = m[[j]]$xbar
			# replace xbar zeroes with min-positive xbar from the same condition
			xbar_j[xbar_j==0] = min(xbar_j[xbar_j>0])
			xbar = cbind(xbar, xbar_j)
		}
		min_xbar_value = apply(xbar, MARGIN=2, FUN=min)
		
		# calculate model-based STN
		model_stn = (xbar[,2] - xbar[,1]) / (calcs(best_model,xbar[,1]) + calcs(best_model,xbar[,2]))
		if (!all(is.finite(model_stn))) { stop("model_stn contains non-finite values") }
		
		# calculate null distribution of model_stn using the baseline (i.e. best fit) condition
		if (best_idx==1) {
			# condition A is the baseline
			pValue = getSTNdistrib(A = read_dataFile()$A, nA = dim(read_dataFile()$A)[2], nB = dim(read_dataFile()$B)[2], iter = read_numIter(), xbar0_replace = min_xbar_value[best_idx], best_model = best_model) %>% getPvals(model_stn, read_numIter(), .)
		} else {
			# condition B is the baseline, so assign it to A (since A is used as baseline while resampling)
			pValue = getSTNdistrib(A = read_dataFile()$B, nA = dim(read_dataFile()$B)[2], nB = dim(read_dataFile()$A)[2], iter = read_numIter(), xbar0_replace = min_xbar_value[best_idx], best_model = best_model) %>% getPvals(model_stn, read_numIter(), .)
		}
		if (!all(is.finite(pValue))) { stop("pValue contains non-finite values") }
		
		return(list(model_stn = model_stn, pValue = pValue))
	})
	
	getSTNdistrib = function(A, nA, nB, iter, xbar0_replace, best_model) {
		# A: data matrix for the baseline (i.e. best fit) condition
		# nA: number of replicates in the baseline condition
		# nB: number of replicates in the other condition
		withProgress(message="generating resampled STN distribution ...", {
			XbarStar = vector("list", length=nrow(A))
			for (i in 1:nrow(A)) {
				if (i %% 100 == 0) { setProgress(detail=paste0(i," of ",nrow(A)), value=i/nrow(A)) }
				Astar_i = matrix(sample(A[i,], size=nA*iter, replace=TRUE), nrow=iter, ncol=nA)
				Bstar_i = matrix(sample(A[i,], size=nB*iter, replace=TRUE), nrow=iter, ncol=nB)
				XbarStar[[i]] = dplyr::data_frame(xbarA = rowMeans(Astar_i), xbarB = rowMeans(Bstar_i)) %>%
												dplyr::mutate(xbarA = ifelse(xbarA==0,xbar0_replace,xbarA),
																			xbarB = ifelse(xbarB==0,xbar0_replace,xbarB),
																			model_stn_dist = (xbarB - xbarA) / (calcs(best_model,xbarB) + calcs(best_model,xbarA)))
			}
			XbarStar = do.call("bind_rows", XbarStar)
			if (XbarStar %>% dplyr::filter(!is.finite(model_stn_dist)) %>% nrow() > 0) { stop("model_stn_dist contains non-finite values") }
			return(XbarStar %>% dplyr::transmute(model_stn_dist))
		})
	}
	
	getPvals = function(model_stn, iter, model_stn_dist) {
		withProgress(message="calculating p-values ...", {
			pValue = rep(NA, length(model_stn))
			F = ecdf(unlist(model_stn_dist, use.names=FALSE))
			for (i in 1:length(model_stn)) {
				p = F(model_stn[i])
				pValue[i] = 2*min(p,1-p)
			}
			if (!all(is.finite(pValue))) { stop("pValue contains non-finite values") }
			return(pValue)
		})
	}
	
	output$stn_pval_plots = renderPlot({
		if (is.null(calc_stn_pval())) { return(NULL) }
		y = calc_stn_pval()
		par(mfrow=c(1,2), mar=c(5,4,4,2)+0.1)
		plot(y$model_stn, y$pValue, type="p", pch=20, col="gray", xlab="STN", ylab="p-value", main="P-value distribution")
		plot(sort(y$model_stn), type="p", pch=20, col="gray", xlab="protein", ylab="STN", main="STN")
	})
	# , width=plotRendW, height=plotRendH, res=plotRendRes
	
	read_numShow = reactive({
		y = tryCatch({ as.numeric(input$num_show) }, error = function(e) { NA }, warning = function(w) { NA })
		if (!is.finite(y)) { return(NULL) }
		return(y)
	})
	
	read_numDigits = reactive({
		y = tryCatch({ as.numeric(input$num_digits) }, error = function(e) { NA }, warning = function(w) { NA })
		if (!is.finite(y)) { return(NULL) }
		return(y)
	})
	
	output$diff_exp_title = renderUI({
		if (is.null(calc_stn_pval())) { return(NULL) }
		y = as.character(strong("Differentially Expressed Proteins"))
		HTML(y)
	})
	
	output$diff_exp_table = renderTable({
		if (is.null(calc_stn_pval())) { return(NULL) }
		if (is.null(read_numShow())) { return(NULL) }
		if (is.null(read_numDigits())) { return(NULL) }
		
		ind = order(calc_stn_pval()$pValue)[1:min(read_numShow(),length(calc_stn_pval()$pValue))]
		return( data.frame( Name = read_dataFile()$Prot[ind], STN = as.character(round(calc_stn_pval()$model_stn[ind],read_numDigits())), pVal = as.character(round(calc_stn_pval()$pValue[ind],read_numDigits())) ) )
	}, include.rownames=FALSE)
	
	output$download_button = renderUI({
		if (is.null(calc_stn_pval())) { return(NULL) }
		downloadButton("download_table", label = "Download table")
	})

	output$download_table = downloadHandler(
		filename = function() {
			"diff-exp.txt"
		},
		content = function(file) {
			ind = order(calc_stn_pval()$pValue)
			write.table(data.frame( Name = read_dataFile()$Prot[ind], STN = as.character(round(calc_stn_pval()$model_stn[ind],read_numDigits())), pVal = as.character(round(calc_stn_pval()$pValue[ind],read_numDigits())) ), file, append=FALSE, quote=FALSE, sep="\t", eol="\n", row.names=FALSE)
		}
	)
	
	output$help_format = renderUI({
		y = paste(
			"<strong>Format</strong>",
			"<ul>",
				"<li>All the data must be in a tab-delimited text file</li>",
				"<li>Column headers MUST be specified</li>",
				"<li>The first column must contain alpha-numeric protein names (cannot accept just numbers!)</li>",
				"<li>Replicates for the first condition (condition A) must appear after the first column, each replicate in a separate column</li>",
				"<li>Then the second condition (condition B) replicates must follow in succeeding columns</li>",
				"<li>No blank columns are allowed between replicates or between conditions</li>",
				"<li>No blank rows or cells are allowed in the data (please enter 0 for undetected proteins)</li>",
			"</ul>"
		)
		HTML(y)
	})
	
	output$help_settings = renderUI({
		y = paste(
			"<strong>Input settings</strong>",
			"<ul>",
				"<li># reps(A) = this is the number of replicates in condition A, i.e. the first set of replicates you specify in the input data file</li>",
				"<li># reps(B) = this is the number of replicates in condition B, i.e. the second set of replicates</li>",
				"<li># iterations = this specifies the number of iterations for re-sampling</li>",
				"<li>fit type = this specifies the type of model that will be fit, either linear or cubic</li>",
			"</ul>",
			"<strong>Output settings</strong>",
			"<ul>",
				"<li># display = number of top (i.e. most differentially expressed) proteins to display in the results table</li>",
				"<li># digits = number of digits (i.e. precision) for the numbers in the results table</li>",
			"</ul>"
		)
		HTML(y)
	})
	
	output$faq_answers = renderUI({
		y = paste(
			"<strong>How long does the program take to complete the analysis?</strong>", "<br>",
			"It depends on the number of proteins you have in your dataset and the number of iterations you 
			specify. During our tests, we used 10,000 iterations on 2000 proteins and it finished in about 2 
			minutes. Updates will be displayed on the webpage while the analysis is in progress.", "<br>",
			"<br>",
			"<strong>How many iterations should I use for resampling? Are 1000 iterations enough?</strong>", "<br>",
			"This would depend on how accurate you want the results to be, and how much memory is available
			on the computer/server that's running the analysis. A computationally accurate, tenable STN
			distribution is obtained when the product of the number of proteins and the number of iterations 
			is less than 10 million. So if your dataset has 2000 proteins, you could use up to 5000 iterations, 
			and increasing it further would likely cause memory outage or increase runtime substantially.", "<br>",
			"<br>",
			"<strong>The program suddenly stopped / no response / weird error. What should I do?</strong>", "<br>",
			"Could be a case of memory outage. Re-start the app by refreshing your web-browser, reduce the
			number of iterations and re-run your analysis.", "<br>",
			"<br>",
			"<strong>Why don't I get the same p-values when I re-run the app with the same dataset and parameters?</strong>", "<br>",
			"Since a randomization is performed on each protein, the simulated STN distribution will differ 
			slightly, leading to slight differences in the p-values even when the app is re-run with identical
			parameters. But the general trends in the results can be expected to be the same from one run to 
			another, i.e. the proteins that have low p-value (which appear at the top of the output file) should 
			continue to have low p-values when re-run with the same parameters. One way to ensure fairly similar 
			p-values from repeated runs is to increase the number of resampling iterations.", "<br>"
		)
		HTML(y)
	})
	
})
