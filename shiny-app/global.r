library(dplyr)

numericInputMini = function(inputId, label, value, min = NA, max = NA, step = NA) {
	div(style="display:inline-block",
	tags$label(label, `for` = inputId), 
	tags$input(id = inputId, type = "number", value = value, min = min, max = max, step = step, class = "input-mini"))
}

textInputMini = function(inputId, label, value = "") {
	div(style="display:inline-block",
	tags$label(label, `for` = inputId), 
	tags$input(id = inputId, type = "text", value = value, class = "input-mini"))
}

# # plot render settings (server.r)
# plotRendW = 2500
# plotRendH = 1250
# plotRendHL = 2500
# plotRendRes = 300

# # plot display settings (ui.r)
# plotDispW = "800px"
# plotDispH = "400px"
# plotDispHL = "800px"

calcs = function(model, x) {
	exp(predict(model, data.frame(xbar_log=log(x))))
}
