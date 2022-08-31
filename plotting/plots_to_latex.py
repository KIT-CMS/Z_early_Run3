import os

#function that writes code to the end of control-plots-slides.tex to make an automatic slideshow
def plots_to_latex(tag, era, postfix, channel):
    directory = "plots/%s/%s_plots_%s/%s" % (tag, era, postfix, channel)
    plotsPath = "plots/earlyRun3" + channel + "/Run2018_plots_fully_classic/"  + channel + "/"
    latex_slides_file = open("presentation/control-plots-slides.tex", "a+")
    list_of_plots = []
    uniqueVariables = list()

    for plot in os.listdir(directory):
        list_of_plots.append(plot)

    #comparing every element in the list to every other element to group plots together by biggest similarities
    for iFirstPlot in range(0, len(list_of_plots)):
        for iSecondPlot in range(iFirstPlot, len(list_of_plots)):
            fileType = '.png'
            #plot names as a string which will be put into the latex file
            firstPlotName = str(list_of_plots[iFirstPlot])
            secondPlotName = str(list_of_plots[iSecondPlot])

            #groups which plots are on each slide for mm/ee
            if channel == "mm" or channel == "ee":
                baseFirstPlot = firstPlotName.replace("_1", "_2")
                baseSecondPlot = secondPlotName.replace("_1", "_2")

                if "_1" not in firstPlotName and "_2" not in firstPlotName and firstPlotName == secondPlotName: 
                    uniqueVariables.append(firstPlotName)

            #groups which plots are on each slide for mmet/emet
            else:
                baseFirstPlot = firstPlotName.replace("_pos", "_neg")
                baseSecondPlot = secondPlotName.replace("_pos", "_neg") 

                if "_pos" not in firstPlotName and "_neg" not in firstPlotName and firstPlotName == secondPlotName: 
                    uniqueVariables.append(firstPlotName)

            sameBase = True if baseFirstPlot == baseSecondPlot else False
            correctFileType = True if fileType in baseFirstPlot and fileType in baseSecondPlot else False

            if sameBase and correctFileType and firstPlotName != secondPlotName:
                latex_slides_file.writelines(
                    ["\\begin{frame}{\insertsubsection}\n",
                    "   \\begin{columns}\n",
                    "   \column{.5\\textwidth}\n",
                    "   \includegraphics[width=0.9\linewidth]{" + plotsPath + secondPlotName + "}\n",
                    "   \column{.5\\textwidth}\n",
                    "   \includegraphics[width=0.9\linewidth]{" + plotsPath + firstPlotName + "}\n",
                    "\end{columns}\n",
                    "\end{frame}\n",
                    "\n"]
                    )
    for plot in uniqueVariables:
        if fileType in plot:
            latex_slides_file.writelines(
                            ["\\begin{frame}{\insertsubsection}\n",
                            "   \\begin{columns}\n",
                            "   \column{.5\\textwidth}\n",
                            "   \includegraphics[width=0.9\linewidth]{" + plotsPath + plot + "}\n",
                            "\end{columns}\n",
                            "\end{frame}\n",
                            "\n"]
                            )
    latex_slides_file.writelines(
       ["\\appendix\n",
        "\\beginbackup\n",
        "\\begin{frame}\n",
        "\centering\n",
        "\\vspace{1cm}\n",
        "\LARGE{\\textbf{Backup}}\n",
        "\end{frame}\n",
        "\n",
        "\\backupend",
        "\end{document}"
        ]
    )
    latex_slides_file.close()
