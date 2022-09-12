library(shiny)
library(DT)
library(plotrix)
options(warn=-1)

# https://stackoverflow.com/questions/43635846/
BinMean <- function (vec, every, na.rm = TRUE) {
  n <- length(vec)
  x <- .colMeans(vec, every, n %/% every, na.rm)
  r <- n %% every
  if (r) x <- c(x, mean.default(vec[(n - r + 1):n], na.rm = na.rm))
  x
}

config_file=read.table("./interface_config.txt",comment.char="",sep="\t",quote="")
temp_file_dir=config_file[2,2]
username=config_file[4,2]
regions_list=read.table(paste(temp_file_dir,"/../user_data/",config_file[1,2],".pipeline_summary.txt",sep=""),comment.char="",sep="\t",quote="")

# https://stackoverflow.com/questions/63882483
mymodal <- function (..., title = NULL, footer = modalButton("Dismiss"), 
                     size = c("m", "s", "l"), easyClose = FALSE, fade = TRUE, idcss = "") 
{
  size <- match.arg(size)
  cls <- if (fade) 
    "modal fade"
  else "modal"
  div(id = "shiny-modal", class = cls, tabindex = "-1", `data-backdrop` = if (!easyClose) 
    "static", `data-keyboard` = if (!easyClose) 
      "false", div(class = paste("modal-dialog", idcss), class = switch(size, 
                                                                        s = "modal-sm", 
                                                                        m = NULL, 
                                                                        l = "modal-lg"), 
                   div(class = "modal-content", 
                       if (!is.null(title)) 
                         div(class = "modal-header", tags$h4(class = "modal-title", 
                                                             title)
                         ), 
                       div(class = "modal-body", ...), 
                       if (!is.null(footer)) 
                         div(class = "modal-footer", footer))
      ), 
    tags$script("$('#shiny-modal').modal().focus();"))
}

# sort by priority, then by size
regions_name_sort1=regions_list[sort(regions_list$V5-regions_list$V4+1,index.return=T,decreasing=T)$ix,]
regions_name=regions_name_sort1$V17[sort(regions_name_sort1$V18,index.return=T)$ix]
results_out=paste(config_file[3,2],"/SCIP_results.txt",sep="")
if (file.exists(results_out)==F){
  file.create(results_out)
}

ui=fluidPage(
  tags$head(
    tags$style(
      HTML(".shiny-notification{position:fixed;top:calc(2.5%);left:calc(60%);
           width:280px;height:50px;font-size:17px;opacity:1;font-family:'Verdana',sans-serif}"))),
  tags$head(tags$style(".modal-dialog.xl1{width:95%}")),

  titlePanel("SCIP: Suite for CNV Interpretation and Prioritization"),
  
  fluidRow(
    column(6,selectInput(inputId="cnv_name",label="CNV Name (Start Typing to Search)",
           choice=regions_name,width="100%")),
    column(2,actionButton(inputId="save",label="Save",width="100%"),style="margin-top:25px;"),
    column(2,actionButton(inputId="prev_var",label="Previous",width="100%"),style="margin-top:25px;"),
    column(2,actionButton(inputId="next_var",label="Next",width="100%"),style="margin-top:25px"),
  ),
  
  fluidRow(
    column(2,selectInput(inputId="qual_override",label="Variant Quality",
           choices=c("No Override","Passed","Failed"),selected="No Override",width="100%")),
    
    column(4,selectInput(inputId="decision",label="Determination",
                         choices=c("Not Evaluated","Ruled Out - Quality Inadequate / Difficult to Assess","Ruled Out - Population Variation",
                                   "Ruled Out - No Gene of Interest Identified","Ruled Out - Incorrect Boundary, Fully Intronic",
                                   "Ruled Out - Non-intragenic DUP","Ruled Out - Other Reasons",
                                   "Deferred - Recessive Gene - Look for Compound Het SNV", "Already Interpreted - Same Variant Identified by Another Caller",
                                   "Further Review - Potentially Reportable","Further Review - Not Likely Reportable"),selected="Not Evaluated",width="100%")),
    
    column(6,textAreaInput(inputId="comment",label="Notes (Optional)",placeholder="Optional Notes about the Variant",
                       width="100%")),
  ),
  textOutput("last_edit_info"),
  hr(),
  
  h4("Variant Summary"), 
  DT::dataTableOutput("basic_info",width="1250px"),
  DT::dataTableOutput("basic_links",width="250px"),
  br(),
  DT::dataTableOutput("highlights1",width="800px"),
  DT::dataTableOutput("highlights2",width="800px"),
  
  hr(),
  h4("Read Depth and Mapping Quality"),
  fluidRow(
    column(2,selectInput(inputId="binsize",label="Bin Size",choices=c("1 bp"=1,"10 bp"=10,"25 bp"=25,"50 bp"=50,"100 bp"=100,"250 bp"=250,"500 bp"=500,
                                                                      "1 kb"=1000,"5 kb"=5000,"10 kb"=10000,"50 kb"=50000,"100 kb"=99999.999),selected="1000",width="100%")),
    column(2,actionButton(inputId="render_depth",label="Update Read Depth",width="100%"),style="margin-top:25px;"),
    column(3,selectInput(inputId="use_12878",label="NA12878",choices=c("View Sample Only"="sample_only","View Sample with NA12878 Overlaid"="both","View NA12878 Only"="12878_only"),
                         selected="both",width="100%")),
    column(2,actionButton(inputId="render_mq",label="Update Mapping Quality",width="100%"),style="margin-top:25px;"),
  ),
  
  plotOutput('depth',dblclick='click_depth',brush=brushOpts(id="brush_depth",resetOnNew=T),height="250px"),
  plotOutput('mq',dblclick='click_mq',brush=brushOpts(id="brush_mq",resetOnNew=T),height="250px"),
  
  hr(),
  h4("Anomalous Reads"),
  fluidRow(
    column(2,actionButton(inputId="render_reads",label="Load Anomalous Reads",width="100%"),style="margin-top:25px;"),
    column(3,selectInput(inputId="read_type",label="Read Type",choices=c("All (Paired-end and Split-reads)"="all","Paired-end Reads Only"="pe",
                                                                    "Normal Orientation, Small Insert Size Only"="si","Normal Orientation, Large Insert Size Only"="li",
                                                                    "Abnormal Orientation Only"="pe_abnormal_orientation",
                                                                    "Split-reads Only"="sr"),selected="all",width="100%")),
    column(3,sliderInput(inputId="pe_percentile",label="Paired-end Outlier Percentile",value=99.5,min=90,max=99.9,step=0.1)),
  ),
  p("Insert Size Estimates"),
  tableOutput("insert_size_estimates"),
  plotOutput("reads",dblclick='click_reads',brush=brushOpts(id="brush_reads",resetOnNew=T),height="400px"),
  DT::dataTableOutput("reads_table",width="100%"),
  
  hr(),
  h4("External and Internal Variant Databases"),
  fluidRow(
    column(3,selectInput(inputId="gnomad_af",label="gnomAD SV Allele Frequency",choices=c("No Filter"=0,"Above 0.01%"=1e-4,"Above 0.05%"=5e-4,"Above 0.1%"=1e-3,
                                                                                                                "Above 0.5%"=5e-3,"Above 1%"=0.01,"Above 5%"=0.05),selected=0,width="100%")),
    column(2,selectInput(inputId="clinvar_subset",label="ClinVar - Consequence",choices=c("No Filter"="all","P/LP Only"="p","Any non-B/LB"="non-b",
                                                                                          "B/LB Only"="b"),selected="all",width="100%")),
    column(2,selectInput(inputId="clinvar_size",label="ClinVar - Variant Size",choices=c("No Filter"=1e20,"Under 10 Mb"=10e6,"Under 5 Mb"=5e6,"Under 1 Mb"=1e6,
                                                                                         "Under 500 kb"=5e5,"Under 250 kb"=2.5e5),selected=5e6,width="100%")),
  ),
  plotOutput("gnomadsv_clinvar",dblclick='click_gnomadsv_clinvar',brush=brushOpts(id="brush_gnomadsv_clinvar",resetOnNew=T),height="600px"),
  br(),
  tabsetPanel(type="tabs",
              tabPanel("gnomAD SV",DT::dataTableOutput("gnomadsv_table",width="800px")),
              tabPanel("ClinVar",DT::dataTableOutput("clinvar_table",width="100%")),
              tabPanel("Internal Cohort",DT::dataTableOutput("cgc_table",width="800px")),
              selected="gnomAD SV"
  ),
  
  hr(),
  h4("Genomic Neighbourhood"),
  fluidRow(
    column(3,actionButton(inputId="render_genes",label="Plot Genomic Neighbourhood",width="100%"))),
  plotOutput("genes",dblclick='click_genes',brush=brushOpts(id="brush_genes",resetOnNew=T),height="400px"),
  tabsetPanel(type="tabs",
              tabPanel("Genes",DT::dataTableOutput("genes_table",width="100%")),
              tabPanel("ClinGen Dosage Map",DT::dataTableOutput("dosage_table",width="1000px")),
              selected="Genes"
  ),

  hr(),
  h4("Important Notes"),
  p("For Manta INS and BND variants: the Read Depth and Mapping Quality section is not applicable. The Anomalous Reads section is experimental. 
    Users are encouraged to double check those variants in IGV."),
  p("In the External and Internal Variant Databases section, CNVs in the opposite direction of the analyzed variant (e.g., DELs when the analyzed variant is a DUP) are not shown."),
  p("In the Genomic Neighbourhood section, ClinGen HI/TS regions/genes that do not overlap the analyzed variant are not shown. For DELs, TS information are not displayed."),
    
  hr(),
  h4("Version Control"),
  tableOutput("version_control1"),
  tableOutput("version_control2"),
  
  hr(),
  p("Interface Version 0.1.7 Public (20220202)"),
  p("Cardiac Genome Clinic, Ted Rogers Centre for Heart Research"),
  p("Division of Clinical and Metabolic Genetics & The Centre for Applied Genomics, The Hospital for Sick Children. © 2022"),
)

server=function(input,output){
  # reactive values
  ranges=reactiveValues(x=NULL,y1=NULL,name=NULL)
  ranges_default=reactiveValues(x=NULL,y1=NULL)
  load=reactiveValues(depth=0,mq=0,reads=0,gene=0)
  master_info=reactiveValues(chr=NULL,start=NULL,end=NULL,proband=NULL,type=NULL)
  
  latest_interpretation=reactiveValues(qual_override=NULL,decision=NULL,comment=NULL,time="Never",user="None")
  insert_size=reactiveValues(lower=NULL,upper=NULL,name=NULL,percentile=NULL)
  abnormal_reads=reactiveValues(pe_orientation=NULL,pe_small=NULL,pe_large=NULL,se=NULL)
  out_reads=reactiveValues(table=NULL)
  sv1=reactiveValues(gnomad=NULL,clinvar=NULL,cgc=NULL,dosage_hi=NULL,dosage_ts=NULL)
  # end reactive values
  
  observeEvent(input$prev_var,{
    current_var=which(regions_name==input$cnv_name)
    st2=0
    if (current_var>1){
      updateSelectInput(inputId="cnv_name",selected=regions_name[current_var-1])
    } else{
      showNotification("Already at the First Variant",type="error",duration=15)
      st2=1
    }
    
    if ((input$qual_override!="No Override" && input$qual_override!=latest_interpretation$qual_override) || 
        (input$decision!="Not Evaluated" && input$decision!=latest_interpretation$decision) || (input$comment!="" && input$comment!=latest_interpretation$comment)){
      write(paste(as.numeric(Sys.time()),input$cnv_name,input$qual_override,input$decision,input$comment,Sys.time(),username,sep="\t"),
            file=results_out,append=T)
      if (st2==0){
        showNotification("Interpretation Saved",type="message",duration=2)
      } else{
        latest_interpretation$qual_override=input$qual_override
        latest_interpretation$decision=input$decision
        latest_interpretation$comment=input$comment
        latest_interpretation$time=Sys.time()
        latest_interpretation$user=username
        showNotification("Interpretation Saved",type="message",duration=1)
      }
    }
  })
  
  observeEvent(input$next_var,{
    current_var=which(regions_name==input$cnv_name)
    st2=0
    if (current_var<length(regions_name)){
      updateSelectInput(inputId="cnv_name",selected=regions_name[current_var+1])
    } else{
      showNotification("Already at the Last Variant",type="error",duration=15)
      st2=1
    }
    
    if ((input$qual_override!="No Override" && input$qual_override!=latest_interpretation$qual_override) || 
        (input$decision!="Not Evaluated" && input$decision!=latest_interpretation$decision) || (input$comment!="" && input$comment!=latest_interpretation$comment)){
      write(paste(as.numeric(Sys.time()),input$cnv_name,input$qual_override,input$decision,input$comment,Sys.time(),username,sep="\t"),
            file=results_out,append=T)
      if (st2==0){
        showNotification("Interpretation Saved",type="message",duration=2)
      } else{
        latest_interpretation$qual_override=input$qual_override
        latest_interpretation$decision=input$decision
        latest_interpretation$comment=input$comment
        latest_interpretation$time=Sys.time()
        latest_interpretation$user=username
        showNotification("Interpretation Saved",type="message",duration=1)
      }
    }
  })
  
  observeEvent(input$save,{
    if ((input$qual_override!="No Override" && input$qual_override!=latest_interpretation$qual_override) || 
        (input$decision!="Not Evaluated" && input$decision!=latest_interpretation$decision) || (input$comment!="" && input$comment!=latest_interpretation$comment)){
      write(paste(as.numeric(Sys.time()),input$cnv_name,input$qual_override,input$decision,input$comment,Sys.time(),username,sep="\t"),
            file=results_out,append=T)
      latest_interpretation$qual_override=input$qual_override
      latest_interpretation$decision=input$decision
      latest_interpretation$comment=input$comment
      latest_interpretation$time=Sys.time()
      latest_interpretation$user=username
      showNotification("Interpretation Saved",type="message",duration=2)
    } else{
      showNotification("No Change from Last Saved",type="message",duration=2)
    }
  })

  observeEvent(input$cnv_name,{
    info=as.data.frame(strsplit(input$cnv_name,"\\."))[,1]
    master_info$proband=info[1]
    master_info$chr=info[2]
    master_info$start=as.numeric(info[3])
    master_info$end=as.numeric(info[4])
    master_info$type=info[5]

    updateSelectInput(inputId="use_12878",selected="both")
    updateSelectInput(inputId="read_type",selected="all")
    updateSliderInput(inputId="pe_percentile",value=99.5)
    updateSelectInput(inputId="gnomad_af",selected=0)
    updateSelectInput(inputId="clinvar_subset",selected="all")
    updateSelectInput(inputId="clinvar_size",selected=5e6)
    
    load$depth=0
    load$mq=0
    load$reads=0
    load$genes=0
    
    length=as.numeric(master_info$end)-as.numeric(master_info$start)+1
    if (length<=10000){
      updateSelectInput(inputId="binsize",selected="250")
    } else if (length<=50000){
      updateSelectInput(inputId="binsize",selected="500")
    } else if (length<=500000){
      updateSelectInput(inputId="binsize",selected="1000")
    } else if (length<=1000000){
      updateSelectInput(inputId="binsize",selected="5000")
    } else if (length<=2000000){
      updateSelectInput(inputId="binsize",selected="10000")
    } else if (length<=5000000){
      updateSelectInput(inputId="binsize",selected="50000")
    } else{
      updateSelectInput(inputId="binsize",selected="99999.999")
    }
    
    current_var=which(regions_name==input$cnv_name)
    x1=try(read.table(results_out,comment.char="",sep="\t",quote=""),silent=T)
    st1=0
    if(class(x1)=="data.frame"){
      x1_sub=x1[which(x1$V2==regions_name[current_var]),]
      if(length(x1_sub[,1])>0){
        x1_sub_sorted=x1_sub[sort(x1_sub$V1,index.return=T,decreasing=T)$ix,]
        updateSelectInput(inputId="qual_override",selected=x1_sub_sorted[1,3])
        updateSelectInput(inputId="decision",selected=x1_sub_sorted[1,4])
        updateTextAreaInput(inputId="comment",value=x1_sub_sorted[1,5])
        latest_interpretation$qual_override=x1_sub_sorted[1,3]
        latest_interpretation$decision=x1_sub_sorted[1,4]
        latest_interpretation$comment=x1_sub_sorted[1,5]
        latest_interpretation$time=x1_sub_sorted[1,6]
        latest_interpretation$user=x1_sub_sorted[1,7]
        st1=1
      }
    }
    if(st1==0){
      updateSelectInput(inputId="qual_override",selected="No Override")
      updateSelectInput(inputId="decision",selected="Not Evaluated")
      updateTextAreaInput(inputId="comment",value="")
      latest_interpretation$qual_override="No Override"
      latest_interpretation$decision="Not Evaluated"
      latest_interpretation$comment=""
      latest_interpretation$time="Never"
      latest_interpretation$user="None"
    }
    
    if(is.null(ranges$x)==T || ranges$name!=input$cnv_name){
      x2=read.table(paste(temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script05_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="")
      ranges$name=input$cnv_name
      ranges$x=c(min(x2$V1,na.rm=T),max(x2$V1,na.rm=T))
      ranges$y1=c(0,quantile(x2$V2,.95,na.rm=T)*1.5)
      ranges_default$x=c(min(x2$V1,na.rm=T),max(x2$V1,na.rm=T))
      ranges_default$y1=c(0,quantile(x2$V2,.95,na.rm=T)*1.5)
    }
  })
  
  # action buttons for loading specific section of the webpage
  observeEvent(input$render_depth,{
    load$depth=load$depth+1
  })
  
  observeEvent(input$render_mq,{
    load$mq=load$mq+1
  })
  
  observeEvent(input$render_reads,{
    load$reads=load$reads+1
  })
  
  observeEvent(input$render_genes,{
    load$genes=load$genes+1
  })
  ### end action buttons
  
  output$last_edit_info=renderText(paste("Last Interpreted: ",latest_interpretation$time," (",latest_interpretation$user,")",sep=""))
  
  output$basic_info=DT::renderDataTable({
    x3=regions_list[which(regions_list$V17==input$cnv_name),]
    auto_result=as.data.frame(strsplit(as.character(x3[1]),"\\s"))[1,1]
    length=x3[5]-x3[4]+1

    if (x3[18]<99){
      priority_print=paste(x3[18],"  ","<span class='label label-danger'>High</span>",sep="")
    }
    else if (x3[18]==99){
      priority_print=paste(x3[18],"  ","<span class='label label-warning'>Medium</span>",sep="")
    }
    else{
      priority_print=paste(x3[18],"  ","<span class='label label-info'>Low</span>",sep="")
    }
    
    as.data.frame(rbind(c("Quality","Priority","Sample ID","Chr","Start","End","Type","Size (kb)","Depth Ratio","Mapping Quality",
                    "Supporting Pairs","Opposing Pairs","Split-reads"),
                  c(auto_result,priority_print,x3[2],x3[3],x3[4],x3[5],x3[6],round(length*1e3/1e3)/1e3,round(x3[9]*1e3)/1e3,round(x3[10]*1e3)/1e3,x3[13],x3[14],x3[16])))
  },options=list(dom="t",ordering=F,headerCallback=JS("function(thead, data, start, end, display){","$(thead).remove();", "}")),rownames=F,selection='none')
  
  output$basic_links=DT::renderDataTable({
    as.data.frame(cbind(paste("<a href='https://www.deciphergenomics.org/browser#q/grch37:",master_info$chr,":",master_info$start,"-",master_info$end,"/location/grch37' target='_blank'>DECIPHER</a>",sep=""),
                    paste("<a href='https://gnomad.broadinstitute.org/region/",master_info$chr,"-",master_info$start,"-",master_info$end,"?dataset=gnomad_sv_r2_1' target='_blank'>gnomAD SV</a>",sep=""),
                    paste("<a href='http://dgv.tcag.ca/gb2/gbrowse/dgv2_hg19/?name=chr",master_info$chr,"%3A",master_info$start,"-",master_info$end,";search=Search' target='_blank'>DGV</a>",sep="")))
  },options=list(dom="t",ordering=F,headerCallback=JS("function(thead, data, start, end, display){","$(thead).remove();", "}")),rownames=F,escape=F,selection='none')
  
  output$highlights1=DT::renderDataTable({
    x18=read.table(paste(temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script08_file2.txt.gz",sep=""),comment.char="",sep="\t",quote="")
    if (length(x18$V1[which(x18$V2==1)])>=1){
      datatable(rbind("Positive Information:",as.data.frame(x18$V1[which(x18$V2==1)])),class="compact",selection='none',
                options=list(dom="t",ordering=F,headerCallback=JS("function(thead, data, start, end, display){","$(thead).remove();", "}")),rownames=F,escape=F) %>% DT::formatStyle(1,color="#ef6548",fontSize="110%")
    } else{
      datatable(as.data.frame("No Positive Information"),class="compact",selection='none',
                options=list(dom="t",ordering=F,headerCallback=JS("function(thead, data, start, end, display){","$(thead).remove();", "}")),rownames=F) %>% DT::formatStyle(1,color="#ef6548",fontSize="110%")
    }
  }) 
  
  output$highlights2=DT::renderDataTable({
    x18=read.table(paste(temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script08_file2.txt.gz",sep=""),comment.char="",sep="\t",quote="")
    if (length(x18$V1[which(x18$V2==2)])>=1){
      datatable(rbind("Negative Information:",as.data.frame(x18$V1[which(x18$V2==2)])),class="compact",selection='none',
                options=list(dom="t",ordering=F,headerCallback=JS("function(thead, data, start, end, display){","$(thead).remove();", "}")),rownames=F) %>% DT::formatStyle(1,color="#02818a",fontSize="110%")
    } else{
      datatable(as.data.frame("No Negative Information"),class="compact",selection='none',
                options=list(dom="t",ordering=F,headerCallback=JS("function(thead, data, start, end, display){","$(thead).remove();", "}")),rownames=F) %>% DT::formatStyle(1,color="#02818a",fontSize="110%")  
    }
  })
  
  output$depth=renderPlot({
    par(mgp=c(2.1,1,0),mai=c(0.4,0.62,0.5,0.1))
    length=as.numeric(master_info$end)-as.numeric(master_info$start)+1
    x6=read.table(paste(temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script05_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="")
    
    plot(0,0,col=rgb(0,0,0,0),
         xlab="",xlim=ranges$x,ylab="Read Depth",ylim=ranges$y1,xaxs="i",yaxs="i",main="Read Depth")
    if (load$depth==0){
      par(new=T)
      plot(0,0,col=rgb(0,0,0,0),xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="")
      text(0.5,0.5,"Click the \"Update Read Depth\" button above to show read depth.\n
           For larger variants, setting bin size to 5 kb or larger (before clicking the button) is recommended for faster performance.",cex=1.5)
    } else if (master_info$type=="BND" || master_info$type=="INS"){
      plot(0,0,col=rgb(0,0,0,0),
           xlab="",xlim=ranges$x,ylab="Read Depth",ylim=ranges$y1,xaxs="i",yaxs="i",main="Read Depth")
      par(new=T)
      plot(0,0,col=rgb(0,0,0,0),xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="")
      text(0.5,0.5,"This section is not applicable for INS or BND type variants.",cex=1.5)
      box()
    } else if (load$depth>0){
      # speed up using BinMean and plotH
      # using isolate so the plot will need to be manually updated
      x6_sub=x6[intersect(which(x6$V1>=ranges$x[1]-as.numeric(isolate(input$binsize))*2),which(x6$V1<=ranges$x[2]+as.numeric(isolate(input$binsize))*2)),]
      
      if (input$use_12878=="both" || input$use_12878=="12878_only"){
        x7=read.table(paste(temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script04_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="")
        norm_rd_12878=median(x7$V2[c(which(x7$V1<=master_info$start),which(x7$V1>=master_info$end))],na.rm=T)/median(x6$V2[c(which(x6$V1<=master_info$start),which(x6$V1>=master_info$end))],na.rm=T)
        x7_sub=x7[intersect(which(x7$V1>=ranges$x[1]-as.numeric(isolate(input$binsize))*2),which(x7$V1<=ranges$x[2]+as.numeric(isolate(input$binsize))*2)),]
        par(new=T)
        plot(0,0,col=rgb(0,0,0,0),
             xlab="",ylab="",xlim=ranges$x,ylim=ranges$y1,xaxs="i",yaxs="i",axes=F)
        
        # speed up using BinMean and plotH
        x7_pos=BinMean(x7_sub$V1,every=as.numeric(isolate(input$binsize)))
        x7_val=BinMean(x7_sub$V2,every=as.numeric(isolate(input$binsize)))/norm_rd_12878
        median_out_x7=median(x7$V2[c(which(x7$V1<=master_info$start),which(x7$V1>=master_info$end))],na.rm=T)/norm_rd_12878
        median_in_x7=median(x7$V2[intersect(which(x7$V1>master_info$start),which(x7$V1<master_info$end))],na.rm=T)/norm_rd_12878
        par(new=T)
        plotH(x7_pos,x7_val,width=as.numeric(isolate(input$binsize)),col="#ff7f0080",xlab="",ylab="",xlim=ranges$x,ylim=ranges$y1,xaxs="i",yaxs="i",axes=F,border=NA)
        lines(x=c(ranges_default$x[1],master_info$start),y=c(median_out_x7,median_out_x7),col="#d95f0e",lty=2)
        lines(x=c(master_info$end,ranges_default$x[2]),y=c(median_out_x7,median_out_x7),col="#d95f0e",lty=2)
        lines(x=c(master_info$start,master_info$end),y=c(median_in_x7,median_in_x7),col="#d95f0e",lty=2)        
        # for(i in seq(1,length(x7$V1),as.numeric(input$binsize))){
        #   end=i+as.numeric(input$binsize)-1
        #   rect(x7$V1[i]-0.5,0,x7$V1[end]+0.5,mean(x7$V2[i:end],na.rm=T)/norm_rd_12878,border=NA,col="#ff7f0080")
        # }
        if (input$use_12878=="both"){
          legend("topleft",legend=c(master_info$proband,"NA12878"),fill=c("#96969680","#ff7f0080"),cex=0.8,bg="white",border=NA)
        }
      }
      
      if(input$use_12878=="both" || input$use_12878=="sample_only"){
        x6_pos=BinMean(x6_sub$V1,every=as.numeric(isolate(input$binsize)))
        x6_val=BinMean(x6_sub$V2,every=as.numeric(isolate(input$binsize)))
        median_out_x6=median(x6$V2[c(which(x6$V1<=master_info$start),which(x6$V1>=master_info$end))],na.rm=T)
        median_in_x6=median(x6$V2[intersect(which(x6$V1>master_info$start),which(x6$V1<master_info$end))],na.rm=T)
        par(new=T)
        plotH(x6_pos,x6_val,width=as.numeric(isolate(input$binsize)),col="#96969680",xlab="",ylab="",xlim=ranges$x,ylim=ranges$y1,xaxs="i",yaxs="i",axes=F,border=NA)
        lines(x=c(ranges_default$x[1],master_info$start),y=c(median_out_x6,median_out_x6),col="#636363",lty=2)
        lines(x=c(master_info$end,ranges_default$x[2]),y=c(median_out_x6,median_out_x6),col="#636363",lty=2)
        lines(x=c(master_info$start,master_info$end),y=c(median_in_x6,median_in_x6),col="#636363",lty=2)
      }
      
      # for(i in seq(1,length(x6$V1),as.numeric(input$binsize))){
      #   end=i+as.numeric(input$binsize)-1
      #   rect(x6$V1[i]-0.5,0,x6$V1[end]+0.5,mean(x6$V2[i:end],na.rm=T),border=NA,col="#969696")
      # }
      
      rect(master_info$start,-1000,master_info$end,1000,col="#8dd3c733",border="#8dd3c74D")
    }
  })
  
  output$mq=renderPlot({
    par(mgp=c(2.1,1,0),mai=c(0.4,0.62,0.5,0.1))

    plot(0,0,col=rgb(0,0,0,0),
         xlab="",xlim=ranges$x,ylab="Mapping Quality",ylim=c(0,65),xaxs="i",yaxs="i",main="Mapping Quality")
    if (load$mq==0){
      par(new=T)
      plot(0,0,col=rgb(0,0,0,0),xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="")
      text(0.5,0.5,"Click the \"Update Mapping Quality\" button above to show mapping quality.\n
           For larger variants, setting bin size to 5 kb or larger (before clicking the button) is recommended for faster performance.",cex=1.5)
    } else if (master_info$type=="BND" || master_info$type=="INS"){
      plot(0,0,col=rgb(0,0,0,0),
           xlab="",xlim=ranges$x,ylab="Mapping Quality",ylim=c(0,65),xaxs="i",yaxs="i",main="Mapping Quality")
      par(new=T)
      plot(0,0,col=rgb(0,0,0,0),xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="")
      text(0.5,0.5,"This section is not applicable for INS or BND type variants.",cex=1.5)
      box()
    } else if (load$mq>0){
      if (input$use_12878=="both" || input$use_12878=="12878_only"){
        x7=read.table(paste(temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script04_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="")
        x7_sub=x7[intersect(which(x7$V1>=ranges$x[1]-as.numeric(isolate(input$binsize))*2),which(x7$V1<=ranges$x[2]+as.numeric(isolate(input$binsize))*2)),]
        
        # speed up using BinMean and plotH
        x7_pos=BinMean(x7_sub$V1,every=as.numeric(isolate(input$binsize)))
        x7_val=BinMean(x7_sub$V3,every=as.numeric(isolate(input$binsize)))
        par(new=T)
        plotH(x7_pos,x7_val,width=as.numeric(isolate(input$binsize)),col="#ff7f0080",xlab="",ylab="",xlim=ranges$x,ylim=c(0,65),xaxs="i",yaxs="i",axes=F,border=NA)
        if (input$use_12878=="both"){
          legend("topleft",legend=c(master_info$proband,"NA12878"),fill=c("#96969680","#ff7f0080"),cex=0.8,bg="white",border=NA)
        }
      }
      
      if(input$use_12878=="both" || input$use_12878=="sample_only"){
        x6=read.table(paste(temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script05_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="")
        x6_sub=x6[intersect(which(x6$V1>=ranges$x[1]-as.numeric(isolate(input$binsize))*2),which(x6$V1<=ranges$x[2]+as.numeric(isolate(input$binsize))*2)),]
        
        # speed up using BinMean and plotH
        x6_pos=BinMean(x6_sub$V1,every=as.numeric(isolate(input$binsize)))
        x6_val=BinMean(x6_sub$V3,every=as.numeric(isolate(input$binsize)))
        par(new=T)
        plotH(x6_pos,x6_val,width=as.numeric(isolate(input$binsize)),col="#96969680",xlab="",ylab="",xlim=ranges$x,ylim=c(0,65),xaxs="i",yaxs="i",axes=F,border=NA)
      }
      rect(master_info$start,-1000,master_info$end,1000,col="#8dd3c733",border="#8dd3c74D")
    }
  })
  
  observeEvent(input$click_depth,{
    brush=input$brush_depth
    if (!is.null(brush)){
      ranges$x=c(brush$xmin,brush$xmax)
      ranges$y1=c(brush$ymin,brush$ymax)
    }else{
      ranges$x=ranges_default$x
      ranges$y1=ranges_default$y1
    }
  })
  
  observeEvent(input$click_mq,{
    brush=input$brush_mq
    if (!is.null(brush)){
      ranges$x=c(brush$xmin,brush$xmax)
    }else{
      ranges$x=ranges_default$x
    }
  })
  
  output$insert_size_estimates=renderTable({
    if (load$reads>0){
      x8=read.table(paste(temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script05_file2.txt.gz",sep=""),comment.char="",sep="\t",quote="")
      # added to remove read pairs/split-reads that have one part completely outside the extended CNV region
      # x8=x8[intersect(which(x8$V5>=ranges_default$x[1]),which(x8$V6<=ranges_default$x[2])),]
      
      all_size=x8$V8[which(x8$V11=="FR_inward")]
      if (is.null(insert_size$lower)==T || insert_size$name!=input$cnv_name || insert_size$percentile!=input$pe_percentile){
        insert_size$lower=quantile(all_size,1-input$pe_percentile/100)
        insert_size$higher=quantile(all_size,input$pe_percentile/100)
        insert_size$name=input$cnv_name
        insert_size$percentile=input$pe_percentile
      }
      
      abnormal_reads$pe_orientation=x8[intersect(which(x8$V11!="FR_inward"),which(x8$V11!="SR")),]
      abnormal_reads$pe_large=x8[intersect(which(x8$V11=="FR_inward"),which(x8$V8>=insert_size$higher)),]
      abnormal_reads$pe_small=x8[intersect(which(x8$V11=="FR_inward"),which(x8$V8<=insert_size$lower)),]
      abnormal_reads$sr=x8[which(x8$V11=="SR"),]
      as.data.frame(rbind(c("Lower Bound",round(insert_size$lower*1e2)/1e2),c("Upper Bound",round(insert_size$higher*1e2)/1e2)))
    } else{
      as.data.frame(c("Click the \"Load Anomalous Reads\" button above to load this section."))
    }
  },colnames=F,spacing="s",digits=2,striped=T,bordered=T,align="c")
  
  output$reads=renderPlot({
    par(mgp=c(2.1,1,0),mai=c(0.4,0.62,0.5,0.1))
    if (load$reads==0){
      plot(0,0,col=rgb(0,0,0,0),xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="",main="Anomalous Reads")
      box()
      text(0.5,0.5,"Click the \"Load Anomalous Reads\" button above to load this section.\n\nExperimental for INS and BND type variants.",cex=2)
    }
    
    if (load$reads>0){
      table_out=NULL
      if (input$read_type=="all"){
        table_out=rbind(abnormal_reads$pe_orientation,abnormal_reads$pe_small,abnormal_reads$pe_large,abnormal_reads$sr)
      } else if (input$read_type=="pe"){
        table_out=rbind(abnormal_reads$pe_orientation,abnormal_reads$pe_small,abnormal_reads$pe_large)
      } else if (input$read_type=="pe_abnormal_orientation"){
        table_out=abnormal_reads$pe_orientation
      } else if (input$read_type=="si"){
        table_out=abnormal_reads$pe_small
      } else if (input$read_type=="li"){
        table_out=abnormal_reads$pe_large
      } else if (input$read_type=="sr"){
        table_out=abnormal_reads$sr
      }
      colnames(table_out)=c("Read_Name","Flag","Chr","Left_Start","Left_End","Right_Start","Right_End","Size","MQ","CIGAR","Read_Type")
      # table_out_interval=table_out[union(intersect(which(table_out$Left_Start>=ranges$x[1]),which(table_out$Right_End<=ranges$x[2])),
      #                                    intersect(which(table_out$Left_Start>=ranges$x[1]),which(table_out$Right_End<=ranges$x[2]))),]
      # table_out_interval=table_out[intersect(which(table_out$Left_Start<=ranges$x[2]),which(table_out$Right_End>=ranges$x[1])),]
      
      table_out_interval=table_out[union(intersect(which(table_out$Left_End>=ranges$x[1]),which(table_out$Left_Start<=ranges$x[2])),
                                         intersect(which(table_out$Right_End>=ranges$x[1]),which(table_out$Right_Start<=ranges$x[2]))),]
      table_out_interval=table_out_interval[sort(table_out_interval$Left_Start,index.return=T,decreasing=T)$ix,]
      table_out_interval$Read_Type[which(table_out_interval$Read_Type=="SR")]="Split-read"
      out_reads$table=table_out_interval
      
      if (length(table_out_interval$Chr)>=1){
        plot(0,0,col=rgb(0,0,0,0),xlim=ranges$x,ylim=c(0,(length(table_out_interval$Chr)+1)),axes=F,xlab="",ylab="",xaxs="i",yaxs="i",main="Anomalous Reads")
        box()
        axis(1)
  
        for (i in 1:length(table_out_interval$Chr)){
          if (length(table_out_interval$Chr)<=50){
            read_name=unlist(strsplit(table_out_interval$Read_Name[i],"\\:"))
            text(table_out_interval$Right_End[i],i+0.1,read_name[length(read_name)],pos=4,cex=0.8)
          }
          
          if(table_out_interval$Read_Type[i]=="FR_inward"){
            if (table_out_interval$Left_End[i]-25>=table_out_interval$Left_Start[i]){
              rect(table_out_interval$Left_Start[i],i-0.25,table_out_interval$Left_End[i]-25,i+0.25,border=NA,col="#e41a1cBF")
            }
            if (table_out_interval$Right_Start[i]+25<=table_out_interval$Right_End[i]){
              rect(table_out_interval$Right_Start[i]+25,i-0.25,table_out_interval$Right_End[i],i+0.25,border=NA,col="#e41a1cBF")
            }
            polygon(x=c(table_out_interval$Left_End[i]-25,table_out_interval$Left_End[i],table_out_interval$Left_End[i]-25),
                    y=c(i+0.25,i,i-0.25),border=NA,col="#e41a1cBF")
            polygon(x=c(table_out_interval$Right_Start[i]+25,table_out_interval$Right_Start[i],table_out_interval$Right_Start[i]+25),
                    y=c(i+0.25,i,i-0.25),border=NA,col="#e41a1cBF")
            lines(x=c(table_out_interval$Left_End[i],table_out_interval$Right_Start[i]),y=c(i,i),col="#e41a1cBF")
  
          } else if (table_out_interval$Read_Type[i]=="FR_outward"){
            if (table_out_interval$Left_Start[i]+25<=table_out_interval$Left_End[i]){
              rect(table_out_interval$Left_Start[i]+25,i-0.25,table_out_interval$Left_End[i],i+0.25,border=NA,col="#377eb8BF")
            }
            if (table_out_interval$Right_End[i]-25>=table_out_interval$Right_Start[i]){
              rect(table_out_interval$Right_Start[i],i-0.25,table_out_interval$Right_End[i]-25,i+0.25,border=NA,col="#377eb8BF")
            }
            polygon(x=c(table_out_interval$Right_End[i]-25,table_out_interval$Right_End[i],table_out_interval$Right_End[i]-25),
                    y=c(i+0.25,i,i-0.25),border=NA,col="#377eb8BF")
            polygon(x=c(table_out_interval$Left_Start[i]+25,table_out_interval$Left_Start[i],table_out_interval$Left_Start[i]+25),
                    y=c(i+0.25,i,i-0.25),border=NA,col="#377eb8BF")
            lines(x=c(table_out_interval$Left_End[i],table_out_interval$Right_Start[i]),y=c(i,i),col="#377eb8BF")
  
          } else if (table_out_interval$Read_Type[i]=="FF"){
            if (table_out_interval$Left_End[i]-25>=table_out_interval$Left_Start[i]){
              rect(table_out_interval$Left_Start[i],i-0.25,table_out_interval$Left_End[i]-25,i+0.25,border=NA,col="#7bccc4BF")
            }
            if (table_out_interval$Right_End[i]-25>=table_out_interval$Right_Start[i]){
              rect(table_out_interval$Right_Start[i],i-0.25,table_out_interval$Right_End[i]-25,i+0.25,border=NA,col="#7bccc4BF")
            }
            polygon(x=c(table_out_interval$Left_End[i]-25,table_out_interval$Left_End[i],table_out_interval$Left_End[i]-25),
                    y=c(i+0.25,i,i-0.25),border=NA,col="#7bccc4BF")
            polygon(x=c(table_out_interval$Right_End[i]-25,table_out_interval$Right_End[i],table_out_interval$Right_End[i]-25),
                    y=c(i+0.25,i,i-0.25),border=NA,col="#7bccc4BF")
            lines(x=c(table_out_interval$Left_End[i],table_out_interval$Right_Start[i]),y=c(i,i),col="#7bccc4BF")
  
          } else if (table_out_interval$Read_Type[i]=="RR"){ #RR
            if (table_out_interval$Left_Start[i]+25<=table_out_interval$Left_End[i]){
              rect(table_out_interval$Left_Start[i]+25,i-0.25,table_out_interval$Left_End[i],i+0.25,border=NA,col="#7bccc4BF")
            }
            if (table_out_interval$Right_Start[i]+25<=table_out_interval$Right_End[i]){
              rect(table_out_interval$Right_Start[i]+25,i-0.25,table_out_interval$Right_End[i],i+0.25,border=NA,col="#7bccc4BF")
            }
            polygon(x=c(table_out_interval$Left_Start[i]+25,table_out_interval$Left_Start[i],table_out_interval$Left_Start[i]+25),
                    y=c(i+0.25,i,i-0.25),border=NA,col="#7bccc4BF")
            polygon(x=c(table_out_interval$Right_Start[i]+25,table_out_interval$Right_Start[i],table_out_interval$Right_Start[i]+25),
                    y=c(i+0.25,i,i-0.25),border=NA,col="#7bccc4BF")
            lines(x=c(table_out_interval$Left_End[i],table_out_interval$Right_Start[i]),y=c(i,i),col="#7bccc4BF")
  
          } else if (table_out_interval$Read_Type[i]=="Split-read"){
            rect(table_out_interval$Left_Start[i],i-0.25,table_out_interval$Left_End[i],i+0.25,border=NA,col="#984ea3")
            rect(table_out_interval$Right_Start[i],i-0.25,table_out_interval$Right_End[i],i+0.25,border=NA,col="#984ea3")
            lines(x=c(table_out_interval$Left_End[i],table_out_interval$Right_Start[i]),y=c(i,i),col="#984ea3")
          }
        }
        
        rect(master_info$start,-100,master_info$end,length(table_out_interval$Chr)*2.1,col="#8dd3c71A",border="#8dd3c74D")
        legend("bottomleft",legend=c("-->  <-- [deletion]","<--  --> [tandem dup or translocation]","-->  --> or <--  <-- [inversion]","Split-reads"),
               fill=c("#e41a1c","#377eb8","#7bccc4","#984ea3"),border=NA,cex=0.8,bg="white")
      } else{
        plot(0,0,col=rgb(0,0,0,0),xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="")
        box()
        text(0.5,0.5,"No reads were found based on the filter(s).",cex=2)
      }
    }
  })

  observeEvent(input$click_reads,{
    brush=input$brush_reads
    if (!is.null(brush)){
      ranges$x=c(brush$xmin,brush$xmax)
    }else{
      ranges$x=ranges_default$x
    }
  })
    
  output$reads_table=DT::renderDataTable({
    if (load$reads>0){
      out_reads$table[sort(out_reads$table$Left_Start,index.return=T,decreasing=F)$ix,]
    }
  },rownames=F,selection='none')
  
  output$gnomadsv_clinvar=renderPlot({
    sv1$gnomad=NULL
    sv1$clinvar=NULL
    sv1$cgc=NULL
    
    x9=read.table(paste(temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script06_file1.txt.gz",sep=""),comment.char="",sep="\t",fill=T)
    x9_gnomad=x9[which(x9$V1=="gnomAD_SV"),]
    
    if (length(as.data.frame(x9_gnomad)[,1])>0){
      x9_gnomad$V4=as.numeric(as.character(x9_gnomad$V4))
      x9_gnomad$V5=as.numeric(as.character(x9_gnomad$V5))
      x9_gnomad$V7=as.numeric(as.character(x9_gnomad$V7))
      gnomad_af_filt=x9_gnomad[c(which(x9_gnomad$V7>=as.numeric(input$gnomad_af)),which(is.na(x9_gnomad$V7)==T)),]
      gnomad_interval=gnomad_af_filt[intersect(which(gnomad_af_filt$V5>=ranges$x[1]),which(gnomad_af_filt$V4<=ranges$x[2])),]
      sv1$gnomad=gnomad_interval
    }
    
    x9_clinvar=x9[which(x9$V1=="ClinVar"),]
    
    if (length(as.data.frame(x9_clinvar)[,1])>0){
      x9_clinvar$V3=as.numeric(as.character(x9_clinvar$V3))
      x9_clinvar$V4=as.numeric(as.character(x9_clinvar$V4))
      clinvar_len=x9_clinvar$V4-x9_clinvar$V3+1
      clinvar_len_filt=x9_clinvar[which(clinvar_len<=as.numeric(input$clinvar_size)),]
      
      if (input$clinvar_subset=="all"){
        clinvar_consequence_filt=clinvar_len_filt
      } else if (input$clinvar_subset=="non-b"){
        clinvar_consequence_filt=clinvar_len_filt[which(clinvar_len_filt$V11!="B"),]
      } else if (input$clinvar_subset=="b"){
        clinvar_consequence_filt=clinvar_len_filt[which(clinvar_len_filt$V11=="B"),]
      } else if (input$clinvar_subset=="p"){
        clinvar_consequence_filt=clinvar_len_filt[which(clinvar_len_filt$V11=="P"),]
      }
      clinvar_interval=clinvar_consequence_filt[intersect(which(clinvar_consequence_filt$V4>=ranges$x[1]),which(clinvar_consequence_filt$V3<=ranges$x[2])),]
      sv1$clinvar=clinvar_interval
    }
    
    x9=read.table(paste(temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script06_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="",fill=T)
    x9_cgc=x9[which(x9$V1=="cohort_CNV"),]
    
    if (length(as.data.frame(x9_cgc)[,1])>0){
      x9_cgc$V3=as.numeric(as.character(x9_cgc$V3))
      x9_cgc$V4=as.numeric(as.character(x9_cgc$V4))
      cgc_interval=x9_cgc[intersect(which(x9_cgc$V4>=ranges$x[1]),which(x9_cgc$V3<=ranges$x[2])),]
      cgc_interval_sorted=cgc_interval[sort(as.character(cgc_interval[,6]),index.return=T)$ix,]
      sv1$cgc=cgc_interval_sorted
    }
    
    par(mgp=c(2.1,1,0),mai=c(0.4,0.62,0.5,0.1))
    plot(0,0,col=rgb(0,0,0,0),xlim=ranges$x,ylim=c(-20,10),axes=F,xlab="",ylab="",xaxs="i",yaxs="i",main="External and Internal Variant Databases")
    box()
    axis(1)
    abline(h=0)
    abline(h=-10)
    
    if(length(sv1$gnomad[,1])>0){
      step=10/(length(sv1$gnomad[,1])+1)
      for(i in 1:length(sv1$gnomad[,1])){
        level=10-step*i
        rect(sv1$gnomad[i,4],level-0.25*step,sv1$gnomad[i,5],level+0.25*step,border="#969696",col="#969696")
        if(length(sv1$gnomad[,1])<=25){
          if (sv1$gnomad[i,5]>=ranges$x[2]){
            text(ranges$x[2],level,paste(sv1$gnomad[i,9],", AF=",sv1$gnomad[i,7],sep=""),cex=0.8,pos=2)
          } else{
            text(sv1$gnomad[i,5],level,paste(sv1$gnomad[i,9],", AF=",sv1$gnomad[i,7],sep=""),cex=0.8,pos=4)
          }
        }
      }
    } else{
      text((ranges$x[1]+ranges$x[2])/2,5,"No gnomAD-SV variants were found based on the filter(s).",cex=2)
    }
    
    if (length(sv1$clinvar[,1])>0){
      step=10/(length(sv1$clinvar[,1])+1)
      for(i in 1:length(sv1$clinvar[,1])){
        level=0-step*i
        col_clinvar="#969696"
        if (sv1$clinvar[i,11]=="P"){
          col_clinvar="#f781bf"
        } else if (sv1$clinvar[i,11]=="B"){
          col_clinvar="#e6f5c9"
        }
        
        rect(sv1$clinvar[i,3],level-0.25*step,sv1$clinvar[i,4],level+0.25*step,border=col_clinvar,col=col_clinvar)
        if(length(sv1$clinvar[,1])<=25){
          clinvar_name1=unlist(strsplit(unlist(strsplit(sv1$clinvar[i,10],"<"))[2],">"))[2]
          if (sv1$clinvar[i,4]>=ranges$x[2]){
            text(ranges$x[2],level,clinvar_name1,cex=0.8,pos=2)
          } else{
            text(sv1$clinvar[i,4],level,clinvar_name1,cex=0.8,pos=4)
          }
        }
      }
    } else{
      text((ranges$x[1]+ranges$x[2])/2,-5,"No ClinVar variants were found based on the filter(s).",cex=2)
    }
    
    if (length(sv1$cgc[,1])>0){
      step=10/(length(sv1$cgc[,1])+1)
      for(i in 1:length(sv1$cgc[,1])){
        level=-10-step*i
        col_cgc="#969696"
        if (sv1$cgc[i,8]=="1"){
          col_cgc="#fb8072"
        }
        
        rect(sv1$cgc[i,3],level-0.25*step,sv1$cgc[i,4],level+0.25*step,border=col_cgc,col=col_cgc)
        if(length(sv1$cgc[,1])<=25){
          if (sv1$cgc[i,4]>=ranges$x[2]){
            text(ranges$x[2],level,sv1$cgc[i,6],cex=0.7,pos=2)
          } else{
            text(sv1$cgc[i,4],level,sv1$cgc[i,6],cex=0.7,pos=4)
          }
        }
      }
    } else{
      text((ranges$x[1]+ranges$x[2])/2,-15,"No variants in other internal samples were found based on the filter(s).",cex=2)
    }

    text(ranges$x[1],0.8,"gnomAD-SV",pos=4)
    text(ranges$x[1],-0.8,"ClinVar",pos=4)
    text(ranges$x[1],-10.8,paste("Internal Cohort (Total=",length(sv1$cgc[,1]),")",sep=""),pos=4)
    
    legend("bottomleft",legend=c("ClinVar P/LP","ClinVar B/LB","Same Family"),fill=c("#f781bf","#e6f5c9","#fb8072"),cex=0.8,bg="white",border=NA)
    rect(master_info$start,-1000,master_info$end,1000,col="#8dd3c71A",border="#8dd3c74D")
  })

  observeEvent(input$click_gnomadsv_clinvar,{
    brush=input$brush_gnomadsv_clinvar
    if (!is.null(brush)){
      ranges$x=c(brush$xmin,brush$xmax)
    }else{
      ranges$x=ranges_default$x
    }
  })
  
  output$gnomadsv_table=DT::renderDataTable({
    if(is.null(sv1$gnomad)==F){
      x10=sv1$gnomad[,c(3:8,2)]
      colnames(x10)=c("Chr","Start","End","Type","Popmax AF","QC Flag","ID")
      x10
    } else{
      datatable(as.data.frame("No gnomAD-SV variants were found based on the filter(s)."),selection='none',options=list(dom="t",ordering=F,
                                                                                                       headerCallback=JS("function(thead, data, start, end, display){","$(thead).remove();", "}")),rownames=F)
    }
  },escape=F,rownames=F,selection='none')
  
  output$clinvar_table=DT::renderDataTable({
    if(is.null(sv1$clinvar)==F){
      x11=sv1$clinvar[,c(2:4,6,5,7:10)]
      colnames(x11)=c("Chr","Start","End","Type","Interpretation","Condition","Allele Origin","Genes","ID")
      x11
    } else{
      datatable(as.data.frame("No ClinVar variants were found based on the filter(s)."),selection='none',options=list(dom="t",ordering=F,
                                                                                                       headerCallback=JS("function(thead, data, start, end, display){","$(thead).remove();", "}")),rownames=F)
    }
  },escape=F,rownames=F,selection='none')
  
  output$cgc_table=DT::renderDataTable({
    if(is.null(sv1$cgc)==F){
      x12=sv1$cgc[,c(2:7)]
      colnames(x12)=c("Chr","Start","End","Type","Sample","Algorithm")
      x12
    } else{
      datatable(as.data.frame("No variants in other CGC samples were found based on the filter(s)."),selection='none',options=list(dom="t",ordering=F,
                                                                                                     headerCallback=JS("function(thead, data, start, end, display){","$(thead).remove();", "}")),rownames=F)
    }
  },escape=F,rownames=F,selection='none')
  
  output$genes_table=DT::renderDataTable({
    x13=read.table(paste(temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script07_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="")
    x13=x13[,c(1,19,5:17)]
    colnames(x13)=c("Name","Strand","Overlap","HI","TS","pLI","LOEUF","pRec","pNull","Expression","GO","OMIM","GenCC","Links","Search")
    
    datatable(x13,escape=F,class="cell-border stripe",rownames=NULL,options=list(pageLength=25),selection=list(mode="single", target="cell")) %>%
      DT::formatStyle(columns=12:13,fontSize='90%') %>%
      DT::formatStyle(columns=c(6:10,14,15),"text-align"='center') %>%
      DT::formatStyle('pLI',backgroundColor=styleInterval(c(0.5,0.9),c("white","#ffffb2","#fcae91"))) %>%
      DT::formatStyle('LOEUF',backgroundColor=styleInterval(c(0.35),c("#fcae91","white"))) %>%
      DT::formatStyle('pRec',backgroundColor=styleInterval(c(0.9),c("white","#bdc9e1"))) %>%
      DT::formatStyle(columns=c(15),width='100px')
  })

  observeEvent(input$genes_table_cells_selected,{
    x13_v1=read.table(paste(temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script07_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="")
    coord=unlist(input$genes_table_cells_selected)
    if (length(coord)==2){
      if (coord[2]==1-1 || coord[2]==2-1 || coord[2]==3-1){
        showModal(mymodal(title="Transcript Information",easyClose=T,
                          HTML(x13_v1$V20[as.numeric(coord[1])]),idcss="xl1"))
      }
    }
  })
    
  output$dosage_table=DT::renderDataTable({
    x9=read.table(paste(temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script06_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="",fill=T)
    x9_dosage=x9[which(x9$V1=="ClinGen_Region"),]
    hi_region_index=c(which(x9_dosage$V6=="Sufficient (3)"),which(x9_dosage$V6=="Emerging (2)"),which(x9_dosage$V6=="Little (1)"),
                      which(x9_dosage$V6=="Recessive (30)"),which(x9_dosage$V6=="Unlikely (40)"))
    ts_region_index=c(which(x9_dosage$V7=="Sufficient (3)"),which(x9_dosage$V7=="Emerging (2)"),which(x9_dosage$V7=="Little (1)"),
                      which(x9_dosage$V7=="Recessive (30)"),which(x9_dosage$V7=="Unlikely (40)"))

    ds_region_index=union(hi_region_index,ts_region_index)
    if (master_info$type=="DEL"){
      ds_region_index=hi_region_index
    }
    if (master_info$type=="DEL"){
      ds_table=x9_dosage[ds_region_index,2:6]
      colnames(ds_table)=c("Name","Chr","Start","End","Haploinsufficiency")
    } else{
      ds_table=x9_dosage[ds_region_index,2:7]
      colnames(ds_table)=c("Name","Chr","Start","End","Haploinsufficiency","Triplosensitivity")      
    }
    ds_table
  },escape=F,rownames=F,options=list(pageLength=25),selection='none')
  
  output$genes=renderPlot({
    par(mgp=c(2.1,1,0),mai=c(0.4,0.62,0.5,0.1))
    
    plot(0,0,col=rgb(0,0,0,0),
         xlab="",xlim=ranges$x,ylab="",ylim=c(0,3.5),xaxs="i",yaxs="i",axes=F,main="Genomic Neighbourhood")
    box()
    
    if (load$genes==0){
      text((as.numeric(ranges$x[1])+as.numeric(ranges$x[2]))/2,1.75,"Click the \"Plot Genomic Neighbourhood\" button above to load this plot.",cex=2)
    }
    if (load$genes>0){
      x9=read.table(paste(temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script06_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="",fill=T)
      x9_dosage=x9[which(x9$V1=="ClinGen_Region"),]
      if (length(as.data.frame(x9_dosage)[,1])>0){
        hi_region_index=c(which(x9_dosage$V6=="Sufficient (3)"),which(x9_dosage$V6=="Emerging (2)"),which(x9_dosage$V6=="Little (1)"),
                          which(x9_dosage$V6=="Recessive (30)"),which(x9_dosage$V6=="Unlikely (40)"))
        ts_region_index=c(which(x9_dosage$V7=="Sufficient (3)"),which(x9_dosage$V7=="Emerging (2)"),which(x9_dosage$V7=="Little (1)"),
                          which(x9_dosage$V7=="Recessive (30)"),which(x9_dosage$V7=="Unlikely (40)"))
        ds_region_index=union(hi_region_index,ts_region_index)
        if (master_info$type=="DEL"){
          ds_region_index=hi_region_index
        }
        sv1$dosage_hi=x9_dosage[hi_region_index,]
        sv1$dosage_ts=x9_dosage[ts_region_index,]
      }
  
      axis(1)
      axis(2,at=seq(0,1,0.2))
      abline(h=c(1,2,2.5,3))
      
      if (length(sv1$dosage_hi[,1])>0){
        for (i in 1:length(sv1$dosage_hi[,1])){
          col_hi="#969696"
          if (sv1$dosage_hi[i,6]=="Sufficient (3)"){
            col_hi="#f03b20"
          } else if (sv1$dosage_hi[i,6]=="Emerging (2)"){
            col_hi="#fd8d3c"
          } else if (sv1$dosage_hi[i,6]=="Little (1)"){
            col_hi="#fecc5c"
          }
          
          step_hi=0.5/(length(sv1$dosage_hi[,1])+1)
          rect(sv1$dosage_hi[i,4],2.5+(i-0.25)*step_hi,sv1$dosage_hi[i,5],2.5+(i+0.25)*step_hi,border=NA,col=col_hi)
          hi_name1=unlist(strsplit(unlist(strsplit(sv1$dosage_hi[i,2],"<"))[2],">"))[2]
          if (sv1$dosage_hi[i,5]>=ranges$x[2]){
            if(sv1$dosage_hi[i,4]<=ranges$x[2]){      
              text(ranges$x[2],2.5+i*step_hi,hi_name1,cex=0.7,pos=2)
            }
          } else{
            text(sv1$dosage_hi[i,5],2.5+i*step_hi,hi_name1,cex=0.7,pos=4)
          }
        }
      } else{
        text((ranges$x[1]+ranges$x[2])/2,2.75,"This variant does not overlap with ClinGen HI regions.")
      }
      
      if (master_info$type=="DEL"){
        text((ranges$x[1]+ranges$x[2])/2,3.25,"TS not applicable because the variant is a DEL.")
      } else{
        if (length(sv1$dosage_ts[,1])>0){
          for (i in 1:length(sv1$dosage_ts[,1])){
            col_ts="#969696"
            if (sv1$dosage_ts[i,7]=="Sufficient (3)"){
              col_ts="#f03b20"
            } else if (sv1$dosage_ts[i,7]=="Emerging (2)"){
              col_ts="#fd8d3c"
            } else if (sv1$dosage_ts[i,7]=="Little (1)"){
              col_ts="#fecc5c"
            }
          
            step_ts=0.5/(length(sv1$dosage_ts[,1])+1)
            rect(sv1$dosage_ts[i,4],3+(i-0.25)*step_ts,sv1$dosage_ts[i,5],3+(i+0.25)*step_ts,border=NA,col=col_ts)
            ts_name1=unlist(strsplit(unlist(strsplit(sv1$dosage_ts[i,2],"<"))[2],">"))[2]
            if (sv1$dosage_ts[i,5]>=ranges$x[2]){
              if(sv1$dosage_ts[i,4]<=ranges$x[2]){ 
                text(ranges$x[2],3+i*step_ts,ts_name1,cex=0.7,pos=2)
              }
            } else{
              text(sv1$dosage_ts[i,5],3+i*step_ts,ts_name1,cex=0.7,pos=4)
            }
          }
        } else{
          text((ranges$x[1]+ranges$x[2])/2,3.25,"This variant does not overlap with ClinGen TS regions.")
        }
      }
      
      x15=read.table(paste(temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script05_file3.txt.gz",sep=""),comment.char="",sep="\t",quote="")
      y_step=1/(max(x15$V2)+2)
      for (i in 1:length(x15$V1)){
        level=1+y_step*(x15$V2[i]+1)
        pos=unlist(strsplit(x15[i,3],"\\:"))
        for (j in seq(1,length(pos),2)){
          rect(pos[j],level-y_step*0.25,pos[j+1],level+y_step*0.25,border="black",col="black")
          if (j!=1){
            lines(x=c(pos[j-1],pos[j]),y=c(level,level))
          }
          if (ranges$x[2]-ranges$x[1]<=5e4){
            if ((j+1)%%4==0){
              text((as.numeric(pos[j])+as.numeric(pos[j+1]))/2,level+y_step*0.25,x15$V1[i],pos=3,cex=0.7)
            } else{
              text((as.numeric(pos[j])+as.numeric(pos[j+1]))/2,level-y_step*0.25,x15$V1[i],pos=1,cex=0.7)
            }
          }
        }
      }
      
      x14=read.table(paste(temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script05_file1.txt.gz",sep=""))
      x14sub1=x14$V1[is.na(x14$V4)==F]
      x14sub4=x14$V4[is.na(x14$V4)==F]
      
      # speed up using plotH
      # for(i in 1:length(x14sub1)){
      #   if (x14sub4[i]>=0.2){
      #     rect(x14sub1[i]-0.5,0,x14sub1[i]+0.5,x14sub4[i],border="#756bb1",col="#756bb1")
      #   } else{
      #     rect(x14sub1[i]-0.5,0,x14sub1[i]+0.5,x14sub4[i],border="#969696",col="#969696")        
      #   }
      # }
      par(new=T)
      plotH(x14sub1,x14sub4,width=1,col="#756bb1",xlab="",ylab="",xlim=ranges$x,ylim=c(0,3.5),xaxs="i",yaxs="i",axes=F,border="#756bb1")

      x17=read.table(paste(temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".script07_file1.txt.gz",sep=""),comment.char="",sep="\t",quote="")
      gnomad_pli_gene=c()
      gnomad_pli_start=c()
      gnomad_pli_end=c()
      gnomad_pli_cat=c()
      for (i in 1:length(x17$V1)){
        if(is.na(x17$V8[i])==F){
          if (x17$V8[i]>=0.9 || x17$V9[i]<=0.35){
            gnomad_pli_gene=c(gnomad_pli_gene,unlist(strsplit(unlist(strsplit(x17$V1[i],"<"))[2],">"))[2])
            gnomad_pli_start=c(gnomad_pli_start,x17$V3[i])
            gnomad_pli_end=c(gnomad_pli_end,x17$V4[i])
            gnomad_pli_cat=c(gnomad_pli_cat,1)
          } else if (x17$V8[i]>=0.5){
            gnomad_pli_gene=c(gnomad_pli_gene,unlist(strsplit(unlist(strsplit(x17$V1[i],"<"))[2],">"))[2])
            gnomad_pli_start=c(gnomad_pli_start,x17$V3[i])
            gnomad_pli_end=c(gnomad_pli_end,x17$V4[i])
            gnomad_pli_cat=c(gnomad_pli_cat,2)
          } else if (x17$V8[i]>=0.1){
            gnomad_pli_gene=c(gnomad_pli_gene,unlist(strsplit(unlist(strsplit(x17$V1[i],"<"))[2],">"))[2])
            gnomad_pli_start=c(gnomad_pli_start,x17$V3[i])
            gnomad_pli_end=c(gnomad_pli_end,x17$V4[i])
            gnomad_pli_cat=c(gnomad_pli_cat,3)
          }
        }
      }

      if (length(gnomad_pli_cat)>0){
        gnomad_pli_color=rep("#ffffcc",length(gnomad_pli_cat))
        gnomad_pli_color[which(gnomad_pli_cat==2)]="#a1dab4"
        gnomad_pli_color[which(gnomad_pli_cat==1)]="#2c7fb8"
        
        x17_step=0.5/(length(gnomad_pli_cat)+1)
        for (i in 1:length(gnomad_pli_cat)){
          rect(gnomad_pli_start[i],2+(i-0.25)*x17_step,gnomad_pli_end[i],2+(i+0.25)*x17_step,border=gnomad_pli_color[i],col=gnomad_pli_color[i])
          if (gnomad_pli_end[i]>=ranges$x[2]){
            if(gnomad_pli_start[i]<=ranges$x[2]){      
              text(ranges$x[2],2+i*x17_step,gnomad_pli_gene[i],cex=0.7,pos=2)
            }
          } else{
            text(gnomad_pli_end[i],2+i*x17_step,gnomad_pli_gene[i],cex=0.7,pos=4)
          }
        }
      } else{
        text((as.numeric(ranges$x[1])+as.numeric(ranges$x[2]))/2,2.25,"This variant does not overlap with genes with pLI >= 0.1.")
      }
            
      text(ranges$x[1],3.4,"ClinGen Triplosensitivity Map",pos=4)
      text(ranges$x[1],2.9,"ClinGen Haploinsufficiency Map",pos=4)
      text(ranges$x[1],2.4,"gnomAD HI Constraint",pos=4)
      text(ranges$x[1],1.9,"Genes",pos=4)
      text(ranges$x[1],0.9,"pext",pos=4)
      legend("topright",legend=c("Sufficient","Emerging","Little","Recessive/Unlikely","pLI 0.9+ / LOEUF 0.35-","pLI 0.5+","pLI 0.1+"),
             fill=c("#f03b20","#fd8d3c","#fecc5c","#969696","#2c7fb8","#a1dab4","#ffffcc"),cex=0.8,bg="white",border=NA)
      # legend("bottomleft",legend=c("pext >= 0.2"),fill=c("#756bb1"),cex=0.8,bg="white",border=NA) # currently all pext bars (regardless of level) are in purple
      rect(master_info$start,-1000,master_info$end,1000,col="#8dd3c71A",border="#8dd3c74D")
    }
  })
  
  observeEvent(input$click_genes,{
    brush=input$brush_genes
    if (!is.null(brush)){
      ranges$x=c(brush$xmin,brush$xmax)
    }else{
      ranges$x=ranges_default$x
    }
  })
  
  output$version_control1=renderTable({
    x16=read.table(paste(temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".log.txt.gz",sep=""),fill=T,comment.char="",sep="\t",quote="",head=F)
    x16[1:3,1:2]
  },colnames=F,spacing="s",striped=T,bordered=T)

  output$version_control2=renderTable({
    x16=read.table(paste(temp_file_dir,"/",master_info$proband,"/",input$cnv_name,".log.txt.gz",sep=""),fill=T,comment.char="",sep="\t",quote="",head=F)
    x16[4:length(x16$V1),1:3]
  },colnames=F,spacing="s",striped=T,bordered=T)
}

shinyApp(ui=ui,server=server)
