median_depth_in=median(x$V4[which(x$V8==1)],na.rm=T)
median_depth_out=median(x$V4[which(x$V8==0)],na.rm=T)
median_mq_in=median(x$V5[which(x$V8==1)],na.rm=T)
median_mq_out=median(x$V5[which(x$V8==0)],na.rm=T)

ct_paired_end_support=0
ct_paired_end_opposite=0
ct_paired_end_non_supporting=0
if (rb-lb<=999){
}else{
  if (is.null(dim(y))==FALSE){
    for (i in 1:length(y$V1)){
      if (abs(y$V1[i]-lb)<=2500 && abs(y$V2[i]-rb)<=2500){
        if ((type=="DUP" && y$V3[i]=="FR_outward") || (type=="DEL" && y$V3[i]=="FR_inward")){
          ct_paired_end_support=ct_paired_end_support+1
        }
        if ((type=="DUP" && y$V3[i]!="FR_outward") || (type=="DEL" && y$V3[i]!="FR_inward")){
          ct_paired_end_opposite=ct_paired_end_opposite+1
        }
      }
    }
  }
  ct_paired_end_non_supporting=length(y$V1)-ct_paired_end_support-ct_paired_end_opposite
}

ct_split_read=0
if (is.null(dim(z))==FALSE){
  for (i in 1:length(z$V1)){
    if (abs(z$V1[i]-lb)<=2500 && abs(z$V2[i]-rb)<=2500){
      ct_split_read=ct_split_read+1
    }
  }
}

if (type!="DUP" && type!="DEL"){
  ct_split_read=0
  ct_paired_end_support=0
  ct_paired_end_opposite=0
  ct_paired_end_non_supporting=0
}
