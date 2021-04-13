## rqtl
## a function to extract marker information (markers and corresponding positions) from either draws or probability
getmap <- function(x) attributes(x$draws)$map
getmap.prob<-function(x) attributes(x$prob)$map


#for Arabidopsis
# extracting chromosome names that are unrepeated when extracting by r/qtl
getChr<-function(x){
      tem<-NULL
     snames<-names(x)
     for (i in 1:length(names(x))){
       # dealing with a chromosome name with quote. i.e. '1' (regular type:unrepeated rows) 
       myText = paste(c('rep(names(',deparse(substitute(x)),')[',toString(i),'],length(x$`',snames[i],'`))'), collapse='')
       pos<-paste(c('x$`',snames[i],'`'),collapse = '')
       tem<-rbind(tem,cbind(eval(parse(text = myText)),eval(parse(text=pos))))
     
     }
     return(tem)
     
   }
## for switchgrass 
## concatenate chromosome names and corresponding marker positions   
getChrom<-function(x){
   tem<-NULL
  snames<-names(x)
  for (i in 1:length(names(x))){
    # dealing with a chromosome name with quote. i.e. '1K'(switchgrass data has repeated values)
    myText = paste(c('rep(names(',deparse(substitute(x)),')[',toString(i),'],length(x$`',snames[i],'`[1,]))'), collapse='')
    pos<-paste(c('x$`',snames[i],'`[1,]'),collapse = '')
    tem<-rbind(tem,cbind(eval(parse(text = myText)),eval(parse(text=pos))))
  
  }
  return(tem)
  
}
  

