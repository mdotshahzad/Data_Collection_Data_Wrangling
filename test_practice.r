#This is a practice file for data wrangling


file.exists('data')

file.exists('test_practice.r')


########  Chech Directory ###########

if (!file.exists('data')) {
  dir.create('data')
}
setwd('./data')
getwd()


setwd('../')
getwd()


########  Download File ###########
download.file() #download a file from internet... hum kisi b qisam ka data dl kar skty hain
#This help with reproduceablity with recent stats and data
#url,destfile,method are imp parameters

fileurl <- "https://www.stats.govt.nz/assets/Uploads/Annual-enterprise-survey/Annual-enterprise-survey-2020-financial-year-provisional/Download-data/annual-enterprise-survey-2020-financial-year-provisional-csv.csv"
download.file(fileurl, destfile = "camera_survey.csv", method = "curl")

setwd('./data')
getwd()

list.files('./data')
cameradata <- read.csv('camera_survey.csv')

str(cameradata)

cdata <- read.table('camera_survey.csv', sep = ',', header = TRUE)
getwd()
list.files('./data')


View(cdata)

"""
#parameters of read.csv
quote()
na.string
nrows
skip
"""

#read.csv2 is a alternative but fast.

# library is data.table() and use fread() function which is very fast

syntax:
  library(data.table) 
  data  = fread("camera_survey.csv")

  fread("A,B\n1,2\n3,4")
  fread("
This is perhaps a banner line or two or ten.
A,B
1,2
3,4
")
  # Detects whether column names are present automatically :
  fread("
1,2
3,4
")



  # Numerical precision :
  
  DT = fread("A\n1.010203040506070809010203040506\n")
  # TODO: add numerals=c("allow.loss", "warn.loss", "no.loss") from base::read.table, +"use.Rmpfr"
  typeof(DT$A)=="double"   # currently "allow.loss" with no option
  
  
  DT = fread("A\n1.46761e-313\n")   # read as 'numeric'
  DT[,sprintf("%.15E",A)]   # beyond what double precision can store accurately to 15 digits
  # For greater accuracy use colClasses to read as character, then package Rmpfr.
  
  
  
  getwd()
list.files()  
newwritexlx <- writexl::write_xlsx(cdata)  
list.files()

#######   WEB SCRAPPING ######

# 1: Read data from web pages
# 2: Read data from the specific section
# 3: Read Required table only

library(rvest)

theurl <- "http://en.wikipedia.org/wiki/Brazil_national_football_team"
theurl

class(theurl)

file <- read_html(theurl)

# Read table from webpage

tables <- html_nodes(file, "table") 

tables

print(tables)

tab1 <- html_table(tables[4], fill = TRUE)
print(tab1)


library(RSQLite)
 date("mtcars")
library(datasets)
 
 data <- datasets::mtcars
data 
 
names_car <- rownames(data)

names_car

rownames(data) <- c() #it will clear the row name

data
head(data)

#######################################
######### Create a Datebase ###########
#######################################


conn <- dbConnect(RSQLite::SQLite(), "datadatebase.db")



dbWriteTable(conn, "carTable", data )

dbListTables(conn)

head(data)

dbListFields(conn, "carTable")
dbReadTable(conn, "carTable")


dbGetQuery(conn, "select * from carTable where carb = '4'")
dbGetQuery(conn, "select * from carTable limit 6")

hpc <- dbGetQuery(conn, "select cyl, avg(hp) as ahp, avg(carb) as acar from carTable
group by cyl
order by ahp")


hpc
data
# Parametersed query
mpg <- 18
cyl <- 6

Result <- dbGetQuery(conn, "select * from carTable 
                     where mpg >= ? and cyl = ? ", params = c(mpg,cyl))

Result


#dbExecute() this will use to insert update delete data 

dbExecute(conn, "delete from carTable where cyl = '6'")

dbGetQuery(conn, "select * from carTable")


res <- dbSendQuery(conn, "select * from carTable where cyl = 4")
res
dbFetch(res)

while(!dbHasCompleted(res)) {
  chunk <- dbFetch(res, n=2)
  print(nrow(chunk))
}


dbClearResult(res)
dbDisconnect(conn)



#######################################
###### Creating new Datebase ##########
#######################################


lob <- datasets::Loblolly
lob

rownames <- rownames(lob)

rownames


rownames(lob) <- c()

lob


conn <- dbConnect(RSQLite::SQLite(), "lob.db")




dbWriteTable(conn, "lobTable", lob )

##########################################################
##########################################################
##########################################################


library(data.table)

DF = data.frame(x=rnorm(9), y=rep(c('A','b','c'), each=3),z=rnorm(9) )


head(DF, 3) 
                

Dt = data.table(x=rnorm(9), y=rep(c('A','b','c'), each=3),z=rnorm(9) )


head(Dt, 3) 

tables()
Dt
Dt[,list(mean(x),median(y) ,sum(z))]
Dt[,table(y)]
#Adding new column in data.table library

Dt
Dt[,w:=z^2]
Dt
Dt[,a:=x>0]
Dt
Dt[10:=1,]




##############################################
############  Megre Data Table################
##############################################






##############################################
############  ordering Table  ################
##############################################

Dt
View(Dt)

Dt[order(Dt$x),] #Acsending
Dt[order(-Dt$x),] #decsending


########################
# data frame with missing data
#Recode/Impute Missing Values




df <- data.frame(col1 = c(1:3, NA), col2 = c("this", NA,"is", "text"),
                 col3 = c(TRUE, FALSE, TRUE, TRUE), col4 = c(2.5, 4.2, 3.2, NA),
                 stringsAsFactors= FALSE)
df



df$col4[is.na(df$col4)]<- mean(df$col4, na.rm = TRUE)
df

##############################################
############       Tidyverse  ################
##############################################
"""
Its collection of libraries that use in 

data cleansing
transformation
exploratery analysis 
statistical analysis
forcasting.



It makes data cleansing and preperation very easy for data scientist
========================

dplyr is a data exploratery and transformation library
• select(): select specific columns from your dataset
• filter(): filter out rows that meet a certain condition(s)
• group_by(): group together observations. 
• summarise(): generate summary results
• arrange(): arrange the columns in ascending or descending order 
• join(): perform joins: left, right, full and inner 
• mutate(): add a new column to existing data





"""


msleep <- read.csv('msleep-200908-125135.csv') 
head(msleep)

library(dplyr)


sleeptime <- select(msleep, ï..name,sleep_total)

msleep

sleeptime

head(sleeptime)


filter(msleep, order == 'Carnivora')

groupss <- group_by(msleep, genus)
groupss


levels(groupss$genus)

mtcars

hmtcard <- head(mtcars)
hmtcard

nc_hmtcar <- mutate(hmtcard, mpg_cyl = mean(cyl))
nc_hmtcar


sdnc_hmtcar <- mutate(hmtcard, mpg_cyl = sd(cyl))
sdnc_hmtcar


sdnc_hmtcar <- mutate(hmtcard, mpg_cyl = 0)
sdnc_hmtcar

head(mutate(mtcars, disp_l = disp / 61.0237), 3)


x <- 10
ifelse(x > 11, "x is greater than 9", "x is not greater than 9")


section <- c("MATH111", "MATH111", "ENG111")
grade <- c(78, 93, 56)
student <- c("David", "Kristina", "Mycroft")
gradebook <- data.frame(section, grade, student)
gradebook


mutate(gradebook, Result = ifelse(grade>60, "Pass", "Fail"))



mutate(gradebook, letter = ifelse(grade %in% 60:69, "D",
                                  ifelse(grade %in% 70:79, "C",
                                         ifelse(grade %in% 80:89, "B",
                                                ifelse(grade %in% 90:99, "A", "F")))))


grepl("MATH", gradebook$section)



mutate(gradebook, department = ifelse(grepl("MATH", section), "Math Department",
                                      ifelse(grepl("ENG", section), "English Department", "Other")))



1
2
3
4
5
title <- c('Data Smart','Orientalism','False Impressions','The Age of Wrath','Making Software')
author <- c('Foreman, John','Said, Edward','Archer, Jeffery','Eraly, Abraham','Oram, Andy')
height <- c('235','197','177','238','235')
year <- c('2010','2011','2012','1999','1998')
bookDF <-  data.frame(title, author, height, year)




bookDF



sortedBookDF <- bookDF[order(title),]
sortedBookDF


sortedBookDF <- bookDF[order(bookDF$year, decreasing = TRUE),]
sortedBookDF


sortedBookDF <- bookDF[order(bookDF$year, bookDF$height),]
sortedBookDF
library(plyr)




orderedBookDF <- arrange(bookDF, title)
orderedBookDF


orderedBookDF <- arrange(bookDF, desc(title))
orderedBookDF


1
2
3
4
5
6
7
id <- c('S1002','S1003','S1005','S1008','S1011','S1015')
age <- c(58, 67,64, 34, 30, 37)
gender <- c('female','male','male','male','male','male')
height <- c(61, 67, 68, 71, 69, 59)
weight <- c(256, 119, 183, 190, 191, 170)
# Create a data frame
subjectDfrm <-  data.frame(id, age, gender, weight)
subjectDfrm



sortedSubjectDfrm <- arrange(subjectDfrm, id)
sortedSubjectDfrm


subjectDfrm %>% arrange(desc(age)) %>% head



sortedSubjectDfrm <- arrange(subjectDfrm, id, age, weight)
sortedSubjectDfrm

# Or with the pipe operator:


subjectDfrm %>% arrange(id, age, weight) %>% head



id <- c('S1022','S1023','S1024','S1025','S1026','S1027')
glyhb <- c(4.64, 4.63, 7.72, 4.81, 4.84, 3.94)
chol <- c(132, 228, 228, 181, 249, 248)
frame <- c('large','large','medium','small','small','medium')
# Create a data frame
testDfrm <-  data.frame(id, glyhb, chol, frame)
testDfrm



install.packages("doBy")
library(doBy)

orderBy(~chol, data = testDfrm) #Ascending order 
orderBy(~-chol, data=testDfrm) #descending order 


orderBy(~-chol+id, data=testDfrm)
orderBy(~-chol-id, data=testDfrm)



head(mtcars)

mtcars %>% group_by(cyl) 

mtcars %>% group_by(cyl) %>% summarise(mean = mean(disp))

a <- c(1,2,3
      )
as.character(a)
z<-list(1,2,c("a","b"))
y<- unlist(z)
y
a<- c(3,3,6.5,8)
b<- c(7,2,5.5,10)
a<b
c<- a<b
c


as.data.frame()

unclass(as.Date("1971-01-01"))
library(wordcloud)

 camdata <- read.csv("camera_survey.csv")
camdata
names(camdata)


camdata$Variable_name

wd <- camdata$Variable_name
wd

library(wordcloud)

library(tm)




id <- c('S1022','S1023','S1024','S1025','S1026','S1027')
glyhb <- c(4.64, 4.63, 7.72, 4.81, 4.84, 3.94)
chol <- c(132, 228, 228, 181, 249, 248)
frame <- c('large','large','medium','small','small','medium')


sentence_data = data.frame(id,glyhb,chol,frame)

sentence_data

library(tm)

sentence_data$frame

wordcloud(sentence_data)

wordcloud(sentence_data$frame,
          scale=c(2,0.5), random.order=FALSE,
          rot.per=0.35, use.r.layout=FALSE,
          colors=brewer.pal(8, "Dark2"))



"""

## Section 1: Manipulating Data

### Question 1 
**The dataset is loaded for you. Perform the following tasks:**  
1. use the USPersonalExpenditure dataset and store in a variable.  
2. Compute summary statistics: mean, standard deviation (sd) across 
time for each of the indicators.  
3. Create a new column that contains average values for each of the indicators.  
"""

df<-USPersonalExpenditure
msleep_data<-as.data.frame(df)
msleep_data
  
library(dplyr)
library(tidyr)

sd<-summarise_at(msleep_data, vars("1940", "1945", "1955", "1960"), funs(mean, sd), na.rm = T)
sd

mutate(msleep_data,average=rowMeans(data))


mutate(msleep_data, mean_col = summarise_at(msleep_data, vars("1940", "1945", "1955", "1960"), funs(mean), na.rm = T) )


vardata <- datasets::USPersonalExpenditure 
vardata  <- as.data.frame(vardata)
vardata

dplyr::summarize()



summarize(vardata, avg_1940 = mean(vardata$`1940`))

summarise_at(vardata, vars(vardata$`1940`,vardata$`1945`,vardata$`1950`,
                           vardata$`1955`,vardata$`1960`), funs(mean, sd), na.rm = T)


summarise_at(vardata, vars("1940","1945","1950","1955","1960"), funs(mean, sd), na.rm = T)


mutate(vardata, mean_col = summarise_at(vardata, vars("1940","1945","1950","1955","1960"), funs(mean), na.rm = T))


"""
### Question 2 
**download the data from the available URL:**  
  1. Create a new column containing the average bodywt (body weight) of
  each genus and of each order.  
2. Create a dataframe containing the average REM sleep time for each order.  
3. How many missing values are there (in total) and per column?  
  4. How would you like to impute missing values? 
  Write Justification.  <i> Hint: Overall Mean/media/mode vs. Groupby 
  Mean/media/mode?

url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/msleep_ggplot2.csv"
filename <- "msleep_ggplot2.csv"
  

"""

library(tidyr)

unclass(url)

sleepdata <- read.csv(url)
sleepdata

head(sleepdata)
colnames(sleepdata)


mutate(sleepdata, average_bodywt = bodywt*2,na.rm=T) %>%
  select(genus,order,average_bodywt) %>% head()


q2 <- data.frame(sleepdata %>% group_by(genus,order) %>% 
                   mutate(rem_sleep = mean(sleep_rem, na.rm = T))
                          %>% select(order,genus,sleep_rem,rem_sleep))
q2

colSums(is.na(sleepdata))

"""

### Question 1
**Use the above dataset and perform the following tasks using any library
from tidyverse:**  
1. Filter results to print average REM sleep and average total sleep  
for those animals who are carnivores and then for those who are primates.  
2. Use the order column and "spread" it across the rest of the observations.  


"""

Results <- sleepdata %>% select(name, order, sleep_total, sleep_rem) %>% filter(order %in% c('Carnivora', 'Primates')) %>% group_by(order) %>% mutate(Avg_REM_Sleep= mean(sleep_rem, na.rm = TRUE),  
           Total_Avg_Sleep= mean(sleep_total, na.rm = TRUE))
Results



spread(sleepdata,order,name)
spread(sleepdata,order,genus)
spread(sleepdata,order,vore)
spread(sleepdata,order,conservation)
spread(sleepdata,order,sleep_total)
spread(sleepdata,order,sleep_rem)
spread(sleepdata,order,sleep_cycle)
spread(sleepdata,order,awake)
spread(sleepdata,order,brainwt)
spread(sleepdata,order,bodywt)



for (name in names(file)){
  spread(sleepdata,order,name)
}

for (name in names(sleepdata)){
  print(paste(name))
}
