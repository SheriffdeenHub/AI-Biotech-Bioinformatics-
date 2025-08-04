#Creating a working directory
setwd("C:\\Users\\sheri\\Documents\\R\\AI & Biotech\\Bioinformatics")

#Confirm the working directory
getwd()

#Create and arrange subfolders 
dir.create("Data")
dir.create("Scripts")
dir.create("Results")

#Load the patient_info data
Data <- read.csv(file.choose())
Data

#Save the csv file in Data directory of the project directory
save(Data, file = "Data/patient_info.csv")

#View Data in spreadsheet form
View(Data)

#Check structure of Data
str(Data)

#Convert 'patient_id' column to factor
Data$patient_id <- as.factor(Data$patient_id)
str(Data)

#Convert 'gender' column to factor
Data$gender <- as.factor(Data$gender)
str(Data)

#Convert 'gender' factor to numeric
Data$gender <- ifelse(Data$gender == "Female", 1, 0)
str(Data)

#COnvert gender numeric to factor
Data$gender <- as.factor(Data$gender)
str(Data)

#Convert 'diagnosis' column to factor
Data$diagnosis <- as.factor(Data$diagnosis)
str(Data)

#Reorder 'diagnosis' factor levels
Data$diagnosis <- factor(Data$diagnosis, levels = c("Cancer","Normal"))
levels(Data$diagnosis)

str(Data)

#Convert 'diagnosis' factor to numeric
Data$diagnosis <- ifelse(Data$diagnosis == "Cancer", 1, 0)
str(Data)

unique(Data$diagnosis)

#Convert 'diagnosis' numeric to factor
Data$diagnosis <- as.factor(Data$diagnosis)
str(Data)

#Convert 'smoker' column to factor
Data$smoker <- as.factor(Data$smoker)
str(Data)

#Save file as csv in results folder
write.csv(Data, file = "Results/Cleaned_Patient_Info.csv")

#Save the entire R workspace
save.image(file = "Hammed_Sheriffdeen_Class_1b_Assignment.RData")


