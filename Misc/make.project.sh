#!/bin/sh

##########################
#Shell script for initiating a project repository structure.
##########################
#SYNTAX:
#sh make.project.sh <project_name>
#with:
#<project_name> being the name of the project...
##########################
#----
#guillert(at)tcd.ie - 2018/05/23
##########################

## Input values
PROJECT_NAME=$1

## Making the project tree
mkdir ${PROJECT_NAME}
#|
#\- Analysis/
    mkdir ${PROJECT_NAME}/Analysis
#|    
#\- Data/
    mkdir ${PROJECT_NAME}/Data
#   |    
#   \- Data/Raw
        mkdir ${PROJECT_NAME}/Data/Raw
#   |    
#   \- Data/Processed
        mkdir ${PROJECT_NAME}/Data/Processed
#   |    
#   \- Data/Results
        mkdir ${PROJECT_NAME}/Data/Results
#|    
#\- Functions/
    mkdir ${PROJECT_NAME}/Functions
#|    
#\- Manuscript/
    mkdir ${PROJECT_NAME}/Manuscript
#   |    
#   \- Manuscript/Figures
        mkdir ${PROJECT_NAME}/Manuscript/Figures
#   |    
#   \- Manuscript/Tables
        mkdir ${PROJECT_NAME}/Manuscript/Tables

## Adding the general files

## README
echo "# ${PROJECT_NAME}" > ${PROJECT_NAME}/README.md

## .gitignore
echo "## Compiled ignore" > ${PROJECT_NAME}/.gitignore
echo "*.pdf" > ${PROJECT_NAME}/.gitignore
echo "*.html" > ${PROJECT_NAME}/.gitignore
echo "" > ${PROJECT_NAME}/.gitignore
echo "## R ignore" > ${PROJECT_NAME}/.gitignore
echo ".Rhistory" >> ${PROJECT_NAME}/.gitignore
echo "*.Rda" >> ${PROJECT_NAME}/.gitignore
echo "*.rda" >> ${PROJECT_NAME}/.gitignore
echo "*.RData" >> ${PROJECT_NAME}/.gitignore
echo "*.Rapp.history" >> ${PROJECT_NAME}/.gitignore

## Initialising the project on git
git init ${PROJECT_NAME}
