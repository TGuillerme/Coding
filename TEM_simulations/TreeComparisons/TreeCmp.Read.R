##########################
#Chain TreeCmp Read
##########################
#Combine the data of .Cmp files from TEM_TreeCmp.sh into R by chains
#v.0.2
#Update: improved input format
#Update: output is now a 'TreeCmp' list class
##########################
#SYNTAX :
#<chain> the chain name to read in
#<suffix> the file chain suffix (default =.Cmp)
#<sep> the separator between the chain name and the constant part of the file name (default=_)
#<header> whether the tables to read have a header or not (default=TRUE)
#<header> whether to be verbose or not (default=FALSE) - useful for large data sets handling
##########################
#----
#guillert(at)tcd.ie - 20/03/2014
##########################
#Requirements:
#-R 3
##########################

TreeCmp.Read<-function(chain, suffix='.Cmp', sep='_', header=TRUE, verbose=FALSE) {

#DATA INPUT
    #chain
    if(class(chain) !='character') {
        stop('No file has been found with the given chain name')
    } else {
        chain.list<-list.files(pattern=chain)
        if(length(chain.list) == 0) {
            stop('No file has been found with the given chain name')
        }
    }

    #suffix
    if(class(suffix) !='character') {
        stop('Suffix not found')
    } else {
        names.list<-unlist(strsplit(chain.list, suffix))
        if(length(names.list) == 0) {
            stop('Suffix not found')
        }
    }

    #separator
    if(class(sep) !='character') {
        stop('Separator not found')
    } else {
        elements.list<-unlist(strsplit(names.list, sep))
        unique.list<-unique(elements.list[-c(grep(chain, elements.list))])
        if(length(unique.list) == 0) {
            stop('Separator not found')
        }
    }

    #header
    if(class(header) != 'logical') {
        stop('Header must be logical')
    }

    #verbose
    if(class(verbose) != 'logical') {
        stop('Verbose must be logical')
    }

#FUNCTION

FUN.TreeCmp.Read<-function(chain, suffix, sep, header, verbose) {
    #Combine the data frames according to unique.list

        #Creates the list of empty tables
        combined.list<-NULL
        for (i in 1:length(unique.list)){
            combined.list[[i]]<-data.frame(NULL)
        }

        #Renaming
        names(combined.list)<-unique.list

        #Creates the list of files to put in the list
        for (j in 1:length(unique.list)) {
            table.list<-grep(unique.list[j], chain.list)

            #Building the list
            for (i in 1:length(table.list)){
                combined.list[[j]]<-rbind(combined.list[[j]], read.table(chain.list[table.list[i]], header=header))
            }

            #Be verbose
            if (verbose == TRUE) {
                cat(format(Sys.time(), '%H:%M:%S'), '-', unique.list[j], 'combined', '\n')
            } 
        }

        return(combined.list)
    }

#COMBINING THE DATA SETS

    combined.list<-FUN.TreeCmp.Read(chain, suffix, sep, header, verbose)
    class(combined.list)<-'TreeCmp'
    return(combined.list)    

#End
}