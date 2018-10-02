#' @title shufflR
#'
#' @description A function for assigning poster to students (with an awesome pun in the name!)
#'
#' @param n_students the number students
#' @param n_posters the number posters per students
#' 
#' @examples
#' ## The basics:
#' shufflR(36, 5)
#'
#' ## Your names
#' vector_of_names <- letters
#' 
#' ## The matrix
#' matrix <- shufflR(26, 5)
#' 
#' ## The backup of the rownames
#' rownames_backup <- rownames(matrix)
#' 
#' ## The loop (in reverse order to avoid having problem with replace the 11th student by "student1student1")
#' for(name in rev(1:length(vector_of_names))) {
#'     ## Replacing in the matrix
#'     matrix <- apply(matrix, 2, function(x) gsub(name, vector_of_names[name], x))
#'     ## Replacing in the row names
#'     rownames_backup <- gsub(name, vector_of_names[name], rownames_backup)
#' }
#' 
#' ## Adding the rownames
#' rownames(matrix) <- rownames_backup
#' 
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

shufflR <- function(n_students, n_posters) {

    ## Create a neat an ordinated matrix
    nice_matrix <- matrix(1:n_students, n_students, n_posters, byrow = TRUE)

    ## Randomise the names of the students in the column
    rownames(nice_matrix) <- sample(1:n_students, n_students)

    ## A function for selecting the self assigned students
    find.self.assigned <- function(matrix) {
        ## Creating a vector of rows to shuffle
        rows_to_shuffle <- NULL

        for(i in 1:nrow(matrix)) {
            ## Checking whether the student "ID" is present in the row
            if(as.numeric(rownames(matrix)[i]) %in% matrix[i,]) {
                ## Adding the wrong row to the vector
                rows_to_shuffle <- c(rows_to_shuffle, i)
            }
        }

        return(rows_to_shuffle)
    }

    ## A function for reshuffling the wrong rows
    reshuffle.self.assigned <- function(matrix, shuffle) {
        ## Apply the reshuffling if there is more than one student to reshuffle
        if(length(shuffle) != 1){
            ## Save the row names (to keep them fixed)
            rownames_save <- rownames(matrix)
            ## Shuffle!
            matrix[shuffle, ] <- matrix[sample(shuffle, length(shuffle)), ]
            ## Keep the same row names as before
            rownames(matrix) <- rownames_save

        } else {
            ## Simply swap the row to reshuffle
            
            ## Save the row names (to keep them fixed)
            rownames_save <- rownames(matrix)

            ## If it is the last row, replace by the one above
            if(shuffle == nrow(matrix)) {
                matrix[c(shuffle - 1, shuffle), ] <- matrix[c(shuffle, shuffle - 1), ]
            } else {
                ## Else replace by the one below
                matrix[c(shuffle, shuffle + 1), ] <- matrix[c(shuffle + 1, shuffle), ]
            }
            ## Keep the same row names as before
            rownames(matrix) <- rownames_save
        }

        return(matrix)
    }

    ## With these two functions we can then use a while loop: while the rows to reshuffle is not null,
    ## Do the same find.self.assigned and reshuffle.self.assigned functions:

    ## Checking if any student is self assigned
    rows_to_shuffle <- find.self.assigned(nice_matrix)

    while(!is.null(rows_to_shuffle)) {
        ## Shuffle the wrong rows
        nice_matrix <- reshuffle.self.assigned(nice_matrix, rows_to_shuffle)

        ## Re-check for self.assigned
        rows_to_shuffle <- find.self.assigned(nice_matrix)

        ## If that was NULL, it'll exit the loop
    }

    return(nice_matrix)
}
