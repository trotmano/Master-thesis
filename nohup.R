##Submitting jobs automatically

print(paste0("kill ", Sys.getpid())) # manual for killing the job if necessary

#########################################
# add-on
command_2 = paste("qstat -u bq_sbaz research") # !!! very important: no typo! used to get your current number of jobs on the cluster
limit = 300 # max number of submitted jobs in the queue (for batch reasonable = 30, for research reasonable = 300)
start_it = 1
len_mz = 6
len_allPairs = 8*7
# len_cond = 12*10
len_cond = 6*11
end_it = len_allPairs * len_mz * len_cond + 1 # jobs 1 to 100 are executed if end_it = 101
# end_it = 3*4*3+1 # jobs 1 to 100 are executed if end_it = 101

# where do you want to save individualized R scripts & corresponding shell files (you need it once, after that it is very much trash)?
prefix=paste0("")
suffix="_Ori.R" # suffix of your individualized R scripts
suffix_shell="_Ori.sh" # suffix of your individualized shell files
ori = paste0("") # name of the original file that serves as template for individualized scripts

# create the individual job files

create.jobs <- function(seed_par){
  
  startvar<-matrix(data=NA, nrow=1,ncol=1)
  startvar[1,1]<-paste("ind <- c(",seed_par,")",sep="")
  
  write.table(startvar,
              paste(prefix, 
                    seed_par,
                    suffix, sep=""),
              row.names=FALSE,
              col.names=FALSE,
              quote=FALSE)
  
}

append.files<- function(seed_par){
  file.append( paste(prefix, 
                     seed_par,
                     suffix, sep=""),
               ori)
  
}

create.submission <- function(seed_par){
  globvar <- matrix(data=NA,nrow=3,ncol=1)
  globvar[1,1] <- "#PBS -j oe"
  globvar[2,1] <- "#PBS -o /dev/null"
  globvar[3,1] <- paste0("Rscript --verbose ",prefix,
                         seed_par,
                         suffix, " > ",prefix, 
                         seed_par,
                         suffix,"out 2>&1")
  write.table(globvar,
              paste0(prefix, seed_par, suffix_shell),
              row.names=FALSE,
              col.names=FALSE,
              quote=FALSE)
}

#submit jobs
submit.jobs <- function(seed_par){
  file <- paste0(prefix, seed_par, suffix_shell)
  print("Submitting")
  command <- paste("qsub -l nodes=1:avail,cput=04:59:00,mem=2gb ", file, sep = "") # for batch make it larger than 4:59:00; for research make it larger than 1:00:00 and smaller than 4:59:00
  print(command)
  system(command)
}

# Main function to create corefiles, mainfiles, jobfiles, submissionfiles and finally to submit the jobs

create.and.submit <- function(seed_par){
  create.jobs(seed_par)
  append.files(seed_par)
  create.submission(seed_par)
  submit.jobs(seed_par)
}

while(TRUE) {
  if (length(system(command_2,intern=TRUE))==0) {
    count = limit
  }
  else {
    count = limit-(length(system(command_2,intern=TRUE))-5) # -5 bc first five lines are not listed jobs
  }
  
  if (count > 0) {
    for(i in start_it:(start_it-1+count)){ 
      if (i >= end_it) {break}
      create.and.submit(seed_par=i)
    }
  }
  Sys.sleep(30) # in seconds
  start_it = i+1
  if (i >= end_it) {break}
}
