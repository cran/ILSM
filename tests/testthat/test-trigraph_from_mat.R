test_that("Determine the type of matrices of inputing", {
   m1<-matrix(1:10,5,2)
   rownames(m1)<-paste0("species",seq=1:5)
   m2<-matrix(1:12,4,3)
   m3<-matrix(1:12,6,2)
   rownames(m3)<-c(paste0("species",seq=1:5),NA)
   m4<-matrix(1:10,5,2)
   rownames(m4)<-c(paste0("species",seq=1:4),NA)
   m5<-matrix(1:15,5,3)
   rownames(m5)<-c(paste0("species",seq=1:4),NA)

   m7<-matrix(1:18,6,3)
   rownames(m7)<-c(paste0("species",seq=c(1,3,2,5,4)),NA)
   m8<-matrix(1:18,6,3)
   rownames(m8)<-c(paste0("species",seq=c(1,3,2,4)),NA,NA)
   expect_message(trigraph_from_mat(m1,m2),
                  "Warning! Row IDs were set as rownames for matching connector nodes since no rownames are provided for the matrices")
   expect_error(trigraph_from_mat(m1,m7),
                "Please make sure the two matrices have appropriate row names. NA is not accepted.")
   expect_error(trigraph_from_mat(m4,m5),
                "Please make sure the two matrices have appropriate row names. NA is not accepted.")
   expect_error(trigraph_from_mat(m3,m7),
                "Please make sure the two matrices have appropriate row names. NA is not accepted.")
   expect_error(trigraph_from_mat(m3,m8),
                "Please make sure the two matrices have appropriate row names. NA is not accepted.")
})


test_that("Make sure the function is implemented", {
   set.seed(1)
   m1<-matrix(sample(c(rep(1,7),rep(0,3))),5,2)
   rownames(m1)<-paste0("species",seq=1:5)
   m2<-matrix(sample(c(rep(1,9),rep(0,6))),5,3)
   rownames(m2)<-c(paste0("species",seq=c(1,3,2,5,4)))
   m3<-matrix(sample(c(rep(1,8),rep(0,4))),6,2)
   m4<-matrix(sample(c(rep(1,9),rep(0,6))),5,3)
   N<-trigraph_from_mat(m1,m2)
   M<-trigraph_from_mat(m3,m4)
   m<-trigraph_from_mat(m4,m3)
   expect_identical(class(N),
                    "igraph")
   expect_identical(length(N),
                    10L)
   expect_identical(length(M),
                    11L)
   expect_true(length(M)==length(m))
   expect_length(icmotif_count(N)[,2],
                    48L)
   expect_true(icmotif_count(M)[1,2]==icmotif_count(m)[1,2])
})
