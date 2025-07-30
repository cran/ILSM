test_that("Build a large network", {
  expect_error(build_toy_net(4,2,3,0.5),
               "Error: please make sure N_a>=3, N_b>=3 and N_c>=3!!!")
   expect_error(build_toy_net(4,2,3,0.5,output_matrices=TRUE),
                "Error: please make sure N_a>=3, N_b>=3 and N_c>=3!!!")
   expect_error(build_toy_net(2,4,3,0.5),
                "Error: please make sure N_a>=3, N_b>=3 and N_c>=3!!!")
   expect_error(build_toy_net(4,3,2,0.5),
                "Error: please make sure N_a>=3, N_b>=3 and N_c>=3!!!")
   expect_error(build_toy_net(2,2,2,0.5),
                "Error: please make sure N_a>=3, N_b>=3 and N_c>=3!!!")
   expect_equal(length(build_toy_net(5,5,5,0.5)),
                15)
   set.seed(1)
   N<-build_toy_net(5,5,5,0.5,output_matrices = TRUE)
   expect_named(N,c( "network","supraadjacency_matrix","subnet_P",
                     "subnet_Q" ))
})

test_that("Ensure the accuracy of the conversion into a matrix", {

   set.seed(3)
   M<-build_toy_net(5,5,5,0.5)
   set.seed(3)
   N<-build_toy_net(5,5,5,0.5,TRUE)
   expect_identical(class(M),"igraph")
   expect_identical(class(N[[1]]),"igraph")
   expect_identical(class(N[[2]]),c("matrix","array"))
   expect_identical(class(N[[3]]),c("matrix","array"))
   expect_identical(class(N[[4]]),c("matrix","array"))
   expect_true(sum(!(as.matrix(M)==as.matrix(N[[1]])))==0)
   expect_identical(sum(!(as.matrix(N[[1]])==N[[2]]))==0,TRUE)
   expect_identical(sum(rownames(N[[3]])==rownames(N[[4]])),5L)
   expect_identical(sum(sort(rownames(N[[2]]))==sort(c(rownames(N[[3]]),
                                                       colnames(N[[3]]),
                                                       colnames(N[[4]])))),15L)
})

