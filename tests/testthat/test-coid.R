test_that("Determine the type of parameters", {


   expect_error(coid(c(1:10)),
                "Please check the type of 'network.or.subnet_mat1' or 'subnet_mat2'")
   expect_error(coid(c("a", "b", "c")),
                "Please check the type of 'network.or.subnet_mat1' or 'subnet_mat2'")
   expect_error(coid(c(T, F, F, T, NA)),
                "Please check the type of 'network.or.subnet_mat1' or 'subnet_mat2'")
   expect_error(coid(matrix(1:10,2,5),c("a", "b", "c")),
                "Please check the type of 'network.or.subnet_mat1' or 'subnet_mat2'")
   m1<-matrix(1:10,5,2)
   rownames(m1)<-paste0("species",seq=1:5)
   m3<-matrix(1:12,6,2)
   rownames(m3)<-c(paste0("species",seq=1:5),NA)
   m4<-matrix(1:10,5,2)
   rownames(m4)<-c(paste0("species",seq=1:4),NA)
   m5<-matrix(1:15,5,3)
   rownames(m5)<-c(paste0("species",seq=1:4),NA)
   m6<-matrix(1:18,6,3)
   rownames(m6)<-paste0("species",seq=6:11)
   m7<-matrix(1:18,6,3)
   rownames(m7)<-c(paste0("species",seq=c(1,3,2,5,4)),NA)
   expect_error(coid(m1,m7),
                "Please make sure matrices either have no row names or have full row names. No NA!")
   expect_error(coid(m1,m6),
                "No connectors existed.")
   expect_error(coid(m4,m5),
                "Please make sure matrices either have no row names or have full row names. No NA!")
   expect_error(coid(m3,m7),
                "Please make sure matrices either have no row names or have full row names. No NA!")
})


test_that("Input a proper network data", {
   ma<-Multi_motif("all")
   for(i in c(26,31)){
      expect_true(class(coid(ma[[i]]))=="numeric")
   }
   set.seed(1)
   MAT <- build_toy_net(11,22,21,0.2,output_matrices=TRUE)
   expect_error(coid(t(MAT[[3]]),t(MAT[[4]])), "No connectors existed.")

   MA<-build_toy_net(5,3,3,0.9)
   expect_identical(class(coid(MA)), "numeric")

   m8<-matrix(1:15,5,3)
   m9<-matrix(1:16,4,4)
   expect_true(!is.na(coid(m8,m9,TRUE)), TRUE)

   m10<-matrix((1:16)/2,4,4)
   m11<-matrix(1:10,5,2)
   expect_true(is.numeric(coid(m10,m11,TRUE)), TRUE)

})
