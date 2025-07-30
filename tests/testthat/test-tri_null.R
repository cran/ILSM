test_that("Determine the type of null model", {


   ma<-Multi_motif("all")
   for(i in 23:31){
      expect_error(tri_null(ma[[i]],null_type = "a"),
                   "Wrong input for null_type")
   }
   set.seed(1)
   MA<-build_toy_net(5,3,3,0.9)
   expect_true(class(tri_null(MA,null_type="sauve")[[1]])=="igraph")

})





test_that("Make sure the function is implemented", {
   m1<-matrix(sample(c(rep(1,9),rep(0,1))),5,2)
   rownames(m1)<-paste0("species",seq=1:5)
   m2<-matrix(sample(c(rep(1,13),rep(0,2))),5,3)
   rownames(m2)<-c(paste0("species",seq=c(1,3,2,5,4)))
   N<-trigraph_from_mat(m1,m2)
   expect_identical(class(N),
                    c("igraph"))
   expect_identical(length(tri_null(N,null_type="sauve")),
                    100L)
   set.seed(2)
   M<-V(tri_null(N,null_type="sub_P",sub_method ="r00")[[1]])$level
   set.seed(2)
   M1<-V(tri_null(N,null_type="both_sub",sub_method ="r1")[[1]])$level
   expect_identical(M,M1)
   expect_identical(length(V(tri_null(N,null_type="sub_P",sub_method ="r2")[[1]])$name),10L)
   expect_identical(length(V(tri_null(N,null_N  = 4,null_type="sub_Q",sub_method ="r00")[[1]])$name),10L)
   expect_identical(length(V(tri_null(N,null_type="sauve")[[1]])$name),10L)
   expect_identical(length(tri_null(N,null_N  = 4,null_type="sauve")),
                    4L)
})
