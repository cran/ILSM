test_that("Determine whether a network can be adapted into an interconnected network", {
   m1<-matrix(c(0,0,0,1,0,0),3,2)
   m2<-matrix(c(0,0,0,2,0,0),2,3)
   N<-trigraph_from_mat(m1,m2)
   expect_error(icmotif_count(N),
                "No connector nodes.")
})
