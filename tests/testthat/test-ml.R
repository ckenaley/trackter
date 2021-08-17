test_that("free.ml.ang works fine", {
  
  cont <- read.csv(system.file("extdata", "cont.csv", package = "trackter"))[,3:4]
  cont <- cont[seq(1,nrow(cont),length.out = 500),]
  fml <- free.ml.ang(as.matrix(cont),smooth.n = 5)
  
  expect_error(free.ml.ang(cont,smooth.n = 5),"must be a matrix")
  
  expect_error(free.ml.ang(as.matrix(cont),smooth.n = 5,red=10),"must be numeric and 0-1")
  
  expect_named(fml,c("ml","cont.sm","cont.sides"))
  
  expect_error(free.ml.ang(as.matrix(cont),smooth.n = 5,red=10,dens=2),"both 'red' and 'dens' are not NULL")

  
})

test_that("free.ml.hull works fine", {
  
  cont <- read.csv(system.file("extdata", "cont.csv", package = "trackter"))[,3:4]
  cont <- cont[seq(1,nrow(cont),length.out = 500),]
  fml <- free.ml.hull(as.matrix(cont),smooth.n = 5)
  
  expect_error(free.ml.hull(cont,smooth.n = 5),"must be a matrix")
  
  expect_error(free.ml.hull(as.matrix(cont),smooth.n = 5,red=10),"must be numeric and 0-1")
  
  expect_named(fml,c("ml","cont.sm","cont.sides"))
  
  expect_error(free.ml.hull(as.matrix(cont),smooth.n = 5,red=10,dens=2),"both 'red' and 'dens' are not NULL")
  
  
})

test_that("free.ml.del works fine", {
  
  cont <- read.csv(system.file("extdata", "cont.csv", package = "trackter"))[,3:4]
  cont <- cont[seq(1,nrow(cont),length.out = 200),]
  fml <- free.ml.del(as.matrix(cont),smooth.n = 5,red=0.5)
  
  expect_error(free.ml.del(cont,smooth.n = 5),"must be a matrix")
  
  expect_error(free.ml.del(as.matrix(cont),smooth.n = 5,red=10),"must be numeric and 0-1")
  
  expect_named(fml,c("ml","cont.sm","cont.sides"))
  
  fml2 <- free.ml.del(as.matrix(cont),smooth.n = 5,dens=1)
  
  expect_true(nrow(fml$ml)!=nrow(fml2$ml))
  
})
