
context("kin functions")

test_that("kin.simple works fine", {
  y <- EBImage::readImage(system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter"))
  t <-tempdir()

  ti <-paste0(tempdir(),"/images")
  tp <- paste0(tempdir(),"/processed_images")
  dir.create(ti)
  dir.create(tp)
  
  EBImage::writeImage(y,paste0(ti,"/sunfish001.jpg"),type = "jpeg")
  
  invisible(capture.output( kin.y <- kin.simple(image.dir = ti,save = TRUE,out.dir =tp)))
  
  
  expect_length(list.files(tp),1)
  
  expect_is(kin.y,"list")
  expect_named(kin.y,c("kin.dat", "midline","cont","cont.sm","mid.pred","all.classes","dim"))
  expect_true(kin.y$kin.dat$head.pval!=0)
  expect_type(kin.y$midline$roi,type = "character")
  expect_type(kin.y$cont$x,type = "integer")
  expect_true(kin.y$all.classes$size>0)
  
  
  expect_error(invisible(capture.output( kin.y <- kin.simple(image.dir = ti,save = TRUE))),"'out.dir' not specified")
  
 dir.create(paste0(t,"/test_images2"))
  
  expect_error(invisible(capture.output( kin.simple(image.dir = paste0(t,"/test_images2"),save=FALSE))),"no images in image.dir")
  
  unlink(paste0(t,"/test_images2"),recursive = TRUE)
  
  expect_error(invisible(capture.output( kin.simple(image.dir = "foo",out.dir=tp,save=TRUE))),"does not exist")
  expect_error(invisible(capture.output( kin.simple(image.dir = ti,out.dir="foo",save=TRUE))),"does not exist")
  
  expect_error(invisible(capture.output(kin.simple(image.dir =ti ,save=TRUE))),"not specified")
  
  expect_error(invisible(capture.output( kin.simple(image.dir =ti,frames=2,save=FALSE))),"out of range")
  
  expect_error(invisible(capture.output( kin.simple(image.dir =ti ,ml.smooth=list(0.5,"foo"),save=FALSE))),"'ml.smooth' must")
  
  
  unlink(ti,recursive = TRUE)
  unlink(tp,recursive = TRUE)

  
})

test_that("kin.search works fine", {
  
  y <- EBImage::readImage(system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter"))
  t <-tempdir()
  ti <-paste0(tempdir(),"/images")
  tp <- paste0(tempdir(),"/processed_images")
  dir.create(ti)
  dir.create(tp)
  
  EBImage::writeImage(y,paste0(ti,"/sunfish001.jpg"),type = "jpeg")

  invisible(capture.output( kin.y <- kin.search(image.dir =ti,save = TRUE,out.dir =tp)))
  
 
  expect_length(list.files(tp),1)
  
  expect_is(kin.y,"list")
  expect_named(kin.y,c("kin.dat", "midline","cont","cont.sm","mid.pred","all.classes","dim"))
  expect_true(kin.y$all.classes$size>0)
  expect_type(kin.y$midline$roi,type = "character")
  expect_type(kin.y$cont$x,type = "integer")
  

  expect_error(invisible(capture.output( kin.y <- kin.search(image.dir = ti,save = TRUE))),"'out.dir' not specified")
  
  
  expect_error(invisible(capture.output( kin.search(image.dir = "foo",out.dir=tp,save=TRUE))),"does not exist")
  expect_error(invisible(capture.output( kin.search(image.dir = ti,out.dir="foo",save=TRUE))),"does not exist")
  
  expect_error(invisible(capture.output(kin.search(image.dir =ti ,save=TRUE))),"not specified")
  
  expect_error(invisible(capture.output( kin.search(image.dir = ti,frames=2,save=FALSE))),"out of range")
  
  expect_error(invisible(capture.output( kin.search(image.dir =ti , thr="foo",save=FALSE))),"must be set to")
  expect_error(invisible(capture.output( kin.search(image.dir =ti ,ml.smooth=list("foo",0.6),save=FALSE))),"'ml.smooth' must")
  
  expect_error(invisible(capture.output( kin.search(image.dir =ti ,search.for="foo",save=FALSE))),"must be set to")
  
  expect_error(invisible(capture.output( kin.search(image.dir =ti,save=TRUE,out.qual=1.1,paste0(t,"/test_images")))),"'out.qual' must be >=0 and <=1")
  
  
  dir.create(paste0(t,"/test_images2"))
  
  expect_error(invisible(capture.output( kin.search(image.dir = paste0(t,"/test_images2"),save=FALSE))),"no images in image.dir")
  
  unlink(paste0(t,"/test_images2"),recursive = TRUE)
  
  unlink(ti,recursive = TRUE)
  unlink(tp,recursive = TRUE)

  
})


test_that("kin.free works fine", {
  
  y <- list.files(system.file("extdata/img", package = "trackter"),full.names = TRUE)
  y <- y[grepl("lamp",y)]
  ti <-paste0(tempdir(),"/images")
  tp <- paste0(tempdir(),"/processed_images")
  dir.create(ti)
  dir.create(tp)
  
  file.copy(y,paste0(ti,"/",basename(y)))
  
  invisible(capture.output( kin.y <- kin.free(image.dir = ti,save = TRUE,out.dir =tp,red=0.5)))
  
  expect_length(list.files(tp),2)
  expect_length(kin.y$all.classes$size,2)
  expect_is(kin.y,"list")
  expect_named(kin.y,c('kin.dat', 'midline', 'cont', 'cont.sm', 'all.classes', 'mid.pred', 'dim'))
  expect_true(kin.y$all.classes$size[1]>0)
  expect_type(kin.y$midline$roi,type = "character")
  expect_type(kin.y$cont$x,type = "integer")

  
  expect_error(invisible(capture.output( kin.free(image.dir = ti,out.dir=tp,save = FALSE))),"To save processed images")
  
  
  expect_error(invisible(capture.output( kin.free(image.dir ="foo",out.dir=tp,save=TRUE))),"does not exist")
  
  expect_error(invisible(capture.output( kin.free(image.dir = ti,out.dir="foo",save=TRUE))),"does not exist")
  
  expect_error(invisible(capture.output(kin.free(image.dir =ti,save=TRUE))),"not specified")
  
  expect_error(invisible(capture.output( kin.free(image.dir = ti,frames=1:3,save=FALSE))),"out of range")
  
  expect_error(invisible(capture.output( kin.free(image.dir = ti,frames=1,save=FALSE))),"number of frames must be")
  
  expect_error(invisible(capture.output( kin.free(image.dir =ti , thr="foo",save=FALSE))),"must be set to")
  expect_error(invisible(capture.output( kin.free(image.dir =ti,ml.smooth = list(1,"foo"),save=FALSE))),"must contain 'loess'")
  
  expect_error(invisible(capture.output( kin.free(image.dir =ti,ml.smooth = list(1,2),save=FALSE))),"'ml.smooth' must be a list of length 2")
  
  expect_error(invisible(capture.output( kin.free(image.dir =ti ,search.for="foo",save=FALSE))),"must be set to")
  
  expect_error(invisible(capture.output( kin.free(image.dir=ti,save=TRUE,out.qual=1.1,out.dir=tp))),"'out.qual' must be >=0 and <=1")
  
  
  
  expect_error(invisible(capture.output( kin.free(image.dir=ti,save=FALSE,ml.smooth=list("spline",1.1)))),"'smooth' must <1")
  

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  
  .Platform$OS.type 
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    cor.n <- 2L
  } else {
    # use all cores in devtools::test()
    cor.n <- parallel::detectCores()
  }
  
  n.cor <- parallel::detectCores()
  if (.Platform$OS.type == "windows" && n.cor > 1) {
    skip("mc.cores > 1 is not supported on Windows.")
  }
  

  invisible(capture.output( kin.yp <- kin.free(image.dir = ti,save = TRUE,out.dir =tp,red=0.5,par=TRUE,cores.n = cor.n)))
  
  expect_identical(kin.y,kin.yp)
  

  
  dir.create(paste0(tempdir(),"/test_images2"))
  
  
  expect_error(invisible(capture.output( kin.free(image.dir = paste0(tempdir(),"/test_images2"),save=FALSE))),"no images in image.dir")
  
  unlink(paste0(tempdir(),"/test_images2"),recursive = TRUE)
  
  unlink(ti,recursive = TRUE)
  unlink(tp,recursive = TRUE)
  
  
})



test_that("fin.kin works fine", {
  
  y <- EBImage::readImage(system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter"))
  t <-tempdir()
  
  ti <-paste0(tempdir(),"/images")
  tp <- paste0(tempdir(),"/processed_images")
  dir.create(ti)
  dir.create(tp)
  
  EBImage::writeImage(y,paste0(ti,"/sunfish001.jpg"),type = "jpeg")
  
  invisible(capture.output( kin.y <- kin.simple(image.dir = ti,save = TRUE,out.dir =tp)))
  
  fin.pos <- c(0.2,0.55)
  fin.y <- fin.kin(kin=kin.y,fin.pos = fin.pos,smooth.n=1,ml.smooth=0.3)
  
  expect_is(fin.y,"list")
  expect_named(fin.y,c("cont","fin","fin.pts","comp","midline","amp","bl"))
  expect_type(fin.y$cont$y,type = "double")
  expect_type(fin.y$comp$y,type = "double")
  expect_type(fin.y$fin.pts$y,type = "double")
  expect_type(fin.y$fin$y,type = "double")
  expect_type(fin.y$bl$bl,type = "double")
  expect_type(fin.y$amp$amp2[1],type = "double")
  
  expect_error(fin.kin(kin.y$cont))
  expect_error(fin.kin(kin.y,fin.pos=0.1))
  expect_error(fin.kin(kin.y,fin.pos=NULL))
})




