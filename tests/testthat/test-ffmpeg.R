context("ffmpeg functions")

test_that("vid.to.images works", {
  
  
  v <- system.file("extdata/vid", "sunfish_BCF.avi", package = "trackter")
  file.copy(v,tempdir())
  vid.to.images(vid.path = paste0(tempdir(),"/sunfish_BCF.avi")  )
  expect_true(dir.exists(paste0(tempdir(),"/images")))
  expect_true(length(list.files(paste0(tempdir(),"/images")))==2)
  expect_error(vid.to.images(vid.path = NULL,silent=T))
  expect_error(vid.to.images(vid.path = paste0(tempdir(),"/foo.avi"),silent=T))

  unlink(paste0(tempdir(),"/images"),recursive = T)
  unlink(paste0(tempdir(),"/sunfish_BCF.avi"))
})

test_that("images.to.video works", {
  if(dir.exists(paste0(tempdir(),"/sunfish"))) unlink(paste0(tempdir(),"/sunfish"),recursive = T)
  dir.create(paste0(tempdir(),"/sunfish"))
 v <- system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter")
 
file.copy(v,paste0(tempdir(),"/sunfish/img_001.jpg"))
file.copy(v,paste0(tempdir(),"/sunfish/img_002.jpg"))
  
images.to.video(image.dir = paste0(tempdir(),"/sunfish"),vid.name = paste0(tempdir(),"/test.mp4"), silent = T)
 

  
  expect_true(file.exists(paste0(tempdir(),"/test.mp4")))
  expect_true(file.size(paste0(tempdir(),"/test.mp4"))>10)
  
  
  unlink(paste0(tempdir(),"/test.mp4"))

  unlink(paste0(tempdir(),"/sunfish"),recursive = T)
  
 })

test_that("vid.to.images2 works", {
  v <- system.file("extdata/vid", "sunfish_BCF.avi", package = "trackter")
  file.copy(v,tempdir())
  vid.to.images2(vid.path = paste0(tempdir(),"/sunfish_BCF.avi"))  
  expect_true(dir.exists(paste0(tempdir(),"/images")))
  expect_true(length(list.files(paste0(tempdir(),"/images")))==2)
  img1 <- EBImage::readImage(paste0(tempdir(),"/images/",list.files(paste0(tempdir(),"/images"))[1]))
  expect_error(vid.to.images2(vid.path = NULL))
  expect_error(vid.to.images2(vid.path = paste0(tempdir(),"/foo.avi")))
  
  vid.to.images2(vid.path = paste0(tempdir(),"/sunfish_BCF.avi"),filt = " -vf scale=200:-1 ") 
  img2 <- EBImage::readImage(paste0(tempdir(),"/images/",list.files(paste0(tempdir(),"/images"))[1]))
  expect_true(dim(img1)[1]>dim(img2)[1]) #images with scaling filter are smaller
  
  unlink(paste0(tempdir(),"/images"),recursive = T)
  unlink(paste0(tempdir(),"/sunfish_BCF.avi"))
  
})

test_that("images.to.video2 works", {
  
  if(dir.exists(paste0(tempdir(),"/sunfish"))) unlink(paste0(tempdir(),"/sunfish"),recursive = T)
  dir.create(paste0(tempdir(),"/sunfish"))
  v <- system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter")
  
  file.copy(v,paste0(tempdir(),"/sunfish/img_001.jpg"))
  file.copy(v,paste0(tempdir(),"/sunfish/img_002.jpg"))
  
  images.to.video2(image.dir = paste0(tempdir(),"/sunfish"),vid.name = paste0(tempdir(),"/test"), silent = T)
  
  
  expect_true(file.exists(paste0(tempdir(),"/test.avi")))
  expect_true(file.size(paste0(tempdir(),"/test.avi"))>10)
  
  f.s1 <- file.size(paste0(tempdir(),"/test.avi"))
  images.to.video2(image.dir = paste0(tempdir(),"/sunfish"),vid.name = paste0(tempdir(),"/test"),silent = T,raw = F,overwrite = T)
  f.s2 <- file.size(paste0(tempdir(),"/test_red.mp4"))
  
  expect_true(f.s1>f.s2) #is compressed file <raw file
  
  images.to.video2(image.dir = paste0(tempdir(),"/sunfish"),vid.name = paste0(tempdir(),"/test"),silent = F,raw = F,filt = " -s 120x80 ",overwrite = T) #reduce scale
  f.s3 <- file.size(paste0(tempdir(),"/test_red.mp4"))
  
  expect_true(f.s1>f.s3)
  expect_true(f.s2>f.s3)
  
  expect_message(images.to.video2(paste0(tempdir(),"/sunfish"),vid.name = paste0(tempdir(),"/test"),silent = T,raw = T,overwrite = T),"video saved to")
  
  unlink(paste0(tempdir(),"/test.mp4"))
  unlink(paste0(tempdir(),"/sunfish"),recursive = T)
  
  
})
