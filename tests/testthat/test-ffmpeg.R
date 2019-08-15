context("ffmpeg functions")

test_that("vid.to.images works", {
  
  
  v <- system.file("extdata/vid", "sunfish_BCF.avi", package = "trackter")
  file.copy(v,getwd())
  vid.to.images(vid.path = "sunfish_BCF.avi")  
  expect_true(dir.exists("images"))
  expect_true(length(list.files("images"))==2)
  expect_error(vid.to.images(vid.path = NULL,silent=T))
  expect_error(vid.to.images(vid.path = "foo.avi",silent=T))

  unlink("images",recursive = T)
  unlink("sunfish_BCF.avi")
})

test_that("images.to.video works", {
  if(dir.exists("sunfish")) unlink("sunfish",recursive = T)
  dir.create("sunfish")
 v <- system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter")
file.copy(v,paste0(getwd(),"/sunfish/img_001.jpg"))
file.copy(v,paste0(getwd(),"/sunfish/img_002.jpg"))
  
 images.to.video(image.dir = paste0(getwd(),"/sunfish"),vid.name = "test.mp4",silent = T)  
  
  expect_true(file.exists("test.mp4"))
  expect_true(file.size("test.mp4")>10)
  unlink("test.mp4")
  unlink("sunfisf",recursive=T)
  
 })

test_that("vid.to.images2 works", {
  
  
  v <- system.file("extdata/vid", "sunfish_BCF.avi", package = "trackter")
  file.copy(v,getwd())
  vid.to.images2(vid.path = "sunfish_BCF.avi")  
  expect_true(dir.exists("images"))
  expect_true(length(list.files("images"))==2)
  img1 <- EBImage::readImage(paste0("images/",list.files("images")[1]))
  
 

  expect_error(vid.to.images2(vid.path = NULL))
  expect_error(vid.to.images2(vid.path = "foo.avi"))
  
  vid.to.images2(vid.path = "sunfish_BCF.avi",filt = " -vf scale=200:-1 ") 
  img2 <- EBImage::readImage(paste0("images/",list.files("images")[1]))
  expect_true(dim(img1)[1]>dim(img2)[1]) #images with scaling filter are smaller
  
  unlink("images",recursive = T)
  unlink("sunfish_BCF.avi")
  
})

test_that("images.to.video2 works", {
  if(dir.exists("sunfish")) unlink("sunfish",recursive = T)
  dir.create("sunfish")
  v <- system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter")
  file.copy(v,paste0(getwd(),"/sunfish/img_001.jpg"))
  file.copy(v,paste0(getwd(),"/sunfish/img_002.jpg"))
  
  images.to.video2(image.dir = paste0(getwd(),"/sunfish"),vid.name = "test",silent = T,raw = T)
  
  expect_true(file.exists("test.avi"))
  expect_true(file.size("test.avi")>10)
  f.s1 <- file.size("test.avi")
  images.to.video2(image.dir = paste0(getwd(),"/sunfish"),vid.name = "test",silent = T,raw = F)
  f.s2 <- file.size("test_red.mp4")
  
  expect_true(f.s1>f.s2) #is compressed file <raw file
  
  images.to.video2(image.dir = paste0(getwd(),"/sunfish"),vid.name = "test",silent = T,raw = F,filt = " -s 120x80 ") #reduce scale
  f.s3 <- file.size("test_red.mp4")

  expect_true(f.s1>f.s3)
  expect_true(f.s2>f.s3)
  
  expect_message(  images.to.video2(image.dir = paste0(getwd(),"/sunfish"),vid.name = "test",silent = T,raw = T),"video saved to")
  
  unlink("test_red.mp4")
  unlink("test.avi")
  unlink("sunfish",recursive=T)
  
})
