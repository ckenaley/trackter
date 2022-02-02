context("ffmpeg functions")



test_that("vid.to.images works", {
  
  v <- system.file("extdata/vid", "sunfish_BCF.avi", package = "trackter")

  dir.create(paste0(tempdir(),"/test_images"))
  
  vid.to.images(vid.path = v,out.dir =paste0(tempdir(),"/test_images"))
  
  expect_true(length(list.files(paste0(tempdir(),"/test_images")))==2)
  expect_error(vid.to.images(vid.path = paste0(tempdir(),"/sunfish_BCF.avi"),out.dir = NULL),"'out.dir' not specified")
  expect_error(vid.to.images(vid.path = paste0(tempdir(),"/sunfish_BCF.avi"),out.dir =paste0(tempdir(),"/foo") ),"does not exist")
  expect_error(vid.to.images(vid.path = paste0(tempdir(),"/foo.avi"),out.dir=paste0(tempdir(),"/test_images")),"does not exist")

  unlink(paste0(tempdir(),"/test_images"),recursive = TRUE)
  unlink(paste0(tempdir(),"/sunfish_BCF.avi"))
})


test_that("images.to.video works", {

  if(dir.exists(paste0(tempdir(),"/sunfish"))) unlink(paste0(tempdir(),"/sunfish"),recursive = TRUE)
  
  dir.create(paste0(tempdir(),"/sunfish"))
  v <- system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter")
  
  file.copy(v,paste0(tempdir(),"/sunfish/img_001.jpg"))
  file.copy(v,paste0(tempdir(),"/sunfish/img_002.jpg"))
  
  images.to.video(image.dir = paste0(tempdir(),"/sunfish"),vid.name = "test.mp4",out.dir=tempdir())
  
  
  
  expect_true(file.exists(paste0(tempdir(),"/test.mp4")))
  expect_true(file.size(paste0(tempdir(),"/test.mp4"))>10)
  
  
  expect_error(images.to.video(image.dir = paste0(tempdir(),"/sunfish"),out.dir = NULL),"not specified")
  
  expect_error(images.to.video(image.dir = paste0(tempdir(),"/sunfish"),out.dir = "/foo"),"does not exist")
  
  expect_error(images.to.video(image.dir = paste0(tempdir(),"/foo"),out.dir = tempdir()),"does not exist")
  
  expect_error(images.to.video(image.dir = paste0(tempdir(),"/sunfish"),vid.name = "test.mp4",out.dir=tempdir(),silent = TRUE,overwrite = FALSE),"'overwrite' must be 'TRUE'")
  
  unlink(paste0(tempdir(),"/test.mp4"))
  
  unlink(paste0(tempdir(),"/sunfish"),recursive = TRUE)
  
 })
