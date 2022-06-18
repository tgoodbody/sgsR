set.seed(2022)

e <- existing[sample(nrow(existing),20),]

o <- sample_clhs(mraster = mraster, existing = e, nSamp = 50, details = TRUE)

test_that("Input classes", {
  expect_error(sample_clhs(mraster = "mraster", existing = e, nSamp = 5), "'mraster' must be type SpatRaster.")
  expect_error(sample_clhs(mraster = mraster, nSamp = "A", existing = e), "'nSamp' must be type numeric.")
  expect_error(sample_clhs(mraster = mraster, existing = "existing", nSamp = 5), "'existing' must be a data.frame or sf object.")
  expect_error(sample_clhs(mraster = mraster, existing = e, iter = "A", nSamp = 50), "'iter' must be type numeric.")
  expect_error(sample_clhs(mraster = mraster, nSamp = 100, existing = e, plot = 1), "'plot' must be type logical.")
  expect_error(sample_clhs(mraster = mraster, nSamp = 100, existing = e, details = "A"), "'details' must be type logical.")
  expect_error(sample_clhs(mraster = mraster, nSamp = 100, existing = e, cost = TRUE), "'cost' must be either type numeric or character.")
  expect_error(sample_clhs(mraster = mraster, nSamp = 100, existing = e, cost = 13), "'cost' index doest not exist within 'mraster'.")
  expect_error(sample_clhs(mraster = mraster, nSamp = 100, existing = data.frame(x = c(1,2,3), y = c(1,2,3))), "'existing' only extracts NA values. Ensure that 'existing' overlaps with 'mraster'.")
  expect_error(sample_clhs(mraster = mraster, nSamp = 10, filename = 56), "'filename' must be a file path character string.")
  expect_error(sample_clhs(mraster = mraster, nSamp = 10, filename = file.path(tempdir(), "temp.shp"), overwrite = "A"), "'overwrite' must be type logical.")
})

test_that("Total outputs", {
  expect_equal(nrow(o$samples), 50L)
  expect_equal(ncol(o$samples), 5L)
  expect_equal(length(o$clhs), 50L)
  expect_s3_class(o$samples,"sf")
  
  expect_message(sample_clhs(mraster = mraster, existing = e, nSamp = 50, cost = 1), "Using `zq90` as sampling constraint.")
  expect_message(sample_clhs(mraster = mraster, existing = e, nSamp = 50, cost = "pzabove2"), "Using `pzabove2` as sampling constraint.")
  
})

test_that("Messages", {
  expect_message(sample_clhs(mraster = mraster, nSamp = 320, existing = existing.df.n.xy.lc), "Column coordinates names for 'existing' are lowercase - converting to uppercase.")
  expect_message(sample_clhs(mraster = mraster, nSamp = 320, existing = existingna), "16 samples are located where metric values are NA.")
  expect_message(sample_clhs(mraster = mraster, nSamp = 20, access = access, buff_inner = 50, buff_outer = 200), "An access layer has been provided. An internal buffer of 50 m and an external buffer of 200 m have been applied.")
  expect_message(sample_clhs(mraster = mraster, nSamp = 20, access = access, buff_outer = 200), "An access layer has been provided. An external buffer of 200 m have been applied.")
  expect_message(sample_clhs(mraster = mraster, nSamp = 20, filename = file.path(tempdir(), "temp.shp") , overwrite = TRUE), "Output samples written to disc.")
})  


test_that("Errors", {
  expect_error(sample_clhs(mraster = mraster, nSamp = 320, cost = "A"),"No layer named 'A' exists in 'mraster'.")
  expect_error(sample_clhs(mraster = mraster, existing = existing, nSamp = 20),"nSamp must be > than number of existing samples.")
  expect_error(sample_clhs(mraster = mraster, nSamp = 100, existing = data.frame(r = NULL, g = NULL)), "'existing' must have columns named 'X' and 'Y'.")
})  

test_that("df input", {
  expect_equal(nrow(sample_clhs(mraster = mraster, existing = existing.df.n.xy, nSamp = 320)),320L)
})  

