test_that("edge removal renumbers only the receiving model", {
  s <- data.frame(From_Emulator = c("Global", "u", "Global", "u"),
    To_Emulator = c("a", "a", "b", "b"), From_Output = 1L, To_Input = c(1L, 2L, 1L, 2L))
  z <- dgpsi:::.prunlgp_structure(s, list(a = 1L, b = 1:2))
  expect_false(any(z$From_Emulator == "u" & z$To_Emulator == "a"))
  expect_true(any(z$From_Emulator == "u" & z$To_Emulator == "b"))
  expect_identical(z$To_Input[z$To_Emulator == "b"], 1:2)
})

test_that("grouped deletion also removes the corresponding origin input", {
  s <- data.frame(From_Emulator = c("Global", "Global", "Global", "u"),
    To_Emulator = c("a0", "a0", "a1", "a1"), From_Output = c(1L,2L,1L,1L),
    To_Input = c(1L,2L,1L,2L))
  z <- dgpsi:::.prunlgp_structure(s, list(a0 = 1L, a1 = 1L))
  expect_equal(nrow(z), 2)
  expect_true(all(z$From_Output == 1L))
})

# Native integration tests are opt-in so CRAN checks never install Python.
native_fixture <- function() {
  skip_if(Sys.getenv("DGPSI_RUN_PRUNLGP_TESTS") != "true")
  init_py(dgpsi_ver = Sys.getenv("DGPSI_TEST_BACKEND", "dev"), verb = FALSE)
  set_thread_num(1L); set_seed(41L)
  x <- matrix(seq(-1, 1, length.out = 18), ncol = 1)
  u <- sin(2*x)
  components <- list(
    gp(x, u, id = "u", name = "sexp", verb = FALSE),
    gp(cbind(x, u), cos(x), id = "a", name = "sexp", verb = FALSE),
    gp(cbind(x, u), u^2, id = "b", name = "sexp", verb = FALSE))
  s <- data.frame(From_Emulator = c("Global","Global","u","Global","u"),
    To_Emulator = c("u","a","a","b","b"), From_Output = 1L,
    To_Input = c(1L,1L,2L,1L,2L))
  list(models = components, network = lgp(s, components, verb = FALSE),
    x = matrix(seq(-.93,.93,length.out=11), ncol=1),
    refit = function(model, keep, seed) {
      dgpsi::set_seed(seed)
      dgpsi::gp(model$data$X[,keep,drop=FALSE],model$data$Y,
        id=model$id,name="sexp",verb=FALSE)
    })
}

test_that("native search preserves the baseline, all outputs and reproducibility", {
  f <- native_fixture(); y <- list(a=cos(f$x),b=sin(2*f$x)^2)
  set.seed(621); rng <- .Random.seed
  a <- prunlgp(f$network,f$models,f$x,y,f$refit,seed=42,max_evals=12,verb=FALSE)
  expect_identical(.Random.seed,rng)
  b <- prunlgp(f$network,f$models,f$x,y,f$refit,seed=42,max_evals=12,verb=FALSE)
  expect_equal(a$score,b$score,tolerance=1e-12)
  expect_equal(a$structure,b$structure)
  expect_lte(a$score,a$initial_score+1e-10)
  expect_setequal(a$model$specs$metadata$Emulator,c("u","a","b"))
  expect_equal(sum(a$history$selected),1)
  expect_equal(a$score,min(a$history$loss[a$history$admissible]),tolerance=1e-10)
  expect_s3_class(predict(a,f$x,full_layer=TRUE),"lgp")
  c <- prunlgp(f$network,f$models,f$x,y,f$refit,seed=42,protect=0,verb=FALSE)
  expect_true(all(c$rmse <= c$initial_rmse+1e-10))
})

test_that("zero internal links predict with native GP moments", {
  f <- native_fixture()
  models <- lapply(f$models,function(m) f$refit(m,1L,5L))
  s <- data.frame(From_Emulator="Global",To_Emulator=c("u","a","b"),From_Output=1L,To_Input=1L)
  m <- lgp(s,models,verb=FALSE)
  a <- predict(m,f$x,full_layer=TRUE)$results
  for(i in seq_along(models)) {
    b <- predict(models[[i]],f$x)$results
    expect_equal(a$mean[[1]][[models[[i]]$id]],b$mean,tolerance=1e-9)
    expect_equal(a$var[[1]][[models[[i]]$id]],b$var,tolerance=1e-9)
  }
  expect_length(predict(m,f$x,full_layer=FALSE)$results$mean,3)
  draw <- predict(m,f$x,method="sampling",sample_size=4,full_layer=TRUE)$results
  expect_equal(dim(draw[[1]]$u$output1),c(nrow(f$x),4L))
})

test_that("failed refits are reported and cannot replace the initial model", {
  f <- native_fixture()
  bad <- function(model,keep,seed) stop("deliberate refit failure")
  z <- prunlgp(f$network,f$models,f$x,list(a=cos(f$x)),bad,max_evals=4,verb=FALSE)
  expect_equal(z$score,z$initial_score)
  expect_true(any(grepl("deliberate refit failure",z$history$error)))
  expect_lte(nrow(z$history),4)
  expect_error(prunlgp(f$network,f$models,f$x,list(a=cos(f$x)),bad,scale=0),"scale")
})

test_that("multi-output components are rejected explicitly", {
  f <- native_fixture()
  f$models[[2]]$data$Y <- cbind(1:18,1:18)
  expect_error(prunlgp(f$network,f$models,f$x,list(a=cos(f$x)),f$refit),"scalar")
})

test_that("parallel candidates match serial candidates", {
  f <- native_fixture(); y <- list(a=cos(f$x),b=sin(2*f$x)^2)
  a <- prunlgp(f$network,f$models,f$x,y,f$refit,seed=43,verb=FALSE)
  b <- prunlgp(f$network,f$models,f$x,y,f$refit,seed=43,workers=2,
    worker_init=function() dgpsi::init_py(dgpsi_ver="dev",verb=FALSE),verb=FALSE)
  expect_equal(a$score,b$score,tolerance=1e-10)
  expect_equal(a$structure,b$structure)
})

test_that("retained upstream predictions affect only connected receivers", {
  f <- native_fixture()
  a <- f$refit(f$models[[2]],1L,7L)
  original <- list(f$models[[1]],a,f$models[[3]])
  moved <- original
  moved[[1]] <- gp(original[[1]]$data$X,original[[1]]$data$Y+.2,
    id="u",name="sexp",verb=FALSE)
  s <- f$network$specs$struc
  s <- s[!(s$From_Emulator=="u" & s$To_Emulator=="a"),]
  p <- predict(lgp(s,original,verb=FALSE),f$x,full_layer=TRUE)$results
  q <- predict(lgp(s,moved,verb=FALSE),f$x,full_layer=TRUE)$results
  pm <- do.call(c,unname(p$mean));qm<-do.call(c,unname(q$mean))
  expect_equal(pm$a,qm$a,tolerance=1e-12)
  expect_gt(max(abs(pm$b-qm$b)),1e-5)
})

test_that("DGP moments and samples survive a one-layer serialisation", {
  f <- native_fixture();set_seed(61)
  m <- dgp(f$models[[1]]$data$X,f$models[[1]]$data$Y,
    id="u",name="sexp",depth=2,node=2,N=3,B=3,verb=FALSE)
  s <- data.frame(From_Emulator="Global",To_Emulator="u",From_Output=1L,To_Input=1L)
  network <- lgp(s,list(m),B=3,verb=FALSE)
  restored <- deserialize(serialize(network,light=FALSE))
  p<-predict(network,f$x,full_layer=TRUE)$results
  q<-predict(restored,f$x,full_layer=TRUE)$results
  expect_equal(p$mean,q$mean,tolerance=1e-12)
  expect_equal(p$var,q$var,tolerance=1e-12)
  expect_true(all(is.finite(q$var[[1]]$u) & q$var[[1]]$u>=0))
  sampled<-predict(restored,f$x,method="sampling",sample_size=4,full_layer=FALSE)$results
  expect_equal(dim(sampled$u$output1),c(nrow(f$x),12L))
})

test_that("fitting groups remove an input from every copy", {
  f<-native_fixture()
  a0<-set_id(f$models[[2]],"a0");a1<-set_id(a0,"a1")
  s<-data.frame(From_Emulator=c("Global","Global","Global","Global","u"),
    To_Emulator=c("u","a0","a0","a1","a1"),
    From_Output=c(1L,1L,2L,1L,1L),To_Input=c(1L,1L,2L,1L,2L))
  ms<-list(f$models[[1]],a0,a1);net<-lgp(s,ms,verb=FALSE)
  result<-prunlgp(net,ms,cbind(f$x,sin(2*f$x)),
    list(a0=cos(f$x),a1=cos(f$x)),f$refit,groups=c(u="u",a0="a",a1="a"),
    objective=function(rmse,scale,weights) 0,verb=FALSE)
  expect_identical(result$kept_inputs$a0,1L)
  expect_identical(result$kept_inputs$a1,1L)
  expect_true(all(result$structure$From_Emulator=="Global"))
  expect_equal(nrow(result$structure),3L)
})

test_that("refits cannot silently change training targets or retained columns", {
  f<-native_fixture()
  wrong<-function(model,keep,seed) dgpsi::gp(model$data$X[,keep,drop=FALSE],
    model$data$Y+1,id=model$id,verb=FALSE)
  z<-prunlgp(f$network,f$models,f$x,list(a=cos(f$x)),wrong,verb=FALSE)
  expect_equal(z$score,z$initial_score)
  expect_true(any(grepl("preserve the ID/response",z$history$error)))
})
