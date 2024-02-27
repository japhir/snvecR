test_that("test_solution() works", {
  ZB18a_head <- structure(list(t = c(0, -146100),
                               aa = c(0.999999570867702, 1.00000267419287),
                               ee = c(0.0167054504495442, 0.016854305837952),
                               inc = c(7.15495578299901, 7.14593294487681),
                               lph = c(27.323464691619, 26.1208300138006),
                               lan = c(179.999248773572, -179.58527791053),
                               arp = c(-152.675784081952, -154.293892075669),
                               mna = c(-2.45433424316846, 1.26740886389077)),
                          row.names = c(NA, -2L),
                          class = c("tbl_df", "tbl", "data.frame"))

  expect_error(prepare_solution(data.frame(x = 1, y = 2), quiet = TRUE))
  expect_equal(colnames(prepare_solution(ZB18a_head, quiet = TRUE)),
               c("t", "t_kyr", #"age",
                 "aa", "ee", "inc", "lph", "lan", "arp", "mna",
                 "lphu", "lanu",
                 "hh", "kk", "pp", "qq", "cc", "dd",
                 "nnx", "nny", "nnz"))
})
