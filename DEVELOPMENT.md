# Testing with the `vdiffr` package

`vdiffr::expect_doppelganger()` statements create a snapshot of plots the first time they are run. 
Subsequent runs will compare the current plot against the snapshot.
If the current plot is different from the snapshot, the test will fail.
If you intentionally change the plot, you can use `testthat::snapshot_review()` to overwrite the snapshot with the current plot.
