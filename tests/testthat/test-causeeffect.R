extract_text_labels <- function(plot) {
  built <- ggplot2::ggplot_build(plot)
  text_layers <- built$data[vapply(built$data, function(df) "label" %in% names(df), logical(1))]

  unlist(lapply(text_layers, function(df) as.character(df$label)), use.names = FALSE)
}

test_that("causeEffectDiagram renders each odd-count branch title and item once", {
  cMan <- c("Reservations staffs undertrained", "Front office staff unfocused")
  cMethods <- c(
    "No credit card blocking",
    "No advance payments\nfor confirmed bookings",
    "Reconfirmation of guest departures",
    "Unpredicted top account\nVIP bookings"
  )
  cMachines <- c(
    "System allows reservation for overbookings due to wait list",
    "Internet issues",
    "OTA portal issues",
    "System errors"
  )

  plot <- causeEffectDiagram(
    cause = list(Man = cMan, Methods = cMethods, Machines = cMachines),
    effect = "Overbooking leading to guest bouncing"
  )

  expect_s3_class(plot, "ggplot")
  expect_no_error(ggplot2::ggplot_build(plot))

  labels <- extract_text_labels(plot)
  branch_counts <- vapply(names(list(Man = cMan, Methods = cMethods, Machines = cMachines)), function(label) {
    sum(labels == label)
  }, integer(1))

  expect_equal(unname(branch_counts), c(1L, 1L, 1L))
  expect_true(all(cMachines %in% labels))
  expect_true(all(cMethods %in% labels))
  expect_true(all(cMan %in% labels))
})

test_that("causeEffectDiagram builds the documentation example (even case)", {
  plot <- causeEffectDiagram(
    cause = list(
      Measurements = c("Micrometers", "Microscopes", "Inspectors"),
      Materials = c("Alloys", "Lubricants", "Suppliers"),
      Personnel = c("Shifts", "Supervisors", "Training", "Operators"),
      Environment = c("Condensation", "Moisture"),
      Methods = c("Brake", "Engager", "Angle"),
      Machines = c("Speed", "Lathes", "Bits", "Sockets")
    ),
    effect = "Surface Flaws"
  )

  expect_s3_class(plot, "ggplot")
  expect_no_error(ggplot2::ggplot_build(plot))
})
