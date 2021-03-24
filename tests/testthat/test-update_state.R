context("state update")

## With no births and deaths, and a single patch
## and with a fixed seed, check that the numbers work out.
test_that("update state with no births and deaths and one patch", {

    state_in <- state(
        s_patches = 1000,
        e_patches = 0,
        i_patches = 10,
        r_patches = 0,
        birth_rates = 0,
        death_rates = 0,
        transmission_rates = 0.2,
        infection_rates = 0.1,
        recovery_rates = 0.1,
        movement_rate = matrix(1, nrow = 1, ncol = 1)
    )
    set.seed(42)
    state_out <- update_state(state_in, 1)
    patch <- state_out[["patches"]][[1]]

    expect_equal(patch$susceptible, 996)
    expect_equal(patch$exposed, 4)
    expect_equal(patch$infected, 8)
    expect_equal(patch$recovered, 2)
})

# Temporarily mask test: the expected answers are not correct
# test_that("update state with no births and deaths and 2 patches", {
# 
#     state_in <- state(
#         s_patches = c(1000, 1000),
#         e_patches = c(0, 0),
#         i_patches = c(10, 10),
#         r_patches = c(0, 0),
#         birth_rates = c(0, 0),
#         death_rates = c(0, 0),
#         transmission_rates = c(0.2, 0.2),
#         infection_rates = c(0.1, 0.1),
#         recovery_rates = c(0.1, 0.1),
#         movement_rate = matrix(
#             c(1.61, 0.22, 0.36, 1.2), nrow = 2, ncol = 2, byrow = TRUE
#         )
#     )
#     set.seed(42)
#     state_out <- update_state(state_in, 1, movement_type = "rate")
# 
#     patch1 <- state_out[["patches"]][[1]]
# 
#     expect_equal(patch1$susceptible, 1124)
#     expect_equal(patch1$exposed, 2)
#     expect_equal(patch1$infected, 7)
#     expect_equal(patch1$recovered, 1)
# 
# 
#     patch2 <- state_out[["patches"]][[2]]
# 
#     expect_equal(patch2$susceptible, 872)
#     expect_equal(patch2$exposed, 2)
#     expect_equal(patch2$infected, 10)
#     expect_equal(patch2$recovered, 2)
# 
# 
# })
