context("state update")

## With no births and deaths, and a single patch
## and with a fixed seed, check that the numbers work out.
test_that("update state works with no births and deaths", {

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

    state_out <- update_state(state_in, 1, seed = 42)
    patch <- state_out[["patches"]][[1]]

    expect_equal(patch$susceptible, 996)
    expect_equal(patch$exposed, 4)
    expect_equal(patch$infected, 8)
    expect_equal(patch$recovered, 2)
})
