# Test properties

ASPECT's test parameter files can contain
certain special markers that indicate how this test should be executed.
Placing such a marker at any point in the test parameter file will indicate to
the test system what to do with this test. The following properties are
currently supported:

-   "`# MPI: x`": Indicates that this test should be executed with
    `x` MPI ranks. This is useful to test certain features in parallel.

-   "`# EXPECT FAILURE`": Indicates that this test should fail.
    The test passes if the run fails, and the test fails, if the run passes.
    This is useful to check if error checks (e.g., assertions) are correctly
    executed.

-   "`# DEPENDS-ON: x`": Indicates that this test depends on the
    completion of another test. The test will only start after the other test
    has finished (successful or not). If this test is executed without the
    other test being executed as well, this test will fail. This feature is
    useful if a test reuses output of another test.
