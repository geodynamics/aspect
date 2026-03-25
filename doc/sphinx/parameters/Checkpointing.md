(parameters:Checkpointing)=
# Checkpointing


## **Subsection:** Checkpointing


::::{dropdown} __Parameter:__ {ref}`Number of checkpoints to keep<parameters:Checkpointing/Number_20of_20checkpoints_20to_20keep>`
:name: parameters:Checkpointing/Number_20of_20checkpoints_20to_20keep
**Default value:** 3

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** The number of checkpoint slots to rotate through in the output/restart directory. Units: None.
::::

::::{dropdown} __Parameter:__ {ref}`Resume checkpoint<parameters:Checkpointing/Resume_20checkpoint>`
:name: parameters:Checkpointing/Resume_20checkpoint
**Default value:** 0

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** If nonzero, resume from the checkpoint with this slot id inside output/restart. This overrides Resume time and the last good checkpoint selection. Units: None.
::::

::::{dropdown} __Parameter:__ {ref}`Resume time<parameters:Checkpointing/Resume_20time>`
:name: parameters:Checkpointing/Resume_20time
**Default value:** -1

**Pattern:** [Double -1...MAX_DOUBLE (inclusive)]

**Documentation:** If non-negative, resume from the checkpoint whose saved model time is closest to this value. This option is ignored when Resume checkpoint is set. Units: Years if the &rsquo;Use years instead of seconds&rsquo; parameter is set; seconds otherwise.
::::

::::{dropdown} __Parameter:__ {ref}`Steps between checkpoint<parameters:Checkpointing/Steps_20between_20checkpoint>`
:name: parameters:Checkpointing/Steps_20between_20checkpoint
**Default value:** 0

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The number of timesteps between performing checkpoints. If 0 and time between checkpoint is not specified, checkpointing will not be performed. Units: None.
::::

::::{dropdown} __Parameter:__ {ref}`Time between checkpoint<parameters:Checkpointing/Time_20between_20checkpoint>`
:name: parameters:Checkpointing/Time_20between_20checkpoint
**Default value:** 0

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The wall time between performing checkpoints. If 0, will use the checkpoint step frequency instead. Units: Seconds.
::::
