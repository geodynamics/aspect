(parameters:Checkpointing)=
# Checkpointing


## **Subsection:** Checkpointing


(parameters:Checkpointing/Steps_20between_20checkpoint)=
### __Parameter name:__ Steps between checkpoint
**Default value:** 0

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The number of timesteps between performing checkpoints. If 0 and time between checkpoint is not specified, checkpointing will not be performed. Units: None.

(parameters:Checkpointing/Time_20between_20checkpoint)=
### __Parameter name:__ Time between checkpoint
**Default value:** 0

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** The wall time between performing checkpoints. If 0, will use the checkpoint step frequency instead. Units: Seconds.
